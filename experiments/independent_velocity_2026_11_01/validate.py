"""Independent dense checks for the independent-velocity block solvers.

References assemble the full linearized systems directly; they do not reuse
production increment formulas. Requires NumPy only for these algebraic checks.
"""
from __future__ import annotations

import argparse
from dataclasses import replace
import importlib
import json
from pathlib import Path

import numpy as np


class Checks:
    def __init__(self):
        self.count = 0
        self.maxima = {}

    def close(self, name, actual, expected, tolerance=3e-9):
        self.count += 1
        actual, expected = np.asarray(actual), np.asarray(expected)
        error = float(np.linalg.norm(actual-expected)/max(1.0, np.linalg.norm(expected)))
        self.maxima[name] = max(error, self.maxima.get(name, 0.0))
        assert np.isfinite(error) and error <= tolerance, (name, error, tolerance)

    def true(self, name, condition):
        self.count += 1
        assert condition, name


def finite_difference(function, point):
    value = np.asarray(function(point))
    derivative = np.empty((value.size, len(point)))
    for column in range(len(point)):
        delta = 2e-6*max(1.0, abs(point[column]))
        plus, minus = point.copy(), point.copy()
        plus[column] += delta
        minus[column] -= delta
        derivative[:, column] = (function(plus)-function(minus))/(2*delta)
    return derivative


def independent(point, xn, vn, h, mass, stiffness, length, gravity,
                epsilon, target, floor, beta, rigid=False):
    x, v, lam, rho, eta = point[:2], point[2:4], point[4], point[5], point[6]
    radius = np.linalg.norm(x)
    normal = x/radius
    constraint = radius-length
    normal_derivative = (np.eye(2)-np.outer(normal, normal))/radius
    alpha = 0.0 if rigid else 1.0/stiffness
    elastic = 0.0 if rigid else 0.5*stiffness*constraint**2
    gradient = -mass*gravity if rigid else stiffness*constraint*normal-mass*gravity
    hessian = (np.zeros((2, 2)) if rigid else stiffness*(
        np.outer(normal, normal)+constraint*normal_derivative))
    potential = elastic-mass*gravity@x
    physical = 0.5*mass*v@v+potential
    mismatch = v-(x-xn)/h
    residual = np.r_[
        mass*(x-xn-h*vn)/h**2-lam*normal-mass*gravity
        -mass*mismatch/(beta*h)+eta*gradient,
        mass*mismatch/beta+eta*mass*v,
        constraint+alpha*lam,
        epsilon*rho+eta*(target-floor),
        physical-(1-rho)*target-rho*floor,
    ]
    identity = np.eye(2)
    metric = np.block([
        [(1+1/beta)*mass*identity/h**2, -mass*identity/(beta*h)],
        [-mass*identity/(beta*h), mass*identity/beta],
    ])
    inverse = np.block([
        [h*h*identity/mass, h*identity/mass],
        [h*identity/mass, (1+beta)*identity/mass],
    ])
    a = np.r_[normal, 0.0, 0.0]
    gradient_energy = np.r_[gradient, mass*v]
    matrix = np.zeros((7, 7))
    matrix[:4, :4] = metric
    matrix[:4, 4] = -a
    matrix[4, :4] = a
    matrix[4, 4] = alpha
    matrix[:4, 6] = gradient_energy
    matrix[6, :4] = gradient_energy
    matrix[5, 5] = epsilon
    matrix[5, 6] = matrix[6, 5] = target-floor
    exact = matrix.copy()
    exact[:2, :2] += -lam*normal_derivative+eta*hessian
    exact[2:4, 2:4] += eta*mass*identity
    return dict(residual=residual, matrix=matrix, exact=exact, metric=metric,
                inverse=inverse, constraint=constraint, normal=normal, a=a,
                gradient=gradient, hessian=hessian, gradient_energy=gradient_energy,
                potential=potential, energy=physical, alpha=alpha, gap=target-floor,
                mismatch=mismatch)


def dualized_coupling(point, pi, xn, vn, h, mass, stiffness, length, gravity,
                      epsilon, target, floor, beta, rigid=False):
    """Explicit pi equations, including their current consistency residual."""
    base = independent(point, xn, vn, h, mass, stiffness, length, gravity,
                       epsilon, target, floor, 1.0 if beta == 0 else beta, rigid)
    x, v, lam, rho, eta = point[:2], point[2:4], point[4], point[5], point[6]
    n, gradient = base['normal'], base['gradient']
    residual = np.r_[
        mass*(x-xn-h*vn)/h**2-lam*n-mass*gravity-pi/h+eta*gradient,
        pi+eta*mass*v,
        v-(x-xn)/h-beta*pi/mass,
        base['residual'][4:],
    ]
    matrix = np.zeros((9, 9))
    matrix[:2, :2] = mass*np.eye(2)/h**2
    matrix[:2, 4:6] = matrix[4:6, :2] = -np.eye(2)/h
    matrix[2:4, 4:6] = matrix[4:6, 2:4] = np.eye(2)
    matrix[4:6, 4:6] = -beta*np.eye(2)/mass
    matrix[:2, 6] = -n
    matrix[6, :2] = n
    matrix[6, 6] = base['alpha']
    matrix[:4, 8] = matrix[8, :4] = base['gradient_energy']
    matrix[7, 7] = epsilon
    matrix[7, 8] = matrix[8, 7] = target-floor
    return matrix, residual


def cases(count=48):
    rng = np.random.default_rng(20261103)
    for index in range(count):
        rigid = index % 4 == 0
        mass = 10**rng.uniform(-0.3, 0.3)
        stiffness = 0.0 if rigid else 10**rng.uniform(1.3, 2.7)
        length = rng.uniform(0.7, 1.3)
        h = rng.uniform(0.03, 0.18)
        epsilon = 10**rng.uniform(-1, 2)
        beta = (.1, 1.0, 10.0)[index % 3]
        angle = rng.uniform(-np.pi, np.pi)
        x = rng.uniform(.8, 1.3)*length*np.array([np.cos(angle), np.sin(angle)])
        xn = x+rng.normal(size=2)*0.04
        vn = rng.normal(size=2)
        v = rng.normal(size=2)
        gravity = np.array([0.0, -2.1]) if index % 2 else np.zeros(2)
        target, floor = rng.uniform(2, 20), rng.uniform(-2, 0)
        point = np.r_[x, v, rng.normal()*10, rng.normal()*.2, rng.normal()*.3]
        args = (xn, vn, h, mass, stiffness, length, gravity, epsilon, target, floor, beta, rigid)
        yield point, args


def independent_checks(checks):
    max_factor = 0.0
    for point, args in cases():
        value = independent(point, *args)
        residual, matrix = value['residual'], value['matrix']
        checks.close('coupling metric inverse', value['metric']@value['inverse'], np.eye(4), 3e-12)
        checks.true('coupling metric positive definite', np.min(np.linalg.eigvalsh(value['metric'])) > 0)
        derivative = finite_difference(lambda q: independent(q, *args)['residual'], point)
        checks.close('full residual finite differences including dropped curvature',
                     derivative, value['exact'], 4e-8)
        checks.close('physical energy finite differences',
                     finite_difference(lambda q: np.array([independent(np.r_[q, point[4:]], *args)['energy']]), point[:4])[0],
                     value['gradient_energy'], 4e-8)
        xn, vn, h, mass, stiffness, length, gravity, epsilon, target, floor, beta, rigid = args
        def primal_objective(q):
            data = independent(np.r_[q, point[4:]], *args)
            inertia = 0.5*mass*np.dot(q[:2]-xn-h*vn, q[:2]-xn-h*vn)/h**2
            coupling = 0.5*mass*np.dot(data['mismatch'], data['mismatch'])/beta
            return np.array([inertia+data['potential']+coupling])
        expected_gradient = value['residual'][:4].copy()-point[6]*value['gradient_energy']
        expected_gradient[:2] += (point[4]+(0.0 if rigid else stiffness*value['constraint']))*value['normal']
        checks.close('original implicit plus coupling objective gradient finite differences',
                     finite_difference(primal_objective, point[:4])[0], expected_gradient, 4e-8)
        joint = np.linalg.solve(matrix, -residual)
        a, g, inverse, alpha = value['a'], value['gradient_energy'], value['inverse'], value['alpha']
        epsilon = args[7]
        diagonal_l = alpha+a@inverse@a
        diagonal_e = g@inverse@g+value['gap']**2/epsilon
        cross = a@inverse@g
        schur = np.array([[diagonal_l, -cross], [-cross, diagonal_e]])
        rhs = np.array([-residual[4]+a@inverse@residual[:4],
                        residual[6]-g@inverse@residual[:4]-value['gap']*residual[5]/epsilon])
        checks.close('independent multiplier Schur solve', np.linalg.solve(schur, rhs), joint[[4, 6]], 2e-8)
        factor = cross**2/(diagonal_l*diagonal_e)
        max_factor = max(max_factor, factor)
        checks.true('frozen GS contraction below one', 0 <= factor < 1)
        for order in ('EL', 'LE'):
            total = np.zeros(7)
            for sweep in range(10000):
                for block in order:
                    ids = [0, 1, 2, 3, 4] if block == 'L' else [0, 1, 2, 3, 5, 6]
                    current = residual+matrix@total
                    increment = np.linalg.solve(matrix[np.ix_(ids, ids)], -current[ids])
                    total[ids] += increment
                    checks.close('frozen block rows after immediate primal push',
                                 (residual+matrix@total)[ids], np.zeros(len(ids)), 3e-8)
                if np.linalg.norm(total-joint)/max(1, np.linalg.norm(joint)) < 3e-9:
                    break
            checks.close('frozen '+order+' versus joint solution', total, joint, 3e-8)

        # Arbitrary pi deliberately violates its constitutive relation. The
        # full consistency residual is required for elimination equivalence.
        pi = np.array([0.37, -1.2])
        dual_matrix, dual_residual = dualized_coupling(point, pi, *args)
        dual_delta = np.linalg.solve(dual_matrix, -dual_residual)
        checks.close('primal versus explicit dual coupling with residual',
                     dual_delta[[0, 1, 2, 3, 6, 7, 8]], joint, 2e-8)
        for block in ('L', 'E'):
            primal_ids = [0, 1, 2, 3, 4] if block == 'L' else [0, 1, 2, 3, 5, 6]
            dual_ids = [0, 1, 2, 3, 4, 5, 6] if block == 'L' else [0, 1, 2, 3, 4, 5, 7, 8]
            primal_delta = np.linalg.solve(matrix[np.ix_(primal_ids, primal_ids)], -residual[primal_ids])
            block_dual = np.linalg.solve(dual_matrix[np.ix_(dual_ids, dual_ids)], -dual_residual[dual_ids])
            checks.close('primal versus explicit dual '+block+' block',
                         np.r_[block_dual[:4], block_dual[6:]], primal_delta, 2e-8)
    return dict(arbitrary_states=48, max_frozen_GS_factor=max_factor)


def production_checks(checks):
    module = importlib.import_module('independent_solver')
    for point, args in cases():
        xn, vn, h, mass, stiffness, length, gravity, epsilon, target, floor, beta, rigid = args
        scene = module.Scene('arbitrary', mass, stiffness, length, gravity, xn, vn, rigid)
        ref = independent(point, *args)
        evaluation = module.evaluate_state(scene, xn, vn, h, point[:2], point[2:4],
                                           *point[4:], epsilon, beta, target=target, kstar=floor)
        checks.close('production analytical residuals',
                     np.r_[evaluation.r_x, evaluation.r_v, evaluation.r_lambda,
                           evaluation.r_rho, evaluation.r_eta], ref['residual'])
        checks.true('arbitrary position stationarity is nonzero', np.linalg.norm(evaluation.r_x) > 0)
        checks.true('arbitrary velocity stationarity is nonzero', np.linalg.norm(evaluation.r_v) > 0)
        checks.close('production physical independent-velocity energy', evaluation.physical_energy, ref['energy'])
        checks.close('production physical energy gradient',
                     np.r_[evaluation.energy_gradient_x, evaluation.energy_gradient_v], ref['gradient_energy'])
        checks.close('production compliant velocity mismatch', evaluation.coupling_residual, ref['mismatch'])
        checks.close('production full coupling Hessian', module.primal_matrix(scene, h, beta), ref['metric'])
        checks.close('production inverse application',
                     np.r_[module.inverse_primal(scene, h, beta, evaluation.r_x, evaluation.r_v)],
                     np.linalg.solve(ref['metric'], ref['residual'][:4]), 2e-10)
        for block in ('L', 'E', 'J'):
            ids = ([0, 1, 2, 3, 4] if block == 'L' else
                   [0, 1, 2, 3, 5, 6] if block == 'E' else list(range(7)))
            expected = np.zeros(7)
            expected[ids] = np.linalg.solve(ref['matrix'][np.ix_(ids, ids)], -ref['residual'][ids])
            delta = (module.material_increment(scene, h, beta, evaluation) if block == 'L' else
                     module.energy_increment(scene, h, epsilon, beta, evaluation) if block == 'E' else
                     module.joint_increment(scene, h, epsilon, beta, evaluation))
            actual = np.r_[delta.dx, delta.dv, delta.dlam, delta.drho, delta.deta]
            checks.close('production '+block+' versus independently assembled dense block',
                         actual, expected, 2e-8)
            checks.close('production '+block+' restores included linearized rows',
                         (ref['residual']+ref['matrix']@actual)[ids], np.zeros(len(ids)), 5e-8)

    replay_visits, distinct_velocities = 0, 0
    for scene in module.scenes():
        for beta in (.1, 1.0, 10.0):
            for order in ('EL', 'LE', 'JOINT'):
                for sweeps in (1, 2, 4):
                    h, epsilon = 1/120, 100.0
                    result = module.block_step(scene, scene.x0, scene.v0, h, epsilon,
                                               beta_v=beta, order=order, sweeps=sweeps, trace=True)
                    checks.true('finite replayed step accepted', result.accepted)
                    blocks = 'J' if order == 'JOINT' else order
                    checks.true('fixed budget visit count', result.visits == len(blocks)*sweeps)
                    checks.true('all post-visit states retained', len(result.history) == result.visits)
                    point = np.r_[scene.x0+h*scene.v0, scene.v0, 0.0, 0.0, 0.0]
                    target, floor = module.energy(scene, scene.x0, scene.v0), module.energy_floor(scene)
                    args = (scene.x0, scene.v0, h, scene.m, scene.k, scene.l, scene.g,
                            epsilon, target, floor, beta, scene.rigid)
                    for visit, snapshot in enumerate(result.history):
                        block = blocks[visit % len(blocks)]
                        checks.true('sequential block order', snapshot['block'] == block)
                        checks.true('sequential sweep index', snapshot['sweep'] == visit//len(blocks))
                        ref = independent(point, *args)
                        ids = ([0, 1, 2, 3, 4] if block == 'L' else
                               [0, 1, 2, 3, 5, 6] if block == 'E' else list(range(7)))
                        point[ids] += np.linalg.solve(ref['matrix'][np.ix_(ids, ids)], -ref['residual'][ids])
                        checks.close('refreshed dense replay of x v and accumulated scalars',
                                     np.r_[snapshot['x'], snapshot['v'], snapshot['lam'],
                                           snapshot['rho'], snapshot['eta']], point, 5e-8)
                        refreshed = independent(point, *args)
                        checks.close('post-visit physical residual trace',
                                     np.r_[snapshot['r_x'], snapshot['r_v'], snapshot['r_lambda'],
                                           snapshot['r_rho'], snapshot['r_eta']], refreshed['residual'], 5e-8)
                        replay_visits += 1
                    checks.close('returned independent velocity matches last primal iterate',
                                 np.r_[result.x, result.v, result.lam, result.rho, result.eta], point, 5e-8)
                    distinct_velocities += int(np.linalg.norm(result.v-(result.x-scene.x0)/h) > 1e-7)
    checks.true('returned velocity is not overwritten from displacement', distinct_velocities > 0)

    cancellation = module.Scene('position energy-gradient cancellation', 1.0, 100.0,
                                 1.0, np.zeros(2), np.array([1.1, 0.0]), np.zeros(2))
    cancel_x, cancel_v, cancel_h = np.array([1.1, 0.0]), np.array([-1.0, 0.0]), .1
    for beta in (.1, 1.0, 10.0):
        evaluation = module.evaluate_state(cancellation, cancellation.x0, cancellation.v0,
                                           cancel_h, cancel_x, cancel_v, -10., .1, .2,
                                           1e8, beta, target=1.0, kstar=0.0)
        checks.close('exact reconstructed energy-gradient cancellation',
                     cancel_h*evaluation.energy_gradient_x+evaluation.energy_gradient_v,
                     np.zeros(2), 1e-12)
        delta = module.energy_increment(cancellation, cancel_h, 1e8, beta, evaluation)
        checks.close('independent velocity leaves positive energy Schur at cancellation',
                     delta.denominator, beta+1e-8)
        checks.true('energy cancellation produces finite independent-velocity update',
                    np.all(np.isfinite(np.r_[delta.dx, delta.dv, delta.drho, delta.deta])))

    # At feasible current velocity the limiting equations equal the archived
    # position-only system. This reference is assembled here, independently.
    for point, original_args in list(cases(12)):
        xn, vn, h, mass, stiffness, length, gravity, epsilon, target, floor, beta, rigid = original_args
        point = point.copy()
        point[2:4] = (point[:2]-xn)/h
        base = independent(point, *original_args)
        normal = base['normal']
        effective_gradient = base['gradient']+mass*point[2:4]/h
        old_residual = np.r_[base['residual'][:2]+base['residual'][2:4]/h, base['residual'][4:]]
        old_matrix = np.zeros((5, 5))
        old_matrix[:2, :2] = mass*np.eye(2)/h**2
        old_matrix[:2, 2] = -normal
        old_matrix[2, :2] = normal
        old_matrix[2, 2] = base['alpha']
        old_matrix[:2, 4] = old_matrix[4, :2] = effective_gradient
        old_matrix[3, 3] = epsilon
        old_matrix[3, 4] = old_matrix[4, 3] = target-floor
        scene = module.Scene('beta limit', mass, stiffness, length, gravity, xn, vn, rigid)
        for beta in (1.0, 1e-3, 1e-6, 1e-9):
            evaluation = module.evaluate_state(scene, xn, vn, h, point[:2], point[2:4],
                                               *point[4:], epsilon, beta, target=target, kstar=floor)
            material = module.material_increment(scene, h, beta, evaluation)
            old_ids = [0, 1, 2]
            old_delta = np.linalg.solve(old_matrix[np.ix_(old_ids, old_ids)], -old_residual[old_ids])
            checks.close('material position response equals position-only equations',
                         np.r_[material.dx, material.dlam], old_delta, 3e-10)
            if beta == 1e-9:
                energy = module.energy_increment(scene, h, epsilon, beta, evaluation)
                old_ids = [0, 1, 3, 4]
                old_delta = np.linalg.solve(old_matrix[np.ix_(old_ids, old_ids)], -old_residual[old_ids])
                checks.close('energy beta tends to zero recovers position-only block',
                             np.r_[energy.dx, energy.drho, energy.deta], old_delta, 2e-7)
                checks.close('energy beta tends to zero recovers derived velocity increment',
                             energy.dv, energy.dx/h, 2e-7)
        hard_matrix, hard_residual = dualized_coupling(point, np.array([.4, -.7]),
                                                      *original_args[:10], 0.0, rigid)
        hard_delta = np.linalg.solve(hard_matrix, -hard_residual)
        old_delta = np.linalg.solve(old_matrix, -old_residual)
        checks.close('exact hard dual coupling equals position-only joint system',
                     hard_delta[[0, 1, 6, 7, 8]], old_delta, 2e-8)
        checks.close('exact hard dual velocity increment', hard_delta[2:4], hard_delta[:2]/h, 2e-8)

    scene = module.Scene('zero energy', 1.0, 100.0, 1.0, np.zeros(2), np.array([1.0, 0.0]), np.zeros(2))
    state = module.evaluate_state(scene, scene.x0, scene.v0, 1/60, scene.x0,
                                  scene.v0, 0.0, 0.0, 0.0, 100.0, 1.0)
    delta = module.energy_increment(scene, 1/60, 100.0, 1.0, state)
    checks.true('compatible zero energy row flagged', delta.compatible_rank_defect)
    checks.close('compatible zero energy row stays at rest', np.r_[delta.dx, delta.dv], np.zeros(4))
    try:
        module.energy_increment(scene, 1/60, 100.0, 1.0, replace(state, r_eta=1.0))
    except FloatingPointError as error:
        checks.true('incompatible zero energy row fails explicitly', str(error) == 'incompatible_zero_energy_row')
    else:
        raise AssertionError('incompatible zero energy row silently accepted')
    return dict(production_arbitrary_states=48, refreshed_replay_visits=replay_visits,
                returned_velocities_distinct_from_displacement=distinct_velocities,
                beta_limit_states=12)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--independent-only', action='store_true')
    parser.add_argument('--write-report', action='store_true')
    args = parser.parse_args()
    checks = Checks()
    details = independent_checks(checks)
    if not args.independent_only:
        details.update(production_checks(checks))
    result = dict(status='passed', checks=checks.count, details=details,
                  largest_relative_errors=checks.maxima)
    print(json.dumps(result, indent=2, sort_keys=True))
    if args.write_report:
        Path(__file__).with_name('VALIDATION.md').write_text(
            '''# Independent validation

The audit independently assembles the dense linearized equations. It does not
reuse the production closed-form block increments as the reference.

## Coverage

- Forty-eight arbitrary spring/rod states at beta_v 0.1, 1, and 10, with nonzero
  position and velocity stationarity residuals and independently carried material,
  energy, and relaxation variables.
- Finite differences of the original implicit objective plus velocity penalty,
  physical mechanical energy, and full residual Jacobian. The latter includes the
  material/energy curvature omitted by the numerical approximation.
- Full positive coupling Hessian and its analytic inverse; material, energy, and
  joint production increments versus independently assembled dense solves.
- Frozen EL and LE iterations versus the same joint linear solution, retaining
  position and velocity stationarity and pushing both primal variables immediately.
- Explicit velocity-coupling multiplier equations, with arbitrary inconsistent
  current multipliers and the full consistency residual. Elimination reproduces
  the primal formulation, including beta_v=1.
- Six hundred thirty refreshed nonlinear visits across all six scenes, both
  orders and JOINT, all three positive beta values, and one/two/four sweeps.
  Every trace state checks refresh timing and accumulation. Returned velocity is
  compared with the final independent velocity iterate, not overwritten from
  displacement.
- Twelve beta-limit states: material position updates equal the old position-only
  equations; the energy block approaches them as beta tends to zero. The explicit
  dual formulation at beta exactly zero reproduces the hard position-only joint
  system. The primal beta-zero penalty remains undefined and is not used.
- Exact cancellation of the reconstructed position energy gradient still leaves
  the beta_v-weighted independent-velocity contribution in the energy denominator.
  Compatible and incompatible zero energy rows are handled explicitly.

## Interpretation

These checks validate the stated adapted objective, curvature approximation,
block equations, and their implementation. They do not establish stable or
physically plausible long-time motion. Frozen GS convergence is distinct from
convergence of the refreshed nonlinear iteration. Independent velocity can also
permit displacement inconsistent with physical velocity; that requires separate
trajectory diagnostics even when actual mechanical energy is accurate.

Run with a NumPy environment:

```bash
OPENBLAS_NUM_THREADS=1 python3 validate.py --write-report
```

## Numerical receipt

```json
'''+json.dumps(result, indent=2, sort_keys=True)+'\n```\n')


if __name__ == '__main__':
    main()
