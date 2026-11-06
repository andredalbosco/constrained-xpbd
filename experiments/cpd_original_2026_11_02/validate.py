"""Independent equation-21 CPD audit for planar spring networks.

The reference assembles PD curvature, objective gradients and full conservation
rows directly, then compares paper Schur elimination with a dense KKT solve.
"""
from __future__ import annotations

import argparse
from dataclasses import dataclass
import importlib
import json
from pathlib import Path

import numpy as np


@dataclass
class Case:
    x0: np.ndarray
    v0: np.ndarray
    mass: np.ndarray
    free: np.ndarray
    edges: list
    gravity: np.ndarray
    linear: bool = True
    angular: bool = True
    anchor_origin: np.ndarray | None = None
    forces: np.ndarray | None = None


def elastic(case, x, projections=None):
    energy = 0.0
    gradient = np.zeros_like(x)
    for index, (a, b, stiffness, rest) in enumerate(case.edges):
        delta = x[a]-x[b]
        radius = np.linalg.norm(delta)
        if projections is None:
            force = stiffness*(radius-rest)*delta/radius
            energy += .5*stiffness*(radius-rest)**2
        else:
            force = stiffness*(delta-projections[index])
            energy += .5*stiffness*np.dot(delta-projections[index], delta-projections[index])
        gradient[a] += force
        gradient[b] -= force
    return energy, gradient


def preparation(case, h, epsilon):
    indices = np.flatnonzero(case.free)
    masses = case.mass[indices]
    x0, v0 = case.x0[indices], case.v0[indices]
    total_mass = masses.sum()
    origin = (np.sum(masses[:, None]*x0, axis=0)/total_mass
              if case.anchor_origin is None else case.anchor_origin.copy())
    external = masses[:, None]*case.gravity if case.forces is None else case.forces[indices]
    prediction = case.x0.copy()
    prediction[indices] += h*v0+h*h*external/masses[:, None]
    target_p = np.sum(masses[:, None]*v0+h*external, axis=0)
    target_l = np.sum(np.cross(np.c_[x0-origin, np.zeros(len(x0))],
                               np.c_[masses[:, None]*v0+h*external, np.zeros(len(x0))])[:, 2])
    target_h = elastic(case, case.x0)[0]+.5*np.sum(masses[:, None]*v0*v0)+h*np.sum(external*v0)
    inertia = np.sum(masses*np.sum((x0-origin)**2, axis=1))
    kinetic_reference = .5*np.dot(target_p, target_p)/total_mass if case.linear else 0.0
    if case.angular:
        if inertia == 0 and target_l != 0:
            raise ValueError('angular target outside zero inertia range')
        kinetic_reference += 0.0 if inertia == 0 else .5*target_l**2/inertia
    position_matrix = np.diag(np.repeat(masses, 2)/h**2)
    mapping = {int(index): 2*i for i, index in enumerate(indices)}
    for a, b, stiffness, rest in case.edges:
        for endpoint in (a, b):
            if endpoint in mapping:
                start = mapping[endpoint]
                position_matrix[start:start+2, start:start+2] += stiffness*np.eye(2)
        if a in mapping and b in mapping:
            ia, ib = mapping[a], mapping[b]
            position_matrix[ia:ia+2, ib:ib+2] -= stiffness*np.eye(2)
            position_matrix[ib:ib+2, ia:ia+2] -= stiffness*np.eye(2)
    matrix = np.zeros((2*len(indices)+1, 2*len(indices)+1))
    matrix[:-1, :-1], matrix[-1, -1] = position_matrix, epsilon
    return dict(indices=indices, masses=masses, origin=origin, prediction=prediction,
                target_p=target_p, target_l=float(target_l), target_h=float(target_h),
                kinetic_reference=float(kinetic_reference), inertia=float(inertia),
                metric=matrix, h=h, epsilon=epsilon)


def evaluate(case, prepared, point):
    ids, mass = prepared['indices'], prepared['masses']
    h, epsilon = prepared['h'], prepared['epsilon']
    x = case.x0.copy()
    x[ids] = point[:-1].reshape(-1, 2)
    rho = point[-1]
    displacement = x[ids]-case.x0[ids]
    velocity = displacement/h
    potential, gradient_u = elastic(case, x)
    gradient = np.r_[(mass[:, None]*(x[ids]-prepared['prediction'][ids])/h**2
                     +gradient_u[ids]).ravel(), epsilon*rho]
    objective = .5*np.sum(mass[:, None]*(x[ids]-prepared['prediction'][ids])**2)/h**2+potential+.5*epsilon*rho*rho
    momentum = np.sum(mass[:, None]*velocity, axis=0)
    relative = x[ids]-prepared['origin']
    angular = float(np.sum(mass*(relative[:, 0]*velocity[:, 1]-relative[:, 1]*velocity[:, 0])))
    internal = .5*np.sum(mass[:, None]*velocity*velocity)+potential
    row_count = 2*int(case.linear)+int(case.angular)+1
    jacobian = np.zeros((row_count, len(point)))
    constraints = np.zeros(row_count)
    row = 0
    if case.linear:
        constraints[row:row+2] = momentum-prepared['target_p']
        for particle, particle_mass in enumerate(mass):
            jacobian[row:row+2, 2*particle:2*particle+2] = particle_mass*np.eye(2)/h
        row += 2
    if case.angular:
        constraints[row] = angular-prepared['target_l']
        old_relative = case.x0[ids]-prepared['origin']
        jacobian[row, :-1] = (mass[:, None]/h*np.column_stack([-old_relative[:, 1], old_relative[:, 0]])).ravel()
        row += 1
    constraints[row] = internal-(1-rho)*prepared['target_h']-rho*prepared['kinetic_reference']
    jacobian[row, :-1] = (mass[:, None]*displacement/h**2+gradient_u[ids]).ravel()
    jacobian[row, -1] = prepared['target_h']-prepared['kinetic_reference']
    return dict(objective=float(objective), gradient=gradient, constraints=constraints,
                jacobian=jacobian, internal_energy=float(internal), momentum=momentum,
                angular=angular, positions=x, velocity=velocity)


def finite_difference(function, point):
    value = np.asarray(function(point))
    derivative = np.empty((value.size, len(point)))
    for column in range(len(point)):
        step = 2e-6*max(1.0, abs(point[column]))
        plus, minus = point.copy(), point.copy()
        plus[column] += step
        minus[column] -= step
        derivative[:, column] = (function(plus)-function(minus))/(2*step)
    return derivative


class Checks:
    def __init__(self):
        self.count, self.maxima = 0, {}

    def close(self, name, actual, expected, tolerance=3e-9):
        self.count += 1
        error = float(np.linalg.norm(np.asarray(actual)-np.asarray(expected))/max(1.0, np.linalg.norm(expected)))
        self.maxima[name] = max(error, self.maxima.get(name, 0.0))
        assert np.isfinite(error) and error <= tolerance, (name, error, tolerance)

    def true(self, name, condition):
        self.count += 1
        assert condition, name


def cases():
    rng = np.random.default_rng(20261102)
    for index in range(48):
        anchored = index % 4 == 0
        if anchored:
            x0 = np.array([[0.0, 0.0], [1.1, .4]])+rng.normal(size=(2, 2))*.05
            x0[0] = 0.0
            free = np.array([False, True])
            edges = [(0, 1, 10**rng.uniform(1, 3), rng.uniform(.8, 1.2))]
            linear, angular = False, bool(index % 8)
        else:
            n = 2 if index % 3 else 3
            x0 = np.array([[-.7, -.1], [.6, .3], [.2, 1.0]])[:n]+rng.normal(size=(n, 2))*.08
            free = np.ones(n, bool)
            edges = [(0, 1, 10**rng.uniform(1, 3), rng.uniform(.8, 1.2))]
            if n == 3:
                edges += [(1, 2, 50.0, .9), (2, 0, 70.0, 1.1)]
            linear = angular = True
        v0 = rng.normal(size=x0.shape)*.7
        v0[~free] = 0.0
        mass = 10**rng.uniform(-.3, .3, len(x0))
        gravity = np.array([0.0, -2.7]) if index % 3 else np.zeros(2)
        forces = (mass[:, None]*gravity+rng.normal(size=x0.shape)*.8
                  if index % 5 == 1 else None)
        case = Case(x0, v0, mass, free, edges, gravity, linear, angular,
                    np.zeros(2) if anchored else None, forces)
        h, epsilon = rng.uniform(.02, .1), 10**rng.uniform(1, 6)
        prepared = preparation(case, h, epsilon)
        x = x0[free]+h*v0[free]+rng.normal(size=(free.sum(), 2))*.04
        point = np.r_[x.ravel(), rng.normal()*.1]
        yield case, prepared, point


def independent_checks(checks):
    for case, prepared, point in cases():
        value = evaluate(case, prepared, point)
        jacobian, gradient = value['jacobian'], value['gradient']
        metric = prepared['metric']
        checks.close('full conservation Jacobian finite differences',
                     finite_difference(lambda q: evaluate(case, prepared, q)['constraints'], point), jacobian, 4e-8)
        checks.close('implicit plus relaxation objective gradient finite differences',
                     finite_difference(lambda q: np.array([evaluate(case, prepared, q)['objective']]), point)[0], gradient, 4e-8)
        inverse_gradient = np.linalg.solve(metric, gradient)
        inverse_rows = np.linalg.solve(metric, jacobian.T)
        schur = jacobian@inverse_rows
        multipliers = np.linalg.solve(schur, value['constraints']-jacobian@inverse_gradient)
        increment = -inverse_gradient-inverse_rows@multipliers
        kkt = np.block([[metric, jacobian.T], [jacobian, np.zeros((len(jacobian), len(jacobian)))]])
        dense = np.linalg.solve(kkt, -np.r_[gradient, value['constraints']])
        checks.close('paper equation21 versus full dense KKT position and relaxation', increment, dense[:len(point)], 2e-8)
        checks.close('paper equation21 fresh multiplier versus dense KKT', multipliers, dense[len(point):], 2e-8)
        checks.close('linearized conservation satisfied jointly', value['constraints']+jacobian@increment,
                     np.zeros(len(jacobian)), 3e-8)
        checks.close('quadratic subproblem stationarity satisfied', gradient+metric@increment+jacobian.T@multipliers,
                     np.zeros(len(point)), 3e-8)
        checks.true('independent planar conservation rows', np.linalg.matrix_rank(jacobian) == len(jacobian))

        # Freeze local projections to recover PD curvature, not GN curvature.
        x = value['positions']
        projection = [(rest*(x[a]-x[b])/np.linalg.norm(x[a]-x[b]))
                      for a, b, stiffness, rest in case.edges]
        ids, masses, h = prepared['indices'], prepared['masses'], prepared['h']
        def frozen_gradient(q):
            trial = case.x0.copy()
            trial[ids] = q[:-1].reshape(-1, 2)
            gu = elastic(case, trial, projection)[1]
            return np.r_[(masses[:, None]*(trial[ids]-prepared['prediction'][ids])/h**2+gu[ids]).ravel(), prepared['epsilon']*q[-1]]
        checks.close('frozen local projection gives full PD curvature',
                     finite_difference(frozen_gradient, point), metric, 4e-8)
        shift = np.array([17.2, -8.4])
        translated = Case(case.x0+shift, case.v0, case.mass, case.free, case.edges,
                          case.gravity, case.linear, case.angular,
                          None if case.anchor_origin is None else case.anchor_origin+shift, case.forces)
        translated_prepared = preparation(translated, h, prepared['epsilon'])
        translated_point = point.copy()
        translated_point[:-1] = (point[:-1].reshape(-1, 2)+shift).ravel()
        translated_value = evaluate(translated, translated_prepared, translated_point)
        checks.close('consistent origin translation preserves conservation rows',
                     translated_value['constraints'], value['constraints'], 2e-10)
        checks.close('consistent origin translation preserves conservation Jacobian',
                     translated_value['jacobian'], jacobian, 2e-10)
        checks.close('consistent origin translation preserves objective gradient',
                     translated_value['gradient'], gradient, 2e-10)
        # The kinetic reference must equal the constrained least-kinetic norm,
        # not merely a sum of formulas with inconsistent origins.
        count_momentum = 2*int(case.linear)+int(case.angular)
        if count_momentum:
            velocity_rows = jacobian[:count_momentum, :-1]*h
            target = np.r_[prepared['target_p'] if case.linear else [],
                           [prepared['target_l']] if case.angular else []]
            inverse_mass = np.diag(1/np.repeat(masses, 2))
            minimal = inverse_mass@velocity_rows.T@np.linalg.solve(
                velocity_rows@inverse_mass@velocity_rows.T, target)
            minimum_kinetic = .5*minimal@np.diag(np.repeat(masses, 2))@minimal
            checks.close('inertia kinetic reference versus constrained least kinetic energy',
                         prepared['kinetic_reference'], minimum_kinetic, 2e-10)

    # Redundant but compatible rows: one free particle in pure translation.
    translation = Case(np.array([[.3, -.2]]), np.array([[.6, -.4]]), np.ones(1),
                       np.ones(1, bool), [], np.zeros(2))
    prep = preparation(translation, .03, 1e4)
    value = evaluate(translation, prep, np.r_[prep['prediction'].ravel(), 0.0])
    checks.true('single particle full conservation rows are dependent',
                np.linalg.matrix_rank(value['jacobian']) < len(value['constraints']))
    checks.close('redundant momentum and energy rows can be exactly compatible',
                 value['constraints'], np.zeros(4), 1e-12)

    # A rotating, unstretched spring at x_n has zero energy gradient. Nonzero
    # incoming kinetic energy makes the linearized energy equation inconsistent.
    rotating = Case(np.array([[-.5, 0.0], [.5, 0.0]]),
                    np.array([[0.0, -.4], [0.0, .4]]), np.ones(2),
                    np.ones(2, bool), [(0, 1, 100.0, 1.0)], np.zeros(2))
    prep = preparation(rotating, .03, 1e4)
    value = evaluate(rotating, prep, np.r_[rotating.x0.ravel(), 0.0])
    checks.close('zero energy gradient at unstretched incoming position', value['jacobian'][-1],
                 np.zeros(5), 1e-12)
    checks.true('zero energy row has incompatible nonzero residual', value['constraints'][-1] != 0.0)

    # Full free-body P/L do not remove the relaxed stationary release branch.
    releasing = Case(np.array([[-.6, 0.0], [.6, 0.0]]), np.zeros((2, 2)), np.ones(2),
                     np.ones(2, bool), [(0, 1, 100.0, 1.0)], np.zeros(2))
    h, epsilon = .03, 1e4
    prep = preparation(releasing, h, epsilon)
    initial = prep['target_h']
    rho = initial/epsilon
    target = initial-initial**2/epsilon
    old_radius, stiffness, rest = 1.2, 100.0, 1.0
    roots = np.roots([.25/h**2+.5*stiffness,
                      -.5*old_radius/h**2-stiffness*rest,
                      .25*old_radius**2/h**2+.5*stiffness*rest**2-target])
    new_radius = min(roots, key=lambda r: abs(r-old_radius))
    point = np.r_[-new_radius/2, 0.0, new_radius/2, 0.0, rho]
    value = evaluate(releasing, prep, point)
    multipliers = np.array([0.0, 0.0, 0.0, -1.0])
    checks.close('full P L E relaxed release branch feasible', value['constraints'], np.zeros(4), 2e-10)
    checks.close('full P L E relaxed release branch stationary',
                 value['gradient']+value['jacobian'].T@multipliers, np.zeros(5), 2e-10)
    checks.true('finite epsilon relaxation need not vanish for a feasible strict target', rho > 0)
    return dict(random_network_states=48, explicit_rank_cases=2, full_free_release_counterexample=True)


def production_checks(checks):
    module = importlib.import_module('cpd_solver')
    for case, prepared, point in cases():
        indices = prepared['indices']
        mapping = {int(index): i for i, index in enumerate(indices)}
        edges = []
        for a, b, stiffness, rest in case.edges:
            if a in mapping and b in mapping:
                edges.append(module.Edge(mapping[a], mapping[b], rest, stiffness))
            else:
                free = a if a in mapping else b
                edges.append(module.Edge(mapping[free], -1, rest, stiffness))
        forces = case.mass[:, None]*case.gravity if case.forces is None else case.forces
        scene = module.Scene('arbitrary', case.mass[indices], case.x0[indices], case.v0[indices],
                             forces[indices], edges, 'free' if case.linear else 'anchored')
        # The optional angular-row ablation intentionally retains the original
        # kinetic reference, a declared root-solver experimental convention.
        ref_case = Case(case.x0, case.v0, case.mass, case.free, case.edges,
                        case.gravity, case.linear, True, case.anchor_origin, case.forces)
        actual_prepared = preparation(ref_case, prepared['h'], prepared['epsilon'])
        ref = evaluate(case, actual_prepared, point)
        h, epsilon = prepared['h'], prepared['epsilon']
        targets = module.make_targets(scene, scene.x0, scene.v0, h, angular=case.angular)
        checks.close('production incoming momentum target', targets.linear, actual_prepared['target_p'])
        checks.close('production incoming angular and torque target', targets.angular, actual_prepared['target_l'])
        checks.close('production incoming internal energy plus power target', targets.energy, actual_prepared['target_h'])
        checks.close('production inertia kinetic reference', targets.kstar, actual_prepared['kinetic_reference'])
        checks.close('production incoming COM or anchor origin', targets.origin, actual_prepared['origin'])
        checks.close('production PD metric not rank-one GN metric',
                     module.pd_matrix(scene, h), actual_prepared['metric'][:-1, :-1])
        value = module.evaluate(scene, scene.x0, scene.v0, h,
                                point[:-1].reshape(scene.x0.shape), point[-1], epsilon, targets)
        checks.close('production true elastic and inertial objective gradient', value.objective_gradient, ref['gradient'])
        checks.close('production true P L energy constraint values', value.constraints, ref['constraints'])
        checks.close('production analytic P L energy gradients', value.constraint_gradients, ref['jacobian'])
        checks.close('production kinetic plus elastic constrained energy', value.internal_energy, ref['internal_energy'])
        physical = ref['internal_energy']-np.sum(scene.forces*point[:-1].reshape(scene.x0.shape))
        checks.close('production mechanical energy diagnostic includes force potential', value.physical_energy, physical)
        dense = np.linalg.solve(np.block([
            [actual_prepared['metric'], ref['jacobian'].T],
            [ref['jacobian'], np.zeros((len(ref['jacobian']), len(ref['jacobian'])))]
        ]), -np.r_[ref['gradient'], ref['constraints']])
        update = module.increment(scene, h, epsilon, value)
        checks.close('production equation21 increment versus independent dense KKT',
                     np.r_[update.dx.ravel(), update.drho], dense[:len(point)], 5e-8)
        checks.close('production equation21 fresh multipliers versus independent dense KKT',
                     update.multipliers, dense[len(point):], 5e-8)
        checks.true('production subproblem reports small linear residual', update.linear_residual < 1e-8)

    replayed = 0
    fresh_distinct_from_accumulated = 0
    trace_errors = {}
    all_scenes = module.scenes()
    for scene in all_scenes:
        h, epsilon = 1/60, 1e4
        # Include a fixed origin vertex in the independent network only when
        # needed by anchored edges; the production unknowns remain free vertices.
        anchor = any(edge.b < 0 for edge in scene.edges)
        n = len(scene.masses)
        old = np.vstack([scene.x0, np.zeros(2)]) if anchor else scene.x0.copy()
        velocity = np.vstack([scene.v0, np.zeros(2)]) if anchor else scene.v0.copy()
        mass = np.r_[scene.masses, 1.0] if anchor else scene.masses.copy()
        free = np.r_[np.ones(n, bool), False] if anchor else np.ones(n, bool)
        forces = np.vstack([scene.forces, np.zeros(2)]) if anchor else scene.forces.copy()
        edges = [(e.a, e.b if e.b >= 0 else n, 0.0 if e.hard else e.stiffness, e.rest)
                 for e in scene.edges]
        case = Case(old, velocity, mass, free, edges, np.zeros(2),
                    scene.family == 'free', True,
                    None if scene.family == 'free' else np.zeros(2), forces)
        prepared = preparation(case, h, epsilon)
        targets = module.make_targets(scene, scene.x0, scene.v0, h)

        def reference(point):
            data = evaluate(case, prepared, point)
            if any(edge.hard for edge in scene.edges):
                rows, values = [], []
                for edge in scene.edges:
                    if edge.hard:
                        x = point[:-1].reshape(scene.x0.shape)
                        delta = x[edge.a]-(x[edge.b] if edge.b >= 0 else np.zeros(2))
                        radius = np.linalg.norm(delta)
                        row = np.zeros(len(point))
                        row[2*edge.a:2*edge.a+2] = delta/radius
                        if edge.b >= 0:
                            row[2*edge.b:2*edge.b+2] -= delta/radius
                        rows.append(row)
                        values.append(radius-edge.rest)
                data['jacobian'] = np.vstack([data['jacobian'][:-1], rows, data['jacobian'][-1:]])
                data['constraints'] = np.r_[data['constraints'][:-1], values, data['constraints'][-1]]
            return data

        point = np.r_[prepared['prediction'][free].ravel(), 0.0]
        check_point = point.copy()
        check_point[:-1] += .003*np.arange(1, len(point))
        check = reference(check_point)
        actual = module.evaluate(scene, scene.x0, scene.v0, h,
                                 check_point[:-1].reshape(scene.x0.shape), check_point[-1], epsilon, targets)
        checks.close('all scene analytic gradients including exact rod extension', actual.constraint_gradients,
                     check['jacobian'])
        checks.close('all scene constraints including exact rod extension', actual.constraints, check['constraints'])
        checks.close('all scene conservation finite differences including hard geometry',
                     finite_difference(lambda q: reference(q)['constraints'], check_point), check['jacobian'], 4e-8)
        result = module.step(scene, scene.x0, scene.v0, h, epsilon, policy='fixed', iterations=4, trace=True)
        checks.true('four-step dense replay accepted', result.accepted)
        checks.true('fixed global iteration budget retained', result.iterations == 4 and len(result.history) == 4)
        accumulated = np.zeros(len(result.multipliers))
        largest_state_error = largest_multiplier_error = 0.0
        for snapshot in result.history:
            data = reference(point)
            jacobian = data['jacobian']
            matrix = np.block([[prepared['metric'], jacobian.T],
                               [jacobian, np.zeros((len(jacobian), len(jacobian)))]])
            change = np.linalg.solve(matrix, -np.r_[data['gradient'], data['constraints']])
            point += change[:len(point)]
            new_multiplier = change[len(point):]
            state_actual = np.r_[np.array(snapshot['x']).ravel(), snapshot['rho']]
            largest_state_error = max(largest_state_error, float(np.linalg.norm(state_actual-point)/max(1.0, np.linalg.norm(point))))
            largest_multiplier_error = max(largest_multiplier_error, float(np.linalg.norm(np.asarray(snapshot['multipliers'])-new_multiplier)/max(1.0, np.linalg.norm(new_multiplier))))
            checks.close('global step trace position and relaxation versus dense replay',
                         np.r_[np.array(snapshot['x']).ravel(), snapshot['rho']], point, 2e-6)
            checks.close('global step trace NEW multiplier versus dense replay',
                         snapshot['multipliers'], new_multiplier, 2e-6)
            fresh_distinct_from_accumulated += int(np.linalg.norm(accumulated+new_multiplier-new_multiplier) > 1e-6)
            accumulated += new_multiplier
            refreshed = reference(point)
            checks.close('refreshed global trace true constraints', snapshot['constraints'], refreshed['constraints'], 3e-6)
            replayed += 1
        trace_errors[scene.name] = dict(state=largest_state_error, multiplier=largest_multiplier_error,
                                       max_scaled_schur_condition=result.condition_max,
                                       max_linear_residual=result.max_linear_residual)
        checks.close('returned velocity reconstructed from physical displacement', result.v, (result.x-scene.x0)/h)
    checks.true('fresh multipliers distinguished from accumulated values', fresh_distinct_from_accumulated > 0)

    # Free COM-zero pairs use reduced mass one, so their relative coordinates
    # must reproduce anchored unit-mass equations with the same angular row.
    pairs = [('spring_release', 'free_release'), ('spring_orbit', 'free_orbit'),
             ('spring_gravity', 'free_opposite_load'), ('spring_stiff', 'free_stiff')]
    lookup = {scene.name: scene for scene in all_scenes}
    for anchored_name, free_name in pairs:
        anchored, free_pair = lookup[anchored_name], lookup[free_name]
        for h in (1/30, 1/120):
            a = module.step(anchored, anchored.x0, anchored.v0, h, 1e4, policy='fixed', iterations=4)
            b = module.step(free_pair, free_pair.x0, free_pair.v0, h, 1e4, policy='fixed', iterations=4)
            checks.true('anchored and free relative coordinate comparison accepted', a.accepted and b.accepted)
            checks.close('free pair relative position versus anchored position', b.x[0]-b.x[1], a.x[0], 2e-7)
            checks.close('free pair relative velocity versus anchored velocity', b.v[0]-b.v[1], a.v[0], 2e-7)
            checks.close('free pair energy relaxation versus anchored relaxation', b.rho, a.rho, 2e-7)

    # A scalar axial finite spring gives a closed-form backward Euler PD step.
    scalar = lookup['spring_release']
    h = .03
    pd = module.step(scalar, scalar.x0, scalar.v0, h, 1e4, policy='fixed', iterations=1, method='PD')
    expected_x = (scalar.x0[0, 0]/h**2+100.0)/(1/h**2+100.0)
    checks.true('analytic PD finite spring accepted', pd.accepted)
    checks.close('one axial PD global solve equals analytic backward Euler', pd.x, [[expected_x, 0.0]], 1e-12)
    checks.true('PD convergence uses its own objective stationarity', pd.stationary)
    precheck = module.step(scalar, scalar.x0, scalar.v0, h, 1e8,
                           policy='paper', initial_precheck=True)
    checks.true('initial conservation precheck exposes zero-iteration release', precheck.accepted and precheck.iterations == 0)
    checks.true('zero-iteration conservation does not imply objective stationarity', not precheck.stationary)
    return dict(production_random_states=48, global_trace_visits=replayed,
                free_pair_anchored_comparisons=8, hard_rod_extension_scenes=2,
                trace_error_by_scene=trace_errors)


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
            '''# Independent CPD validation

The independent reference assembles planar spring-network PD curvature, the true
objective gradient, incoming momentum/power targets, and conservation rows. It
compares paper equation 21 with a dense KKT solve, rather than duplicating the
production Schur increment formula.

## Coverage

- Forty-eight arbitrary two/three-particle and anchored network states, nonzero
  velocities, unequal masses, gravity and nonuniform external forces, off-solution
  positions and relaxation values, and explicit angular-row ablations.
- Finite differences of the objective gradient and full P/L/energy Jacobian;
  frozen projection curvature equals the full PD matrix rather than GN normal
  outer products. Translations use a consistently translated incoming origin.
- Kinetic reference checked against an independent minimum-kinetic-energy solve
  with the enabled incoming momentum rows. The production angular-row ablation
  deliberately keeps the original kinetic reference unchanged.
- Production target values, gradients, internal versus mechanical energy, PD
  metric, global position/relaxation increment, and newly solved conservation
  multipliers compared with independent assembly and dense KKT solutions.
- Fifty-two refreshed global iterations across all thirteen scenes, including
  analytic hard-geometry gradients for two explicitly separate rod extensions.
  The trace checks newly solved multipliers, not their cumulative sum.
- Eight free-pair versus anchored comparisons: reduced mass one, equal relative
  spring state and force, and equal angular target. Relative positions, velocities,
  and energy relaxation agree to roundoff.
- One axial PD global solve agrees with analytic backward Euler. PD convergence
  uses objective stationarity, independently of inactive CPD constraints.
- Compatible dependent conservation rows and an incompatible zero energy row are
  exposed analytically. No reference step silently substitutes a pseudoinverse.
- A full-free-body release counterexample satisfies P/L/relaxed-energy feasibility
  and KKT stationarity with a nonzero relaxation value despite a feasible strict
  target. An initial conservation-only precheck can return zero iterations while
  objective stationarity remains unsatisfied.

## Numerical limitations

Random production increments and multipliers match dense KKT to approximately
5.1e-14. Finite-spring trace agreement is below 1e-10. The hard rod extensions
are substantially less well conditioned: the kick example reaches scaled Schur
condition 1.28e10 and frozen linear residual 1.67e-7. Its four-step dense replay
differs by 8.57e-7 in state and 1.14e-6 in multiplier, within the explicitly
declared 2e-6 trace comparison tolerance. The receipt exposes these values rather
than presenting that extension as a well-conditioned original-CPD test.

These are equation/implementation checks, not proof of acceptable motion. The
paper conservation stop, full KKT stationarity, physical energy under external
loads, and visual motion remain separate diagnostics. The source cross-check in
CPP_PARITY.md uses the prior unchanged local C++ baseline, not authors' code.

```bash
OPENBLAS_NUM_THREADS=1 python3 validate.py --write-report
```

## Numerical receipt

```json
'''+json.dumps(result, indent=2, sort_keys=True)+'\n```\n')


if __name__ == '__main__':
    main()
