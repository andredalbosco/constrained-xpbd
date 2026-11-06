"""Audit primal constrained-objective curvature at selected saved steps.

This tests the original objective in (x, rho), with spring energy retained,
against the tangent of its one relaxed physical-energy equality. It does not
test the solver's approximate matrix or the independently introduced lambda.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np

import solver


ROOT = Path(__file__).resolve().parent
CONFIGURATIONS = [
    ('spring_orbit', 'EL', 60, 16, 1e8, None),
    ('spring_orbit', 'LE', 60, 16, 1e8, None),
    ('spring_orbit', 'JOINT', 60, 16, 1e8, None),
    ('spring_release', 'JOINT', 60, 16, 1e8, None),
    ('spring_release', 'LE', 60, 1, 1e4, 1),
]


def potential_hessian(scene, x):
    radius = float(np.linalg.norm(x))
    if radius <= 1e-14*scene.l:
        raise FloatingPointError('undefined_distance_normal')
    normal = x/radius
    radial = np.outer(normal, normal)
    return scene.k*(radial+(radius-scene.l)/radius*(np.eye(len(x))-radial))


def audit_step(scene, previous, current, h, epsilon):
    incoming_x = np.array([float(previous['x']), float(previous['y'])])
    incoming_v = np.array([float(previous['vx']), float(previous['vy'])])
    x = np.array([float(current['x']), float(current['y'])])
    lam, rho, eta = (float(current[key]) for key in ('lambda', 'rho', 'eta'))
    evaluation = solver.evaluate_state(scene, incoming_x, incoming_v, h, x,
                                       lam, rho, eta, epsilon)
    reconstructed_v = (x-incoming_x)/h
    saved_v = np.array([float(current['vx']), float(current['vy'])])
    assert np.allclose(reconstructed_v, saved_v, rtol=0, atol=1e-11)

    # Scale position coordinates by ell, while rho is already dimensionless.
    scale = np.diag([scene.l, scene.l, 1.])
    primal_hessian = np.zeros((3, 3))
    primal_hessian[:2, :2] = (1+eta)*(scene.m*np.eye(2)/(h*h)+potential_hessian(scene, x))
    primal_hessian[2, 2] = epsilon
    constraint_gradient = np.r_[evaluation.energy_gradient, evaluation.target-evaluation.kstar]
    scaled_gradient = scale.T@constraint_gradient
    scaled_hessian = scale.T@primal_hessian@scale
    _, singular_values, vh = np.linalg.svd(scaled_gradient[None, :], full_matrices=True)
    largest_singular_value = float(singular_values[0])
    rank = int(largest_singular_value > 0.)
    tangent = vh[rank:].T
    reduced_hessian = tangent.T@scaled_hessian@tangent
    reduced_hessian = .5*(reduced_hessian+reduced_hessian.T)
    eigenvalues = np.linalg.eigvalsh(reduced_hessian)
    relative_threshold = 1e-8
    absolute_threshold = relative_threshold*float(np.max(np.abs(eigenvalues)))
    positive = int(np.count_nonzero(eigenvalues > absolute_threshold))
    negative = int(np.count_nonzero(eigenvalues < -absolute_threshold))
    nearzero = len(eigenvalues)-positive-negative
    classification = ('negative_curvature' if negative else
                      'nearzero_curvature' if nearzero else 'positive_curvature')
    scaled_residual = solver.scaled_residual(scene, h, epsilon, evaluation)
    constitutive_consistent = abs(evaluation.r_lambda) <= 1e-9*scene.l
    primal_r_x = (scene.m*(x-incoming_x-h*incoming_v)/(h*h)
                  +evaluation.potential_gradient+eta*evaluation.energy_gradient)
    stationarity_eligible = bool(int(current['converged']) and scaled_residual <= 1e-9
                                and constitutive_consistent and rank == 1)
    return dict(
        step=int(current['step']), time=float(current['time']),
        incoming_x=incoming_x.tolist(), incoming_v=incoming_v.tolist(),
        x=x.tolist(), rho=rho, eta=eta, lam=lam,
        energy_target=evaluation.target, energy_reference=evaluation.kstar,
        row_converged=bool(int(current['converged'])),
        recomputed_scaled_residual=scaled_residual,
        r_x=evaluation.r_x.tolist(), r_lambda=evaluation.r_lambda,
        original_primal_r_x=primal_r_x.tolist(),
        r_rho=evaluation.r_rho, r_eta=evaluation.r_eta,
        constitutive_consistent=bool(constitutive_consistent),
        first_order_eligible=stationarity_eligible,
        constraint_rank=rank,
        constraint_gradient=constraint_gradient.tolist(),
        primal_hessian=primal_hessian.tolist(),
        position_coordinate_scale=scene.l,
        scaled_tangent_basis=tangent.tolist(),
        tangent_gradient_defect=float(np.linalg.norm(scaled_gradient@tangent)),
        reduced_hessian=reduced_hessian.tolist(),
        eigenvalues=eigenvalues.tolist(),
        relative_sign_threshold=relative_threshold,
        absolute_sign_threshold=absolute_threshold,
        positive_eigenvalues=positive, negative_eigenvalues=negative,
        nearzero_eigenvalues=nearzero, classification=classification,
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--results', type=Path, default=ROOT/'results')
    parser.add_argument('--output', type=Path, default=ROOT/'diagnostics'/'curvature_audit.json')
    args = parser.parse_args()
    scene_by_name = {scene.name: scene for scene in solver.scenes()}
    reports = []
    for name, order, hz, sweeps, epsilon_ratio, maximum_steps in CONFIGURATIONS:
        scene = scene_by_name[name]
        assert not scene.rigid
        path = args.results/'trajectories'/f'{name}_{order}_hz{hz}_s{sweeps}_eps{epsilon_ratio:g}.csv'
        with path.open(newline='') as stream:
            rows = list(csv.DictReader(stream))
        epsilon = epsilon_ratio*(solver.energy(scene, scene.x0, scene.v0)-solver.energy_floor(scene))
        steps = []
        for previous, current in zip(rows, rows[1:]):
            if maximum_steps is not None and len(steps) >= maximum_steps:
                break
            if not int(previous['accepted']) or not int(current['accepted']):
                break
            steps.append(audit_step(scene, previous, current, 1/hz, epsilon))
        counts = {key: sum(row['classification'] == key for row in steps)
                  for key in ('positive_curvature', 'negative_curvature', 'nearzero_curvature')}
        eligible_counts = {key: sum(row['first_order_eligible'] and row['classification'] == key for row in steps)
                           for key in counts}
        report = dict(
            scene=name, order=order, hz=hz, sweeps=sweeps,
            epsilon_ratio=epsilon_ratio, epsilon=epsilon,
            source=path.name, source_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
            audited_steps=len(steps), classification_counts=counts,
            first_order_eligible_steps=sum(row['first_order_eligible'] for row in steps),
            first_order_eligible_classification_counts=eligible_counts,
            max_material_residual=max((abs(row['r_lambda']) for row in steps), default=0.),
            max_recomputed_scaled_residual=max((row['recomputed_scaled_residual'] for row in steps), default=0.),
            steps=steps,
        )
        reports.append(report)
        print(name, order, 'audited', len(steps), counts,
              'first-order-eligible', report['first_order_eligible_steps'],
              'max |r_lambda|', report['max_material_residual'])
    data = dict(
        method='Exact primal Lagrangian Hessian in (x,rho), restricted to one physical-energy equality tangent',
        hessian_xx='(1+eta)*(m/h^2*I + full Hessian U)', hessian_rho_rho='epsilon',
        tangent='Nullspace from SVD of [ell*energy_position_gradient, E_star-K_star]',
        sign_threshold='1e-8 times largest absolute reduced eigenvalue, per step',
        qualification='Negative tangent curvature excludes a local minimum at an exact regular stationary feasible point; approximate residual roots require numerical qualification. Positive curvature is a local second-order check, not global optimality or trajectory accuracy. Nearzero curvature is inconclusive.',
        reports=reports,
    )
    args.output.parent.mkdir(exist_ok=True)
    args.output.write_text(json.dumps(data, indent=2, allow_nan=False)+'\n')
    print('Wrote', args.output)


if __name__ == '__main__':
    main()
