"""Strictly infeasible CPD endpoint targets: standalone one-step fixtures.

These are outside the primary 13-scene grid. Both particles are free; there are
no pins, hard rows, prescribed loads, or changes to the production solver.
"""
from __future__ import annotations

import dataclasses
import json
import math
from pathlib import Path

import numpy as np
from scipy.optimize import brentq, minimize_scalar

import cpd_solver


ROOT = Path(__file__).resolve().parent


def native(value):
    if dataclasses.is_dataclass(value):
        return native(dataclasses.asdict(value))
    if isinstance(value, np.ndarray):
        return native(value.tolist())
    if isinstance(value, np.generic):
        return native(value.item())
    if isinstance(value, float) and not math.isfinite(value):
        return str(value)
    if isinstance(value, dict):
        return {key: native(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [native(item) for item in value]
    return value


def certify_minimum_gap(h, speed, stiffness):
    # P=0 and L=speed imply relative chord velocity (radial_velocity, speed).
    def radius(radial_velocity):
        return math.hypot(1+h*radial_velocity, h*speed)

    def gap(radial_velocity):
        return .5*radial_velocity**2+.5*stiffness*(radius(radial_velocity)-1)**2

    def derivative(radial_velocity):
        current_radius = radius(radial_velocity)
        return radial_velocity+stiffness*(current_radius-1)*h*(1+h*radial_velocity)/current_radius

    # Any global minimum has gap <= gap(0); the kinetic term confines it here.
    radial_bound = math.sqrt(2*gap(0.))
    lower, upper = -radial_bound, radial_bound
    assert 1+h*lower > 0.
    minimum_radius_bound = math.hypot(1+h*lower, h*speed)
    second_derivative_lower_bound = 1+stiffness*h*h*(1-(h*speed)**2/minimum_radius_bound**3)
    assert second_derivative_lower_bound > 0.
    assert derivative(lower) < 0. < derivative(upper)
    minimizer = brentq(derivative, lower, upper, xtol=1e-15, rtol=1e-14)
    # Independent bounded scalar minimization on the globally certified interval.
    bounded = minimize_scalar(gap, bounds=(lower, upper), method='bounded',
                              options={'xatol': 1e-14})
    assert bounded.success
    assert abs(bounded.fun-gap(minimizer)) < 1e-12
    # Also inspect its two halves, preventing reliance on one initial bracket.
    half_minima = [minimize_scalar(gap, bounds=bounds, method='bounded',
                                  options={'xatol': 1e-14})
                   for bounds in ((lower, 0.), (0., upper))]
    return dict(
        radial_velocity_minimizer=minimizer,
        minimum_internal_energy_gap=gap(minimizer),
        zero_radial_velocity_gap=gap(0.),
        global_search_interval=[lower, upper],
        interval_justification='Any global minimizer has .5*radial_velocity^2 <= gap(0)',
        analytic_second_derivative_lower_bound=second_derivative_lower_bound,
        uniqueness_certificate='Strict convexity on interval containing every global minimizer',
        bounded_minimizer=bounded.x, bounded_minimum_gap=bounded.fun,
        half_interval_minima=[dict(success=result.success, x=result.x, gap=result.fun) for result in half_minima],
    )


def fixture(name, h, speed, stiffness):
    scene = cpd_solver.Scene(
        name=name, masses=np.full(2, 2.),
        x0=np.array([[.5, 0.], [-.5, 0.]]),
        v0=np.array([[0., speed/2], [0., -speed/2]]),
        forces=np.zeros((2, 2)), edges=[cpd_solver.Edge(0, 1, 1., stiffness)],
        family='free', energy_floor=0., characteristic_length=1.,
    )
    targets = cpd_solver.make_targets(scene, scene.x0, scene.v0, h)
    assert targets.energy == targets.kstar
    cases = []
    for ratio in (1e4, 1e8):
        epsilon = ratio*targets.energy
        for policy in ('paper', 'fixed', 'stationary'):
            result = cpd_solver.step(scene, scene.x0, scene.v0, h, epsilon,
                                     policy=policy, iterations=16, max_iterations=200,
                                     tol=1e-4, trace=True)
            cases.append(dict(
                epsilon_ratio=ratio, epsilon=epsilon, policy=policy,
                accepted=result.accepted, converged=result.converged,
                stationary=result.stationary, reason=result.reason,
                iterations=result.iterations, x=result.x, v=result.v,
                rho=result.rho, multipliers=result.multipliers,
                measures=result.measures, rank_min=result.rank_min,
                condition_max=result.condition_max,
                max_linear_residual=result.max_linear_residual,
                max_block_displacement=result.max_block_displacement,
                final_evaluation=result.evaluation, trace=result.history,
            ))
    for first, second in zip(cases[:3], cases[3:]):
        assert first['policy'] == second['policy']
        for field in ('x', 'v', 'rho', 'multipliers'):
            assert np.array_equal(first[field], second[field]), (name, field)
    return dict(name=name, h=h, speed=speed, stiffness=stiffness, scene=scene,
                targets=targets, relaxation_energy_derivative=targets.energy-targets.kstar,
                certified_gap=certify_minimum_gap(h, speed, stiffness), cases=cases)


def main():
    data = dict(
        qualification='Standalone exact-feasibility fixtures outside the primary campaign; unchanged validated Eq.21 production solver',
        fixtures=[fixture('rest_tangent_strong', 1/30, 2., 1000.),
                  fixture('rest_tangent_small', 1/60, .8, 100.)],
    )
    path = ROOT/'feasibility_diagnostic_results.json'
    path.write_text(json.dumps(native(data), indent=2, allow_nan=False)+'\n')
    for report in data['fixtures']:
        print(report['name'], 'exact minimum energy gap', report['certified_gap']['minimum_internal_energy_gap'])
        for case in report['cases']:
            print(' ', case['epsilon_ratio'], case['policy'], 'accepted', case['accepted'],
                  'converged', case['converged'], 'reason', case['reason'],
                  'iterations', case['iterations'],
                  'raw residual', case['measures'].get('raw_constraints'),
                  'rank', case['rank_min'])
    print('Wrote', path)


if __name__ == '__main__':
    main()
