"""Compare the known position-only EL failure from exactly the same incoming state."""
from __future__ import annotations

import argparse
import csv
import importlib.util
import json
import math
from pathlib import Path
import sys

import numpy as np

import independent_solver


ROOT = Path(__file__).resolve().parent
DEFAULT_OLD_SOURCE = Path('/home/blacknoir/Projects/constrained-xpbd/experiments/block_orders_2026_11_01/solver.py')


def native(value):
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return native(value.item())
    if isinstance(value, float) and not math.isfinite(value):
        return str(value)
    if isinstance(value, dict):
        return {key: native(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [native(item) for item in value]
    return value


def load_old_solver(path):
    specification = importlib.util.spec_from_file_location('archived_position_solver_diagnostic', path)
    module = importlib.util.module_from_spec(specification)
    sys.modules[specification.name] = module
    specification.loader.exec_module(module)
    return module


def one_sweep(module, beta_v, x0, v0, h, epsilon):
    scene = module.scenes()[0]
    x, v = (x0+h*v0).copy(), v0.copy()
    lam = rho = eta = 0.
    trace = []
    independent = beta_v is not None

    def evaluate():
        if independent:
            return module.evaluate_state(scene, x0, v0, h, x, v, lam, rho, eta, epsilon, beta_v)
        return module.evaluate_state(scene, x0, v0, h, x, lam, rho, eta, epsilon)

    for block in ('E', 'L'):
        evaluation = evaluate()
        before = dict(x=x.copy(), v=v.copy(), lam=lam, rho=rho, eta=eta,
                      evaluation=vars(evaluation).copy())
        if independent:
            if block == 'E':
                change = module.energy_increment(scene, h, epsilon, beta_v, evaluation)
            else:
                change = module.material_increment(scene, h, beta_v, evaluation)
            rx_response, rv_response = module.inverse_primal(scene, h, beta_v, evaluation.r_x, evaluation.r_v)
            position_from_rx = -h*h*evaluation.r_x/scene.m
            position_from_rv = -h*evaluation.r_v/scene.m
            combined = h*evaluation.potential_gradient+scene.m*v
            denominator = (combined@combined/scene.m+beta_v*scene.m*(v@v)
                           +(evaluation.target-evaluation.kstar)**2/epsilon)
            numerator = (evaluation.r_eta-evaluation.potential_gradient@rx_response
                         -(scene.m*v)@rv_response
                         -(evaluation.target-evaluation.kstar)*evaluation.r_rho/epsilon)
        else:
            if block == 'E':
                change = module.energy_increment(scene, h, epsilon, evaluation)
            else:
                change = module.material_increment(scene, h, evaluation)
            position_from_rx = -h*h*evaluation.r_x/scene.m
            position_from_rv = np.zeros_like(x)
            denominator = h*h/scene.m*(evaluation.energy_gradient@evaluation.energy_gradient)+(evaluation.target-evaluation.kstar)**2/epsilon
            numerator = (evaluation.r_eta-h*h/scene.m*(evaluation.energy_gradient@evaluation.r_x)
                         -(evaluation.target-evaluation.kstar)*evaluation.r_rho/epsilon)
        x = x+change.dx
        if independent:
            v = v+change.dv
        else:
            v = (x-x0)/h
        lam += change.dlam
        rho += change.drho
        eta += change.deta
        after = evaluate()
        trace.append(dict(
            block=block, before=before, increment=vars(change).copy(),
            current_energy_numerator=float(numerator),
            current_energy_denominator=float(denominator),
            position_response_to_r_x=position_from_rx,
            position_response_to_r_v=position_from_rv,
            after=dict(x=x.copy(), v=v.copy(), lam=lam, rho=rho, eta=eta,
                       evaluation=vars(after).copy()),
        ))
        if (np.linalg.norm(x) > 1e4*scene.l or np.linalg.norm(v) > 1e6*scene.l
                or np.linalg.norm((x-x0)/h) > 1e6*scene.l):
            raise FloatingPointError('runaway_guard')

    if independent:
        public = module.block_step(scene, x0, v0, h, epsilon, beta_v, 'EL', 1, trace=True)
        control = module.block_step(scene, x0, v0, h, epsilon, beta_v, 'EL', 16, trace=True)
    else:
        public = module.block_step(scene, x0, v0, h, epsilon, 'EL', 1, trace=True)
        control = module.block_step(scene, x0, v0, h, epsilon, 'EL', 16, trace=True)
    assert public.accepted
    assert np.array_equal(public.x, x) and np.array_equal(public.v, v)
    assert public.lam == lam and public.eta == eta and public.rho == rho

    first = trace[0]
    change = first['increment']
    after = first['after']['evaluation']
    if independent:
        expected_rx = change['deta']*scene.k*change['dx']
        expected_rv = change['deta']*scene.m*change['dv']
        expected_energy = .5*scene.k*(change['dx']@change['dx'])+.5*scene.m*(change['dv']@change['dv'])
        remainder = dict(
            exact_r_x_remainder=expected_rx, exact_r_v_remainder=expected_rv,
            exact_energy_remainder=float(expected_energy),
            r_x_identity_error=float(np.linalg.norm(after['r_x']-expected_rx)),
            r_v_identity_error=float(np.linalg.norm(after['r_v']-expected_rv)),
            energy_identity_error=abs(after['r_eta']-expected_energy),
        )
    else:
        expected_rx = change['deta']*(scene.m/(h*h)+scene.k)*change['dx']
        expected_energy = .5*(scene.m/(h*h)+scene.k)*(change['dx']@change['dx'])
        remainder = dict(
            exact_r_x_remainder=expected_rx, exact_energy_remainder=float(expected_energy),
            r_x_identity_error=float(np.linalg.norm(after['r_x']-expected_rx)),
            energy_identity_error=abs(after['r_eta']-expected_energy),
        )
    return dict(
        formulation='independent_velocity' if independent else 'position_only_derived_velocity',
        beta_v=beta_v,
        production_api_parity='Bitwise identical one-sweep endpoint state and multipliers',
        extra_cold_energy_denominator_from_independent_velocity=(beta_v*scene.m*float(v0@v0) if independent else 0.),
        one_EL_sweep=dict(x=x, v=v, lam=lam, rho=rho, eta=eta,
                          physical_energy=module.energy(scene, x, v), trace=trace),
        first_energy_visit_exact_remainders=remainder,
        sixteen_EL_sweeps=dict(x=control.x, v=control.v, lam=control.lam,
                              rho=control.rho, eta=control.eta, accepted=control.accepted,
                              converged=control.converged, reason=control.reason,
                              physical_energy=module.energy(scene, control.x, control.v),
                              max_internal_position_increment=control.max_block_displacement,
                              evaluation=(vars(control.evaluation) if control.evaluation is not None else None),
                              trace=control.history),
    )


def saved_factor_screen():
    summary_path = ROOT/'results'/'summary.csv'
    if not summary_path.exists():
        return dict(status='No saved campaign CSV available')
    with summary_path.open(newline='') as stream:
        summaries = list(csv.DictReader(stream))
    selected = [row for row in summaries if float(row['beta_v']) == 1.
                and int(row['hz']) == 60 and int(row['sweeps']) in (4, 16)
                and row['order'] in ('EL', 'LE', 'JOINT')]
    records = []
    for summary in selected:
        path = ROOT/'results'/'trajectories'/(summary['run_id']+'.csv')
        with path.open(newline='') as stream:
            rows = [row for row in csv.DictReader(stream) if int(row['accepted']) and int(row['step']) > 0]
        if not rows:
            continue
        factor = np.array([1+float(row['eta']) for row in rows])
        speed = np.array([float(row['speed']) for row in rows])
        chord_speed = np.array([float(row['chord_speed']) for row in rows])
        times = np.array([float(row['time']) for row in rows])
        tail = times >= times[-1]-1.
        velocity = np.array([[float(row['vx']), float(row['vy'])] for row in rows])
        chord_velocity = np.array([[float(row['chord_vx']), float(row['chord_vy'])] for row in rows])
        identity_defect = np.linalg.norm(chord_velocity-factor[:, None]*velocity, axis=1)
        records.append(dict(
            run_id=summary['run_id'], completed=bool(int(summary['completed'])),
            minimum_1_plus_beta_eta=float(factor.min()), maximum_1_plus_beta_eta=float(factor.max()),
            fraction_negative_factor=float(np.mean(factor < 0)),
            fraction_abs_factor_below_0_1=float(np.mean(np.abs(factor) < .1)),
            tail_stalled=bool(int(summary['tail_stalled'])),
            max_accepted_energy_error=float(summary['max_accepted_rel_energy']),
            max_velocity_chord_mismatch=float(summary['max_velocity_chord_mismatch']),
            reversed_velocity_fraction=float(summary['reversed_velocity_fraction']),
            max_residual=float(summary['max_residual']),
            tail_max_speed=float(np.max(speed[tail])),
            tail_max_chord_speed=float(np.max(chord_speed[tail])),
            tail_minimum_factor=float(np.min(factor[tail])),
            tail_maximum_factor=float(np.max(factor[tail])),
            tail_max_velocity_stationarity_residual=max(float(row['r_v_norm']) for row, use in zip(rows, tail) if use),
            max_factor_identity_defect=float(np.max(identity_defect)),
            qualification='Factor identity holds at stationarity; large-residual iterates need not follow it',
        ))
    return dict(status='Read-only snapshot of available completed summary rows; no new campaign solves',
                selected_run_count=len(records), runs=records)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--old-source', type=Path, default=DEFAULT_OLD_SOURCE)
    parser.add_argument('--output', type=Path, default=ROOT/'diagnostics.json')
    args = parser.parse_args()
    old = load_old_solver(args.old_source)
    x0 = np.array([1.2432432432432432, 0.])
    v0 = np.array([-0.4054054054054079, 0.])
    h, epsilon = 1/60, 31250.
    cases = [one_sweep(old, None, x0, v0, h, epsilon)]
    cases += [one_sweep(independent_solver, beta, x0, v0, h, epsilon) for beta in (.1, 1., 10.)]
    data = dict(
        incoming_x=x0, incoming_v=v0, h=h, epsilon=epsilon,
        energy_target=old.energy(old.scenes()[0], x0, v0),
        initialization='x=x_n+h*v_n, independent v=v_n, lambda=eta=rho=0',
        old_source=str(args.old_source),
        qualification='Exactly identical incoming state; one EL sweep and bounded 16-sweep control only; no global stability proof',
        cases=cases, saved_beta1_factor_screen=saved_factor_screen(),
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(native(data), indent=2, allow_nan=False)+'\n')
    for case in cases:
        energy_visit = case['one_EL_sweep']['trace'][0]
        material_visit = case['one_EL_sweep']['trace'][1]
        print(case['formulation'], case['beta_v'],
              'E denominator', energy_visit['current_energy_denominator'],
              'deta', energy_visit['increment']['deta'],
              'E dx', np.linalg.norm(energy_visit['increment']['dx']),
              'L dx', np.linalg.norm(material_visit['increment']['dx']),
              'one-sweep energy', case['one_EL_sweep']['physical_energy'],
              '16-sweep energy', case['sixteen_EL_sweeps']['physical_energy'])
    print('Wrote', args.output)


if __name__ == '__main__':
    main()
