#!/usr/bin/env python3
"""Exact nonlinear E -> unchanged XPBD L, with explicit row traces.

See SCALAR_MATH.md for the objective, root policy, and validation boundaries.
Only Python's standard library is needed by the experiment itself.
"""
import argparse
import csv
import math
from pathlib import Path
import sys
import time

EPS = sys.float_info.epsilon


def energy(q, v, m, k):
    return 0.5 * (m * v * v + k * q * q)


def endpoint_energy(q, old_q, h, m, k):
    return energy(q, (q - old_q) / h, m, k)


def energy_roots(old_q, h, target, m, k):
    sigma = k * h * h / m
    center = old_q / (1.0 + sigma)
    minimum = 0.5 * k * old_q * old_q / (1.0 + sigma)
    excess = target - minimum
    tolerance = 64 * EPS * max(abs(target), abs(minimum), sys.float_info.min)
    if excess < -tolerance:
        return None
    radius = math.sqrt(2 * max(0.0, excess) / (m / (h * h) + k))
    return center - radius, center + radius


def nearest_root(roots, current, predictor):
    distances = [abs(r - current) for r in roots]
    tolerance = 64 * EPS * max(abs(current), abs(predictor), *(abs(r) for r in roots), sys.float_info.min)
    ambiguous = abs(distances[0] - distances[1]) <= tolerance
    if ambiguous:
        index = min(range(2), key=lambda i: (abs(roots[i] - predictor), roots[i]))
    else:
        index = min(range(2), key=lambda i: distances[i])
    return roots[index], index, ambiguous


def objective_root(roots, old_q, old_v, h, m, k):
    predictor = old_q + h * old_v
    values = [0.5 * m * ((r - predictor) / h) ** 2 + 0.5 * k * r * r for r in roots]
    tolerance = 64 * EPS * max(*(abs(f) for f in values), sys.float_info.min)
    ambiguous = abs(values[0] - values[1]) <= tolerance
    if ambiguous:
        index = min(range(2), key=lambda i: (abs(roots[i] - predictor), roots[i]))
    else:
        index = min(range(2), key=lambda i: values[i])
    return roots[index], index, ambiguous


def material_visit(q, lam, h, m, k):
    alpha = 1.0 / (k * h * h)
    delta_lam = -(q + alpha * lam) / (1.0 / m + alpha)
    return q + delta_lam / m, lam + delta_lam


def validate():
    checks = 0
    for h in [1 / 60, 1 / 120, 1 / 480]:
        for old_q in [-.3, 0.0, .25]:
            for old_v in [-.8, 0.0, .8]:
                target = energy(old_q, old_v, 1, 100)
                roots = energy_roots(old_q, h, target, 1, 100)
                assert roots is not None
                for root in roots:
                    assert abs(endpoint_energy(root, old_q, h, 1, 100) - target) < 2e-12 * max(1, target)
                    checks += 1
                chosen, _, _ = objective_root(roots, old_q, old_v, h, 1, 100)
                fs = [0.5 * ((r - old_q - h * old_v) / h) ** 2 + 50 * r * r for r in roots]
                fc = 0.5 * ((chosen - old_q - h * old_v) / h) ** 2 + 50 * chosen * chosen
                assert fc <= min(fs) + 2e-12 * max(1, target)
                checks += 1
        q, lam = material_visit(.24, -.0001, h, 1, 100)
        assert abs(q + lam / (100 * h * h)) < 1e-14
        checks += 1
        roots = energy_roots(.25, h, 3.125, 1, 100)
        selected, _, tie = objective_root(roots, .25, 0, h, 1, 100)
        assert tie and abs(selected - .25) < 1e-15
        checks += 1
        q, v, target = .25, 0, 3.125
        for _ in range(2):
            old_q = q
            q = q + h * v
            residual = endpoint_energy(q, old_q, h, 1, 100) - target
            gradient = (q - old_q) / (h * h) + 100 * q
            q -= residual / gradient
            q, _ = material_visit(q, 0, h, 1, 100)
            v = (q - old_q) / h
        sigma = 100 * h * h
        expected = -.25 * (.5 + sigma) / (1 + sigma) ** 2
        assert abs(q - expected) < 2e-8
        checks += 1
    assert energy_roots(.25, 1 / 480, 0, 1, 100) is None
    checks += 1
    return checks


def run(case, initial_v, hz, sweeps, method, target_mode, duration, states_writer, visits_writer):
    m, k, initial_q = 1.0, 100.0, .25
    h = 1.0 / hz
    omega = math.sqrt(k / m)
    initial_energy = energy(initial_q, initial_v, m, k)
    amplitude = math.sqrt(2 * initial_energy / k)
    q, v = initial_q, initial_v
    metadata = dict(case=case, method=method, target_mode=target_mode, hz=hz, sweeps=sweeps)
    stats = dict(status='complete', completed_steps=0, accepted_time=0.0, attempted_step=0,
                 max_energy_relative=0.0, final_energy_relative=0.0,
                 max_post_E_target_relative=0.0, max_post_L_target_relative=0.0,
                 max_E_displacement=0.0, max_L_displacement=0.0,
                 no_root=0, ambiguous_roots=0, minus_roots=0, plus_roots=0,
                 within_step_branch_changes=0, zero_gradient_skips=0,
                 row_visits=0, path_length=0.0, frozen_steps=0)
    square_q_error = square_v_error = 0.0
    begin = time.perf_counter()
    total_steps = round(duration * hz)
    for step in range(1, total_steps + 1):
        old_q, old_v = q, v
        predictor = old_q + h * old_v
        target = initial_energy if target_mode == 'initial' else energy(old_q, old_v, m, k)
        lam, eta = 0.0, 0.0
        q = predictor
        branch_before = None
        success = True
        for sweep in range(1 if method == 'objective_exact' else sweeps):
            before_E = q
            ambiguous = False
            branch = ''
            if method in ('energy_exact', 'objective_exact'):
                roots = energy_roots(old_q, h, target, m, k)
                if roots is None:
                    stats['no_root'] += 1
                    stats['status'] = 'no_real_energy_root'
                    success = False
                    break
                if method == 'objective_exact':
                    q, branch, ambiguous = objective_root(roots, old_q, old_v, h, m, k)
                else:
                    q, branch, ambiguous = nearest_root(roots, q, predictor)
                stats['ambiguous_roots'] += int(ambiguous)
                stats['minus_roots' if branch == 0 else 'plus_roots'] += 1
                if branch_before is not None and branch != branch_before:
                    stats['within_step_branch_changes'] += 1
                branch_before = branch
            else:
                residual = endpoint_energy(q, old_q, h, m, k) - target
                kinetic_gradient = m * (q - old_q) / (h * h)
                spring_gradient = k * q
                gradient = kinetic_gradient + spring_gradient
                tolerance = 64 * EPS * max(abs(kinetic_gradient) + abs(spring_gradient), 1.0)
                tolerance += 4 * EPS * m / (h * h) * (abs(q) + abs(old_q))
                if abs(gradient) <= tolerance:
                    stats['zero_gradient_skips'] += 1
                else:
                    delta_eta = -residual / (gradient * gradient / m)
                    q += gradient * delta_eta / m
                    eta += delta_eta
            stats['row_visits'] += 1
            post_E = q
            energy_E = endpoint_energy(q, old_q, h, m, k)
            stats['max_E_displacement'] = max(stats['max_E_displacement'], abs(q - before_E))
            stats['max_post_E_target_relative'] = max(stats['max_post_E_target_relative'], abs(energy_E - target) / initial_energy)
            before_L = q
            if method != 'objective_exact':
                q, lam = material_visit(q, lam, h, m, k)
                stats['row_visits'] += 1
            energy_L = endpoint_energy(q, old_q, h, m, k)
            stats['max_L_displacement'] = max(stats['max_L_displacement'], abs(q - before_L))
            stats['max_post_L_target_relative'] = max(stats['max_post_L_target_relative'], abs(energy_L - target) / initial_energy)
            if step <= 8:
                visits_writer.writerow({**metadata, 'step': step, 'sweep': sweep + 1,
                    'old_q': old_q, 'old_v': old_v, 'predictor': predictor, 'target': target,
                    'before_E': before_E, 'post_E': post_E, 'post_L': q,
                    'energy_E': energy_E, 'energy_L': energy_L, 'lambda': lam,
                    'eta': eta if method == 'energy_linear' else '', 'branch': branch,
                    'ambiguous': int(ambiguous)})
            if not all(math.isfinite(a) for a in (q, lam, eta, energy_L)) or abs(q) > 1e4 or abs((q - old_q) / h) > 1e6:
                stats['status'] = 'runaway_or_nonfinite'
                success = False
                break
        stats['attempted_step'] = step
        if not success:
            q, v = old_q, old_v
            break
        v = (q - old_q) / h
        t = step * h
        ref_q = initial_q * math.cos(omega * t) + initial_v / omega * math.sin(omega * t)
        ref_v = -omega * initial_q * math.sin(omega * t) + initial_v * math.cos(omega * t)
        physical_energy = energy(q, v, m, k)
        relative_energy = (physical_energy - initial_energy) / initial_energy
        stats['max_energy_relative'] = max(stats['max_energy_relative'], abs(relative_energy))
        stats['final_energy_relative'] = relative_energy
        stats['completed_steps'] = step
        stats['accepted_time'] = t
        stats['path_length'] += abs(q - old_q)
        stats['frozen_steps'] += int(abs(q - old_q) <= 1e-12 * amplitude)
        square_q_error += (q - ref_q) ** 2
        square_v_error += (v - ref_v) ** 2
        states_writer.writerow({**metadata, 'step': step, 'time': t, 'q': q, 'v': v,
            'energy': physical_energy, 'target': target, 'reference_q': ref_q,
            'reference_v': ref_v, 'lambda': lam,
            'material_residual': q + lam / (k * h * h) if method != 'objective_exact' else '',
            'last_E_energy': energy_E, 'last_L_energy': energy_L})
    count = stats['completed_steps']
    stats['position_rms'] = math.sqrt(square_q_error / count) if count else math.nan
    stats['velocity_rms'] = math.sqrt(square_v_error / count) if count else math.nan
    stats['position_rms_over_amplitude'] = stats['position_rms'] / amplitude
    stats['velocity_rms_over_velocity_amplitude'] = stats['velocity_rms'] / (omega * amplitude)
    stats['final_q'], stats['final_v'] = q, v
    stats['cpu_seconds'] = time.perf_counter() - begin
    return {**metadata, **stats}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, default=Path(__file__).parent / 'scalar_results')
    parser.add_argument('--duration', type=float, default=5.0)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    checks = validate()
    metadata = ['case', 'method', 'target_mode', 'hz', 'sweeps']
    states_fields = metadata + ['step', 'time', 'q', 'v', 'energy', 'target', 'reference_q', 'reference_v', 'lambda', 'material_residual', 'last_E_energy', 'last_L_energy']
    visits_fields = metadata + ['step', 'sweep', 'old_q', 'old_v', 'predictor', 'target', 'before_E', 'post_E', 'post_L', 'energy_E', 'energy_L', 'lambda', 'eta', 'branch', 'ambiguous']
    results = []
    with (args.output / 'states.csv').open('w', newline='') as sf, (args.output / 'visits_first8.csv').open('w', newline='') as vf:
        sw, vw = csv.DictWriter(sf, states_fields), csv.DictWriter(vf, visits_fields)
        sw.writeheader()
        vw.writeheader()
        for case, initial_v in [('release', 0.0), ('mixed', .8)]:
            for target_mode in ['initial', 'previous_step']:
                for hz in [60, 120, 240, 480]:
                    for method in ['energy_linear', 'energy_exact', 'objective_exact']:
                        for sweeps in ([1] if method == 'objective_exact' else [1, 4, 16]):
                            results.append(run(case, initial_v, hz, sweeps, method, target_mode, args.duration, sw, vw))
    with (args.output / 'summary.csv').open('w', newline='') as f:
        writer = csv.DictWriter(f, list(results[0]))
        writer.writeheader()
        writer.writerows(results)
    report = ['# Scalar exact-energy experiment results', '',
              f'{checks} analytic validation checks passed. {len(results)} runs, requested duration {args.duration:g} s.', '',
              'Every E visit is followed immediately by the unchanged compliant XPBD L visit, with accumulated material multiplier. The exact-objective control solves the original scalar constrained objective directly.', '',
              '## 480 Hz results', '',
              '| Case | Target | Method | Sweeps | Status/time | Max energy error | Position RMS | Path length | Ambiguous roots |',
              '|---|---|---|---:|---|---:|---:|---:|---:|']
    for r in results:
        if r['hz'] == 480:
            report.append(f"| {r['case']} | {r['target_mode']} | {r['method']} | {r['sweeps']} | {r['status']}, {r['accepted_time']:.5g}s | {r['max_energy_relative']:.6g} | {r['position_rms']:.6g} | {r['path_length']:.6g} | {r['ambiguous_roots']} |")
    report += ['', 'Energy errors are fractions of initial physical energy, not percentages. Failed runs report their accepted prefix and must not be ranked as complete trajectories.', '',
               '`summary.csv` contains every run; `states.csv` contains every accepted timestep; `visits_first8.csv` contains every row pair for the first eight timesteps of each run, including energy immediately after E and after L. Exact projections do not have a claimed KKT-consistent accumulated energy multiplier.', '']
    (args.output / 'RESULTS.md').write_text('\n'.join(report))
    print(f'validation_checks={checks} runs={len(results)} output={args.output}')
    for method in ['energy_linear', 'energy_exact', 'objective_exact']:
        rows = [r for r in results if r['method'] == method]
        print(method, 'complete', sum(r['status'] == 'complete' for r in rows), '/', len(rows), 'no_root', sum(r['no_root'] for r in rows))


if __name__ == '__main__':
    main()
