"""Reproduce the energy/motion reassessment from existing campaign CSVs.

No simulator, solver, or continuous-reference code is imported or executed.
Requires numpy. Default input is the archived 468-run campaign; --source permits
using a relocated copy of that results directory. Outputs default beside this
script. All motion thresholds are descriptive reporting cutoffs.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import csv
import json
import math
from pathlib import Path

import numpy as np


DEFAULT_SOURCE = Path('/home/blacknoir/Projects/constrained-xpbd/experiments/block_orders_2026_11_01/results')
BANDS = ('energy_within_0p1pct', 'energy_within_1pct', 'energy_within_5pct')
LABELS = {
    'spring_release': 'Spring release', 'spring_orbit': 'Spring orbit',
    'spring_gravity': 'Gravity spring', 'pendulum_release': 'Pendulum release',
    'pendulum_kick': 'Pendulum kick', 'spring_stiff': 'Stiff spring',
}


def write_csv(output, name, rows):
    with (output/name).open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def sign_turns(values, threshold):
    filtered = values[np.abs(values) > threshold]
    return int(np.count_nonzero(filtered[:-1]*filtered[1:] < 0))


def analyze(source):
    summary = list(csv.DictReader((source/'summary.csv').open()))
    metadata = json.loads((source/'metadata.json').read_text())
    scenes = {scene['name']: scene for scene in metadata['scenes']}
    rows = []
    for original in summary:
        scene = scenes[original['case']]
        length, mass = scene['l'], scene['m']
        data = list(csv.DictReader(
            (source/'trajectories'/(original['run_id']+'.csv')).open()))
        kept = [record for record in data if record['accepted'] == '1']
        position = np.array([[float(r['x']), float(r['y'])] for r in kept])
        velocity = np.array([[float(r['vx']), float(r['vy'])] for r in kept])
        times = np.array([float(r['time']) for r in kept])
        speeds = np.linalg.norm(velocity, axis=1)
        radii = np.linalg.norm(position, axis=1)
        steps = np.linalg.norm(np.diff(position, axis=0), axis=1)
        speed_bound = math.sqrt(2*float(original['energy_scale'])/mass)
        tail = times >= max(0, times[-1]-1.0)
        tail_steps = (steps[times[1:] >= max(0, times[-1]-1.0)]
                      if len(steps) else np.array([]))
        tail_path, path = float(tail_steps.sum()), float(steps.sum())
        complete = int(original['completed'])
        category = (
            'guard_failure' if not complete else
            'near_stationary' if path < 1e-4*length else
            'tail_near_stationary' if tail_path < 1e-4*length else
            'large_endpoint_jump' if len(steps) and steps.max() > .5*length else
            'motion_present')
        angles = np.unwrap(np.arctan2(position[:, 0], -position[:, 1]))
        angular = mass*(position[:, 0]*velocity[:, 1]-position[:, 1]*velocity[:, 0])
        radial = np.einsum('ij,ij->i', position, velocity)/radii
        reversals = ((np.einsum('ij,ij->i', velocity[:-1], velocity[1:]) < 0)
                     & (speeds[:-1] > .1*speed_bound)
                     & (speeds[1:] > .1*speed_bound))
        row = dict(original)
        row.update(
            actual_energy_error_percent=100*float(original['max_rel_energy']),
            motion_category=category, accepted_path=path, tail_1s_path=tail_path,
            x_extent=float(np.ptp(position[:, 0])), y_extent=float(np.ptp(position[:, 1])),
            max_distance_from_start=float(np.max(np.linalg.norm(position-position[0], axis=1))),
            min_radius=float(radii.min()), max_radius=float(radii.max()),
            angle_extent=float(np.ptp(angles)), net_angle=float(angles[-1]-angles[0]),
            max_endpoint_speed=float(speeds.max()), energy_based_initial_speed_bound=speed_bound,
            max_endpoint_step=float(steps.max()) if len(steps) else 0,
            max_endpoint_step_over_length=float(steps.max()/length) if len(steps) else 0,
            tail_1s_mean_speed=float(speeds[tail].mean()),
            fast_velocity_reversal_count=int(np.count_nonzero(reversals)),
            angular_velocity_reversal_count=sign_turns(angular, 1e-6*mass*length*speed_bound),
            radial_velocity_reversal_count=sign_turns(radial, 1e-6*speed_bound))
        for name, tolerance in zip(BANDS, (.001, .01, .05)):
            row[name] = int(complete and float(original['max_rel_energy']) <= tolerance)
        rows.append(row)

    by_scene_order = defaultdict(list)
    for row in rows:
        by_scene_order[row['case'], row['order']].append(row)
    bands, strongest = [], []
    for (case, order), group in by_scene_order.items():
        result = dict(case=case, order=order, runs=len(group),
                      completed=sum(int(r['completed']) for r in group),
                      guard_failures=sum(r['motion_category'] == 'guard_failure' for r in group),
                      near_stationary=sum(r['motion_category'] == 'near_stationary' for r in group),
                      tail_near_stationary=sum(r['motion_category'] == 'tail_near_stationary' for r in group))
        for name in BANDS:
            result[name] = sum(r[name] for r in group)
            result[name+'_motion_present'] = sum(
                r[name] and r['motion_category'] == 'motion_present' for r in group)
        bands.append(result)
        for selection, pool in (
                ('best_energy_all', group),
                ('best_energy_motion_present', [r for r in group if r['motion_category'] == 'motion_present'])):
            if not pool:
                continue
            best = min(pool, key=lambda r: float(r['max_rel_energy']))
            keys = ['run_id', 'hz', 'sweeps', 'epsilon_ratio', 'actual_energy_error_percent',
                    'motion_category', 'accepted_path', 'tail_1s_path', 'x_extent', 'y_extent',
                    'angle_extent', 'max_endpoint_speed', 'max_endpoint_step',
                    'fast_velocity_reversal_count', 'max_geometry_error', 'max_tangency_error', 'completed']
            strongest.append(dict(case=case, order=order, selection=selection,
                                  **{key: best[key] for key in keys}))
    grouped_settings = defaultdict(list)
    for row in rows:
        grouped_settings[row['order'], row['hz'], row['sweeps'], row['epsilon_ratio']].append(row)
    settings = []
    for (order, hz, sweeps, epsilon), group in grouped_settings.items():
        result = dict(order=order, hz=hz, sweeps=sweeps, epsilon_ratio=epsilon,
                      scenes=len(group), completed=sum(int(r['completed']) for r in group),
                      motion_present=sum(r['motion_category'] == 'motion_present' for r in group),
                      worst_actual_energy_percent=max(r['actual_energy_error_percent'] for r in group),
                      worst_nonrelease_energy_percent=max(
                          r['actual_energy_error_percent'] for r in group if r['case'] != 'spring_release'))
        for name in BANDS:
            result[name] = sum(r[name] for r in group)
            result[name+'_motion_present'] = sum(
                r[name] and r['motion_category'] == 'motion_present' for r in group)
        for row in group:
            result[row['case']+'_energy_percent'] = row['actual_energy_error_percent']
            result[row['case']+'_motion_category'] = row['motion_category']
        settings.append(result)
    return rows, bands, strongest, settings


def render_audit(rows, bands, settings):
    def selected(case, order='LE', hz=120, sweeps=16, epsilon=1e4):
        return next(row for row in rows if row['case'] == case and row['order'] == order
                    and int(row['hz']) == hz and int(row['sweeps']) == sweeps
                    and float(row['epsilon_ratio']) == epsilon)
    band_lookup = {(r['case'], r['order']): r for r in bands}
    energy_band_table = '\n'.join(
        '| '+label+' | '+' | '.join(' / '.join(str(band_lookup[case, order][name])
                                               for name in BANDS) for order in ('EL', 'LE'))+' |'
        for case, label in LABELS.items())
    shared_settings_table = '\n'.join(
        '| '+label+' | '+' | '.join(f"{selected(case, hz=hz, sweeps=sweeps)['actual_energy_error_percent']:.3f}"
                                   for hz, sweeps in ((120, 1), (120, 4), (60, 16), (120, 16)))+' |'
        for case, label in LABELS.items())
    motion_table = '\n'.join(
        f"| {label} | {selected(case)['accepted_path']:.3f} m | "
        f"{selected(case)['tail_1s_path']:.3f} m | {selected(case)['max_endpoint_step']:.4f} m | "
        f"{selected(case)['angle_extent']:.3f} rad |" for case, label in LABELS.items())
    # Text is intentionally specific to the original six-scene campaign. Check
    # its central claims before emitting that interpretation for a relocated input.
    all_scene_candidates = [r for r in settings if r['energy_within_1pct_motion_present'] == len(LABELS)]
    assert len(all_scene_candidates) == 1
    assert (all_scene_candidates[0]['order'], int(all_scene_candidates[0]['hz']),
            int(all_scene_candidates[0]['sweeps']), float(all_scene_candidates[0]['epsilon_ratio'])) == ('LE', 120, 16, 1e4)
    assert all(selected(case)['fast_velocity_reversal_count'] == 0 for case in LABELS)
    counts = Counter(r['motion_category'] for r in rows)
    assert counts['tail_near_stationary'] == 9
    assert counts['large_endpoint_jump'] == 8
    assert counts['guard_failure'] == 5
    assert selected('spring_stiff')['radial_velocity_reversal_count'] == 11
    return AUDIT_TEMPLATE.format(
        runs=len(rows), energy_band_table=energy_band_table,
        shared_settings_table=shared_settings_table, motion_table=motion_table,
        all_band_counts='/'.join(str(sum(r[name] for r in rows)) for name in BANDS),
        release_radius_min=selected('spring_release')['min_radius'],
        release_radius_max=selected('spring_release')['max_radius'],
        release_max_speed=selected('spring_release')['max_endpoint_speed'],
        release_path=selected('spring_release')['accepted_path'],
        el_orbit_energy=selected('spring_orbit', order='EL', sweeps=1)['actual_energy_error_percent'],
        el_gravity_energy=selected('spring_gravity', order='EL', sweeps=1)['actual_energy_error_percent'],
        el_release_energy=selected('pendulum_release', order='EL', sweeps=1)['actual_energy_error_percent'],
        el_kick_energy=selected('pendulum_kick', order='EL', sweeps=1)['actual_energy_error_percent'],
        el_stiff_energy=selected('spring_stiff', order='EL', sweeps=1)['actual_energy_error_percent'],
        orbit_angular_error=float(selected('spring_orbit')['max_angular_error']),
        stiff_angular_error=float(selected('spring_stiff')['max_angular_error']),
        stiff_angle_extent=selected('spring_stiff')['angle_extent'],
        near_stationary=counts['near_stationary'],
        tail_near_stationary=counts['tail_near_stationary'],
        large_endpoint_jump=counts['large_endpoint_jump'], guard_failures=counts['guard_failure'])


AUDIT_TEMPLATE = "# Actual energy and motion reassessment\n\nThis reassessment reads the completed {runs}-run campaign and its saved\ntrajectories. It runs no solver. Ranking uses actual mechanical-energy error and\ndescribes motion separately; RMS, phase/frequency agreement, residual convergence,\nand stationary curvature are not rejection criteria here.\n\nThe reported energy error is the largest absolute change from initial mechanical\nenergy divided by initial energy above the potential floor, exactly as in the\ncampaign. A failed candidate remains included in the maximum, and its guard\nfailure remains visible. The 0.1%, 1%, and 5% bands below are reporting bands,\nnot a tolerance specified by the user.\n\n## Energy band counts\n\nEach EL or LE cell covers 24 settings. Counts require a complete five-second run\nbut can include nearly stationary motion.\n\n| Scene | EL: within 0.1% / 1% / 5% | LE: within 0.1% / 1% / 5% |\n|---|---:|---:|\n{energy_band_table}\n\nAcross all orders, including JOINT and XPBD controls, {all_band_counts} runs meet these\nthree energy bands. Energy-only counts need the motion qualification below.\n\n## Practical shared settings\n\nThese are maximum actual energy errors, in percent, over the same five seconds.\n\n| Scene | LE, 120 Hz, 1 sweep, ratio 1e4 | LE, 120 Hz, 4 sweeps, ratio 1e4 | LE, 60 Hz, 16 sweeps, ratio 1e4 | LE, 120 Hz, 16 sweeps, ratio 1e4 |\n|---|---:|---:|---:|---:|\n{shared_settings_table}\n\nLE at 120 Hz, 16 sweeps, ratio 1e4 is the only tested shared setting with all six\nscenes complete and below the 1% energy band. Its spring release keeps oscillating:\nradius spans {release_radius_min:.4f} to {release_radius_max:.4f} m, maximum speed is {release_max_speed:.3f} m/s, and total travel is\n{release_path:.3f} m. Release remains the fragile case when lowering sweeps or timestep rate.\n\nAt one sweep, LE is substantially better than EL for the five nonrelease scenes.\nAt 120 Hz and ratio 1e4, EL's errors on those five are {el_orbit_energy:.2f}%, {el_gravity_energy:.2f}%, {el_release_energy:.2f}%,\n{el_kick_energy:.2f}%, and {el_stiff_energy:.2f}%; LE's are shown in the table. This does not make the one-sweep\nmethod successful for spring release.\n\n## What the saved motion says\n\nFor the shared LE/120 Hz/16-sweep/1e4 setting:\n\n| Scene | Travel over 5 s | Travel in last 1 s | Maximum endpoint step | Angular extent |\n|---|---:|---:|---:|---:|\n{motion_table}\n\nThese trajectories have sustained motion and small endpoint displacements. They\nhave no consecutive velocity reversals for which both speeds exceed 10% of the\ninitial energy-based speed bound. This is a descriptive diagnostic; it does not\ncertify plausible mechanics or reject legitimate turning points.\n\nAngular transfer deserves separate scrutiny: maximum angular-momentum change is\n{orbit_angular_error:.3f} for the ordinary orbit and {stiff_angular_error:.3f} for the stiff spring, both starting at\n0.88. The stiff spring progresses through {stiff_angle_extent:.3f} radians and has eleven radial\nvelocity reversals. Its changing balance of radial oscillation and circulation\nis a mechanical effect beyond a scalar energy measurement. The tables preserve\nit for visual assessment rather than declaring every moving run plausible.\n\nVery small energy errors can accompany freezing. For example, JOINT spring\nrelease at 60 Hz/one sweep/1e8 has only 0.000300% energy error but travels\n3.75e-7 m. EL pendulum kick at 120 Hz/16 sweeps/1e8 has 0.001675% error, travels\n0.159 m initially, then only 8.50e-8 m during the last second. Such energy scores\ncannot be interpreted as successful dynamics.\n\n## Machine-readable tables\n\n- `per_run_energy_and_motion.csv`: all {runs} original summary records plus travel,\n  extent, speed, endpoint-step, tail-motion, and reversal diagnostics. Original\n  RMS remains available but was not used for selection.\n- `energy_bands_by_scene_order.csv`: band counts, guard failures, and coarse\n  stationary flags for every scene/order, including JOINT and XPBD controls.\n- `strongest_by_scene_order.csv`: lowest-energy setting, and lowest-energy setting\n  without a coarse stationary/jump flag, for each scene/order. The latter is\n  named `motion_present`; it is not an automatic plausibility verdict.\n- `cross_scene_settings.csv`: all six scenes aligned by shared solver settings.\n\nThe coarse flags are explicit descriptors: less than 1e-4 rest lengths of total\nor last-second travel, or an endpoint step larger than half a rest length.\nThere are {near_stationary} near-stationary runs, nine stationary-tail runs, eight runs with a\nlarge endpoint step, and five guard failures. These flags aid inspection; their\ncutoffs are not user-specified acceptance criteria. All five guard failures are\nretained in the tables.\n\n## Reproduce without solver execution\n\n```bash\npython3 audit_energy_motion.py\n```\n\nUse `--source /path/to/existing/results` for a relocated copy of the campaign,\nand `--output /path/to/report/directory` to write the five outputs elsewhere.\nThe script imports NumPy and reads saved CSV/JSON files; it never imports or\nexecutes the numerical solver. It checks the report's central campaign-specific\nclaims before rendering the interpretation.\n"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, default=DEFAULT_SOURCE,
                        help='Existing results directory containing metadata, summary and trajectories')
    parser.add_argument('--output', type=Path, default=Path(__file__).resolve().parent)
    args = parser.parse_args()
    rows, bands, strongest, settings = analyze(args.source)
    report = render_audit(rows, bands, settings)
    args.output.mkdir(parents=True, exist_ok=True)
    for name, data in (
            ('per_run_energy_and_motion.csv', rows),
            ('energy_bands_by_scene_order.csv', bands),
            ('strongest_by_scene_order.csv', strongest),
            ('cross_scene_settings.csv', settings)):
        write_csv(args.output, name, data)
    (args.output/'AUDIT.md').write_text(report)
    print(json.dumps(dict(source=str(args.source.resolve()), output=str(args.output.resolve()),
                          runs=len(rows), files=5, solver_executed=False), indent=2))


if __name__ == '__main__':
    main()
