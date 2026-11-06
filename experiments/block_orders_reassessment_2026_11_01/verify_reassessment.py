"""Independently recompute physical metrics for the selected recorded runs.

Reads recorded states only; does not import or execute the numerical solver.
The scene constants below are the six campaign definitions. No reference
positions, residual convergence, or stationary curvature enter these checks.
"""
import argparse
import csv
import json
import math
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, default=Path(
        '/home/blacknoir/Projects/constrained-xpbd/experiments/block_orders_2026_11_01/results'))
    parser.add_argument('--output', type=Path, default=Path(__file__).resolve().parent)
    args = parser.parse_args()
    summaries = list(csv.DictReader((args.source / 'summary.csv').open()))
    constants = {
        'spring_release': (100., 0., False),
        'spring_orbit': (100., 0., False),
        'spring_gravity': (100., 9.81, False),
        'pendulum_release': (0., 9.81, True),
        'pendulum_kick': (0., 9.81, True),
        'spring_stiff': (1000., 0., False),
    }
    selected = []
    errors = dict(energy=0., normalized_energy=0., energy_maximum=0.,
                  chord_velocity=0., path_length=0., rod_length=0.)
    checked_states = 0
    for summary in summaries:
        hz = int(summary['hz'])
        ratio = float(summary['epsilon_ratio'])
        if (summary['order'] not in ('EL', 'LE', 'JOINT')
                or int(summary['sweeps']) != 16
                or (hz, ratio) not in ((120, 1e4), (60, 1e8))):
            continue
        rows = list(csv.DictReader((args.source / 'trajectories' /
                                   (summary['run_id'] + '.csv')).open()))
        k, gravity, rigid = constants[summary['case']]
        energy0 = float(summary['initial_energy'])
        scale = float(summary['energy_scale'])
        max_energy = 0.
        max_geometry = 0.
        path = 0.
        previous = None
        accepted_times = []
        accepted_speeds = []
        for row in rows:
            x, y, vx, vy = (float(row[name]) for name in ('x', 'y', 'vx', 'vy'))
            radius = math.hypot(x, y)
            physical = .5 * (vx * vx + vy * vy) + gravity * y
            if not rigid:
                physical += .5 * k * (radius - 1.) ** 2
            rel = (physical - energy0) / scale
            errors['energy'] = max(errors['energy'], abs(physical - float(row['energy'])))
            errors['normalized_energy'] = max(errors['normalized_energy'], abs(rel - float(row['rel_energy'])))
            max_energy = max(max_energy, abs(rel))
            if rigid:
                max_geometry = max(max_geometry, abs(radius - 1.))
            accepted = row['accepted'] == '1'
            if previous is not None:
                chord = ((x - previous[0]) * hz, (y - previous[1]) * hz)
                errors['chord_velocity'] = max(errors['chord_velocity'],
                                              math.hypot(vx - chord[0], vy - chord[1]))
                if accepted:
                    path += math.hypot(x - previous[0], y - previous[1])
            if accepted:
                previous = (x, y)
                accepted_times.append(float(row['time']))
                accepted_speeds.append(math.hypot(vx, vy))
            checked_states += 1
        errors['energy_maximum'] = max(errors['energy_maximum'],
                                      abs(max_energy - float(summary['max_rel_energy'])))
        errors['path_length'] = max(errors['path_length'],
                                   abs(path - float(summary['path_length'])))
        errors['rod_length'] = max(errors['rod_length'],
                                  abs(max_geometry - float(summary['max_geometry_error'])))
        tail = [v for t, v in zip(accepted_times, accepted_speeds) if t >= 4.]
        selected.append(dict(
            run_id=summary['run_id'], scene=summary['case'], order=summary['order'],
            hz=hz, sweeps=16, epsilon_ratio=ratio,
            completed=int(summary['completed']), max_energy_error_percent=100 * max_energy,
            path_m=path, max_rod_length_error_m=max_geometry,
            tail_speed_max_m_per_s=max(tail, default=math.nan),
            tail_speed_rms_m_per_s=math.sqrt(sum(v*v for v in tail) / len(tail)) if tail else math.nan))
    # Independent arithmetic order can differ slightly for extreme rejected states.
    assert len(selected) == 36, len(selected)
    assert errors['energy'] < 1e-8, errors
    assert errors['normalized_energy'] < 1e-8, errors
    assert errors['energy_maximum'] < 1e-8, errors
    assert errors['chord_velocity'] < 1e-8, errors
    assert errors['path_length'] < 1e-8, errors
    assert errors['rod_length'] < 1e-10, errors
    args.output.mkdir(parents=True, exist_ok=True)
    with (args.output / 'selected_metrics.csv').open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(selected[0]))
        writer.writeheader()
        writer.writerows(selected)
    audit = dict(selected_runs=len(selected), checked_states=checked_states,
                 maximum_absolute_discrepancies=errors,
                 checks='Physical energy, run-initial normalized energy maximum, chord velocity, path, rod length.',
                 source=str(args.source), solver_executed=False)
    (args.output / 'VERIFICATION.json').write_text(json.dumps(audit, indent=2) + '\n')
    print(json.dumps(audit, indent=2))


if __name__ == '__main__':
    main()
