"""Read saved trajectories for structural motion signatures; run no solves."""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np

DEFAULT_SOURCE = Path('/home/blacknoir/Projects/constrained-xpbd/experiments/block_orders_2026_11_01/results')


def energy_floor(scene):
    force = scene.m*float(np.linalg.norm(scene.g))
    return -force*scene.l-(0.0 if scene.rigid else force*force/(2*scene.k))


def motion_metrics(scene, summary, source):
    path = source/'trajectories'/(summary['run_id']+'.csv')
    with path.open(newline='') as stream:
        records = [row for row in csv.DictReader(stream) if int(row['accepted'])]
    columns = ('time', 'x', 'y', 'vx', 'vy', 'energy', 'eta')
    values = np.array([[float(row[key]) for key in columns] for row in records])
    time, position, velocity = values[:, 0], values[:, 1:3], values[:, 3:5]
    initial_energy_above_floor = values[0, 5]-energy_floor(scene)
    accepted_energy_error = (values[:, 5]-values[0, 5])/initial_energy_above_floor
    radius = np.linalg.norm(position, axis=1)
    normal = position/radius[:, None]
    speed = np.linalg.norm(velocity, axis=1)
    h = 1/float(summary['hz'])
    acceleration = np.diff(velocity, axis=0)/h
    chord = np.diff(position, axis=0)
    angular_momentum = scene.m*(position[:, 0]*velocity[:, 1]-position[:, 1]*velocity[:, 0])
    angle = np.unwrap(np.arctan2(position[:, 0], -position[:, 1]))
    angle_increments = np.diff(angle)
    angular_signs = np.sign(angle_increments[np.abs(angle_increments) > 1e-7])
    radial_velocity = np.sum(normal*velocity, axis=1)
    if scene.rigid:
        drive = scene.g[None, :]-np.sum(normal*scene.g, axis=1)[:, None]*normal
    else:
        drive = scene.g[None, :]-(scene.k/scene.m)*(radius-scene.l)[:, None]*normal
    drive_norm = np.linalg.norm(drive, axis=1)
    tail = time >= max(0., time[-1]-1)
    future_min = np.minimum.accumulate(radius[::-1])[::-1]
    future_max = np.maximum.accumulate(radius[::-1])[::-1]
    settle = {}
    for tolerance in (1e-3, 1e-4):
        index = int(np.flatnonzero(future_max-future_min < tolerance*scene.l)[0])
        settle[f'remaining_radius_range_below_{tolerance:g}_ell'] = float(time[index])
    future_max_speed = np.maximum.accumulate(speed[::-1])[::-1]
    stalled = np.flatnonzero((future_max_speed < 1e-4*scene.l) & (drive_norm > 1.))
    # This is an exact algebraic identity for reconstructed chord velocity.
    chord_identity_error = np.sum(position[1:]*velocity[1:], axis=1)-(
        radius[1:]**2-radius[:-1]**2+np.sum(chord*chord, axis=1))/(2*h)
    peak_acceleration_index = int(np.argmax(np.linalg.norm(acceleration, axis=1)))
    return dict(
        run_id=summary['run_id'], completed=bool(int(summary['completed'])),
        accepted_time=float(time[-1]), accepted_steps=len(records)-1,
        max_energy_error_above_floor=float(np.max(np.abs(accepted_energy_error))),
        final_energy_error_above_floor=float(accepted_energy_error[-1]),
        summary_max_energy_error_including_failed_candidate=float(summary['max_rel_energy']),
        radius_min=float(radius.min()), radius_max=float(radius.max()),
        last_second_radius_range=float(np.ptp(radius[tail])),
        speed_min=float(speed.min()), speed_max=float(speed.max()),
        last_second_speed_mean=float(np.mean(speed[tail])),
        last_second_speed_max=float(np.max(speed[tail])),
        last_second_unbalanced_drive_mean=float(np.mean(drive_norm[tail])),
        persistent_stall_below_1em4_ell_per_second_with_drive_gt1_time=(float(time[stalled[0]]) if len(stalled) else None),
        final_angle_from_downward_vertical=float(angle[-1]),
        angular_direction_reversals_above_1em7_rad=int(np.sum(angular_signs[1:] != angular_signs[:-1])),
        initial_angular_momentum=float(angular_momentum[0]),
        final_angular_momentum=float(angular_momentum[-1]),
        angular_momentum_min=float(angular_momentum.min()),
        angular_momentum_max=float(angular_momentum.max()),
        max_endpoint_displacement=float(np.max(np.linalg.norm(chord, axis=1))),
        max_velocity_increment=float(np.max(np.linalg.norm(np.diff(velocity, axis=0), axis=1))),
        max_finite_difference_acceleration=float(np.max(np.linalg.norm(acceleration, axis=1))),
        peak_acceleration_interval=[float(time[peak_acceleration_index]), float(time[peak_acceleration_index+1])],
        velocity_before_peak_acceleration=velocity[peak_acceleration_index].tolist(),
        velocity_after_peak_acceleration=velocity[peak_acceleration_index+1].tolist(),
        max_adjacent_acceleration_change=float(np.max(np.linalg.norm(np.diff(acceleration, axis=0), axis=1))),
        max_radial_velocity=float(np.max(np.abs(radial_velocity))),
        max_rod_length_error=(float(np.max(np.abs(radius-scene.l))) if scene.rigid else None),
        max_chord_radial_identity_error=float(np.max(np.abs(chord_identity_error))),
        max_abs_energy_multiplier=float(np.max(np.abs(values[:, 6]))),
        max_internal_block_displacement=float(summary['max_block_displacement']),
        **settle,
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, default=DEFAULT_SOURCE,
                        help='Existing campaign results directory containing metadata.json and trajectories/')
    parser.add_argument('--output', type=Path, default=Path(__file__).resolve().parent,
                        help='Directory in which to write motion_metrics.json')
    args = parser.parse_args()
    with (args.source/'summary.csv').open(newline='') as stream:
        summaries = list(csv.DictReader(stream))
    metadata = json.loads((args.source/'metadata.json').read_text())
    scene_by_name = {}
    for definition in metadata['scenes']:
        scene = SimpleNamespace(**definition)
        scene.g = np.asarray(scene.g, dtype=float)
        scene_by_name[scene.name] = scene
    selected = [row for row in summaries if float(row['epsilon_ratio']) in (1e4, 1e8)
                and int(row['hz']) in (60, 120) and int(row['sweeps']) == 16
                and row['order'] in ('EL', 'LE', 'JOINT')]
    data = dict(
        source='Existing saved trajectory CSVs only; no new solves or reference comparisons',
        qualification='Structural trajectory screen, not direct video inspection. No phase, frequency, RMS, or curvature criterion is used.',
        energy_normalization='Absolute physical energy change divided by run-initial energy above the shared potential floor',
        tail='Last one second of the accepted trajectory prefix',
        radial_settle='Earliest time after which the remaining saved radius range is below the declared ell-relative tolerance; late end-of-run hits do not establish settling',
        runs=[motion_metrics(scene_by_name[row['case']], row, args.source) for row in selected],
    )
    args.output.mkdir(parents=True, exist_ok=True)
    path = args.output/'motion_metrics.json'
    path.write_text(json.dumps(data, indent=2, allow_nan=False)+'\n')
    print('Wrote', path, 'for', len(selected), 'saved runs')


if __name__ == '__main__':
    main()
