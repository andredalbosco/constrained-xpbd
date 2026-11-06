"""Audit saved independent-velocity metrics without executing the solver.

Only saved CSV/JSON inputs and the independent dense matrix assembly from
validate.py are used. No production energy or increment function is called.
"""
from __future__ import annotations

import argparse
from collections import Counter
import csv
import json
from pathlib import Path

import numpy as np

from validate import independent


ROOT = Path(__file__).resolve().parent


class Audit:
    def __init__(self):
        self.checks = 0
        self.maxima = {}

    def close(self, name, actual, expected, tolerance=3e-8):
        self.checks += np.size(expected)
        actual, expected = np.asarray(actual, float), np.asarray(expected, float)
        nan = np.isnan(actual) & np.isnan(expected)
        error = np.where(nan, 0.0, np.abs(actual-expected)/np.maximum(1.0, np.abs(expected)))
        maximum = float(np.max(error)) if error.size else 0.0
        self.maxima[name] = max(self.maxima.get(name, 0.0), maximum)
        assert np.isfinite(maximum) and maximum <= tolerance, (name, maximum, tolerance)

    def true(self, name, condition):
        self.checks += 1
        assert condition, name


def load_csv(path):
    with path.open(newline='') as stream:
        return list(csv.DictReader(stream))


def physical_values(scene, x, v):
    """Mechanical quantities assembled here, not through production code."""
    radius = np.linalg.norm(x, axis=1)
    elastic = np.zeros(len(x)) if scene['rigid'] else .5*scene['k']*(radius-scene['l'])**2
    potential = elastic-scene['m']*(x@np.array(scene['g']))
    kinetic = .5*scene['m']*np.sum(v*v, axis=1)
    return radius, potential, kinetic, potential+kinetic


def audit_results(source):
    metadata = json.loads((source/'metadata.json').read_text())
    traces = json.loads((source/'first_step_traces.json').read_text())
    summaries = load_csv(source/'summary.csv')
    scenes = {scene['name']: scene for scene in metadata['scenes']}
    audit = Audit()
    counts = Counter()
    trace_visits = rows_checked = 0
    min_abs_velocity_factor = float('inf')
    negative_factor_rows = near_zero_factor_rows = 0

    for summary in summaries:
        run_id = summary['run_id']
        scene = scenes[summary['case']]
        m, k, length = scene['m'], scene['k'], scene['l']
        gravity = np.array(scene['g'])
        rigid = scene['rigid']
        h, beta, epsilon = 1/int(summary['hz']), float(summary['beta_v']), float(summary['epsilon'])
        rows = load_csv(source/'trajectories'/(run_id+'.csv'))
        rows_checked += len(rows)
        columns = {key: np.array([float(r[key]) for r in rows]) for key in rows[0]}
        x = np.column_stack([columns['x'], columns['y']])
        v = np.column_stack([columns['vx'], columns['vy']])
        chord = np.vstack([np.array(scene['v0']), np.diff(x, axis=0)/h])
        radius, potential, kinetic, physical = physical_values(scene, x, v)
        chord_kinetic = .5*m*np.sum(chord*chord, axis=1)
        chord_energy = chord_kinetic+potential
        mismatch = np.linalg.norm(v-chord, axis=1)
        speed, chord_speed = np.linalg.norm(v, axis=1), np.linalg.norm(chord, axis=1)
        radial_speed = np.sum(x*v, axis=1)/radius
        radial_kinetic = .5*m*radial_speed**2
        tangential_kinetic = np.maximum(0.0, kinetic-radial_kinetic)
        angular = m*(x[:, 0]*v[:, 1]-x[:, 1]*v[:, 0])
        force = m*np.linalg.norm(gravity)
        floor = -force*length-(0.0 if rigid else force*force/(2*k))
        scale = max(physical[0]-floor, 1e-12)
        geometry = abs(radius-length) if rigid else np.zeros(len(rows))
        tangency = abs(np.sum(x*v, axis=1))/length if rigid else np.zeros(len(rows))
        chord_tangency = abs(np.sum(x*chord, axis=1))/length if rigid else np.zeros(len(rows))
        displacement = np.r_[0.0, np.linalg.norm(np.diff(x, axis=0), axis=1)]
        velocity_jump = np.r_[0.0, np.linalg.norm(np.diff(v, axis=0), axis=1)]
        chord_velocity_jump = np.r_[0.0, np.linalg.norm(np.diff(chord, axis=0), axis=1)]
        accepted = columns['accepted'].astype(bool)
        rejected_indices = np.flatnonzero(~accepted)
        audit.true('at most one rejected last candidate', len(rejected_indices) <= 1 and
                   (not len(rejected_indices) or rejected_indices[0] == len(rows)-1))
        audit.true('initial state accepted', accepted[0])
        naccepted = int(accepted.sum())-1
        counts['runs'] += 1
        counts['completed'] += int(summary['completed'])
        counts['failed'] += bool(len(rejected_indices))
        counts['completed_with_unconverged_steps'] += int(summary['completed']) and np.any(
            columns['converged'][1:][accepted[1:]] == 0)
        computed = dict(energy=physical, rel_energy=(physical-physical[0])/scale,
                        potential=potential, kinetic=kinetic, chord_energy=chord_energy,
                        chord_rel_energy=(chord_energy-physical[0])/scale,
                        chord_kinetic=chord_kinetic, velocity_chord_mismatch=mismatch,
                        displacement_velocity_mismatch=h*mismatch,
                        pi_x=m*(v[:, 0]-chord[:, 0])/beta,
                        pi_y=m*(v[:, 1]-chord[:, 1])/beta,
                        velocity_chord_dot=np.sum(v*chord, axis=1),
                        speed=speed, chord_speed=chord_speed, radius=radius,
                        radial_kinetic=radial_kinetic, tangential_kinetic=tangential_kinetic,
                        angular_momentum=angular, geometry=geometry, tangency=tangency,
                        chord_tangency=chord_tangency, displacement=displacement,
                        velocity_jump=velocity_jump, chord_velocity_jump=chord_velocity_jump,
                        chord_vx=chord[:, 0], chord_vy=chord[:, 1])
        for field, expected in computed.items():
            audit.close('trajectory '+field, columns[field], expected)
        audit.close('fixed initial energy', float(summary['initial_energy']), physical[0])
        audit.close('potential energy floor', float(summary['energy_floor']), floor)
        audit.close('energy normalization above initial potential floor', float(summary['energy_scale']), scale)
        audit.close('fixed dimensional epsilon', epsilon, float(summary['epsilon_ratio'])*scale)

        # Reconstruct every actual stationarity and constraint residual from
        # consecutive saved states, including the velocity coupling terms.
        normal = x[1:]/radius[1:, None]
        gradient = (-m*gravity+np.zeros_like(normal) if rigid else
                    k*(radius[1:]-length)[:, None]*normal-m*gravity)
        difference = v[1:]-chord[1:]
        lam, rho, eta = columns['lambda'][1:], columns['rho'][1:], columns['eta'][1:]
        rx = (m*(x[1:]-x[:-1]-h*v[:-1])/h**2-lam[:, None]*normal-m*gravity
              -m*difference/(beta*h)+eta[:, None]*gradient)
        rv = m*difference/beta+eta[:, None]*m*v[1:]
        alpha = 0.0 if rigid else 1/k
        rlam = radius[1:]-length+alpha*lam
        rrho = epsilon*rho+eta*(physical[:-1]-floor)
        reta = physical[1:]-(1-rho)*physical[:-1]-rho*floor
        res_scale = np.maximum(np.maximum(abs(physical[:-1]-floor),
                                         m*np.linalg.norm(gravity)*length if rigid else k*length**2), 1e-12)
        residual = np.maximum.reduce([
            np.linalg.norm(rx, axis=1)*h*h/(m*length),
            np.linalg.norm(rv, axis=1)*h/(m*length),
            abs(rlam)/length, abs(rrho)/np.maximum(epsilon, res_scale), abs(reta)/res_scale])
        for field, values in dict(r_x_norm=np.linalg.norm(rx, axis=1),
                                  r_v_norm=np.linalg.norm(rv, axis=1), r_lambda=rlam,
                                  r_rho=rrho, r_eta=reta, residual=residual).items():
            audit.close('reconstructed physical '+field, columns[field][1:], values, 2e-7)

        # A compliant velocity mismatch is allowed. Auditing its recorded value
        # never treats nonzero mismatch or chord-energy error as a solver failure.
        factors = 1+beta*columns['eta'][accepted]
        min_abs_velocity_factor = min(min_abs_velocity_factor, float(np.min(abs(factors))))
        negative_factor_rows += int(np.count_nonzero(factors < 0))
        near_zero_factor_rows += int(np.count_nonzero(abs(factors) < .05))
        tail = accepted & (columns['time'] > max(0.0, naccepted*h-1.0)+1e-10)
        tail[0] = False
        last_accepted = np.flatnonzero(accepted)[-1]
        longest_stall = current_stall = 0.0
        for value in chord_speed[1:][accepted[1:]]:
            current_stall = current_stall+h if value < 1e-4*length else 0.0
            longest_stall = max(longest_stall, current_stall)
        path = float(displacement[accepted].sum())
        velocity_travel = float(h*speed[1:][accepted[1:]].sum())
        expected_summary = dict(
            accepted_time=naccepted*h,
            completed=int(naccepted == round(float(summary['duration'])/h) and not len(rejected_indices)),
            max_rel_energy=float(np.max(abs(computed['rel_energy']))),
            max_accepted_rel_energy=float(np.max(abs(computed['rel_energy'][accepted]))),
            final_rel_energy=computed['rel_energy'][last_accepted],
            max_accepted_chord_rel_energy=float(np.max(abs(computed['chord_rel_energy'][accepted]))),
            final_chord_rel_energy=computed['chord_rel_energy'][last_accepted],
            failed_trial_rel_energy=computed['rel_energy'][-1] if len(rejected_indices) else np.nan,
            path_length=path, velocity_travel=velocity_travel,
            position_velocity_travel_ratio=path/max(velocity_travel, 1e-30),
            tail_path_length=float(displacement[tail].sum()),
            tail_velocity_travel=float(h*speed[tail].sum()),
            tail_stalled=int(bool(tail.sum()) and displacement[tail].sum() < 1e-4*length),
            longest_stall_duration=longest_stall,
            max_accepted_displacement=float(displacement[accepted].max()),
            max_accepted_velocity_jump=float(velocity_jump[accepted].max()),
            max_accepted_chord_velocity_jump=float(chord_velocity_jump[accepted].max()),
            max_velocity_chord_mismatch=float(mismatch[accepted].max()),
            rms_velocity_chord_mismatch=float(np.sqrt(np.mean(mismatch[accepted]**2))),
            max_displacement_velocity_mismatch=float(h*mismatch[accepted].max()),
            reversed_velocity_fraction=float(np.mean(computed['velocity_chord_dot'][1:][accepted[1:]] < -1e-10)) if naccepted else np.nan,
            max_geometry_error=float(geometry[accepted].max()),
            max_tangency_error=float(tangency[accepted].max()),
            max_chord_tangency_error=float(chord_tangency[accepted].max()),
            max_residual=float(columns['residual'][accepted].max()),
            converged_fraction=float(columns['converged'][1:][accepted[1:]].sum()/max(naccepted, 1)),
            max_abs_rho=float(abs(columns['rho'][accepted]).max()),
            max_abs_eta=float(abs(columns['eta'][accepted]).max()),
            initial_angular_momentum=angular[0], final_angular_momentum=angular[last_accepted],
            max_angular_error=float(abs(angular[accepted]-angular[0]).max()),
            angular_conservation_expected=int(not rigid and np.linalg.norm(gravity) == 0),
            min_radius=float(radius[accepted].min()), max_radius=float(radius[accepted].max()),
            final_radial_kinetic=radial_kinetic[last_accepted],
            final_tangential_kinetic=tangential_kinetic[last_accepted])
        for field, expected in expected_summary.items():
            audit.close('summary '+field, float(summary[field]), expected)
        audit.true('explicit failure reason matches rejection', bool(summary['reason']) == bool(len(rejected_indices)))
        if not len(rejected_indices):
            expected_visits = naccepted*int(summary['sweeps'])*(1 if summary['order'] == 'JOINT' else 2)
            audit.close('completed fixed budget visits', int(summary['visits']), expected_visits)
        elif summary['reason'] == 'runaway_guard':
            audit.true('runaway candidate violates an explicit guard',
                       radius[-1] > 1e4*length or speed[-1] > 1e6*length or chord_speed[-1] > 1e6*length)
        audit.true('accepted candidate within explicit guards',
                   np.all(radius[accepted] <= 1e4*length) and
                   np.all(speed[accepted] <= 1e6*length) and
                   np.all(chord_speed[accepted] <= 1e6*length))

        # Dense independent replay of serialized first-step states.
        if run_id in traces:
            point = np.r_[np.array(scene['x0'])+h*np.array(scene['v0']), scene['v0'], 0.0, 0.0, 0.0]
            blocks = 'J' if summary['order'] == 'JOINT' else summary['order']
            args = (np.array(scene['x0']), np.array(scene['v0']), h, m, k, length,
                    gravity, epsilon, physical[0], floor, beta, rigid)
            for visit, snapshot in enumerate(traces[run_id]):
                block = blocks[visit % len(blocks)]
                audit.true('serialized block order', snapshot['block'] == block)
                audit.true('serialized sweep index', snapshot['sweep'] == visit//len(blocks))
                reference = independent(point, *args)
                ids = ([0, 1, 2, 3, 4] if block == 'L' else
                       [0, 1, 2, 3, 5, 6] if block == 'E' else list(range(7)))
                point[ids] += np.linalg.solve(reference['matrix'][np.ix_(ids, ids)], -reference['residual'][ids])
                saved = np.r_[snapshot['x'], snapshot['v'], snapshot['lam'], snapshot['rho'], snapshot['eta']]
                audit.close('serialized first-step dense state replay', saved, point, 2e-6)
                updated = independent(point, *args)
                audit.close('serialized first-step analytical residuals',
                            np.r_[snapshot['r_x'], snapshot['r_v'], snapshot['r_lambda'],
                                  snapshot['r_rho'], snapshot['r_eta']], updated['residual'], 5e-6)
                trace_visits += 1

    if '--trace-first' in metadata['command']:
        audit.true('all requested first-step traces saved', len(traces) == len(summaries))
    return dict(status='passed', checks=int(audit.checks), counts={key: int(value) for key, value in counts.items()},
                trajectory_rows=rows_checked, serialized_trace_visits=trace_visits,
                serialized_trace_runs=len(traces), trace_first_requested='--trace-first' in metadata['command'],
                min_abs_one_plus_beta_eta=min_abs_velocity_factor,
                negative_one_plus_beta_eta_rows=negative_factor_rows,
                abs_one_plus_beta_eta_below_0p05_rows=near_zero_factor_rows,
                largest_scaled_errors=audit.maxima,
                source=str(source.resolve()), numerical_solver_executed=False)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--results', type=Path, default=ROOT/'results')
    parser.add_argument('--write-report', action='store_true')
    args = parser.parse_args()
    result = audit_results(args.results)
    print(json.dumps(result, indent=2, sort_keys=True))
    if args.write_report:
        (ROOT/'METRICS_AUDIT.md').write_text(
            '''# Independent metrics audit

All saved trajectory mechanical quantities and summary statistics are recomputed
from CSV positions and independent velocities. The audit never calls the
production energy or increment functions, and it does not execute a rollout.

## What was checked

- Actual energy uses one half m times returned velocity squared plus the radial
  spring and gravitational potentials. Chord energy instead uses position
  displacement divided by h. Both are checked against their explicitly labeled
  columns; nonzero disagreement is allowed by the compliant formulation.
- Coupling mismatch, its h-scaled displacement measure, implicit coupling
  multipliers, radial/tangential kinetic energy, angular momentum, geometry,
  velocity/chord tangency, endpoint motion, and changes of velocity.
- Every recorded position and velocity stationarity residual and material,
  relaxation, and energy residual, reconstructed from the incoming accepted state.
  Energy targets use incoming actual energy and remain fixed during each solve;
  epsilon uses the initial energy scale and remains fixed throughout a run.
- Summary energy extrema, final energies, accepted-prefix path and tail travel,
  longest stalls, endpoint jumps, velocity disagreement, geometry/tangency,
  angular drift, and convergence fractions. Completed runs use the requested fixed
  visit budget. Completion does not imply residual convergence or plausible motion.
- The explicit guard conditions are checked on accepted states. Rejection and
  failed-candidate summary semantics are audited in the source and reconstructed
  when failed rows exist. This campaign contains no rejected rows, so it supplies
  no runtime failing-candidate example.

## Population conventions

`max_rel_energy` includes a retained rejected candidate when present.
`max_accepted_rel_energy`, most other maxima, and motion statistics cover accepted
states. `final_rel_energy` and `final_chord_rel_energy` describe the last accepted
state; `failed_trial_rel_energy` separately describes a rejected candidate.
This differs from the earlier position-only runner's final-energy convention.

The tail excludes the row exactly one second before the final accepted time.
Stall speed and tail-motion thresholds are descriptive diagnostics, not physical
acceptance requirements. A nonzero compliant velocity mismatch is not itself a
failure, and chord energy is not substituted for actual mechanical energy.

Serialized first-step traces are audited by independent dense block replay when
present. The primary campaign did not request `--trace-first`, so its saved trace
dictionary is empty and this receipt reports zero trace visits. The separate
solver validation covers 630 dense-replayed refreshed visits; this saved-data
audit does not claim those as serialized campaign traces.

The receipt also counts accepted rows with small or negative `1+beta_v*eta`.
These are diagnostics of the velocity stationarity relation, not an additional
acceptance rule.

Run with the experiment's NumPy environment:

```bash
OPENBLAS_NUM_THREADS=1 python3 metrics_audit.py --write-report
```

Use `--results /path/to/saved/results` for another saved campaign. The summary
block-displacement and block-velocity maxima, rank-defect visit counts, and CPU
timings cannot be reconstructed from endpoint CSVs alone; they remain recorded
solver diagnostics rather than independently verified endpoint statistics.

## Numerical receipt

```json
'''+json.dumps(result, indent=2, sort_keys=True)+'\n```\n')


if __name__ == '__main__':
    main()
