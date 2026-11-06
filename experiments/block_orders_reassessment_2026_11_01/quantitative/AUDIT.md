# Actual energy and motion reassessment

This reassessment reads the completed 468-run campaign and its saved
trajectories. It runs no solver. Ranking uses actual mechanical-energy error and
describes motion separately; RMS, phase/frequency agreement, residual convergence,
and stationary curvature are not rejection criteria here.

The reported energy error is the largest absolute change from initial mechanical
energy divided by initial energy above the potential floor, exactly as in the
campaign. A failed candidate remains included in the maximum, and its guard
failure remains visible. The 0.1%, 1%, and 5% bands below are reporting bands,
not a tolerance specified by the user.

## Energy band counts

Each EL or LE cell covers 24 settings. Counts require a complete five-second run
but can include nearly stationary motion.

| Scene | EL: within 0.1% / 1% / 5% | LE: within 0.1% / 1% / 5% |
|---|---:|---:|
| Spring release | 0 / 0 / 0 | 0 / 1 / 1 |
| Spring orbit | 8 / 12 / 12 | 10 / 18 / 18 |
| Gravity spring | 6 / 6 / 10 | 8 / 10 / 18 |
| Pendulum release | 4 / 5 / 10 | 2 / 8 / 14 |
| Pendulum kick | 2 / 3 / 7 | 2 / 9 / 14 |
| Stiff spring | 4 / 6 / 12 | 8 / 8 / 15 |

Across all orders, including JOINT and XPBD controls, 92/130/201 runs meet these
three energy bands. Energy-only counts need the motion qualification below.

## Practical shared settings

These are maximum actual energy errors, in percent, over the same five seconds.

| Scene | LE, 120 Hz, 1 sweep, ratio 1e4 | LE, 120 Hz, 4 sweeps, ratio 1e4 | LE, 60 Hz, 16 sweeps, ratio 1e4 | LE, 120 Hz, 16 sweeps, ratio 1e4 |
|---|---:|---:|---:|---:|
| Spring release | 835.195 | 174.386 | 13.047 | 0.355 |
| Spring orbit | 0.045 | 0.004 | 0.005 | 0.004 |
| Gravity spring | 0.429 | 0.042 | 0.064 | 0.042 |
| Pendulum release | 0.830 | 0.108 | 1.673 | 0.151 |
| Pendulum kick | 0.919 | 0.100 | 0.077 | 0.131 |
| Stiff spring | 3.398 | 0.015 | 0.019 | 0.016 |

LE at 120 Hz, 16 sweeps, ratio 1e4 is the only tested shared setting with all six
scenes complete and below the 1% energy band. Its spring release keeps oscillating:
radius spans 0.7500 to 1.2500 m, maximum speed is 2.500 m/s, and total travel is
7.682 m. Release remains the fragile case when lowering sweeps or timestep rate.

At one sweep, LE is substantially better than EL for the five nonrelease scenes.
At 120 Hz and ratio 1e4, EL's errors on those five are 47.65%, 97.96%, 104.52%,
121.59%, and 69.42%; LE's are shown in the table. This does not make the one-sweep
method successful for spring release.

## What the saved motion says

For the shared LE/120 Hz/16-sweep/1e4 setting:

| Scene | Travel over 5 s | Travel in last 1 s | Maximum endpoint step | Angular extent |
|---|---:|---:|---:|---:|
| Spring release | 7.682 m | 1.575 m | 0.0208 m | 0.000 rad |
| Spring orbit | 6.118 m | 1.275 m | 0.0107 m | 5.784 rad |
| Gravity spring | 6.866 m | 1.312 m | 0.0176 m | 1.251 rad |
| Pendulum release | 6.394 m | 1.272 m | 0.0167 m | 1.300 rad |
| Pendulum kick | 5.349 m | 1.077 m | 0.0144 m | 1.118 rad |
| Stiff spring | 16.084 m | 3.272 m | 0.0272 m | 15.757 rad |

These trajectories have sustained motion and small endpoint displacements. They
have no consecutive velocity reversals for which both speeds exceed 10% of the
initial energy-based speed bound. This is a descriptive diagnostic; it does not
certify plausible mechanics or reject legitimate turning points.

Angular transfer deserves separate scrutiny: maximum angular-momentum change is
0.404 for the ordinary orbit and 2.399 for the stiff spring, both starting at
0.88. The stiff spring progresses through 15.757 radians and has eleven radial
velocity reversals. Its changing balance of radial oscillation and circulation
is a mechanical effect beyond a scalar energy measurement. The tables preserve
it for visual assessment rather than declaring every moving run plausible.

Very small energy errors can accompany freezing. For example, JOINT spring
release at 60 Hz/one sweep/1e8 has only 0.000300% energy error but travels
3.75e-7 m. EL pendulum kick at 120 Hz/16 sweeps/1e8 has 0.001675% error, travels
0.159 m initially, then only 8.50e-8 m during the last second. Such energy scores
cannot be interpreted as successful dynamics.

## Machine-readable tables

- `per_run_energy_and_motion.csv`: all 468 original summary records plus travel,
  extent, speed, endpoint-step, tail-motion, and reversal diagnostics. Original
  RMS remains available but was not used for selection.
- `energy_bands_by_scene_order.csv`: band counts, guard failures, and coarse
  stationary flags for every scene/order, including JOINT and XPBD controls.
- `strongest_by_scene_order.csv`: lowest-energy setting, and lowest-energy setting
  without a coarse stationary/jump flag, for each scene/order. The latter is
  named `motion_present`; it is not an automatic plausibility verdict.
- `cross_scene_settings.csv`: all six scenes aligned by shared solver settings.

The coarse flags are explicit descriptors: less than 1e-4 rest lengths of total
or last-second travel, or an endpoint step larger than half a rest length.
There are 16 near-stationary runs, nine stationary-tail runs, eight runs with a
large endpoint step, and five guard failures. These flags aid inspection; their
cutoffs are not user-specified acceptance criteria. All five guard failures are
retained in the tables.

## Reproduce without solver execution

```bash
python3 audit_energy_motion.py
```

Use `--source /path/to/existing/results` for a relocated copy of the campaign,
and `--output /path/to/report/directory` to write the five outputs elsewhere.
The script imports NumPy and reads saved CSV/JSON files; it never imports or
executes the numerical solver. It checks the report's central campaign-specific
claims before rendering the interpretation.
