# Independent velocity: motion and energy inspection

This is a Python experiment, not native application integration. The numerical
campaign was started after the separate solver-validation pass. The runner
never replaces the returned independent velocity with position difference.

## Campaign and measured diagnostics

The primary grid contains 648 five-second runs: six scenes, EL/LE/JOINT,
beta_v in 0.1/1/10, 60/120 Hz, 1/4/16 sweeps, and energy-penalty ratios
10000/100000000. All 648 reach the requested five seconds without guard
rejection. Forty-four have nearly stationary final-second position travel;
completion is therefore distinct from satisfying a motion criterion.

The 468 prior derived-velocity summaries are copied as an explicitly reused
baseline. Their trajectories are read from the canonical archive, without
rerunning the old solver. No continuous-reference RMS or phase score is used.

Actual energy is `0.5*m*dot(v,v) + U(x)` with the independently returned v.
Chord energy substitutes `(x_new-x_old)/h` only for a diagnostic. The two
energies need not match in the chosen formulation. Both energy-error columns
normalize by initial mechanical energy above minimum potential.

## Selected view

Figures and playback use one shared setting: beta_v=1, 60 Hz, four sweeps,
energy-penalty ratio 100000000. Every figure includes EL, LE, and JOINT.

| Case | LE max actual energy error | LE max chord-energy diagnostic | Max velocity/chord difference |
|---|---:|---:|---:|
| Released spring | 0.09237% | 2.8447% | 0.2475 m/s |
| Orbiting spring | 0.000000249% | 0.8929% | 0.00885 m/s |
| Spring with gravity | 0.000231% | 0.8132% | 0.03199 m/s |
| Released pendulum | 0.01616% | 0.2728% | 0.04893 m/s |
| Kicked pendulum | 0.01600% | 0.3113% | 0.03717 m/s |
| Stiff spring | 0.006367% | 15.8716% | 0.4341 m/s |

LE continues moving in all six selected examples. In sampled frames, both
pendulums swing and the orbit/gravity springs have ongoing smooth-looking
paths. Spring release is substantially more useful than the frozen JOINT
branch, but still has turning-point pulses in the velocity/chord discrepancy
and coupling factor. These pulses should not be described as fully resolved.

The stiff spring loses radial oscillation while spinning up toward a nearly
circular path. Its angular momentum changes substantially. Small actual
energy error does not validate this transfer between motion modes.

JOINT leaves the released spring effectively frozen and the kicked pendulum
nearly stationary after reaching an apex. The released pendulum stays near
its initial pose for several seconds, then moves to another pose and nearly
stalls. Its coupling factor crosses zero. These branches show why a small
actual-energy residual or joint-solve reference is insufficient as a motion
quality test.

The LE selected coupling factor `1 + beta_v*eta` stays positive. Its minimum
is about 0.0463 for the released spring, 0.9894 for the orbiting spring,
0.7227 with gravity, 0.3762 for released pendulum, 0.4742 for kicked
pendulum, and 0.8016 for the stiff spring. A value near zero accompanies
position motion much slower than stored velocity at the turning correction.

## More sweeps and smaller timesteps

At beta_v=1, LE, 60 Hz, 16 sweeps, ratio 1e8, the worst actual energy
error across the six scenes falls to about 0.000705%. The maximum chord
diagnostic remains about 15.887%, with velocity/chord difference about
0.435 m/s. At 120 Hz and 16 sweeps, those worst values are approximately
0.000606%, 4.082%, and 0.167 m/s. These measurements separate conservation
of stored-velocity energy from the independent kinematic approximation.

## Inspection boundary and outputs

All six contact sheets were visually inspected at 0, 0.25, 0.5, 1, 2, 3,
4, and 5 seconds, together with the one-second poster and full energy,
coupling-factor, path, length, and speed traces. This is sampled-frame and
numerical-trace inspection; it is not a claim of watching every video frame
as a continuous animation.

The final movie passes a complete `ffmpeg -v error -i comparison.mp4 -f null -`
decode with exit status zero. It contains 151 H.264 frames, 1800 by 1500,
30 fps, lasting 5.033 seconds. An independently decoded three-second frame
(`decoded_t3.png`) was also visually inspected.

- `comparison.mp4`: synchronized EL/LE/JOINT motion and actual/chord energy.
- `playback_poster.png`: one-second sample of that view.
- `energy_coupling.png`: actual-v energy, chord-energy diagnostic, mismatch,
  and `1 + beta_v*eta`.
- `paths_lengths_speeds.png`: accepted paths, lengths, and actual/chord speeds.
- `*_contacts.png`: six inspectable contact sheets.
- `SELECTED_RESULTS.md`: selected new methods and matching archived baselines.
- `settings.csv`: 108 shared configurations summarized across all six scenes.
- `selected_runs.csv`: exact run IDs used by this view.

The runner records actual/chord energies, reconstructed coupling multiplier,
kinematic differences, geometry/tangency diagnostics, radial/tangential kinetic
energy, angular momentum, jump/stall measures, and stationarity residuals.
These diagnostics support further investigation; they do not establish
robustness for networks, contacts, or arbitrary initial conditions.
