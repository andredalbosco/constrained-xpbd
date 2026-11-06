# Reassessment: plausible motion and actual energy

LE is a useful candidate under the animation criterion. The strongest shared
tested setting is **material then energy, 120 Hz, 16 sweeps, epsilon ratio 10000**.
All six scenes sustain motion and stay below 0.36% maximum mechanical-energy
error over five seconds. This is approximate conservation, not an energy-exact
result or a robust result at cheaper settings.

This reuses the completed [468-run campaign](../block_orders_2026_11_01/README.md).
No solver was rerun or changed. Phase, frequency, reference-position RMS, and
stationary curvature are excluded from the acceptance discussion. Visual review
uses recorded paths, frame sequences, and velocity traces; the playback files
allow inspection of motion in time.

EL visits the coupled energy block before the material block; LE reverses that
order. JOINT denotes the fully coupled solve of the same approximate system.

## Shared-setting comparison

Both orders use 120 Hz, 16 sweeps, and epsilon ratio 10000. Entries are the largest
absolute energy change over the full five seconds, expressed as a percentage of
initial energy above the potential minimum. Gravity potential is included. These
are measurements of actual kinetic plus potential energy, not a relaxed residual.

| Scene | EL maximum energy error | LE maximum energy error | Motion assessment for LE |
|---|---:|---:|---|
| Released spring | 29.128% | 0.355% | Keeps oscillating; brief turning-point pauses and sharp restarts remain |
| Orbiting spring | 0.00429% | 0.00429% | Smooth moving orbit; radial oscillation gradually transfers to rotation |
| Spring with gravity | 0.0423% | 0.0423% | Sustained bouncing and swinging; promising |
| Released pendulum | 0.530% | 0.151% | Sustained swings; promising |
| Kicked pendulum | 0.471% | 0.131% | Sustained swings; promising |
| Stiff orbiting spring | 0.0164% | 0.0164% | Smooth, but radial vibration fades into nearly circular motion |

Maximum LE rod-length errors are 0.0794 mm for the released pendulum and
0.0420 mm for the kicked pendulum, for one-metre rods. No large endpoint jumps or
freezing occur in these six selected LE runs.

The orbit cases are no longer rejected for accumulated phase or frequency
differences. The remaining physical concern is energy redistribution: without
external torque, angular momentum rises from 0.88 to about 1.284 in the ordinary
orbit, and to about 3.278 in the stiff case. In the stiff case, radius variation
falls below 0.1 mm after about 1.93 seconds; speed grows from 0.8 to about
3.245 m/s as radial vibration disappears. Whether this smooth circularization
looks acceptable is a motion-design judgment; low energy error alone does not
settle it.

The released spring has a concrete turning-point hitch. At one inner turning
point its velocity rises from about 0.00026 to 0.41359 m/s in one 1/120-second
step. This repeats at later turns. The resulting finite-difference acceleration
is about twice the spring-force acceleration there. The pause lasts about a
step, so frame playback matters more than the acceleration diagnostic alone.

## What still fails

Increasing epsilon to 1e8 can make actual energy almost constant while producing
non-equilibrium stalls. At 60 Hz and 16 sweeps, JOINT spring release moves only
3.75e-7 m over five seconds, despite the stretched spring. EL and JOINT released
pendulums barely move; kicked pendulums move to an apex and then stall despite
nonzero tangential gravity. These are motion failures under the current criterion.
EL and LE spring-release runs at this setting also hit explicit failure guards.

The usable LE setting is sensitive to budget. With the same epsilon ratio, the
released spring's maximum energy error rises from 0.355% to 174.4% when reducing
120 Hz / 16 sweeps to four sweeps, or to 13.05% at 60 Hz / 16 sweeps. No tested EL
released-spring setting stays within even the 5% reporting band. The focused
campaign therefore supplies a candidate, not a general stable energy-conserving
method.

## Inspect the evidence

- [Main playback: all six scenes, matched EL and LE](visuals/main_motion.mp4)
- [Stronger-energy-setting playback: failures and stalls](visuals/high_penalty_motion.mp4)
- [Paths, spring lengths, and actual energy changes](visuals/main_paths_length_energy.png)
- [Turning-point velocity detail](visuals/release_velocity_detail.png)
- [Visual review and playback details](visuals/VISUAL_REVIEW.md)
- [All-run energy and motion audit](quantitative/AUDIT.md)
- [Detailed motion critique](critical/REPORT.md)
- [Selected-run numerical verification](VERIFICATION.json)

Energy bands of 0.1%, 1%, and 5% in the audit are reporting aids, not tolerances
chosen by the user. A moving run is not automatically classified as plausible.
Five-second evidence does not establish long-term behavior. Multi-pendulums and
constraint networks are outside this focused campaign.

## Reproduce the reassessment

All scripts read the original saved data; none run the integrator. The source
archive is `../block_orders_2026_11_01/results` in the canonical experiments
directory. `SOURCE.sha256` identifies its input CSVs, and `FILES.sha256` identifies
the reassessment artifacts. The source data remain unchanged.

```bash
python3 verify_reassessment.py
OPENBLAS_NUM_THREADS=1 /tmp/xpbd-nonlinear-venv/bin/python quantitative/audit_energy_motion.py
OPENBLAS_NUM_THREADS=1 /tmp/xpbd-nonlinear-venv/bin/python critical/audit_motion.py
```

The selected-run check independently rebuilds kinetic, elastic, and gravity
energy from recorded positions and velocities and verifies their normalization,
returned chord velocities, travel, and rod-length errors. It covers 36 selected
runs and 15,685 states; the largest discrepancy in mechanical energy is below
1e-14 joules. Additional analysis and plotting scripts are in the three
subdirectories.
