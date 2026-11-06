# Validation — 2026-10-23

Baseline XPBD was compiled and tested before active feedback. These are headless numerical checks; no visual/rendering validation was performed.

## Baseline and passive accounting

- The pinned spring matches the analytic implicit-Euler single step and 240-step recurrence, including its known energy decay.
- The hard pendulum swings, preserves its unit length to 1e-12, and keeps its anchor fixed.
- The planar grid moves under gravity, remains planar, and keeps its pins fixed.
- The initially horizontal cloth sags in 3D and keeps its pins fixed.
- An unpinned unequal-mass spring preserves total linear momentum within 1e-10.
- Passive mode has bit-for-bit identical particle positions and velocities to baseline XPBD on all four scenes.
- Independent reconstruction checks the gravity work and every particle and spring work account.

## Ten-second comparison

Reference settings: 2400 timesteps, h=1/240 s, 20 sweeps, active beta=gamma=1000 J. Values below are changes from each scene's initial mechanical energy and include gravitational potential. The reservoir sum has the opposite change for passive/active modes.

| Scene | Initial mechanical energy (J) | XPBD/passive final change (J) | Active final change (J) | Maximum active relative spring error |
|---|---:|---:|---:|---:|
| Spring | 3.125 | -3.076375 | +6.436973 | 0.419872 |
| Pendulum | -7.809582 | -0.771356 | +0.202757 | <1e-12 |
| 2D grid | 30.901500 | -0.030224 | -0.030224 | 0.006832 |
| Cloth | 59.595750 | -20.814073 | -10.004676 | 0.044914 |

All four default active runs stayed finite and within the tested scene deformation bounds. The spring gains too much mechanical energy. Cloth loses less energy in this run, but this does not establish general improvement. The grid is nearly unchanged. No mechanical-conservation assertion is made by the active tests.

Maximum accumulated extended-energy drift was below 3e-13 J at these settings. This verifies the accounting identity; it does not verify mechanical energy conservation. Finite feedback and position stationarity residuals are reported independently.

## Feedback behavior and limitations

- Two symmetric free particles with +0.2 J stored per particle gain about 0.204528 J of total mechanical energy on the first step at beta=1. With -0.2 J per particle, they lose about 0.194472 J. Carrying the remaining reservoir changes the following step.
- A variable-timestep test distinguishes the required previous timestep in lambda/h_previous^2 from using the current timestep. The wrong denominator would change measured work by about 0.443 J in that test.
- A one-sweep cloth solve has a nonzero feedback residual (about 1.65e-6 J) even though its accepted-state work accounts close. The tests explicitly check that these are separate diagnostics.
- Fixed compliant endpoints still update their material and reservoir multipliers without moving.
- Stronger feedback is unstable: beta=gamma=1 caused the cloth to run away in a five-second trial; beta=gamma=10 did so in a ten-second trial. beta=gamma=100 remained bounded for ten seconds but gained about 56 J. The default 1000 is a weak-feedback reference setting, not a stability guarantee.
- No general linear/angular momentum constraints have been added. Pinned scenes include external support reactions, and particle feedback can change momentum.
- No claim is made about arbitrary timesteps, stiffness, iteration counts, long-term bounded reservoirs, collision handling, or convergence to the original optimization objective.

Reproduce the reference comparison with:

```sh
./build/constrained-xpbd --steps 2400 --output results
```

Run the two CTest suites with `ctest --test-dir build --output-on-failure`. The active suite contains actual transfer and motion checks in addition to algebraic bookkeeping checks.

Both suites also passed AddressSanitizer and UndefinedBehaviorSanitizer in a RelWithDebInfo build. LeakSanitizer could not run under the tool environment's ptrace setup, so that rerun used `ASAN_OPTIONS=detect_leaks=0`; leak detection is not claimed.
