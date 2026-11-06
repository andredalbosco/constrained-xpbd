# Sequential work feedback: observed behavior

Measured 2026-10-25 using the existing C++ solver, with the order `lambda_j -> nu_j -> mu_a -> mu_b` and the previous-force work estimate. No production solver or base scene was changed. The method is not mechanically energy conserving: depending on compliance, timestep, and scene it damps motion, pumps energy into motion, or approaches a distorted stationary state. Splitting the material/potential updates avoids some catastrophic paired-update failures, but does not restore correct dynamics.

## Test conditions

The 64 CLI configurations cover the spring, pendulum, 6×6 planar spring grid, and 9×9 cloth for 10 seconds, with 20 sweeps per step and timesteps 1/240 and 1/60 s. Each scene/timestep compares XPBD, sequential feedback at beta=gamma in {1,10,100,1000,1e6}, and paired-local feedback at {1,1000}. Smaller compliance means stronger feedback; compliance has energy units. Reservoirs start at zero and carry between timesteps.

Additional C++ probes record all 32 spring/pendulum trajectories, three kicked-grid trajectories, three cloth motion histories, and eight isolated two-mass spring trajectories. Initial conditions are unchanged except in the explicitly named `grid_kick`, whose free particles receive `vy=-2 sin(pi column/5) sin(pi row/5)`, adding 1.25 J of initial kinetic energy. Every run completed with finite values; some finite results are physically catastrophic. This is a fixed-duration study, not a long-time stability claim.

## Motion at 240 Hz

This table uses h=1/240 s, 20 sweeps, and **beta=gamma=1000**. Energy change is final minus initial mechanical energy, including gravitational potential. A conservative method would have zero change.

| Scene | XPBD energy change (J) | Sequential energy change (J) | Observed sequential behavior |
|---|---:|---:|---|
| Spring | -3.07638 | +6.43714 | Oscillation amplitude grows; final energy is 3.060 times initial. |
| Pendulum | -0.771356 | +0.202757 | Swing grows gently; final energy above the nominal bottom is 10.14% too high. |
| Resting 2D grid | -0.0302241 | -0.0302241 | Small gravity transient dies out, almost identical to XPBD. |
| Kicked 2D grid | -1.28022 | -1.28022 | The kick is dissipated; essentially no motion remains. |
| Cloth | -20.8141 | -10.0034 | More motion remains than in XPBD, but mechanical energy still falls. |

![Motion and energy histories](../../../results/behavior/figures/motion_energy.png)

For the spring, initial extension is 0.25 m; the last measured absolute peak is 0.0329 m with XPBD and 0.4199 m with sequential feedback. For the pendulum, initial angle is 0.65 rad; the last absolute peak is 0.5160 rad with XPBD and 0.6741 rad with sequential feedback. The pendulum's energy reference is `E+mgL`, initially 2.00041794 J, which avoids ratios based on an arbitrary negative gravitational energy zero.

## Strong feedback can freeze distorted geometry

At **beta=gamma=1**, h=1/240 s:

- The spring still oscillates, but ends at 8.027 times its initial energy.
- The pendulum first gains energy, then develops rod-length error after 8.6375 s and rapidly loses swing energy. Final length is 1.004964 m instead of 1 m; it is still moving at 10 s.
- The kicked 2D grid gains **29.1353 J** and then becomes almost stationary. Final elastic energy is **30.4045 J**, compared with 0.03017 J for XPBD. Mean kinetic energy during 8–10 s is only **2.46e-10 J**; peak relative spring-length error is **12.44%**. Both pins remain fixed and the grid remains exactly planar.
- Cloth gains **202.018 J**, reaches **45.67%** peak spring-length error, and nearly freezes. Mean kinetic energy during 8–10 s is **6.81e-5 J**, versus 5.43 J for XPBD and 8.81 J for sequential compliance 1000. Final stored elastic energy is 193.804 J.

![Kicked grid final states](../../../results/behavior/figures/grid_shapes.png)

These are not successful energy-preserving equilibria. The final full position-stationarity residual is about 610 for the kicked grid and 1429 for strong-feedback cloth, despite very little remaining motion. Accepted-state reservoirs offset the mechanical-energy gains.

![Cloth behavior](../../../results/behavior/figures/cloth_behavior.png)

The matched paired-local cloth case at compliance 1 and h=1/240 is worse: final energy is about **1.03e27 J**, while all recorded values remain finite. Sequential avoids that particular blow-up but leaves a different unphysical state. Stronger enforcement is not a monotone improvement: cloth compliance 10 gains 1303 J at the same timestep.

## Timestep sensitivity and convergence

At h=1/60 s with the same compliance 1000, the spring instead dissipates almost everything: its final energy is 9.90e-7 times initial, versus 3.060 times initial at h=1/240. Pendulum excitation falls to 38.74% of initial, rather than growing to 110.14%. Cloth loses 27.3504 J, rather than 10.0034 J. The same compliance therefore does not deliver a reliable conservation behavior across these stepsizes.

Strong coarse-step spring feedback can stop at an incorrect length: compliance 1 ends at length 0.4156745 m, with late RMS speed 1.36e-8 m/s. Compliance 10 stops on the other side of the anchor at length 2.7285639 m. The material and feedback rows are near roundoff there, but the omitted position-stationarity residual is 26.49 and 59.73 N respectively. Increasing confidence based solely on the scalar rows would miss these failures. The strong coarse pendulum instead also leaves its hard material row unsatisfied and ends at length 1.0140526 m.

## Isolated momentum check

The separate rotating unequal-mass pair has no anchors or external forces. At h=1/240, XPBD preserves linear momentum to 1.82e-12 kg m/s. Sequential feedback with compliance 1000 changes it by as much as **2.238 kg m/s**, and angular momentum changes by 1.214 kg m²/s about the fixed origin. Material and feedback rows are near roundoff. Each kinetic update moves one particle, and no equation forces the two corrections to sum to zero. The order alone therefore does not conserve momentum.

See [free-system analysis](../../../results/behavior/free_system/report.md) for all eight cases, center-of-mass angular momentum, and an independent work-bias decomposition.

## Interpretation

There are distinct observed failure modes: continued energy injection into oscillations, excessive damping, constraint violations, and false stationary states containing elastic strain. The sign and magnitude depend on the scene and parameters. Positive compliance softens feedback but does not make its work estimate exact or restore the dropped position equation.

For an affine spring coordinate with an exact old elastic force,

    d = Delta U + f_old dot Delta x = (Delta C)^2/(2 alpha) >= 0.

Thus its signed potential reservoir `S_new=S_old-d` can only decrease; once negative, it cannot repay toward zero through that work rule. Nonlinear distance constraints add geometric and force-history terms, so the sign statement is conditional, not universal. The observed behavior is consistent with this bias and the omitted stationarity residual; these experiments do not isolate one as the sole cause. A correct accounting identity alone is insufficient to validate the dynamics.

## Evidence and reproduction

- [Spring/pendulum report](../../../results/behavior/spring_pendulum/report.md), [64-case components: spring/pendulum CSV](../../../results/behavior/spring_pendulum/summary.csv).
- [Grid/cloth report](../../../results/behavior/grid_cloth/report.md), [grid/cloth plus probes CSV](../../../results/behavior/grid_cloth/combined_summary.csv).
- [Free-pair report](../../../results/behavior/free_system/report.md), [eight-case CSV](../../../results/behavior/free_system/summary.csv).
- [All saved runs and scripts](../../../results/behavior). Each case records commands, settings, diagnostics and state outputs; probes additionally record trajectories or snapshots.

The baseline spring agrees with the analytic backward-Euler trajectory to 3.91e-14 m and 4.30e-13 m/s. Probe and CLI spring/pendulum mechanical energies agree at every step within 1e-12 J. Cloth probes reproduce the CLI final energies within 1e-10 J. The grid/cloth probes check pins and finite state each step, and grid planarity exactly. These checks validate instrumentation and implementation behavior; they do not assert physical conservation.

From the repository root, reproduce a base run with:

```sh
./build/constrained-xpbd --scene grid --mode active --order sequential \
  --dt 0.004166666666666667 --steps 2400 --iterations 20 \
  --beta 1000 --gamma 1000 --output results/grid-sequential
```

Compile the kicked-grid probe and regenerate its strong-feedback case:

```sh
c++ -std=c++17 -O2 -Wall -Wextra -Wpedantic -I. -isystem external/eigen \
  experiments/sequential_work/behavior/grid_kick.cpp -o /tmp/xpbd-grid-kick
/tmp/xpbd-grid-kick sequential 1 results/behavior/grid_cloth/grid_kick_h240_sequential_b1
```

The other probes are `trajectory_probe.cpp` (takes an output directory) and `free_pair.cpp` (writes CSVs in the current directory). Complete sweep and analysis scripts are retained with the result groups. To regenerate the figures with Python and matplotlib installed:

```sh
python3 experiments/sequential_work/behavior/plot_behavior.py results/behavior
```

The solver header SHA-256 for all measurements is `642230240b7d50e56f8c3ea9aabd2da521bd22f1609aeb444e0bfd6224a361c8`. Raw scientific results are under the repository's ignored `results/` directory; probe source and this report are under `experiments/`.
