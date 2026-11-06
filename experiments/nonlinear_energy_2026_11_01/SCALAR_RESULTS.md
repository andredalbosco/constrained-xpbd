# Scalar result: removing the linear energy spike does not recover oscillation

The exact nonlinear E visit completed all 48 scalar E→L runs without runaway or an infeasible energy surface. Energy immediately after every exact E visit agreed with its target to at most 1.37e-14 of initial energy. This removes the specific near-zero-gradient division that caused raw EL's release spike.

It does not fix the physical motion. At 480 Hz over five seconds with the fixed initial energy target:

| Case | E/L sweeps | Maximum final-step energy error | Total travelled distance | Result |
|---|---:|---:|---:|---|
| Release | 1 | 2.104% | 0.499837 m | Reaches the opposite turning point, then stops |
| Release | 4 | 7.08e-12% | 8.85e-15 m | Frozen at the initial position |
| Release | 16 | 0 | 0 | Frozen at the initial position |
| Mixed velocity | 1 | 1.179% | 0.012317 m | Reaches the nearby turning point, then stops |
| Mixed velocity | 4 | 9.95e-11% | 0.012488 m | Reaches the nearby turning point, then stops |
| Mixed velocity | 16 | 6.06e-13% | 0.012488 m | Reaches the nearby turning point, then stops |

The analytical release oscillator travels approximately 7.99 m in those five seconds. Endpoint position RMS errors for these exact-E variants are 0.295–0.322 m, comparable to or greater than their oscillation amplitudes. Energy closure is therefore a poor criterion for selecting this iteration.

This branch problem survives an exact optimization control. Solving the original constrained implicit objective exactly by comparing its two feasible roots freezes release and brings the mixed case to a stop at its turning point. At zero velocity, the objective is constant on the energy-feasible set. Our deterministic nearest-predictor tie selects the zero-motion root. Different tie policies can choose a different branch; they require an additional dynamical criterion that is absent from this objective.

With only one E→L sweep, the following L visit reopens the energy residual. Resetting the energy target to each previous timestep does not repair that error: it accumulates, losing 94.35% of initial energy for release and 80.00% for the mixed case at 480 Hz. Fixed-initial-target and previous-step-target results must remain separate.

## Validation and historical control

`scalar_exact.py` passed 91 analytic checks covering completed-square roots, exact objective minimization over roots, material constitutive closure, the released-spring two-step spike, and infeasible roots. It ran 112 cases: 48 raw linear E/L, 48 exact E/L, and 16 exact constrained-objective controls. Raw linear E/L completed 27/48 runs; all exact methods completed. Failed runs retain accepted-prefix metrics and their failure status.

`scalar_historical_check.cpp` includes the original, unmodified `traditional_cases/solver.hpp` and runs the first two released-spring steps at 60,120,240,480 Hz. All eight states match the independently coded affine scalar formula; maximum position difference is 4.43e-10 m and velocity difference is 2.13e-7 m/s. The slight discrepancy is expected amplification of double-precision representation differences by the near-zero energy gradient. At 480 Hz the original implementation jumps from 3.125 J initially to 16191.4096107 J on step two. The standalone scalar experiment reproduces this event.

After a radial trajectory crosses its anchor, the original distance spring and this affine scalar spring have different potentials. Consequently the later raw scalar rollouts are not claimed to reproduce historical radial rollouts. The first two steps remain on the common positive-position branch and provide the intended independent control.

## Files and reproduction

- `SCALAR_MATH.md`: objective, exact roots, branch policy, E and L visits, multiplier interpretation, degeneracy.
- `scalar_exact.py`: standalone experiment, requiring only the Python standard library.
- `scalar_results/summary.csv`: every run and all status/cost/accuracy diagnostics.
- `scalar_results/states.csv`: every accepted timestep with energy after the last E and L visits.
- `scalar_results/visits_first8.csv`: every E/L pair in the first eight timesteps, including intermediate positions and the accumulated material multiplier.
- `scalar_results/historical_check.csv`: eight original-solver comparisons.

Run `python3 scalar_exact.py` from this directory. The historical control compiles with C++17 and include directories `constrained-xpbd/experiments`, `experiments/gradient_policies`, and `experiments/energy_ordering`. Its source names the verified original header by absolute path. No original repository source was modified.
