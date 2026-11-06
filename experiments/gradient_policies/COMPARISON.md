# Which tested endpoint-energy methods work best?

There is no universal winner. The strongest method depends on the scene, timestep, and whether energy or trajectory error matters more. These results compare the 45 solver configurations in the gradient-policy study, without changing the original endpoint-energy objective. They do not rank the separate midpoint or work-balance formulations.

## Representative results: 60 Hz, five sweeps, ten seconds

All rows below complete the full interval. Position RMS is measured against a refined continuous RK4 trajectory, in meters. Energy error is the maximum absolute error divided by initial energy above minimum potential.

| Case | Method | Max energy error | Position RMS | Interpretation |
|---|---|---:|---:|---|
| Released spring | Fresh EL, gamma 4, drop r_x | 2.884% | 0.05633 | Lowest position error at these settings; sustained oscillation |
| Radial spring | Fresh EL, gamma 4, drop r_x | 2.884% | 0.07318 | Lowest position error at these settings |
| Released spring | Fresh LE, gamma 0.25, drop r_x | 10.380% | 0.08455 | Useful alternative; LE need not freeze |
| Radial spring | Fresh LE, gamma 0.25, drop r_x | 10.380% | 0.07900 | Similar motion accuracy, larger energy dips |
| Released pendulum | Fresh LE, gamma 1, drop r_x | 0.008211% | 0.02965 | Strong energy/motion compromise |
| Released pendulum | Fresh EL, no cutoff, retain r_x | 1.055% | 0.02002 | Lowest position error, but larger energy error |
| Moving pendulum | Fresh LE, gamma 1, drop r_x | 0.002728% | 0.08730 | Very small energy error, visible phase error |
| Moving pendulum | Predictor-frozen LE, gamma 1, drop r_x | 2.069% | 0.03425 | Lowest position error, but larger energy error |
| Tangential spring | Ordinary XPBD | 85.232% | 0.87609 | Least position error here, but severe energy loss |
| Tangential spring | Fresh EL, gamma 4, drop r_x | 0.00000171% | 1.51523 | Energy accuracy does not imply a correct orbit |

EL means energy then material; LE means material then energy. Both multiplier visits use one frozen linearization and the second visit uses the first increment. Positions receive their sum afterward. `gamma` bounds the raw energy-correction displacement relative to `h|v_n|+h^2|f_n|/m`; a candidate exceeding the bound is skipped. It is not an energy compliance or a damping factor.

## Dependence on timestep and sweep count

At 240 Hz and five sweeps:

| Spring method | Released: energy / position RMS | Radial: energy / position RMS |
|---|---:|---:|
| Fresh EL, gamma 4, drop r_x | 0.1733% / 0.06836 | 0.1895% / 0.07459 |
| Fresh LE, gamma 0.25, drop r_x | 0.6914% / 0.02727 | 0.6914% / 0.03165 |

Thus EL gamma 4 is the strongest spring result at the representative 60 Hz setting, but LE gamma 0.25 is more position-accurate at 240 Hz. Reducing the timestep does not improve every metric monotonically.

Five and twenty sweeps give nearly identical spring results for EL gamma 4. One sweep is insufficient: the released spring reaches 224% maximum energy error at 60 Hz. More sweeps can worsen pendulum phase error: fresh LE gamma 1 has released-pendulum RMS 0.02965 at five sweeps and 0.15139 at twenty. These finite-iteration solvers are not monotonically converging to the continuous trajectory.

For the moving pendulum at 240 Hz, predictor-frozen LE gamma 1 has position RMS 0.003407 and energy error 0.08507%; fresh LE gamma 1 has position RMS 0.019975 and energy error 0.000574%. Reusing a gradient can be useful in this case, despite failing or freezing other cases.

## What the results do not establish

- No tested endpoint-energy method gives convincing tangential-spring dynamics. The baseline loses energy; the energy-corrected variants follow the wrong orbit and do not conserve angular momentum.
- Completion means the run did not trigger a numerical safety limit. Frozen, highly inaccurate, or oscillatory failures can still complete. Of 1,350 runs, 178 stop early; these are marked STOP and excluded from completed-run rankings.
- Small energy error alone can reward a nearly motionless solution. The released spring with predictor-frozen EL is an explicit example.
- Retaining r_x is not uniformly bad, nor is dropping it a general conservation guarantee. Retained-r_x unbounded EL is best for released-pendulum position at the representative setting; other retained variants diverge.
- A threshold chosen for one scene is not established as optimal for another. No single weighted score is used to hide these tradeoffs.

## Figures and data

The trajectory figure shows four labeled candidates for each of the five cases, with energy above minimum, x(t), y(t), and distance to the continuous reference. Its legend includes maximum energy error and position RMS. It is a representative comparison, not a plot of every curve on the same axes.

The six exhaustive scorecards cover every configuration at each pair of timestep and sweep count: 45 methods x 5 cases x 6 settings = 1,350 runs. They include both retained and dropped r_x, both orderings, all tested gradient ages, and every cutoff. Failed prefixes are marked STOP, never ranked as low-error solutions.

- [Trajectory image](../../results/gradient_policies/full_results/best_methods_trajectories.png)
- [Trajectory vector PDF](../../results/gradient_policies/full_results/best_methods_trajectories.pdf)
- [All 1,350 runs: six-page scorecard PDF](../../results/gradient_policies/full_results/all_gradient_policy_scorecards.pdf)
- [60 Hz, five sweeps: every configuration](../../results/gradient_policies/full_results/scorecard_h60_sweeps5.png)
- [60 Hz, one sweep](../../results/gradient_policies/full_results/scorecard_h60_sweeps1.png)
- [60 Hz, twenty sweeps](../../results/gradient_policies/full_results/scorecard_h60_sweeps20.png)
- [240 Hz, one sweep](../../results/gradient_policies/full_results/scorecard_h240_sweeps1.png)
- [240 Hz, five sweeps](../../results/gradient_policies/full_results/scorecard_h240_sweeps5.png)
- [240 Hz, twenty sweeps](../../results/gradient_policies/full_results/scorecard_h240_sweeps20.png)
- [Selected numerical values](../../results/gradient_policies/full_results/selected_metrics.csv)
- [All raw summary data](../../results/gradient_policies/summary.csv)
- [Original methods and test details](README.md)

Generate the representative figure with:

```bash
python3 experiments/gradient_policies/plot_best_methods.py results/gradient_policies
python3 experiments/gradient_policies/plot_scorecards.py --summary results/gradient_policies/summary.csv --output results/gradient_policies/full_results
```

The plot independently reconstructs energy errors and position RMS from all 601 stored samples per selected run. No solver equations or simulation data were changed to make these figures.
