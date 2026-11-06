# Measured preconditioning results

720 configurations across 36 physical scenarios. Every run retained a finite 0.4-second trajectory; many steps were deliberately propagated before nonlinear convergence. **Finite completion is not a success criterion.**

## Identity checks

- Analytic directional Jacobians: maximum relative discrepancy 9.25e-09.
- Energy/residual identity away from a root: maximum normalized discrepancy 2.84e-13.
- All 60 converged metric/case combinations agree with dense Newton: maximum combined position/velocity difference 9.15e-14.
- Maximum converged first-step energy defect: 1.68e-12 of initial excitation.

These checks show that the coordinate and row transformations preserve the solved equations. They do not establish convergence for every later step.

## Finite-budget results

Counts use all 36 scenarios. Median errors use only the 28 scenarios with a complete converged-DG reference; the other eight retain explicitly limited reference prefixes. Errors compare the same time discretization, not the continuous ODE.

| Budget | Globalizer | Metric | All steps converged / 36 | Median max modal error, complete reference | Median max energy drift, complete reference |
|---:|---|---|---:|---:|---:|
| 12 | line | none | 10 | 0.000527 | 6.07e-05 |
| 12 | line | mass | 11 | 1.08e-06 | 6.72e-08 |
| 12 | line | diagonal | 10 | 1.62e-06 | 4.72e-08 |
| 12 | line | graph | 13 | 3.75e-07 | 2.89e-08 |
| 12 | line | svd_rows | 17 | 2.71e-08 | 6.5e-09 |
| 12 | trust | none | 10 | 0.00055 | 3.15e-05 |
| 12 | trust | mass | 11 | 1.08e-06 | 6.72e-08 |
| 12 | trust | diagonal | 10 | 1.62e-06 | 4.72e-08 |
| 12 | trust | graph | 13 | 3.75e-07 | 2.89e-08 |
| 12 | trust | svd_rows | 16 | 2.71e-08 | 6.6e-09 |
| 48 | line | none | 19 | 2.43e-08 | 1.95e-09 |
| 48 | line | mass | 19 | 2.53e-08 | 3.37e-09 |
| 48 | line | diagonal | 19 | 2.34e-08 | 2.76e-09 |
| 48 | line | graph | 20 | 3.79e-08 | 1.11e-09 |
| 48 | line | svd_rows | 24 | 1.23e-08 | 1.54e-09 |
| 48 | trust | none | 18 | 2.53e-08 | 3.14e-09 |
| 48 | trust | mass | 18 | 3.39e-08 | 5.32e-09 |
| 48 | trust | diagonal | 18 | 3.4e-08 | 3.33e-09 |
| 48 | trust | graph | 19 | 6.51e-08 | 1.3e-09 |
| 48 | trust | svd_rows | 23 | 1.4e-08 | 1.54e-09 |

![Convergence and solve-work scaling](figures/convergence_scaling.png)

## What improved

Graph and mass scaling improve several contrasting-mass cases under the same Krylov budget. For N=32, mass ratio 100, h=0.04, line search, and budget 12:

| Metric | Maximum modal error | Maximum energy drift | Converged steps / 10 |
|---|---:|---:|---:|
| none | 0.187 | 0.00602 | 0 |
| mass | 0.0367 | 0.00106 | 0 |
| diagonal | 0.0472 | 0.00277 | 0 |
| graph | 0.0233 | 0.000182 | 0 |
| svd_rows | 0.00194 | 7.89e-05 | 6 |

![Finite-budget mode diagnostics](figures/mode_error.png)

Bounded SVD row balancing gives the most fully converged trajectories in this particular sweep. It uses a dense SVD at every nonlinear iteration, so fewer Krylov products must not be presented as a runtime improvement. Graph scaling is also not uniformly beneficial: its isotropic stiffness surrogate can overweight directions whose true tangent stiffness is small.

![Conditioning](figures/conditioning.png)

## Globalization is not a universal repair

Across 360 paired configurations, the trust-radius variant reduced energy drift by more than 10× in 3 pairs and increased it by more than 10× in 13 pairs, ignoring differences below 1e-6. This simple radius-capped Newton variant is not uniformly better than line search.

| Case | Metric / budget | Line-search drift | Trust-radius drift |
|---|---|---:|---:|
| N=8, stiff100, large, h=0.1 | graph / 12 | 0.088 | 0.00316 |
| N=8, stiff100, large, h=0.1 | svd_rows / 48 | 0.234 | 0.013 |
| N=32, uniform, large, h=0.1 | graph / 48 | 0.00217 | 3.96e-05 |

## Failed references and missing convergence

The following dense-Newton references failed to converge within their declared iteration/evaluation limit or stalled in line search. No replacement integrator was substituted. Mode/phase comparisons stop at the listed prefix; energy and residual diagnostics continue for the finite-budget methods.

| N | Regime | Motion | h | Reference prefix | Reason |
|---:|---|---|---:|---:|---|
| 8 | stiff100 | mixed | 0.04 | 6/10 | budget_or_iteration_limit |
| 8 | both100 | mixed | 0.01 | 6/40 | line_search_stall |
| 8 | both100 | mixed | 0.04 | 4/10 | line_search_stall |
| 32 | stiff100 | mixed | 0.04 | 8/10 | budget_or_iteration_limit |
| 32 | stiff100 | large | 0.1 | 2/4 | budget_or_iteration_limit |
| 32 | both100 | mixed | 0.01 | 6/40 | line_search_stall |
| 32 | both100 | mixed | 0.04 | 4/10 | line_search_stall |
| 32 | both100 | large | 0.1 | 3/4 | budget_or_iteration_limit |

![Failed reference prefix](figures/failed_reference_prefix.png)

The combined mass/stiffness contrast at N=32 and h=0.04 remains a failure for every finite-budget variant: none converges at every timestep, and the dense reference stops after four of ten steps. Better early-step mode errors do not validate the later trajectories. The transforms improve linear conditioning but cannot by themselves repair a stalled nonlinear solve or choose its physical branch.

## Isolated setup-inclusive timing

After the other trajectory jobs finished, the timing mode ran seven serial repeats with rotated method order (420 measured trajectories). The table shows median process CPU time for the same N=32, mass100, h=0.04, ten-step scenario; setup is measured separately and is included in total time. These are dense-prototype timings, not a production sparse-solver benchmark.

| Budget | Metric | Total CPU ms | Setup/factorization ms | Converged steps / 10 | Max modal error |
|---:|---|---:|---:|---:|---:|
| 12 | none | 0.261 | 0.001278 | 0 | 0.187 |
| 12 | mass | 0.258 | 0.000873 | 0 | 0.0367 |
| 12 | diagonal | 0.257 | 0.000944 | 0 | 0.0472 |
| 12 | graph | 0.457 | 0.01552 | 0 | 0.0233 |
| 12 | svd_rows | 76.18 | 75.81 | 6 | 0.00194 |
| 48 | none | 0.576 | 0.00082 | 0 | 0.000253 |
| 48 | mass | 0.588 | 0.000815 | 4 | 0.000406 |
| 48 | diagonal | 0.573 | 0.000878 | 0 | 1.97e-05 |
| 48 | graph | 1.23 | 0.0118 | 4 | 9.94e-06 |
| 48 | svd_rows | 83.52 | 83.09 | 10 | 2.73e-09 |

At budget 48 in this case, diagonal scaling has about the same measured cost as identity and roughly 13× lower mode error; graph scaling costs about twice as much and has roughly 25× lower mode error. Dense SVD balancing reaches the tightest residuals but is roughly 145× slower than identity here; almost all its time is setup. Those comparisons apply to the stated scenario and prototype.

![Isolated timing including setup](figures/timing_scaling.png)

Counts of linear products alone therefore overstate the practical advantage of dense spectral balancing. A block approximation or a carefully tested approximate matrix function would be a new experiment; it is not validated by these results.

## Conclusions for the next experiment

1. Keep mass/graph scaling as numerical options: they can improve finite-budget motion without changing the converged DG state.
2. Use bounded SVD balancing diagnostically: isolated timing shows dense setup overwhelms its product-count advantage at N=32.
3. Do not adopt the current trust-radius rule as a default. Its regressions need to be separated into radius scaling, model quality, and inadequate nonlinear budgets.
4. Require convergence-gated acceptance or timestep subdivision for the unresolved contrast cases; these tests deliberately show the damage from accepting unconverged steps.
5. This study measures nonlinear-solve distortion relative to converged DG. Physical phase accuracy versus a refined ODE reference remains a separate question.
