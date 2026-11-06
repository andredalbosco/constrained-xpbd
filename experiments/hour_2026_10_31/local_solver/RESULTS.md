# Numerical local-work results

Generated from retained CSVs. All options solve the same final quarter objective with explicit W, beta=gamma=1, and a frozen start-only midpoint-predicted DG target. Finite unconverged steps are accepted deliberately, counted, and never classified as solved objectives.

20/20 five-environment analytic-vs-AD derivative checks pass.
2/2 larger-chain checks pass (N8 against AD, N32 directional finite differences).
63 reference refinement checks: maximum position difference 2.488e-10, velocity difference 4.783e-08, normalized energy difference 1.497e-10.

## Initialization and globalization

81 one-second trajectories cover the five environments, with nominal and k×10 compliant cases. The rod has one nominal case. Each row below therefore covers nine trajectories.

| Method | All steps converged | Inexact steps | Worst balance | Mean evaluations/step |
| --- | --- | --- | --- | --- |
| AL_BE_512 | 7/9 | 4 | 2.859e-06 | 141.4 |
| AL_midpoint_512 | 8/9 | 2 | 1.563e-07 | 86.29 |
| AL_midpoint_warm_512 | 8/9 | 2 | 1.245e-08 | 84.02 |
| continuation_BE | 5/9 | 6 | 0.06013 | 13.29 |
| kkt_BE | 4/9 | 9 | 0.2319 | 6.243 |
| kkt_midpoint | 9/9 | 0 | 9.116e-09 | 2.524 |
| kkt_midpoint_warm | 9/9 | 0 | 9.116e-09 | 2.524 |
| merit_BE | 5/9 | 5 | 0.06295 | 6.676 |
| residual_LM_BE | 4/9 | 8 | 0.2319 | 6.761 |

The midpoint initializer is the decisive repair in these tests. Continuation, residual damping, and the merit safeguard do not make the BE initializer uniformly reliable. More solver iterations cannot repair a line search that has already stopped. Warm multipliers add little on the small cases because the midpoint predictor already lands close to the relevant branch.

| Case | Regime | RMS x error | Max normalized energy error | Evaluations/step |
| --- | --- | --- | --- | --- |
| chain_mixed | nominal | 0.001568 | 9.54e-10 | 3 |
| chain_mixed | stiff10 | 0.03506 | 6.526e-08 | 2.983 |
| elastic_mixed | nominal | 0.0006995 | 2.174e-08 | 2.917 |
| elastic_mixed | stiff10 | 0.0145 | 1.606e-08 | 2.95 |
| radial_orbit | nominal | 0.0008308 | 1.642e-08 | 2.933 |
| radial_orbit | stiff10 | 0.02839 | 7.455e-09 | 3 |
| rod_kick | nominal | 0.0002762 | 1.146e-08 | 2.933 |
| scalar_release | nominal | 0.0022 | 8.527e-15 | 1 |
| scalar_release | stiff10 | 0.07128 | 4.547e-15 | 1 |

The scalar harmonic predictor is itself the constrained stationary solution, hence one evaluation. It still has the phase error of midpoint mechanics, especially when stiff; zero energy error is not exact dynamics.

## Same-state branch comparison

378 method comparisons use common physical old states at three phases and two time steps. Each objective gap is against the midpoint-initialized joint reference for that exact old state. Warm-only duplicates are omitted.

| Method | KKT converged | Objective gap above1e-6 | Negative reduced curvature | Largest converged objective gap |
| --- | --- | --- | --- | --- |
| AL_BE_512 | 54/54 | 0 | 0 | 3.775e-10 |
| AL_midpoint_512 | 54/54 | 0 | 0 | 3.215e-10 |
| continuation_BE | 53/54 | 1 | 1 | 0.8093 |
| kkt_BE | 53/54 | 0 | 0 | 7.815e-12 |
| kkt_midpoint | 54/54 | 0 | 0 | 0 |
| merit_BE | 54/54 | 0 | 0 | 7.802e-12 |
| residual_LM_BE | 53/54 | 0 | 0 | 7.815e-12 |

A converged KKT point is classified further by reduced curvature. Negative curvature is a stationary saddle, and positive curvature is evidence of a local minimum only. Different old states along separately evolving trajectories are not used for objective-gap claims.

## Shared multipliers on longer chains

72 eight-step trajectories cover N2/8/32, nominal/stiffness×10/alternating-mass×10, fixed evaluation budgets, and warm/cold shared multipliers.

| N | Method | All steps converged | Worst balance | Worst stationarity | Mean evaluations |
| --- | --- | --- | --- | --- | --- |
| 2 | kkt_midpoint_cold | 3/3 | 2.943e-09 | 3.788e-09 | 3 |
| 2 | kkt_midpoint_warm | 3/3 | 2.93e-09 | 3.794e-09 | 3 |
| 2 | AL_cold_128 | 0/3 | 4.506e-06 | 1.615e-05 | 128 |
| 2 | AL_warm_128 | 0/3 | 5.164e-06 | 1.738e-05 | 128 |
| 2 | AL_cold_512 | 1/3 | 5.51e-08 | 1.167e-07 | 286.1 |
| 2 | AL_warm_512 | 1/3 | 5.502e-08 | 1.41e-07 | 282.2 |
| 2 | AL_rhoN_cold | 3/3 | 9.851e-09 | 9.944e-08 | 244.9 |
| 2 | AL_rhoN_warm | 3/3 | 9.917e-09 | 9.915e-08 | 231.8 |
| 8 | kkt_midpoint_cold | 3/3 | 8.794e-09 | 1.272e-08 | 3.792 |
| 8 | kkt_midpoint_warm | 3/3 | 9.588e-09 | 1.306e-08 | 3.75 |
| 8 | AL_cold_128 | 0/3 | 0.0002004 | 0.003589 | 128 |
| 8 | AL_warm_128 | 0/3 | 0.0002971 | 0.003285 | 128 |
| 8 | AL_cold_512 | 0/3 | 1.851e-05 | 0.0001466 | 512 |
| 8 | AL_warm_512 | 0/3 | 2.077e-05 | 0.00016 | 512 |
| 8 | AL_rhoN_cold | 0/3 | 4.097e-06 | 0.004825 | 512 |
| 8 | AL_rhoN_warm | 0/3 | 4.089e-06 | 0.004824 | 512 |
| 32 | kkt_midpoint_cold | 3/3 | 6.515e-09 | 1.995e-08 | 4.583 |
| 32 | kkt_midpoint_warm | 3/3 | 4.373e-09 | 1.933e-08 | 4.583 |
| 32 | AL_cold_128 | 0/3 | 0.04167 | 1.206 | 128 |
| 32 | AL_warm_128 | 0/3 | 0.02795 | 1.019 | 128 |
| 32 | AL_cold_512 | 0/3 | 0.001571 | 0.04968 | 512 |
| 32 | AL_warm_512 | 0/3 | 0.00179 | 0.04718 | 512 |
| 32 | AL_rhoN_cold | 0/3 | 0.00228 | 2.672 | 512 |
| 32 | AL_rhoN_warm | 0/3 | 0.002283 | 2.674 | 512 |

At a fixed evaluation budget the original split partition performs fewer complete passes as N grows. With N32 and128 evaluations, it can update the dual only after a partial three-pass primal cycle. Warm duals and rho=10N do not remove the unresolved primal/dual modes. Rho tuning is numerical; it does not change the constrained objective.

## Grouping position and velocity in local blocks

108 same-state AL solves compare vertex-full and edge-full blocks, rho1/10/100, and256/1024 evaluations. Objective gaps use a converged joint solve of exactly the same objective and old state.

| N | Block | Budget | KKT converged | Median balance | Median stationarity |
| --- | --- | --- | --- | --- | --- |
| 2 | vertex_full | 256 | 6/9 | 8.431e-09 | 5.215e-08 |
| 2 | vertex_full | 1024 | 6/9 | 8.431e-09 | 1.749e-08 |
| 2 | edge_full | 256 | 6/9 | 9.376e-09 | 2.068e-13 |
| 2 | edge_full | 1024 | 6/9 | 9.376e-09 | 9.667e-13 |
| 8 | vertex_full | 256 | 0/9 | 7.601e-07 | 2.418e-06 |
| 8 | vertex_full | 1024 | 0/9 | 6.045e-08 | 2.359e-07 |
| 8 | edge_full | 256 | 0/9 | 7.573e-07 | 2.522e-07 |
| 8 | edge_full | 1024 | 1/9 | 9.374e-08 | 1.145e-07 |
| 32 | vertex_full | 256 | 0/9 | 0.00014 | 0.001039 |
| 32 | vertex_full | 1024 | 0/9 | 1.05e-05 | 8.498e-05 |
| 32 | edge_full | 256 | 0/9 | 3.78e-05 | 0.0004805 |
| 32 | edge_full | 1024 | 0/9 | 2.835e-06 | 1.663e-05 |

## Stable line-search differences

The exact incremental objective/balance difference changes arithmetic only. It avoids cancellation of nearly equal total gravitational/objective values in the Armijo test.

120/120 incremental-difference checks pass across all five small environments and an N8 mass-contrast chain. The 45 additional one-second trajectories retain all baseline CSVs.

| Method | All steps converged | Inexact steps | Mean evaluations/step |
| --- | --- | --- | --- |
| AL_BE_512_stable | 7/9 | 4 | 137.9 |
| AL_edge_full_stable | 9/9 | 0 | 8.004 |
| AL_global_stable | 9/9 | 0 | 6.752 |
| AL_midpoint_512_stable | 8/9 | 2 | 83.57 |
| AL_midpoint_warm_512_stable | 9/9 | 0 | 80.12 |

In one/two-mass systems, an edge-full block can contain the entire primal system. These small-case successes are not evidence that purely local solves scale.

Longer same-state comparisons converge in 13/28 cases with ordinary subtraction and 26/28 with stable differences. Both remaining stable failures are N32 mass-contrast/rho10 cases; fixed budgets still matter.

| N | Regime | Block | Ordinary evaluations | Stable evaluations | Stable balance | Stable stationarity |
| --- | --- | --- | --- | --- | --- | --- |
| 8 | nominal | edge_full | 1524 | 241 | 3.574e-09 | 9.372e-08 |
| 32 | nominal | edge_full | 16384 (inexact) | 3457 | 4.109e-10 | 7.426e-08 |
| 32 | mass10 | edge_full | 16384 (inexact) | 4801 | 6.351e-10 | 9.091e-08 |
| 32 | nominal | global | 79 | 41 | 7.432e-09 | 1.283e-12 |
| 32 | mass10 | global | 108 | 58 | 9.641e-09 | 1.149e-12 |

The global block is a control that solves the full primal matrix on each visit. It is not a local method or a claim of equal factorization cost. The tiny joint KKT solve remains much cheaper on these fixtures.

## Local iteration spectrum

3/3 finite-difference checks validate the derived local full-step AL iteration matrix, with maximum relative error 1.148e-07. These local asymptotic predictions assume the nearby full-step regime; they do not guarantee nonlinear convergence or predict an arbitrary initial error vector.

| N32 block | rho | Spectral radius | Cycles for six decades | Block visits for six decades |
| --- | --- | --- | --- | --- |
| split | 10 | 0.962 | 356.8 | 6.851e+04 |
| split | 100 | 0.9306 | 191.9 | 3.685e+04 |
| edge_full | 100 | 0.7656 | 51.73 | 4966 |
| edge_full | 1000 | 0.9377 | 214.8 | 2.062e+04 |
| global | 100 | 0.724 | 42.78 | 128.3 |
| global | 1000 | 0.2078 | 8.794 | 26.38 |

The slow modes remain even after the line-search arithmetic fix. Increasing rho can improve dual correction while worsening local primal coupling. See [TURNING_AND_LOCALITY.md](TURNING_AND_LOCALITY.md) for the derivation and the scalar turning-point trigger.

## Warm starts after a converged joint first step

24 retained trajectories compare cold duals, ordinary warm duals, and one joint first step followed by warm local AL. All use stable differences, edge-full blocks, and rho100. The joint first step is counted and labeled as a hybrid initialization.

| N | Budget | Dual initialization | All steps converged | Inexact steps | Worst balance | Mean evaluations |
| --- | --- | --- | --- | --- | --- | --- |
| 8 | 512 | cold | 2/2 | 0 | 5.673e-09 | 242.5 |
| 8 | 512 | warm | 2/2 | 0 | 5.673e-09 | 241 |
| 8 | 512 | joint_bootstrap | 2/2 | 0 | 2.852e-09 | 206.8 |
| 8 | 2048 | cold | 2/2 | 0 | 5.673e-09 | 242.5 |
| 8 | 2048 | warm | 2/2 | 0 | 5.673e-09 | 241 |
| 8 | 2048 | joint_bootstrap | 2/2 | 0 | 2.852e-09 | 206.8 |
| 32 | 512 | cold | 0/2 | 16 | 1.421e-05 | 512 |
| 32 | 512 | warm | 0/2 | 16 | 1.421e-05 | 512 |
| 32 | 512 | joint_bootstrap | 0/2 | 14 | 1.542e-05 | 448.6 |
| 32 | 2048 | cold | 0/2 | 16 | 4.307e-07 | 2048 |
| 32 | 2048 | warm | 0/2 | 16 | 4.307e-07 | 2048 |
| 32 | 2048 | joint_bootstrap | 0/2 | 14 | 7.15e-08 | 1793 |

## Single-rod reaction control

All main hour rod runs include the old physical reaction in the quarter predictor. A separate omitted-reaction control also appears second order for this single-degree-of-freedom pendulum because energy, position, and tangency already restrict its motion strongly. It must not be extrapolated to multi-rod chains.

| Old reaction | Hz | Two-second RMS x error | Max energy error | Inexact steps |
| --- | --- | --- | --- | --- |
| omitted | 30 | 0.00331 | 6.234e-09 | 0 |
| omitted | 60 | 0.0008357 | 3.534e-08 | 0 |
| omitted | 120 | 0.0002094 | 2.02e-07 | 0 |
| omitted | 240 | 5.233e-05 | 1.507e-06 | 0 |
| included | 30 | 0.00211 | 5.323e-09 | 0 |
| included | 60 | 0.0005258 | 3.432e-08 | 0 |
| included | 120 | 0.0001312 | 2.019e-07 | 0 |
| included | 240 | 3.273e-05 | 1.507e-06 | 0 |

The fixed per-step residual tolerance accumulates to energy error about 1.5e-6 at 240 Hz. Both the convergence order and that accumulation are retained, rather than treating a green residual as exact conservation.

## Boundaries

- Work-target accuracy, local energy feasibility, stationary-objective convergence, and trajectory accuracy are distinct measurements.
- Some inexact large-chain updates stay closer to the predictor and therefore to the high-resolution trajectory. That does not make them accurate solves of the fixed objective.
- The chain family changes physical lengths, stiffnesses, and total normalization with N. It does not establish replication invariance or a uniform continuum limit.
- Dense analytic derivative storage and dense KKT solves remain in this prototype. Recorded timings are not isolated scaling benchmarks.
- Positive reduced curvature supports local minimality; it never certifies a global optimum.
