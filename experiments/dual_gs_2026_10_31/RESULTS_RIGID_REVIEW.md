# Rigid-pendulum review of the multiplier GS results

This review uses the regenerated, final-engine `results/one_step.csv`, `controls.csv`, and `short_rollouts.csv` read at 20:11 UTC on 2026-10-31, plus the completed `long_rollouts.csv` and its recorded trajectories independently checked after the research pause. Earlier CSVs were superseded after the rigid support/scheduling fixes. The two primary rigid cases are release from angle 0.65 with zero velocity and the same angle with angular speed 1.2. The objective, work prediction and hard constraints are identical across solvers.

## What the tests establish

The independent dense oracle converges on both cases at h=0.02,0.08,0.2 to regular stationary points with positive reduced curvature. Across all 33 primary/compliant-chain reference cases, it needed at most six evaluations and no fallback. This certifies a checked local minimum, not a global minimum of every finite-step problem.

For the single rod, **four schedule labels produce identical numerical results**: `material_storage_kinetic`, `edge`, `edge_reverse`, and `edge_storage_first`. Reversing a one-edge traversal does nothing, and the rod's zero-work row commutes with its length row in these visits. There are only three effective orders in this rigid subset. Do not count the duplicate labels as independent evidence for an ordering.

There are five rows, so a budget multiplier of 4 means 20 actual multiplier visits, and 128 means 640 visits. Tangency is explicitly visited last in the three effective orders. Restored variants also incur their recorded extra row evaluations or stationarity passes; equal visits do not mean equal arithmetic or time.

## Conservation and stationarity separate clearly

The following use `material_storage_kinetic`, omega=1, and 640 visits:

| Case / h | Method | Relative energy defect | Stationarity norm | Scaled distance to dense root |
|---|---|---:|---:|---:|
| Kick / 0.08 | Drop stationarity correction | 9.79e-16 | 1.15e-5 | 1.98e-7 |
| Kick / 0.08 | Local correction | 1.31e-15 | about 2e-16 | 8.80e-16 |
| Kick / 0.08 | Frozen full correction | 9.79e-16 | about 1e-16 | 1.08e-15 |
| Release / 0.20 | Drop correction | 4.22e-11 | 2.42e-4 | 6.61e-5 |
| Release / 0.20 | Local correction | 1.46e-10 | about 4.5e-12 | 1.08e-9 |
| Release / 0.20 | Frozen full correction | 3.60e-10 | 9.44e-12 | 3.83e-9 |

The dropped-correction method can satisfy geometry, work balances and energy while settling at a different point from the objective's stationary root. Extra conservation accuracy is therefore not sufficient evidence that it solves the intended objective.

The corrected variants are not uniformly fast. For release at h=0.02, even 640 visits leave the frozen-full method with a maximum constraint defect 2.16e-6 and root gap 7.70e-5; the local-correction method reaches 1.39e-9 and 4.96e-8. At that turning start, one five-row pass gives a root gap 6.40e-7 for local correction versus 3.13e-4 for frozen-full correction. The whole stationarity restoration is a useful consistency control, but its first finite-budget pass need not preserve the accuracy of the supplied predictor.

Changing order materially changes the dropped method's limiting stationarity. At kick/h=0.08/640 visits, its stationarity ranges from 1.15e-5 for material–storage–kinetic to 1.56e-3 for kinetic–storage–material, while energy remains at roughly floating-point precision. In the same case, both correction variants converge to the same dense root across all effective orders. At small budgets their convergence rates still depend on order.

## Near rest: full rank, severe conditioning, and a separate exact degeneracy

The near-rest control starts at angle 0.00065 and zero velocity. Its dense oracle is **full rank 5** with positive reduced curvature about 1.03871. None of the 144 tested combinations of method, labeled schedule, relaxation and budget reaches the full stopping criterion. This is not an infeasible physical constraint set.

At h=0.08, omega=1, 2,560 visits, material–storage–kinetic gives:

| Method | Scaled length defect | Relative energy defect | Stationarity | Scaled root gap |
|---|---:|---:|---:|---:|
| Drop correction | 3.96e-7 | 8.57e-10 | 6.41e-6 | 8.72e-4 |
| Local correction | 3.96e-7 | 5.04e-7 | 8.07e-14 | 8.66e-4 |
| Frozen full correction | 3.95e-7 | 1.71e-9 | 2.22e-12 | 8.66e-4 |

Correcting stationarity does not remove the remaining feasibility mode. Reducing omega to 0.5 does not repair this control: at the same budget, the local-correction root gap is about 0.0754 and its relative energy defect about 0.00609. The frozen-full correction applies its global restoration undamped while omega damps row pushes, so this relaxation is not a uniform damping of the whole update.

Two conditioning mechanisms must be distinguished:

1. At **exact vertical rest**, the five-row Jacobian has rank four. With the scaled gravity position gradient `b=b_n n`, the dependency is `j_K+j_R+j_G-b_n j_P=0`. The oracle validation confirms rank 4/5 at bottom and top rest and rank 3/4 at zero-gravity rest. Every individual scalar GS denominator is nevertheless at least one. A diagonal denominator check cannot detect the full dependency. The nonvertical rest control remains regular.
2. The **near-rest experiment is not exactly rank deficient**. Its excitation normalization makes the normal gravity-storage gradient grow as `h/sqrt(E)`. At the actual midpoint initializer, independently reconstructing `J B J^T` gives condition number about 5.25e11 at near-rest/h=0.08, versus 104 for nominal release. Scaling each row to unit metric norm still leaves about 1.01e7 versus 58.1. This adds a large numerical conditioning effect to the physical approach toward a singular equilibrium. The remaining length defect is in scaled coordinates, not metres; small absolute motion should not be described as a macroscopic blow-up.

The nominal release also becomes more difficult for this scalar iteration as h decreases: its initial Schur condition number is about 1,261 at h=0.02, 104 at h=0.08 and 90 at h=0.20. This is consistent with the slow small-step turning-point convergence, without implying that larger physical time steps are more accurate.

The control called `rod_release_stiff10` is an alias of the nominal rigid release: spring stiffness is not a physical parameter of this rigid model. Its identical output is a useful check, not an additional stiff-rod experiment.

## Short trajectories: optimization error can cancel time-discretization error

All 148 rigid short trajectories finish their requested horizon. The driver accepts finite inexact steps, so **completion is not full per-step convergence**. The horizons are 1.00 seconds for h=0.02 and 1.04 seconds for h=0.08. Paired physical-reference differences stay below 6.78e-12 in the reported scaled state norm.

For release at h=0.08, material–storage–kinetic:

| Method | Visits per step | Max relative energy drift | Max stationarity | Position RMS in metres |
|---|---:|---:|---:|---:|
| Dense same-objective oracle | Dense control | 6.22e-14 | 2.35e-14 | 0.004412 |
| Local correction | 20 | 0.005179 | 2.07e-5 | 0.003414 |
| Local correction | 80 | 0.0002340 | 6.55e-7 | 0.004360 |
| Frozen full correction | 20 | 0.01278 | 2.25e-5 | 0.001809 |
| Frozen full correction | 80 | 0.0007428 | 1.17e-6 | 0.003739 |

Here more accurate optimization approaches the dense trajectory, but does not monotonically reduce error against the physical trajectory. The low-budget error can partially cancel the objective's finite-step error. Its lower position RMS must not be reported as stronger convergence to the intended constrained objective.

The dense objective itself has nonzero physical trajectory error despite near-exact energy: for kick it has position RMS 0.0003981 at h=0.02 and 0.006448 at h=0.08. Local correction with 80 visits per step gives 0.006439 at h=0.08, nearly the same trajectory accuracy, but energy drift 1.39e-4 and maximum stationarity 9.27e-7. Thus conservation, objective convergence and physical accuracy remain three separate measurements.

The short-rollout summary CSV measures position/velocity error, not an isolated phase error. The crossing-time diagnostic below separates one timing observable; it does not attribute all physical error to phase. The differing horizons also prevent treating the two short-trajectory RMS values as an exact convergence-order calculation.

## Longer trajectories and measured crossing delay

The completed long-run CSV contains 228 trajectories, including 12 dense same-objective controls. Its rigid subset contains **38 completed kick trajectories**, including two dense controls. It tests the same single rod for 5.00 seconds at h=0.02 and 5.04 seconds at h=0.08, using 40 or 160 visits per step. This is a modest horizon extension, not a long-time stability certification. The physical reference refinement gap is at most 2.91e-11 in the recorded scaled state norm. No near-rest or multi-rod rollout is part of this subset.

At h=0.08, omega=1 and 160 visits per step, MSK denotes material–storage–kinetic and SKM denotes storage–kinetic–material:

| Method / order | Max relative energy drift | Max stationarity | Position RMS in metres |
|---|---:|---:|---:|
| Dense same-objective oracle | 5.40e-11 | 4.18e-10 | 0.0288030 |
| Drop correction / MSK | 8.84e-6 | 2.04e-5 | 0.0288047 |
| Local correction / MSK | 8.55e-6 | 2.82e-8 | 0.0288065 |
| Frozen full correction / MSK | 6.54e-4 | 1.09e-7 | 0.0280454 |
| Drop correction / SKM | 7.35e-4 | 1.16e-3 | 0.0286329 |
| Local correction / SKM | 7.90e-4 | 3.02e-8 | 0.0279691 |
| Frozen full correction / SKM | 8.19e-4 | 1.44e-7 | 0.0279618 |

Local correction repairs the stationarity plateau without necessarily producing a visibly better physical trajectory at this step size. For MSK, the dropped and corrected trajectories have almost identical RMS, while their maximum stationarity differs by about 720 times. For SKM, the dropped stationarity is about 1.16e-3 at both budgets; local correction lowers it from 6.67e-6 at 40 visits to 3.02e-8 at 160 visits. The remaining feasibility errors and their accumulated effect are distinct from this stationarity repair. Reversed edge traversal remains exactly equivalent to MSK for the single rod.

At h=0.02, local-correction/MSK/160-visits has position RMS 0.00179032 versus the dense control's 0.00178756; their energy drifts are 7.36e-7 and 3.50e-10. At h=0.08 the corresponding RMS values are 0.0288065 and 0.0288030. The roughly sixteenfold RMS difference between these step sizes is consistent with second-order behavior, but the unequal horizons and only two step sizes do not constitute a certified order study. More visits at fixed h cannot remove the dense objective's physical discretization error.

A separate timing diagnostic can be extracted without rerunning any solver. In `long_trajectories.csv` and `long_references.csv`, detect consecutive samples with opposite signs of the horizontal coordinate `x0`, linearly interpolate the zero within that interval, and pair events by crossing direction. There are five alternating bottom crossings per recorded rigid trajectory. The fifth crossing delays relative to the reference sampled on the same time grid are:

| Method, MSK / 160 visits where applicable | h=0.02 delay | h=0.08 delay |
|---|---:|---:|
| Dense same-objective oracle | 1.78 ms | 28.4 ms |
| Drop correction | 1.78 ms | 28.4 ms |
| Local correction | 1.78 ms | 28.4 ms |
| Frozen full correction | 1.80 ms | 27.9 ms |

For the dense h=0.08 control, the first-to-fifth crossing delay grows from 3.8 ms to 28.4 ms despite nearly exact energy. This directly identifies accumulated timing error in one observable. These are **linear-interpolation estimates from saved samples**, not separately refined event times: for example, the reference's fifth interpolated crossing differs by about 0.116 ms between the two sample grids. The timing estimates should not be reported with more precision than their sampling supports, and the smaller delay of an inexact method is again compatible with error cancellation. Neither the RMS table nor this crossing test establishes global minimum selection or general long-time stability.

## Recommended interpretation and next checks

For this rigid subset, keep refreshed local correction with material–storage–kinetic as a candidate: it preserves the useful small-budget predictor behavior while repairing the dropped method's stationarity plateau. Retain the dropped method as an explicit inexact baseline and the frozen-full method as the Schur consistency control. The current evidence does not establish any of them as a near-rest solution.

Future numerical controls, without changing the objective, are exact elimination of the single `W_r=0` variable/row with its penalty retained as a constant, and a coupled geometry/tangency block. Neither was measured in this result set. They must retain the gravity/kinetic balances and full stationarity checks; a geometry block alone does not remove the cross-row equilibrium dependency. More direct Schur/rank diagnostics should accompany near-rest tests. A plain uniform relaxation factor is not supported as a cure by these controls.

## Snapshot provenance

SHA-256 of the analyzed final-engine CSVs:

- `one_step.csv`: `21ec6b92e5aaa1eb96df584f967a708ab3f69ef6f37c906aa8ac6a8ea3cc2933`
- `controls.csv`: `3dc72513fa3967f781a5db8dde64798cf5f89ace6283f03516dc27114e0f7f08`
- `short_rollouts.csv`: `04e43d4a00c8fe897a05aa082ea4d4b4124e7e06e07805c7d540db645d27df2a`
- `long_rollouts.csv`: `bc15faff07b33e48a339a6a41213175363b906317910e69cea28ce8f177d3544`

- `long_trajectories.csv`: `f9d3d79bd8f6c9d9e8dd8e6ecaf2c4eceda5fbb2ef98f995b091552b10bd2de1`

- `long_references.csv`: `60275c8114a940ed4a317e8efb0f2da487e7eb3d0fc5d1fc3db9e3a11bece847`

All source equations and the exact degeneracy argument are in `RIGID_MATH_REVIEW.md`. The independent oracle validation is in `results/oracle_validation/`.
