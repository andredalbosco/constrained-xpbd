# Local-work results

Generated from the complete CSVs by `python3 analyze.py`. These are research prototypes. The original work-target comparison retains its position objective; a separate final cross experiment also tests the repaired quarter objective.

## Start-only prediction

270 physical state/timestep probes, 4320 per-incidence target comparisons. Maximum normalized 128→256 reference-work refinement: 5.944e-09.

The following is mean absolute per-incidence work error divided by the initial energy above the minimum. It pools the five compliant initial conditions and includes their gravitational work incidences; rigid-rod work is reported separately.

| Scenario | h | frozen_power | acceleration_chord | taylor2 | predicted_midpoint | predicted_dg | implicit_dg | implicit_midpoint_dg |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| nominal | 1/120 | 0.001321 | 0.0008109 | 8.983e-05 | 4.557e-05 | 3.215e-05 | 0.001581 | 1.662e-05 |
| nominal | 1/30 | 0.02123 | 0.01381 | 0.005815 | 0.002922 | 0.002033 | 0.02145 | 0.0011 |
| stiff | 1/120 | 0.2669 | 0.2081 | 0.1106 | 0.057 | 0.04436 | 0.09914 | 0.02412 |
| stiff | 1/30 | 0.7195 | 3.45 | 4.581 | 5.898 | 10.87 | 0.2726 | 0.3974 |

`stiff` means k×100. The start-only midpoint linear solve improves the useful-step-size prediction error and avoids the explicit predictor’s large stiff excursion. The backward-Euler linear predictor is more biased at small h but can be more robust at the largest stiff step. Neither gives an identity at the actual optimized endpoint.

| Rigid-rod target | Max absolute error |
| --- | --- |
| acceleration_chord | 0.74 |
| frozen_power | 3.793e-15 |
| implicit_dg | 4.797e-14 |
| implicit_midpoint_dg | 4.797e-14 |
| predicted_dg | 0.1126 |
| predicted_midpoint | 0.08557 |
| projected_dg | 4.797e-14 |
| taylor2 | 4.334e-15 |

The Taylor cancellation assumes a feasible rod and tangent old velocity. `projected_dg` uses a feasible predicted rod endpoint. The two linear-solve predictors alias this projected construction for the rod.

## Local augmented-Lagrangian solves

All work variables and all vertex multipliers are shared. The methods minimize the same explicit-W objective; only their solve schedules differ. A solved row is not counted as convergence unless the full feasibility and stationarity tolerances pass.

| Method | Evaluation budget | KKT converged | Median feasibility | Median stationarity |
| --- | --- | --- | --- | --- |
| joint_budget | 32 | 24/24 | 2.367e-12 | 1.637e-10 |
| local_AL | 32 | 6/576 | 0.002219 | 0.001726 |
| joint_budget | 128 | 24/24 | 2.367e-12 | 1.637e-10 |
| local_AL | 128 | 128/576 | 0.0001541 | 0.0001697 |
| joint_budget | 512 | 24/24 | 2.367e-12 | 1.637e-10 |
| local_AL | 512 | 265/576 | 9.679e-09 | 9.834e-08 |

| Initial rho | Penalty schedule | KKT converged at512 | Median feasibility |
| --- | --- | --- | --- |
| 0.1 | fixed | 0/72 | 0.006634 |
| 0.1 | progress adaptive | 15/72 | 0.0001681 |
| 1 | fixed | 3/72 | 6.006e-05 |
| 1 | progress adaptive | 13/72 | 1.372e-05 |
| 10 | fixed | 68/72 | 6.673e-09 |
| 10 | progress adaptive | 68/72 | 6.553e-09 |
| 100 | fixed | 49/72 | 3.721e-09 |
| 100 | progress adaptive | 49/72 | 3.721e-09 |

ρ≈10 works well with this particular normalization. These are full-AD evaluation budgets, not a claim of equal FLOPs or distributed speedup. The dense joint reference is much cheaper on the small fixtures.

Different converged branches remain possible:

| Fixture | Objective gap | Local reduced curvature | Joint reduced curvature |
| --- | --- | --- | --- |
| scalar_mixed_implicit_midpoint_dg_soft_midpoint | 0.009369 | 2.171 | 4.661 |
| scalar_mixed_live_dg_soft_midpoint | 0.00775 | 2.156 | 4.637 |
| scalar_mixed_predicted_dg_soft_midpoint | 0.009396 | 2.171 | 4.661 |
| scalar_mixed_taylor2_soft_midpoint | 0.009459 | 2.17 | 4.661 |

Positive reduced curvature confirms that these are different local minima, not merely unfinished residual corrections. The joint reference is a checked local minimum, not a certified global optimum.

## Frozen targets in dynamics

320 one-second runs, statuses: {'completed': 320}. 4 runs contain at least one finite but inexact optimization step.

For nominal stiffness, h=1/60, and soft endpoint kinematics:

| Case | Target | RMS x error | First displacement / reference | Max energy error | Inexact steps |
| --- | --- | --- | --- | --- | --- |
| scalar_release | frozen_power | 0.0253 | 2.088 | 7.677e-08 | 0 |
| scalar_release | predicted_dg | 0.01878 | 2.139 | 8.72e-08 | 0 |
| scalar_release | implicit_midpoint_dg | 0.01877 | 2.139 | 8.78e-08 | 0 |
| scalar_release | live_dg | 0.01882 | 2.198 | 8.688e-08 | 0 |
| radial_orbit | frozen_power | 0.09248 | 1.022 | 2.834e-08 | 0 |
| radial_orbit | predicted_dg | 0.09201 | 1.023 | 2.584e-08 | 0 |
| radial_orbit | implicit_midpoint_dg | 0.0921 | 1.023 | 2.585e-08 | 0 |
| radial_orbit | live_dg | 0.09182 | 1.024 | 2.515e-08 | 0 |
| rod_kick | frozen_power | 0.01464 | 0.9578 | 8.701e-08 | 0 |
| rod_kick | predicted_dg | 0.01457 | 0.9578 | 8.635e-08 | 0 |
| rod_kick | live_dg | 0.01446 | 0.9577 | 8.548e-08 | 0 |
| elastic_mixed | frozen_power | 0.02405 | 0.9937 | 9.057e-08 | 0 |
| elastic_mixed | predicted_dg | 0.02496 | 0.9954 | 8.657e-08 | 0 |
| elastic_mixed | implicit_midpoint_dg | 0.02515 | 0.9953 | 8.628e-08 | 0 |
| elastic_mixed | live_dg | 0.0276 | 0.997 | 8.051e-08 | 0 |
| chain_mixed | frozen_power | 0.04196 | 1.026 | 8.995e-08 | 0 |
| chain_mixed | predicted_dg | 0.0398 | 1.027 | 9.209e-08 | 0 |
| chain_mixed | implicit_midpoint_dg | 0.03999 | 1.027 | 9.18e-08 | 0 |
| chain_mixed | live_dg | 0.04001 | 1.028 | 1.002e-07 | 0 |

More accurate predicted work does not consistently improve motion. For example, the scalar first-step displacement remains biased despite a converged objective, and the elastic/chain soft-midpoint trajectories retain substantial phase error. Work-target quality and mechanical/kinematic consistency are separate questions.

## Local-solver rollouts

88 half-second rollouts at h=1/60. Statuses: {'completed': 46, 'completed_inexact': 42}. Multipliers reset each time step; no trajectory is silently removed for inexact solves.

| Method | Budget | rho | All steps converged | Median max energy error | Median RMS x error |
| --- | --- | --- | --- | --- | --- |
| joint | 128 | 0 | 12/12 | 2.737e-08 | 0.02403 |
| local_AL | 128 | 10 | 0/36 | 4.307e-05 | 0.02669 |
| local_AL | 512 | 10 | 30/36 | 9.854e-08 | 0.0268 |
| local_AL | 512 | 100 | 4/4 | 7.794e-08 | 0.05933 |

Every recorded prefix is available in `results/local_traces/`; summary rows retain unfinished solves and their residuals.

## Corrected objective and start-only work

24 one-second runs at h=1/60, with unchanged explicit-W balances, beta=gamma=1, and soft midpoint velocity penalty. Joint solves get 128 evaluations; local AL gets 512 at fixed rho10. Six derivative checks pass.

| Case | Position objective | Start target | RMS x error | Inexact joint steps | Max balance |
| --- | --- | --- | --- | --- | --- |
| scalar_release | original | frozen_power | 0.03622 | 0 | 9.738e-09 |
| scalar_release | quarter | frozen_power | 0.03011 | 2 | 0.00199 |
| scalar_release | quarter | implicit_midpoint_dg | 0.03044 | 2 | 0.001823 |
| radial_orbit | original | frozen_power | 0.08586 | 0 | 7.921e-09 |
| radial_orbit | quarter | frozen_power | 0.003991 | 0 | 7.755e-09 |
| radial_orbit | quarter | implicit_midpoint_dg | 0.0008308 | 0 | 9.922e-09 |
| chain_mixed | original | frozen_power | 0.2163 | 0 | 5.655e-09 |
| chain_mixed | quarter | frozen_power | 0.008144 | 0 | 9.693e-09 |
| chain_mixed | quarter | implicit_midpoint_dg | 0.001568 | 0 | 9.554e-09 |

The repaired objective makes the more accurate start-only work target useful for radial motion and the chain. Their reported joint trajectories converge at every step; local AL closely matches them, with 4/1 mildly inexact steps for corrected chain frozen-power/predicted-DG. Scalar joint results include two substantial line-search stalls and must not be read as solved-objective accuracy.

The scalar diagnosis replays all four stalled old states without changing their objective or constraints. Increasing the BE-initialized budget from128 to4096 does not help: the same line search stalls after69–112 evaluations. The constraint Jacobian retains minimum singular value0.1887–0.1892, so this is not loss of constraint rank. A start-only linear-midpoint initialization converges in1–3 evaluations; rho10 AL converges in395–453. An independent dense feasible-circle scan agrees with the successful local minimum. The stalled state approaches near-zero reduced curvature while still infeasible; it is not a verified stationary branch. See [CROSS_OBJECTIVE.md](CROSS_OBJECTIVE.md) and `results/scalar_stalls.csv`.

## Limits and failure categories

- Predictor error: explicit Taylor endpoints can overshoot stiff motion; the reference-work error exposes this before optimization.
- Geometric failure: unprojected rod predictions generally assign nonzero rod work; tangent power/Taylor cancellation and feasible predicted DG are distinct fixes.
- Inexact solve: finite trajectory steps with unresolved feasibility/stationarity are counted and retained.
- Local branch: positive-curvature local minima can differ by order and rho even after convergence.
- Physical-model error: converged energy/work balances can coexist with biased displacement, phase, or velocity because the objective/kinematics remain unchanged.
- Scalability is untested: the local algorithm restricts primal variables but uses a dense AD evaluator for this prototype.
