# Measured results: cheap preconditioning and continuation

Local 2×2 factors make the SVD-inspired direction balancing inexpensive. They improve convergence at the larger budget, but tight linear solves can waste a short budget before the nonlinear model is refreshed. Adaptive inexact-Newton forcing helps that tradeoff. For these serial chains, an exact O(N) block Newton solve remains an especially strong control.

Continuation solves more hard discrete steps; it does not correct their under-resolved physical phase. These are separate findings. Full equations, budgets, and verified primary-paper links are in [README.md](README.md).

## Validation and scope

The main sweep contains 798 metric/direct trajectories over 42 scenarios; the curvature graph adds 84. All 798 main trajectories remained finite through 0.4 s, but only 274 converged at every step. Thirteen of the 42 same-timestep direct reference trajectories stop early; 247 main comparison rows therefore have only matched-prefix trajectory errors. Missing continuation/reference intervals are not counted as accurate motion.

All 54 principal metric/root audits pass. Maximum directional Jacobian error is 2.9e-08; block versus dense solve difference 3.84e-15; off-root energy identity defect 2.5e-14; converged-root difference 1.69e-15. Root difference here is mass-weighted position plus h times mass-weighted velocity difference. Additional audits cover 69 Krylov-policy, forcing, cached-factor, and curvature-graph cases.

## Convergence at fixed Krylov budgets

| Metric | Budget 12 | Budget 48 |
|---|---|---|
| none | 4/42 | 15/42 |
| mass | 5/42 | 17/42 |
| diagonal | 4/42 | 15/42 |
| graph | 8/42 | 19/42 |
| vertex_jacobi | 9/42 | 20/42 |
| vertex_rows | 9/42 | 20/42 |
| vertex_right | 9/42 | 21/42 |
| edge_rows | 13/42 | 20/42 |
| edge_right | 14/42 | 23/42 |
| curvature_graph | 15/42 | 23/42 |

The direct block control fully converges in 29/42 cases. Its allowance is larger and its actual cost is measured; it is not credited with an artificial Krylov budget.

![Convergence](figures/convergence.png)

## Isolated setup and solve cost

Five serial repeats with rotated method order, all sibling numerical jobs paused/finished. Values below are N128, mass contrast 100, h = .04, ten steps; iterative budget 48. CPU includes setup.

| Method | CPU ms | Setup ms | Converged steps |
|---|---|---|---|
| diagonal | 0.849 | 0.001 | 0/10 |
| graph | 1.59 | 0.002 | 6/10 |
| vertex_right | 0.904 | 0.356 | 10/10 |
| edge_right | 5.34 | 4.77 | 10/10 |
| direct_block | 0.27 | 0.106 | 10/10 |

For iterative methods, setup is preconditioner construction. For the direct control, that column includes the complete interleaved block factor-and-solve routine. Total CPU is comparable; the setup columns have these different scopes.

At N256, vertex-right Jacobi takes 1.67 ms and diagonal scaling 1.66 ms. Their convergence differs, so time alone would miss the useful result. The direct block control takes 0.493 ms. General meshes need a separate study.

![CPU costs](figures/isolated_cost.png)

## Why short budgets can get worse

For N32 mass100 at h=.04, the first four-product vertex-row linear solve has physical relative residual 0.00872, compared with 0.0325 for diagonal scaling. Nevertheless, the propagated budget 12 trajectory is worse. Moving from four-product to two-product chunks reduces its maximum relative energy defect from 0.173 to 1.66e−4, while allowing more nonlinear model refreshes. It remains unconverged. Physical-residual stopping alone does not change this example. A physical QR minimum improves some errors but adds cost and is not uniformly better.

| Linear forcing / residual policy | Budget 12 fully converged | Budget 48 fully converged |
|---|---|---|
| fixed_transformed | 59/252 | 121/252 |
| fixed_physical | 59/252 | 121/252 |
| loose_physical | 63/252 | 127/252 |
| ew2_physical | 72/252 | 137/252 |

The EW2 example uses the same 120 total products over ten steps, but increases Newton refreshes from 30 to 50 and reduces maximum energy defect to 3.61e−5. Setup work grows. Across 42 cases, adaptive forcing is not uniformly better: the diagonal metric at budget 12 falls from 4 fully converged horizons to 3. This is evidence for distributing work across nonlinear models, not a universal cure or a claim that the resulting inexact states are valid solves.

![Nonlinear budget](figures/nonlinear_budget.png)

## Follow-up isolated cost and cached factors

Seven rotated repeats in a second isolated window. N128 mass100, h=.04, ten steps, budget 12 unless direct.

| Variant | CPU ms | Max energy defect | Converged steps |
|---|---|---|---|
| vertex_rows | 0.462 | 0.0113 | 0/10 |
| vertex_chunk2 | 0.76 | 3.44e-05 | 0/10 |
| vertex_ew2 | 0.698 | 8.3e-06 | 0/10 |
| edge_right | 3.25 | 0.0143 | 1/10 |
| edge_step_cache | 1.27 | 0.025 | 1/10 |
| curvature_graph | 1.25 | 1.08e-09 | 8/10 |
| curvature_step_cache | 0.788 | 0.000275 | 2/10 |
| curvature_ew2 | 1.06 | 1.81e-10 | 10/10 |
| direct_block | 0.259 | 6.31e-13 | 10/10 |

At budget 12, curvature-aware graph scaling plus adaptive forcing converges all ten steps in 1.06 ms in this case, compared with 0.259 ms for the direct control. At budget 48, caching edge factors once per step reduces CPU from 4.99 ms to 1.52 ms while retaining ten converged steps. Caching holds only the metric; every nonlinear residual and Jacobian is still updated. The dedicated cache sweep includes regressions, for example per-step vertex-right caching changes 10/10 to 9/10 convergence on N128 mass100 at budget 48. No unconditional cache recommendation follows.

![Follow-up costs](figures/followup_cost.png)

## Branch continuation and physical phase

| Guess / continuation | Fully converged horizons |
|---|---|
| ballistic | 13/36 |
| acceleration | 3/36 |
| secant | 11/36 |
| continuation | 34/36 |

Continuation agrees with the existing ballistic root on recovered prefixes to about 1e−12 in the nominal contrast 100 examples. Four hard cases find materially different roots on the shared prefix. Two continuation cases still fail: N32 both100, large motion, h=.1 (2/4 steps), and N32 both1000, mixed, h=.01 (23/40). Failed stages are included in cost and work.

All nine physical RK4 references pass their paired-resolution checks: worst position gap 2.27e-07 m and scaled velocity gap 8.88e-06. For N32 both100 at h=.04, continuation finishes 10/10 with about 7.94e−14 energy defect, yet physical RMS velocity error is 3.08 m/s and maximum modal phase error 2.80 rad. Comparing the completed continuation horizon with only the ballistic failure prefix would be misleading.

![Continuation and phase](figures/continuation_and_phase.png)

## Timestep and floating-point limits

A separate 0.02 s probe scales the N8 both100 perturbation to 0.005 of its original amplitude. With the labeled cancellation-aware stopping floor, the last three endpoint phase errors are 0.07537, 0.01888, 0.004724 rad under timestep halving: approximately second-order physical phase convergence. The strict 1e−10 residual test stops early at h=1.25e−4 and produces zero accepted steps at smaller h, due to cancellation in absolute-position subtraction. Both versions are retained. The relaxed probe's largest actual residual is 8.69e−10 and maximum relative energy defect 3.86e−9; the estimated tolerance ceiling is looser than those achieved values.

![Phase refinement and stopping](figures/phase_refinement.png)

## What to carry forward

- Use the 2×2 right block inverse as a cheap iterative control; preserve the direct block control on serial chains.
- Couple preconditioning experiments to nonlinear refresh frequency and forcing tolerances. A fixed tiny linear tolerance can spend the budget in the wrong place.
- Cache expensive patch factors only with explicit convergence checks; keep measuring setup.
- Use fixed-old-state continuation to investigate a stalled discrete root, while checking branch changes and physical references separately.
- Treat phase resolution and floating-point stopping as separate problems from energy preservation.
