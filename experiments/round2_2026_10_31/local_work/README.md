# Start-only work targets and a local augmented-Lagrangian solver

This experiment tests two gaps in the first composition study:

1. Can work predicted entirely at the start be useful beyond `h f_n·v_n`?
2. Can a genuine local primal/dual solve handle the vertex-owned kinetic multipliers more successfully than repeated row corrections?

The previous suite is read-only. Its `common.hpp` and `independent.hpp` are included through a compiler include path. This directory contains the new prediction formulas, an explicit-work-variable formulation, solver implementations, validation, CSVs, and short trajectories. [RESULTS.md](RESULTS.md) is generated from those CSVs; it reports failures and inexact solves as well as successes.

## Reproduce

The prior suite and Eigen directories are build inputs. When both experiment directories are siblings under `experiments/`, their relative paths from this directory are `../../composition_2026_10_31` and `../../../external/eigen`.

```bash
prior_composition_dir=../../composition_2026_10_31
eigen_include_dir=../../../external/eigen
g++ -std=c++17 -O2 -pthread -I"$prior_composition_dir" -I"$eigen_include_dir" study.cpp -o study
g++ -std=c++17 -O2 -pthread -I"$prior_composition_dir" -I"$eigen_include_dir" al_rollouts.cpp -o al_rollouts
g++ -std=c++17 -O2 -pthread -I"$prior_composition_dir" -I"$eigen_include_dir" cross_objective.cpp -o cross_objective
g++ -std=c++17 -O2 -pthread -I"$prior_composition_dir" -I"$eigen_include_dir" scalar_stalls.cpp -o scalar_stalls
./study all results
./al_rollouts results
./cross_objective results
./scalar_stalls results
python3 analyze.py
```

`study` also accepts `validate`, `predictions`, `al`, and `dynamics` as separate modes. Each takes an optional output directory as its second argument. `al_rollouts` takes the output directory as its first argument. Solves use four worker threads; data are deterministic apart from recorded wall time and CSV row order. This is not an isolated timing benchmark.

## Work available at the start

For incidence `(i,j)`, let `f⁰_ij` be the old force, `a⁰_i` the acceleration under all old forces, and `ḟ⁰_ij` the analytic force derivative along the old physical velocity. Spring derivatives include both changing strength and changing direction. The rigid pendulum derivative includes the reaction's change; old rod states are feasible and tangent in the prediction audit.

| Target | Formula or construction |
|---|---|
| `frozen_power` | `h f⁰_ij·v⁰_i` |
| `acceleration_chord` | `f⁰_ij·(h v⁰_i + h² a⁰_i/2)` |
| `taylor2` | `h f⁰_ij·v⁰_i + h²(f⁰_ij·a⁰_i + ḟ⁰_ij·v⁰_i)/2` |
| `predicted_midpoint` | Evaluate `h f_ij(x̂_mid,v̂_mid)·v̂_mid,i` at the old-acceleration midpoint prediction |
| `predicted_dg` | Construct a second-order predicted endpoint, then evaluate averaged-strength discrete-gradient work between old and predicted positions |
| `projected_dg` | For the rod, predict its angle and map back onto the fixed-length circle before evaluating discrete work; aliases `predicted_dg` for compliant springs |
| `implicit_dg` | Use a start-only backward-Euler linear prediction, then evaluate its discrete work |
| `implicit_midpoint_dg` | Use a start-only midpoint linear prediction, then evaluate its discrete work |

The two stiffness-aware predictions use the positive semidefinite old spring matrix

\[
B_n=\sum_j k_j\nabla C_j(x_n)\nabla C_j(x_n)^T.
\]

Their displacements solve, respectively,

\[
(M+h^2B_n)\,\delta\hat x=hMv_n+h^2f_n,
\]

and

\[
(M+\tfrac14h^2B_n)\,\delta\hat x=hMv_n+\tfrac12h^2f_n.
\]

Both matrices and right-hand sides are known initially. They require a small linear solve, so “available at the start” does not mean “free to evaluate.” The midpoint expression is the old-force linearization of midpoint mechanics; the backward-Euler expression is deliberately included to test a more damped predictor. The matrix drops geometric spring curvature, keeping the predictor matrix positive definite. This does not guarantee accurate nonlinear motion. For the rigid rod these two choices alias the projected-angle prediction, because the compliant spring matrix is not applicable.

The estimated work is **frozen during that time step's optimization**. A discrete gradient evaluated between old and predicted positions is potential-exact for those two positions. It is not an identity for the eventual optimized endpoint.

For an ideal fixed rod with tangent old velocity, `f·v=0` and its time derivative `ḟ·v+f·a=0`. Thus both start power and the complete second-order Taylor work vanish. Keeping only the acceleration contribution destroys this cancellation. A feasible projected predicted endpoint provides a separate geometric zero-work guarantee.

The accuracy audit integrates each incidence's physical power with Simpson quadrature along refined RK4 trajectories. It checks 128 versus 256 subdivisions and independently checks summed spring storage and each particle's kinetic work. The cases include scalar release/mixed motion, planar orbit, rigid kick, elastic pendulum, and two-mass chain; five parameter regimes include fast motion, k×100, mass×100, and combined stiffness/mass changes. Three phases and three time steps give 270 physical probes. Gravity has its own work incidence.

## The same objective, with explicit work

The physical objective is the prior XPBD position objective plus a mass-weighted soft velocity relation. After division by the fixed energy scale `E`, the extra work term is

\[
\frac{\gamma}{2}\sum_{ij}
\left(\frac{W_{ij}-\hat W_{ij}}{E}\right)^2.
\]

For a live estimate, `Ŵ` depends on the current position/velocity and its derivatives remain in this objective. For a frozen estimate it is constant. The experiments keep `β=γ=1`. The soft velocity residual is either

\[
v-(x-x_n)/h,
\qquad\text{or}\qquad
v+v_n-2(x-x_n)/h.
\]

The independent unknowns are positions, velocities, and **actual `W/E` values**. The original suite used corrections around an estimated work. `ExplicitProblem` changes coordinates and removes the estimate from the balance equations while preserving its intended penalty in the objective. AD checks verify this change: kinetic balances have no position derivatives, and storage balances have no velocity derivatives.

The dimensionless balance rows are

\[
g^K_i=\frac{K_i(v_i)-K_i(v_i^n)-\sum_j W_{ij}}{E},
\qquad
g^U_j=\frac{U_j(x)-U_j(x_n)+\sum_iW_{ij}}{E}.
\]

There is one kinetic multiplier `μ_i` per vertex, one storage multiplier `ν_j` per spring, and analogous gravity-storage multipliers. The rigid dynamics experiments additionally enforce position feasibility and tangency. `γ=0` is not presented as a new physical formulation: eliminating free work in a connected graph recovers global total-energy feasibility, subject to the same position/velocity conditions. The reported local-work solves use `γ=1`.

## Actual local primal/dual updates

Write all balance rows as `g(z)` and all their shared multipliers as `y`. The solver minimizes the augmented Lagrangian

\[
\mathcal L_\rho(z,y)=F(z)+y^Tg(z)+\frac\rho2\|g(z)\|^2.
\]

Its primal blocks are:

- **Spring block j:** both endpoint positions and the work values delivered by that spring.
- **Vertex block i:** its velocity and every incident work value, including gravity.

These are overlapping primal blocks. The same position/work/multiplier storage is reused across visits. Updating one endpoint includes its effects on all adjacent balances and objective terms. No vertex multiplier is copied into per-spring versions. Each block uses the full augmented-Lagrangian derivatives restricted to those primal variables, an eigenvalue-shifted positive-definite block Hessian, and an Armijo line search. The shift changes the direction, not the objective.

Three primal sweeps at fixed multipliers are followed by one shared dual update `y ← y + ρg`. Forward, reverse, and symmetric block orders are compared. The symmetric order visits the forward list and then its reverse; its extra visits count against the budget. Fixed `ρ` values are 0.1, 1, 10, 100. The adaptive variant doubles `ρ` every five dual sweeps when feasibility has not improved by at least 30%, capped at 1e5. Unscaled `y` is retained when `ρ` changes. This is a progress schedule, not a claimed ADMM dual-residual formula.

This is a nonlinear block augmented-Lagrangian prototype. It is motivated by alternating primal/dual optimization, but it does not satisfy the usual separable convex ADMM setup. No convex convergence theorem or parallel-consensus speedup is claimed. The conceptual reference is [Boyd et al., *Distributed Optimization and Statistical Learning via ADMM* (2011)](https://web.stanford.edu/~boyd/papers/admm_distr_stats.html).

## What is compared and what remains uncertain

The one-step comparison uses 24 fixtures: scalar release, scalar mixed motion, and the two-mass chain; four targets; and both soft kinematics. Every method starts from the same backward-Euler-based guess. The joint benchmark solves the full exact-curvature KKT equations. Both methods use feasibility `<1e-8` and stationarity `<1e-7` as convergence criteria. Budgets are 32, 128, and 512 full AD evaluations. An additional joint reference has a 4096-evaluation cap, usually stopping much earlier. Positive reduced-Hessian eigenvalues check local minimality at converged branches; no global-minimum claim follows.

The budget includes objective/constraint evaluations and backtracking. It does not equate FLOPs: joint and block factorizations differ. The prototype evaluates dense AD arrays even for local blocks. Reported evaluation counts are meaningful solver-work comparisons on these small problems, not evidence of a scalable implementation. Post-solve curvature diagnostics are outside the optimization budget.

The frozen-target dynamics comparison runs 320 one-second trajectories at 60 and 240 Hz, nominal/k×10 stiffness, and both soft kinematics. Rod tangency is enforced. The local-solver comparison additionally runs 88 half-second trajectories at 60 Hz with joint steps and selected local AL configurations, including the order-sensitive scalar mixed case. Dual variables reset each time step. Finite unconverged steps are accepted deliberately to measure their consequences, counted, and labeled; they are not silently treated as solved constrained objectives.

Validation covers analytic force time derivatives, objective and constraint AD gradients, the explicit-W locality identities, rod Taylor cancellation, reference-work refinement, and 64→128 reference trajectory subdivisions. Useful work prediction, converged optimization, and faithful physical motion remain separate criteria. The tables in [RESULTS.md](RESULTS.md) make this separation explicit.

## Composing the corrected position objective with start-only work

`cross_objective.cpp` is a separate 24-run experiment. It compares the original position objective against the quarter predictor / quarter potential objective from the sibling `objectives` experiment. With `a_n=M⁻¹f(x_n)`, the latter replaces the position term by

\[
\frac{1}{2E}\| (x-x_n)/h-v_n-\tfrac14 h a_n\|_M^2
+\frac{1}{4E}U(x).
\]

The mass-weighted soft midpoint velocity term, explicit actual work values, work-target penalty, and per-vertex/per-potential balance constraints remain identical. Both coefficients are one: `β=γ=1`. This comparison isolates the position objective while also crossing it with `frozen_power` and `implicit_midpoint_dg` targets. The corrected objective does not remove the influence of the work penalty or make these local-work systems equivalent to unconstrained trapezoidal mechanics.

The three compliant cases are scalar release, planar orbit, and the two-mass chain. All runs request one second at 60 Hz. Solvers are joint KKT (128 full AD evaluations per step) and forward local AL with fixed ρ=10 (512 evaluations). Both start from the same old BE-based guess, and multipliers reset per time step. Six finite-difference checks validate the modified objective gradient/Hessian and verify that balance values/Jacobians are unchanged. Summaries, all trajectory prefixes, inexact-step counts, and minimum reduced curvature over converged steps are saved separately in `results/objective_work_cross*.csv` and `results/cross_traces/`. See [CROSS_OBJECTIVE.md](CROSS_OBJECTIVE.md) for the measured result.
