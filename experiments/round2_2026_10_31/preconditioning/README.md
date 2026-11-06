# Numerical preconditioning of the same midpoint-DG equations

This experiment implements controlled transfers from Large Steps and Muon-inspired spectral balancing. It changes the numerical solve, not the physical mass, spring energy, force law, or midpoint kinematics. The previous composition suite is read-only; `common.hpp` is a local copy of its model helpers.

## Equations

For one step from `(x_n,v_n)`, eliminate velocity exactly:

\[
v=2(x-x_n)/h-v_n,
\qquad R(x)=M(x-x_n-hv_n)-\frac{h^2}{2}F_{DG}(x_n,x)=0.
\]

An edge has relative vectors `a` at the old state and `b` at the candidate, lengths `r_a,r_b`, rest length `l`, and stiffness `k`. The force on its second endpoint is

\[
f=-\frac{k}{2}\left(1-\frac{2l}{r_a+r_b}\right)(a+b).
\]

This is average reaction strength times the exact distance discrete gradient. Constant gravity is included exactly. The analytic derivative with respect to `b` is

\[
\frac{\partial f}{\partial b}
=-\frac{k}{2}\left[
\left(1-\frac{2l}{r_a+r_b}\right)I
+\frac{2l}{(r_a+r_b)^2}(a+b)\frac{b^T}{r_b}
\right].
\]

The derivative generally is nonsymmetric, so the inner solver is restarted GMRES rather than conjugate gradients on the physical Jacobian. The exact identity

\[
E(x,v)-E(x_n,v_n)=\frac{2}{h^2}R(x)^T(x-x_n)
\]

is tested even away from a root. This distinguishes energy drift caused by an incomplete solve from changing the work formula.

## Five numerical metrics

Each Newton linear model is `A dx = -R`, with `A=∂R/∂x`. We solve the transformed model `L A Q dy = -L R`, then set `dx=Q dy`. Every transform is invertible. Each method evaluates acceptance and stopping using the same original mass-weighted residual, never its own transformed residual.

| Method | Left transform L | Right transform Q |
|---|---|---|
| none | I | I |
| mass | M^(-1/2) | M^(-1/2) |
| diagonal | D^(-1/2) | D^(-1/2) |
| graph | B^(-1) | B^(-T) |
| svd_rows | T M^(-1/2) | M^(-1/2) |

Here `H = M + (h²/4)(incidenceᵀ diag(k) incidence ⊗ I₂)`, `D=diag(H)`, and `H=BBᵀ`. The anchored graph makes H positive definite. H is an isotropic stiffness surrogate, not the exact DG derivative. B is factored once per trajectory; the graph and diagonal metrics stay fixed.

For `svd_rows`, form `A₀=M^(-1/2) A M^(-1/2)=U Σ Vᵀ` and use

\[
T=U\operatorname{diag}\left((\sigma_i^2+\epsilon^2)^{-1/2}\right)U^T,
\qquad \epsilon=.05\operatorname{median}(\sigma).
\]

T is rebuilt at each nonlinear iteration. This is **bounded row balancing**, not Muon applied to positions or forces. There is no momentum history, raw polar force, or arbitrary reshape. All original singular directions remain present; it cannot remove true rank loss.

[Large Steps](https://bnicolet.com/publications/Nicolet2021Large.html) motivates graph-aware coordinate changes without changing the objective. Our Cholesky factor uses physical mass and a spring stiffness surrogate, rather than copying the rendering paper's exact Laplacian parameterization. [Muon's original explanation](https://kellerjordan.github.io/posts/muon/) motivates inspecting the spectrum separately from its scale. The row-balancing transform above is our proposal, not the Muon algorithm.

## Globalization

Both choices use `φ=0.5 ||M^(-1/2)R||²` and the same inexact Newton direction.

- `line`: Armijo backtracking, at most 14 trials.
- `trust`: cap the correction in the physical metric `||M^(1/2)dx||/[h sqrt(2 E_scale)]`; use the ratio of actual to predicted reduction in φ to accept and resize the radius. It starts at radius 2, shrinks by 1/4 below ratio 1/4, and doubles below the cap when the ratio exceeds 3/4. Acceptance requires positive actual decrease and ratio above 1e-4.

A non-descent linear update invokes the same physical-residual steepest-descent fallback in both variants. The trust variant scales the inexact Newton direction; it does not solve a full trust-region quadratic subproblem, and it is not claimed to implement PPO or a complete filter-SQP algorithm. [TRPO's primary paper](https://proceedings.mlr.press/v37/schulman15.pdf) supplies the locality/metric analogy, not a mechanical stability guarantee.

## Cases and budgets

- Serial anchored chains with N=2,8,32 moving masses in 2D; gravity `(0,-9.81)`.
- Rest lengths 1; nominal stiffness 1,000.
- `uniform`: equal masses/stiffnesses.
- `mass100`: alternating masses with ratio 100, normalized to mean mass 1.
- `stiff100`: alternating stiffness 1,000 and 100,000, equal masses.
- `both100`: both contrasts, with stiff edges adjacent to the light masses.
- `mixed`: a smooth angular perturbation plus alternating stretch and velocity components around the hanging equilibrium. h=0.01 and 0.04.
- `large`: four times those perturbations, h=0.1.
- Duration 0.4 seconds, 36 physical scenarios and 720 algorithm configurations.

Finite budgets are **12 or 48 Krylov matrix products per timestep**, shared across nonlinear iterations, with at most 12 nonlinear iterations and 80 solver trial/residual evaluations plus one final diagnostic residual evaluation. Each linear solve gets at most `max(4,budget/4)` products before relinearization. Actual/predicted model products, analytic Jacobian assemblies, residual evaluations, rejected trials, and factorization time are reported separately. These are common upper budgets, not a claim of identical actual work.

An iteration-limited finite state is propagated so that accumulated distortion remains observable. Every such step is counted. Nonfinite/collapsed states or enormous energy/position end the run and preserve its prefix. The CSV distinguishes a completed trajectory from one whose nonlinear steps all converged. Stalled and budget-limited steps are separate counts.

## Validation and references

`identity_tests.txt` checks the analytic Jacobian by directional finite differences, the energy-residual identity, transformed condition numbers, and agreement of all fully converged metric variants with dense Newton on the same first-step problem. These checks precede the finite-budget benchmark.

For every physical scenario, the trajectory reference solves the **same timestep and DG equations** using dense Newton to a normalized residual of 1e-10. Failed reference prefixes are retained in `references.csv`; no reduced-timestep replacement is silently substituted. Position, velocity, and mode errors compare only available matched prefixes. Thus these errors isolate incomplete-solve distortion; they do not measure DG's temporal discretization error against the continuous ODE.

The shared normalized stopping residual is

\[
\frac{\|M^{-1/2}R\|}{h\sqrt{2E_{scale}}},\quad
E_{scale}=E(x_0,v_0)-E_{minimum}.
\]

Finite-step convergence uses 1e-8. Energy drift uses the initial energy and the same excitation scale.

For modal diagnostics, diagonalize the physical small-oscillation stiffness at the hanging equilibrium in mass-whitened coordinates. Record relative modal phase-space error, energy-weighted phase difference, and the top-quarter modal energy difference divided by total reference modal energy. These fixed linearized modes are diagnostic coordinates; their energies are not exact nonlinear invariants. Position and velocity errors are mass RMS per coordinate.

Timing is serial with Eigen's thread count fixed to one. Process CPU and wall time include numerical metric setup and factorization; reference generation is excluded. This small prototype stores dense Jacobians and uses dense SVD/Cholesky, so its N-scaling and timings do not demonstrate large sparse production performance. Other experiments may run concurrently on the workstation; CPU time is the better local cost indicator, and isolated repeated timing would be needed for performance claims.

## Reproduce

```bash
cd /home/blacknoir/Global/constrained-xpbd-round2-2026-10-31/preconditioning
g++ -std=c++17 -O2 -DNDEBUG \
  -I/home/blacknoir/Projects/constrained-xpbd/external/eigen \
  benchmark.cpp -o benchmark
./benchmark --tests
./benchmark > progress.log
MPLCONFIGDIR=/tmp/xpbd-preconditioning-mpl \
  /home/blacknoir/Projects/mandelbrot-diffusion/.venv/bin/python analyze.py
```

`./benchmark --timing` runs seven serial timing repeats and writes `timing.csv`; run it when competing benchmark jobs are idle.

`runs.csv` is the full result table. `selected_traces.csv` retains every step for the N32 mass-contrast and combined-contrast cases. `RESULTS.md` and `figures/` are generated from the actual CSVs.
