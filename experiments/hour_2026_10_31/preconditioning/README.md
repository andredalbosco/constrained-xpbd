# Cheap preconditioning and nonlinear branch continuation

This experiment keeps the midpoint discrete-gradient spring equations fixed. It tests ways to solve those equations with less linear algebra, and separates that question from whether a converged discrete trajectory has accurate physical phase. See [RESULTS.md](RESULTS.md) for the measurements and figures.

## The unchanged physical equations

Endpoint velocity is eliminated exactly:

\[
v=2(x-x_n)/h-v_n,
\qquad
R(x)=M(x-x_n-hv_n)-\frac{h^2}{2}F_{DG}(x_n,x)=0.
\]

For an edge with old relative vector \(a\), new vector \(b\), rest length \(l\), and stiffness \(k\),

\[
f_{DG}=-\frac{k}{2}\left(1-\frac{2l}{\|a\|+\|b\|}\right)(a+b).
\]

Gravity is constant and exact. These equations satisfy the off-root identity

\[
E(x,v)-E(x_n,v_n)=\frac{2}{h^2}R(x)^T(x-x_n).
\]

Thus zero mechanical residual implies exact energy balance in exact arithmetic. Small energy error alone does not imply a small residual, nor accurate phase. Physical masses, forces, and the discrete-gradient formula are identical for every metric.

The chain Jacobian is stored as nonsymmetric 2×2 blocks on its three block diagonals. Matrix-vector products cost O(N). A block elimination control solves the linear Newton equation in O(N), subject to a checked block pivot threshold. This is an unusually favorable structure for a direct solver; conclusions about this baseline do not automatically transfer to general meshes.

## Invertible linear transformations

Each Newton equation is solved as

\[
LJQ\,\Delta y=-LR,\qquad \Delta x=Q\Delta y.
\]

All final stopping and line-search decisions use the original mass-weighted residual

\[
\rho=\frac{\|M^{-1/2}R\|}{h\sqrt{2E_a}},
\qquad E_a=E(x_0,v_0)-E_{minimum}.
\]

| ID | Construction |
|---|---|
| `none` | Identity left and right factors. |
| `mass` | \(L=Q=M^{-1/2}\). |
| `diagonal` | Diagonal of \(H=M+\tfrac14h^2(B^T\operatorname{diag}(k)B\otimes I_2)\), with \(L=Q=\operatorname{diag}(H)^{-1/2}\). |
| `graph` | \(H=CC^T\), \(L=C^{-1},Q=C^{-T}\); scalar tridiagonal graph factors, shared across coordinates. |
| `vertex_rows` | Bounded inverse-square-root balancing of each 2×2 diagonal block of \(A_0=M^{-1/2}JM^{-1/2}\); applied on the left after mass scaling. |
| `vertex_jacobi` | Guarded inverse of each 2×2 diagonal block of \(A_0\), applied on the left. |
| `vertex_right` | The same guarded block inverse, applied on the right; left scaling remains \(M^{-1/2}\). |
| `edge_rows` | Overlapping pairs of vertices: bounded 4×4 patch row factors, summed with symmetric coverage normalization. |
| `edge_right` | The same positive patch operator, applied on the right. |
| `curvature_graph` | Symmetrize each spring's Jacobian contribution, project that 2×2 contribution to positive semidefinite, assemble it with the original mass matrix, and use a 2×2 block Cholesky factor. Rebuilt per nonlinear iteration by default. |
| `direct_block` | Exact block-tridiagonal Newton linear solve; a control, not an iterative preconditioner. |

For a local square matrix \(A=U\Sigma V^T\), the row factor is

\[
T=U\operatorname{diag}\!\left((\sigma_i^2+\epsilon^2)^{-1/2}\right)U^T,
\qquad \epsilon=0.05\max(1,\sigma_{\lfloor d/2\rfloor}).
\]

It is positive definite with gain at most 20 in the mass-scaled coordinates. The guarded Jacobi inverse instead uses \(V\operatorname{diag}(1/\max(\sigma_i,\epsilon))U^T\). The edge operator is \(W\sum_b R_b^TT_bR_bW\), where \(W\) is inverse square root of coverage count; the patches cover all vertices, so the operator is positive definite. These factors remain invertible even when the physical Jacobian is singular. That preserves the root equation but does not make a singular physical Jacobian nonsingular.

`curvature_graph` modifies only a preconditioning matrix. The nonlinear residual and accepted step still use the full original Jacobian and force. Likewise, caching a preconditioner never freezes the physical Jacobian.

## Connection to the papers

[Muon, Keller Jordan's technical writeup](https://kellerjordan.github.io/posts/muon/) motivates balancing singular directions and avoiding a full dense SVD. Here the transfer is confined to invertible linear-system factors. There is no neural-network momentum, arbitrary matrix reshape, or polar replacement of physical forces. Small exact SVDs are used for the bounded local factors; no uncontrolled Newton–Schulz polynomial is used.

[Large Steps in Inverse Rendering of Geometry, Nicolet et al.](https://bnicolet.com/publications/Nicolet2021Large.html) motivates using a graph-based coordinate transformation to communicate among vertices. Our mass-plus-stiffness matrix and curvature-aware graph factor belong to the mechanical linear solve. They are not extra smoothing energies in the physical objective.

[Eisenstat and Walker, Choosing the Forcing Terms in an Inexact Newton Method, 1996](https://users.wpi.edu/~walker/Papers/forcing_terms,SISC_17,1996,16-32.pdf) motivates avoiding an unnecessarily accurate solve of a stale linear model. The tested `ew2` adaptation uses

\[
\eta_0=0.5,\quad
\widetilde\eta_k=0.9(\rho_k/\rho_{k-1})^2,
\]

then raises \(\widetilde\eta_k\) to at least \(0.9\eta_{k-1}^2\) if that safeguard exceeds 0.1, and clamps the result to \([10^{-4},0.5]\). It requests

\[
\|M^{-1/2}(R+J\Delta x)\|\leq
\eta_k\|M^{-1/2}R\|.
\]

This is a specified adaptation of the paper's safeguarded second choice, with explicit finite-budget floors and caps; its general convergence theorem is not asserted for our bounded runs.

## Budgets, restarts, and caches

The main sweep permits 12 or 48 Krylov products per physical step, at most 12 nonlinear iterations, and at most 80 residual/trial evaluations plus a final diagnostic evaluation. Initial chunks are `max(4,budget/4)`. GMRES uses reorthogonalization and a fixed preconditioner within each chunk. Nonlinear iterations reassemble the Jacobian and normally refresh its metric.

The main iterative stopping threshold is \(\rho<10^{-8}\). Full/direct root audits use \(10^{-10}\), 80 nonlinear iterations, and 900 residual evaluations. Direct controls and continuation have their own larger work allowances; their actual cost and work counts are retained.

`krylov_ablation.csv` varies chunks and compares transformed-residual stopping, physical-residual stopping, and minimizing the physical residual over the generated Krylov subspace using QR. The latter is a diagnostic with additional QR cost. `forcing.csv` compares fixed \(10^{-3}\), loose 0.1, and adaptive forcing. `cache.csv` compares fresh factors, updates every second Newton iteration, one factor per physical step, and one factor for the entire trajectory. These are all fixed factors within each individual GMRES call, so flexible GMRES is not required.

Budget exhaustion deliberately propagates the last finite state in the main finite-budget sweep. Completion of that trajectory is not classified as convergence. Physical-reference and branch studies stop at their first unconverged step and retain the prefix.

## Continuation and physical references

Ballistic, acceleration, and secant guesses are compared. The latter two clip only their extra initial-guess displacement at \(2h\sqrt{2E_a}\) in the mass norm; Newton is free to leave this region.

The continuation solve varies the timestep parameter from zero to the desired h while holding the old state fixed. Every intermediate problem has the same anchor \((x_n,v_n)\); its solved state is not advanced as a physical substep. The predictor differentiates \(R(x,h)=0\):

\[
J\frac{dx}{dh}=Mv_n+hF_{DG}.
\]

The initial parameter increment is 1/8 of h, grows by 1.5 up to 1/4, and halves after failure. A minimum increment 1/4096 and 100 stage limit remain explicit failures. The final problem is exactly the original full-h equation. All failed stages and predictor factor solves count toward work and CPU time. Continuation may choose a different coarse discrete root; it is not a phase-accuracy guarantee.

The main coarse reference is direct block Newton at the same timestep; failed prefixes remain marked. A separate physical study uses paired-resolution RK4 references, accepting them only if position gap is below 1e−6 m, velocity gap divided by the available-energy velocity scale is below 1e−5, and relative energy drift is below 1e−6. Position/velocity errors and linearized modal phase are reported separately. The modal basis is the physical Hessian at hanging equilibrium; those modal energies are diagnostics rather than exact invariants of the nonlinear chain.

The small-amplitude refinement probe preserves a strict-tolerance run and a separately labeled run with an estimated cancellation floor:

\[
\rho_{tol}=\max\left(10^{-10},
64\epsilon_{mach}\frac{2\|M^{1/2}x\|+h\|M^{1/2}v\|}{h\sqrt{2E_a}}\right).
\]

This is a heuristic stopping floor for loss of digits in absolute-coordinate subtraction, not a rigorous error bound. It is applied only in `phase_refinement_roundoff.csv`; the strict failures and the actual achieved residuals are also saved. No main-sweep tolerance is silently relaxed.

## Cases, timing, and reproduction

The main study uses N=8,32,128 serial planar spring chains, nominal stiffness 1000, gravity (0,−9.81), alternating mass/stiffness ratios of 100 and 1000, and mixed low/high-frequency perturbations. Larger N changes gravitational tension and modal frequencies as well as problem size. It is not a replicated independent-system experiment. The horizon is 0.4 s and h is 0.01 or 0.04 s. Branch tests additionally use fourfold larger perturbations at h=0.1 s. The timing sweep also includes N=256.

`timing.csv` contains five rotated serial repeats. The root and both sibling agents held their numerical jobs during this window. Setup/factorization CPU time is recorded separately and included in total solver CPU time. Model construction, eigendecomposition for diagnostic modes, and reference generation are outside solver timing. Other sweep CPU columns are diagnostic; use the dedicated timing files for performance comparisons. This is CPU/Eigen prototype evidence, not a GPU or general sparse-mesh benchmark.

For the direct block control, the `setup_us` column includes its complete block factor-and-solve routine; it does not split that tiny interleaved routine further. Iterative metrics record only preconditioner construction in that column. `timing_followup.csv` contains a second isolated window with seven rotated repeats, measuring the added forcing, cache, and curvature options together.

The standalone CMake project requires only `XPBD_ROOT` (default `$HOME/Projects/constrained-xpbd`) or an explicit `XPBD_EIGEN_INCLUDE_DIR`. It has no runtime dependency on the old composition or round-two checkout. Executable targets use the prefix `xpbd_hour_preconditioning_`, followed by the driver filename. Three CTest checks validate principal roots, the curvature graph, and cached factors; their CSVs go into separate build-tree validation directories.

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release \
  -DXPBD_ROOT=/home/blacknoir/Projects/constrained-xpbd
cmake --build build -j2
ctest --test-dir build --output-on-failure
```

```bash
cd /home/blacknoir/Global/constrained-xpbd-hour-2026-10-31/preconditioning
for source in benchmark physical ablation curvature_graph forcing cache refinement timing_followup; do
  g++ -std=c++17 -O2 -DNDEBUG \
    -I/home/blacknoir/Projects/constrained-xpbd/external/eigen \
    "$source.cpp" -o "$source"
done
./benchmark --tests
./benchmark --sweep
./benchmark --branches
./physical
./ablation
./curvature_graph
./forcing
./cache
./refinement
# Run the next commands only in an otherwise quiet numerical window.
./benchmark --timing
./timing_followup
MPLCONFIGDIR=/tmp/xpbd-hour-mpl \
  /home/blacknoir/Projects/mandelbrot-diffusion/.venv/bin/python analyze.py
```
