# Multiplier Gauss–Seidel: derivation, implementation, and measured results

Completed study of the local iteration, begun October 31 and finalized November 1, 2026. The main solver in `constrained-xpbd` is unchanged. This directory contains the mathematical specification, independent C++ experiment code, raw measurements, validation and plots.

**The useful next baseline is frozen-gradient multiplier GS with one stationarity-restoration pass per sweep.** Keep material→storage→kinetic and storage→kinetic→material as explicit ordering controls. Dropping stationarity is a useful comparison, but it can satisfy every material/work row while failing to solve the objective. Retaining stationarity inside every local visit often improves convergence per visit, but costs substantially more and can be unstable. None of these results establishes a robust low-sweep solver for large systems.

## 1. Focus and completed work

The experiment holds the physical objective, work target, normalization and initialization fixed, then changes the numerical update and the order of accumulated multiplier visits. Every visit immediately pushes its supported primal variables. There is one shared kinetic multiplier per moving vertex.

| Stage | What was completed |
|---|---|
| Derive | Full objectives, constraints, gradients/Hessians, stationarity, material lift, approximate Newton system and local scalar updates for all five examples |
| Implement | Three numerical methods, six explicit visit schedules, sparse row evaluation and incident-multiplier gathering |
| Validate | 4,729 equation/update/locality checks; independent dense exact-Hessian oracle and reduced-curvature checks |
| Measure one-step convergence | 4,158 primary configurations; 1,296 stiffness, mass, near-rest and relaxation controls |
| Test propagation | 540 configurations on 3-, 8- and 32-mass chains |
| Compare numerical cost | 432 configurations with matched row-gradient evaluation ceilings |
| Test accumulated motion | 792 local short trajectories plus 22 dense controls; 216 local five-second trajectories plus 12 dense controls |
| Explain instability | Local sweep stability probe, derived separately before implementation |

The short duration is 1–1.04 s and the long duration is 5–5.04 s, because the number of steps is rounded to an integer. A finite inexact step is allowed to advance and is labeled by its residuals; completion does **not** mean convergence or physical accuracy. Failed prefixes remain in the CSVs.

## 2. Objective and the five examples

Positions, velocities and incidence works are independent unknowns. With old physical acceleration `a_n` (including the old rod reaction), fixed excitation scale `E`, and a work target computed entirely from the old state, the objective used here is

\[
\Phi(x,v,W)=
\frac{1}{2E}\left\|\frac{x-x_n}{h}-v_n-\frac h4a_n\right\|_M^2
+\frac{U(x)}{4E}
+\frac{1}{2E}\left\|v+v_n-2\frac{x-x_n}{h}\right\|_M^2
+\frac{1}{2E^2}\|W-\widehat W\|^2.
\]

Thus beta=gamma=1, with a **soft midpoint velocity target**. This is the repaired independent-velocity objective from the preceding experiments. This study fixes that choice; it does not compare all earlier work laws or derive a universal integrator.

At each particle and potential store, respectively,

\[
k_i=\frac{K_i(v_i)-K_i(v_{n,i})-\sum_{j\ni i}W_{ji}}E=0,
\qquad
s_j=\frac{U_j(x)-U_j(x_n)+\sum_{i\in j}W_{ji}}E=0.
\]

Positive `W_ji` delivers energy to the particle. Gravity has separate storage/work leaves. A hard rod has zero stored energy, plus length and endpoint-tangency constraints. Adding a total-energy row would duplicate the sum of the local balance rows and is not done.

| Example | Potential and additional equations | Primary variants |
|---|---|---|
| Affine scalar spring | `U=k(x-l)^2/2`; no absolute-distance substitution | Release and mixed velocity |
| Planar radial spring | `U=k(||x||-l)^2/2` | Radial release and tangential motion |
| Rigid pendulum | Gravity potential; `||x||=l`, `n(x)·v=0`, rod work zero | Release and kick |
| Elastic pendulum | Radial spring plus gravity; no hard length/tangency rows | Mixed motion and vertical release |
| Two-spring pendulum | Two compliant springs, two moving masses and gravity | Mixed motion and larger release |

The N3 extension makes a genuinely shared middle vertex explicit; N8/N32 examine propagation. The chain constructor changes lengths and stiffness with N, so these are also harder physical cases, not a pure complexity benchmark. Full equations and exact scene parameters are in [MATH_EXAMPLES.md](MATH_EXAMPLES.md).

## 3. From the objective to actual multiplier updates

Let `z=(q,u,w)`, with `q=(x-x_n)/x_s`, `u=v/v_s`, `w=W/E`, `v_s=sqrt(E/m_ref)` and `x_s=h*v_s`. Introduce the spring material multiplier **after** formulating the objective. With

\[
\overline C_j=C_j/x_s,\qquad
\widetilde\alpha_j=\frac{4E}{k_jx_s^2},
\]

the saddle term `-lambda_j Cbar_j - alpha_j lambda_j^2/2` eliminates to exactly `U_j/(4E)`. The material row is `Cbar_j+alpha_j lambda_j=0`; the storage row still uses the full spring energy change. Let `f_0` be the objective without quarter spring energy, `y` collect material/kinetic/storage/geometry multipliers, and `D` contain material compliances only. The full equations are

\[
r_z=\nabla f_0-Gy=0,\qquad r_c=b(z)+Dy=0,
\]
\[
\begin{bmatrix}H&-G\\G^T&D\end{bmatrix}
\begin{bmatrix}\delta z\\\delta y\end{bmatrix}
=-\begin{bmatrix}r_z\\r_c\end{bmatrix},\qquad
H=\nabla^2f_0-\sum_a y_a\nabla^2b_a.
\]

Here the columns of `G` are the individual gradients. The numerical approximation replaces `H` by `A=∇²f_0`, explicitly dropping all multiplier curvature. Its inverse is local:

\[
A_i^{-1}=\frac{m_{ref}}{m_i}
\begin{bmatrix}I&2I\\2I&5I\end{bmatrix},\qquad A_w^{-1}=I.
\]

Consequently a kinetic multiplier can move **both x and v**, as well as every incident work variable. This behavior is implemented and tested.

| Method | One visit / sweep |
|---|---|
| `refreshed_drop` | Reevaluate the row; solve `delta y_a=-(b_a+D_a y_a)/(∇b_aᵀ A⁻¹ ∇b_a+D_a)`; immediately push `z += A⁻¹∇b_a delta y_a` |
| `refreshed_local_rx` | Also gather full stationarity restricted to that row's fixed structural support; solve the small restricted KKT equations and update only those variables |
| `frozen_full_rx` | Freeze row gradients for a sweep; apply `-A⁻¹r_z` once in a separate pass; visit multiplier coordinates using the latest accumulated primal displacement and frozen linear row residuals |

The last method is exactly GS on the **frozen approximate Schur system**, expressed by local primal pushes. It uses no dense primal solve. Its separate restoration pass is counted. All methods reset duals to zero each timestep and start from the explicit minimizer of `f_0`.

The six schedules are three separated passes and three edge traversals: material→storage→kinetic, kinetic→storage→material, storage→kinetic→material, edge-interleaved, edge-interleaved with storage before material, and reversed edge traversal with the same within-edge order. Tangency is explicit. Shared kinetic multipliers are revisited, never duplicated. Exact row lists and timing boundaries are in [IMPLEMENTATION.md](IMPLEMENTATION.md); the full derivation is in [MATH_GS_REVIEW.md](MATH_GS_REVIEW.md).

## 4. What the measurements establish

### Row satisfaction does not establish stationarity

For the two-spring mixed case, `h=.08`, material→storage→kinetic, 1,024 visits:

| Method | Largest row residual | Stationarity norm | Relative energy error | Gap to dense root |
|---|---:|---:|---:|---:|
| Drop correction | 2.46e-12 | **2.90e-3** | 5.95e-12 | 3.52e-3 |
| Retain local correction | 2.29e-12 | 5.67e-14 | 5.55e-12 | 7.32e-12 |
| Frozen sweep + restoration | 2.57e-11 | 2.85e-13 | 6.21e-11 | 9.91e-11 |

This is explained by the exact drift identity for the dropped update:

\[
r_z^{new}=r_z^{old}-[G(z^{new})-G(z^{old})]y^{new}.
\]

Once every row is satisfied, the dropped method stops changing the state even if this accumulated stationarity error remains. A scalar mixed case at `h=.2`, storage→kinetic→material, makes the distinction larger: energy error `1.29e-16`, stationarity `0.0966`, and root gap `0.0143`.

### Ordering matters, and equal visits are not equal work

At N3, `h=.08`, 192 visits, local correction gives full residual `7.93e-5` with material→storage→kinetic, versus `6.17e-3` with edge/storage-first and `3.76e-3` with reverse-edge. This is the same formulation and shared multipliers, with different schedules.

However, that local-correction run evaluates **1,426 rows**, versus 192 for the dropped and frozen methods. It gathers neighboring multiplier contributions to stationarity. The matched-evaluation experiment changes the practical ranking:

| Two-spring mixed, h=.08, material→storage→kinetic | Actual visits under 1,000 row evaluations | Full residual |
|---|---:|---:|
| Drop correction | 1,000 | 2.90e-3 |
| Retain local correction | 146 | 2.94e-5 |
| Frozen sweep + restoration | 1,000 | **4.02e-11** |

Frozen restoration has the smallest full residual in **30 of 36 matched case/order groups**, at each of the 1,000 and 3,000 evaluation ceilings. Row evaluations are a useful work proxy, not a full runtime model: restoration arithmetic, writes, sparse gathers and the old-state predictor also cost work. Timings are descriptive, not a hardware performance benchmark. See [RESULTS_COST_REVIEW.md](RESULTS_COST_REVIEW.md).

![Residuals at matched row-gradient evaluation budgets](results/figures/matched_cost.png)

The separate [N3 visit-budget plot](results/figures/chain3_convergence.png) shows why comparing only visit counts favors the more expensive local correction.

### Stability and physical accuracy remain separate requirements

The sweep map was differentiated at the independently solved KKT root. Central differences at two scales were checked against the derived product of individual visit derivatives. For a sweep Jacobian, spectral radius below one indicates local attraction; above one indicates linear instability. The accepted probes agree to roughly `1e-7` in spectral radius.

| Local-correction sweep near the correct root | Material→storage→kinetic | Storage→kinetic→material | Reverse edge |
|---|---:|---:|---:|
| N3, h=.2 | **1.1653** | **1.1653** | **6.1508** |
| N8, h=.08 | 0.9851 | 0.9851 | **3.0216** |

The N3 large-step failures therefore include a locally unstable iteration: improving initialization cannot make that root attractive under those unmodified sweeps. In contrast, the separated N8 sweeps at `h=.08` are attractive near the root even though the standard initialization can fail; basin size and transient growth matter there. N3 at `h=.02/.08` has radii about `.99325/.99126`, explaining slow asymptotic convergence. The dropped method has the expected unit modes along its manifold of row-satisfying fixed points.

Three N8 `h=.2` local probes failed the fixed-point/finite-difference reliability checks; their large numerical spectra are not used as trusted evidence. There are 108 recorded probe rows, representing two difference scales for 54 maps, with 51 maps accepted. See [STABILITY_MATH.md](STABILITY_MATH.md) and [RESULTS_STABILITY.md](RESULTS_STABILITY.md). This analysis covers the two refreshed methods; frozen GS was not included in this stability probe.

- All six local-correction schedules eventually fail the state guard for N3 at `h=.2`. At N32, `h=.08`, all local-correction schedules fail by four visits per distinct row. Frozen restoration also fails in some high-budget N32 runs; its surviving results at that timestep remain inaccurate.
- Near-rest rigid motion is extremely ill-conditioned under the current normalization. The regular near-rest oracle exists, but scalar GS remains slow; exact rest has a genuine rank defect. Positive scalar denominators do not detect either problem. See [RESULTS_RIGID_REVIEW.md](RESULTS_RIGID_REVIEW.md).
- Of 216 local five-second trajectories, 215 finish. The reverse-edge local method on N3 at `h=.08`, 32 visits per distinct row, fails after 2.4 s. Finishing trajectories can still have large energy or motion errors.
- The dense oracle conserves energy accurately but does not eliminate time-discretization phase error. For the affine oscillator it selects implicit midpoint: `omega_h=(2/h) atan(omega*h/2)`. With `omega=10` and `h=.08`, its frequency is about `9.513`, so even an exactly solved step accumulates substantial phase lag over five seconds.

For N3, `h=.08`, material→storage→kinetic and 32 visits per distinct row, the five-second maximum relative energy errors are `0.00280` dropped, `0.00162` local, `0.0405` frozen, and `1.11e-9` dense. Those are matched **visit** budgets, not matched compute. The same runs' position RMS errors are `0.0255`, `0.0262`, `0.0265`, and `0.0262` m. Lower position error in an inexact solve can come from cancellation with discretization error; it is not evidence that the objective was solved more accurately.

![Five-second trajectories and energy, with independent physical reference](results/figures/trajectories.png)

## 5. Validation and boundaries

The 4,729 checks cover sparse derivatives, independent objective/Hessian agreement, material-lift equivalence, local stationarity gathering, untouched coordinates, shared multipliers, restricted KKT solves, frozen Schur-GS equivalence and dropped-residual drift. Maximum errors are `3.37e-9` for gradient differences, `2.03e-10` for Hessian differences, `1.24e-15` for objective equivalence, and about `1e-16` for local-update/Schur equivalence.

All 33 primary oracle cases converge to regular roots with positive reduced curvature. This certifies local minima under the regularity assumptions, not global uniqueness. Exact rigid equilibria are separately classified as rank deficient. The physical reference uses paired RK4 subdivisions; maximum pair discrepancy is `1.23e-8` in short trajectories and `1.01e-7` in long trajectories, far below the reported motion errors.

The frozen work predictor is computed at the start of the timestep using an old-state stiffness solve. It remains a dense setup operation in these small-scene helpers. Only the **iteration visits** have the audited local structure. General graphs, multiple rigid rods, a new work law, a different independent-velocity objective, warm starting, GPU parallelism and a production animation viewer are outside this study.

## 6. What to keep and what to test next

Keep the full residual diagnostics, sparse accumulated multipliers, explicit schedules and independent oracle. Use frozen restoration as the economical reference iteration, with both material→storage→kinetic and storage→kinetic→material orders. Retain local correction as a diagnostic/control rather than assuming it is the best default.

The next numerical experiment should target the slow or unstable coupled modes: derive a small block of related material/storage/kinetic multipliers, then compare it against scalar GS at matched numerical work. Treat near-rest row scaling and geometry/balance coupling separately. A timestep subdivision control should distinguish excessive nonlinear step size from poor propagation. Any stabilization must be named and counted; successful energy projection alone is not an acceptance criterion.

## 7. Reproduce and inspect

Requires C++17/Eigen from the existing `constrained-xpbd` checkout. Plots need Python with NumPy and Matplotlib. No repository solver is overwritten.

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_COMPILER=clang++
cmake --build build -j3
ctest --test-dir build --output-on-failure
./build/experiments one results
./build/experiments controls results
./build/experiments scaling results
./build/experiments trace results
./build/experiments short results
./build/experiments long results
./build/fair_cost results
./build/stability_probe results/stability.csv
MPLCONFIGDIR=/tmp/xpbd-dual-gs-mpl python3 analyze.py results
```

Set `-DXPBD_ROOT=/path/to/constrained-xpbd` when needed. The existing local Python environment with plotting dependencies is `/home/blacknoir/Projects/mandelbrot-diffusion/.venv/bin/python`.

Entry points: [local solver](dual_gs.hpp), [experiment driver](experiments.cpp), [independent oracle](oracle.hpp), [validation](validation.cpp), [cost comparison](fair_cost.cpp), [protocol](EXPERIMENT_PROTOCOL.md), [raw results](results/), and [sweep-stability derivation](STABILITY_MATH.md).
