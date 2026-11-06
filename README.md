# Constrained XPBD: local work reservoirs

A small, headless C++17 implementation of three stages:

1. **xpbd**: ordinary distance-constraint XPBD.
2. **passive**: identical XPBD trajectories, with unmatched work stored after each timestep.
3. **active**: reservoir feedback through compliant local work balances, with selectable paired or sequential multiplier updates.

The active mode is an experiment. It can return or absorb stored energy, but it does **not** guarantee conservation of mechanical energy, momentum, or angular momentum. Strong feedback can destabilize cloth. See [VALIDATION.md](VALIDATION.md) for measured outcomes, including failures.

The [implicit-work experiments](experiments/work_balance/README.md) compare kinetic-exact work, midpoint discrete gradients, and curvature-dropped Gauss–Seidel with and without the position residual. They include spring and pendulum failure cases and a separate derivation of the additional dual variables needed for a variational model.

The [direct-energy coupling study](experiments/energy_coupling/README.md) tests joint material/energy multiplier updates against sequential Gauss–Seidel on a 2D spring and pendulum. It records unstable runs, singular rows, frozen turning points, and converged cases separately.

The [global-energy ordering study](experiments/energy_ordering/README.md) compares six linearized schedules, hard/compliant energy, and retained-position-residual controls across 1,020 spring/pendulum runs. It measures motion against a continuous reference as well as energy, material, and stationarity residuals. The [Jacobi versus GS rerun](experiments/energy_ordering/JACOBI_VS_GS.md) adds independent simultaneous multiplier updates and provides energy/position graphs across an expanded 1,230-run matrix.

The [sequential-work study](experiments/sequential_work/README.md) tests the per-spring order `lambda -> nu -> mu_a -> mu_b`, retaining one shared kinetic multiplier per particle. Separate potential and kinetic compliances control feedback strength. This comparison keeps the previous-force work estimate. The [behavior report](experiments/sequential_work/behavior/README.md) shows measured oscillation growth, damping, and distorted stationary states on the spring, pendulum, kicked 2D grid, and cloth.

The [small-step EL comparison](experiments/small_steps/README.md) uses one energy-then-material sweep per substep throughout. Its 465-run matrix tests relaxation, correction caps, energy compliance, and material-adjusted direction tests on five spring-pendulum motions, recording energy, trajectory, torque balance, and solver work.

## Recent research studies

- [Original CPD comparison](experiments/cpd_original_2026_11_02/README.md):
  paper-equation implementation, momentum controls and motion failures.
- [Independent endpoint velocity](experiments/independent_velocity_2026_11_01/README.md):
  approximate energy conservation with separate velocity variables.
- [Block-order reassessment](experiments/block_orders_reassessment_2026_11_01/README.md):
  motion and energy evidence for material/energy update orderings.

Source, derivations and compact numerical summaries are included. Generated
figures, videos, CSV trajectories, logs and viewers remain local and are listed with
sizes and SHA256 hashes in [the artifact manifest](experiments/local_artifacts.json).
Some report links refer to these locally retained artifacts; regenerate them
with the commands in the corresponding study README. Existing provenance
manifests describe the full local campaigns, not only the committed subset.

## Build and run

```sh
git submodule update --init external/eigen
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
ctest --test-dir build --output-on-failure

# All four scenes and all three methods, ten simulated seconds.
./build/constrained-xpbd --steps 2400 --output results

# Select one comparison. Larger beta/gamma mean weaker feedback.
./build/constrained-xpbd --scene cloth --mode active --beta 1000 --gamma 1000
# Requested per-edge scalar order; beta is kinetic, gamma is potential compliance.
./build/constrained-xpbd --scene spring --mode active --order sequential --beta 1000 --gamma 1000
./build/constrained-xpbd --help
```

Eigen is the only compiled dependency. The existing raylib submodule is unused; this prototype has no viewer. No network access is needed after Eigen is initialized.

The scenes are a pinned scalar spring, a hard-distance pendulum, a planar 6x6 spring grid, and a 9x9 cloth patch. The grids have structural and diagonal distance springs. Pins are stationary. Cloth has no collision handling or bending model.

## Work accounting

The implicit objective is

\[
F(x)=\tfrac12(x-y)^T(M/h^2)(x-y)+\sum_j u_j(x),
\qquad y=x_n+h v_n+h^2g.
\]

For compliant springs, \(u_j=C_j^2/(2\alpha_j)\). A hard constraint has \(\alpha_j=0\), imposes \(C_j=0\), and has no elastic potential assigned to its multiplier.

Use **changes**, not absolute energies, in the work rows:

\[
K_i=\tfrac12m_i(\|v_i^+\|^2-\|v_i^-\|^2),\qquad
U_j=u_j(x)-u_j(x_n),\qquad v_i^+=(x_i-x_{n,i})/h.
\]

Freeze the previous material force for the entire timestep:

\[
W_{ij}=f_{ij}^{-T}(x_i-x_{n,i}),\qquad
f_{ij}^-=\nabla_i C_j(x_n)\lambda_j^-/h_{\rm previous}^2.
\]

On the first timestep there is no history: use the analytical elastic force for compliant springs, and zero for hard constraints. Previous XPBD force estimates need not satisfy the exact constitutive law if their solve was incomplete. Energy accounting always uses the analytical position potential.

With constant gravity, define

\[
b_i=K_i-\sum_jW_{ij}-m_i g^T(x_i-x_{n,i}),\qquad
d_j=U_j+\sum_iW_{ij}.
\]

The signed reservoirs obey the local work balances

\[
b_i+R_i^+-R_i^-=0,\qquad d_j+S_j^+-S_j^-=0.
\]

Their sum conserves \(T+U+V_g+\sum_iR_i+\sum_jS_j\), where \(V_g=-\sum_i m_i g^Tx_i\). Positive reservoir contents represent energy available to return; negative contents represent energy to take back. These are signed error accounts, not nonnegative physical storage elements.

Passive mode only closes these accounts after ordinary XPBD. Active mode adds

\[
\sum_i (R_i^+)^2/(2\beta)+\sum_j(S_j^+)^2/(2\gamma)
\]

to the objective, subject to the same balances. Penalties act on the accumulated contents, not just the latest change. Eliminating reservoirs gives \(R_i^+=R_i^--b_i\), \(S_j^+=S_j^--d_j\), so each penalty encourages repayment through its local work row.

## Active update and approximations

Dualize only the objective's elastic energy; keep the analytical \(U_j(x)\) in the work rows. Eliminating reservoirs yields

\[
\mathcal L=I-\sum_j\left(\lambda_j C_j/h^2+\alpha_j\lambda_j^2/(2h^4)\right)
-\sum_i\mu_i(b_i-R_i^-)-\sum_j\nu_j(d_j-S_j^-)
-\tfrac\beta2\|\mu\|^2-\tfrac\gamma2\|\nu\|^2.
\]

This equivalent formulation avoids dividing by \(1-\nu_j\). It also keeps position energy gradients valid away from the material constitutive manifold.

The default `--order paired` retains the original paired material/work solve followed by one kinetic visit per particle. `--order paired-local` uses the same pair solve followed immediately by the two endpoint kinetic visits. `--order sequential` instead visits the four scalar rows `lambda_j -> nu_j -> mu_a -> mu_b`, relinearizing between updates. Endpoint visits share one accumulated multiplier per particle and include work from all incident springs. Isolated particles are visited once per sweep.

Store scaled work multipliers \(q_\mu=h^2\mu\), \(q_\nu=h^2\nu\). Linearize, replace the position Hessian by the mass matrix, and drop the position stationarity residual \(r_x\), as in the XPBD approximation. With \(g_j=\nabla C_j\), \(z_j=\nabla d_j\), each material/work pair solves

\[
\begin{bmatrix}
g_j^TM^{-1}g_j+\alpha_j/h^2 &g_j^TM^{-1}z_j\\
z_j^TM^{-1}g_j&z_j^TM^{-1}z_j+\gamma/h^2
\end{bmatrix}
\begin{bmatrix}\Delta\lambda_j\\\Delta q_{\nu,j}\end{bmatrix}
=-\begin{bmatrix}C_j+\alpha_j\lambda_j/h^2\\d_j-S_j^-+\gamma q_{\nu,j}/h^2\end{bmatrix}.
\]

Apply \(\Delta x=M^{-1}(g_j\Delta\lambda_j+z_j\Delta q_{\nu,j})\). Then visit each particle balance with the scalar XPBD rule for \(b_i-R_i^-\), gradient \(\nabla b_i\), and compliance \(\beta/h^2\). Update positions and relinearize between visits. A local residual backtracking step limits nonlinear overshoot; it is not a global objective line search. The 2x2 solve uses a positive-term Schur complement to avoid determinant cancellation.

Dropping curvature and dropping \(r_x\) are separate approximations. The method does not claim to solve the full optimization stationarity equations. Diagnostics expose the omitted \(r_x\) and the remaining material/feedback residuals.

For the sequential schedule, each scalar row has update

\[
\Delta q=-\frac{f(x)+c q/h^2}{\nabla f^TM^{-1}\nabla f+c/h^2},
\qquad \Delta x=M^{-1}\nabla f\,\Delta q.
\]

Use \((f,q,c)=(C_j,\lambda_j,\alpha_j)\) for the material row,
\((d_j-S_j^-,q_{\nu,j},\gamma)\) for the potential row, and
\((b_i-R_i^-,q_{\mu,i},\beta)\) for the kinetic row. The implementation accepts finite positive beta/gamma; smaller means stronger feedback, larger means cheaper reservoir storage. These parameters have energy units. Their accumulated multiplier terms make the rows compliant; multiplying a position correction by a gain would be a different operation. As either compliance approaches zero, its row approaches a hard work condition, which need not be compatible with the other rows.

At the accepted positions, close the accounts exactly via \(R^+=R^--b\), \(S^+=S^--d\), including after an incomplete solve. The feedback conditions \(R^+=\beta\mu\), \(S^+=\gamma\nu\) may remain unsatisfied; their discrepancy is reported separately. Exact bookkeeping alone is not evidence of correct dynamics.

With \(E\) particle/spring incidences, a paired sweep costs \(O(E)\) and a particle sweep \(O(N)\). Scratch memory is \(O(E+N)\); the global Schur matrix is never assembled. This implementation prioritizes readable equations over allocation/performance optimization.

## Files and diagnostics

- `xpbd.hpp`: model, baseline, accounting, active local solves, diagnostics.
- `scenes.hpp`: the four reproducible base cases.
- `main.cpp`: command line, comparison summaries, per-step CSV and final particle positions.
- `tests.cpp`: baseline physical checks and passive trajectory parity.
- `active_tests.cpp`: independent work audits, return/absorption, history scaling, finite-iteration residuals, and ten-second scene runs.
- `sequential_tests.cpp`: independent ordered-update calculations, shared particle rows, compliance controls, and accounting checks.

CSV output separates mechanical energy, particle and spring reservoirs, their norm, extended-energy drift, deformation, and solver residuals. Timing includes diagnostics and CSV output when enabled; it is not an optimized solver benchmark.

Local elastic repayment remains an open issue. For affine constraints and an exact old constitutive force, \(d_j=(\Delta C_j)^2/(2\alpha_j)\ge0\), so an independent elastic reservoir cannot repay accumulated negative contents. Reservoir exchange or a different work rule would change this. This implementation preserves the formulation we selected rather than silently replacing it with a global energy correction.

The [spring repair study](experiments/spring_repair/README.md) diagnoses the direct-energy objective's turning-point problem, tests midpoint XPBD, and implements a scalar work-balance spring solve that conserves energy and anchor angular momentum in the tested 2D cases.

The [gradient-policy study](experiments/gradient_policies/README.md) tests skipping oversized energy increments and reusing older gradients within the original endpoint-energy equations. It includes 1,350 spring/pendulum runs and separate checks of energy, motion, and retained versus dropped position residuals.

The [elastic-pendulum study](experiments/elastic_pendulum/README.md) tests fresh EL under gravity before extending to multiple springs. It includes baseline checks, exact vertical-motion validation, and 270 runs; small swings can work, while turning-point freezing and coupled-motion errors remain.
