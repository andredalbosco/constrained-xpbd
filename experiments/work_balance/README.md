# XPBD work balance: GS experiments and diagnostic references

The target is a **Gauss–Seidel solver with linearization and omitted curvature**, not a replacement by a dense Newton solver or an exact trajectory integrator. The main experiment asks whether retaining the position/momentum residual closes the gap between implicit material work and kinetic work. Exact scalar, closed-form pendulum, and Newton-solved radial examples are diagnostic references only.

The results support retaining that residual while dropping curvature. They do **not** establish exact conservation after an arbitrary finite number of GS sweeps, convergence for general cloth, or an equivalent variational dualization of the nonlinear discrete-gradient equations.

## Midpoint equations and local GS update

Let \(x_-,v_-\) be the physical starting state, \(q=(x_++x_-)/2\) the unknown midpoint position, and \(g\) constant gravitational acceleration. Define

\[
y=x_-+\frac h2v_-+\frac{h^2}{4}g,
\qquad x_+=2q-x_-,
\qquad v_+=\frac4h(q-x_-)-v_-.
\]

For compliant elasticity \(U=\tfrac12C^T\alpha^{-1}C\), use

\[
\bar C=\frac{C(x_+)+C(x_-)}2,
\qquad \bar J(x_+-x_-)=C(x_+)-C(x_-).
\]

The compatible residual equations are

\[
\boxed{
r_x=M(q-y)-\bar J^T\lambda,
\qquad r_\lambda=\bar C+\frac{4\alpha}{h^2}\lambda.
}
\]

The midpoint force factor is \(h^2/4\), so the scaled compliance is **\(4\alpha/h^2\)**. Changing the compliance alone without the corresponding reconstruction is a different method.

Since \(D_q\bar C=J_{\rm end}=J(x_+)\), linearizing and dropping the curvature in \(D_qr_x\) gives

\[
\left(J_{\rm end}M^{-1}\bar J^T+\frac{4\alpha}{h^2}\right)\Delta\lambda
=-r_\lambda+J_{\rm end}M^{-1}r_x,
\qquad
\Delta q=M^{-1}(\bar J^T\Delta\lambda-r_x).
\]

For one scalar block the matrix expression is a scalar denominator. For coupled constraints the study solves the reduced linearization with inner GS sweeps, then updates the nonlinear state with backtracking. It uses no constraint Hessians. A further tested approximation replaces \(J_{\rm end}\) by \(\bar J\) in the linear solve while continuing to measure the actual residuals. An arbitrary nonlinear configuration need not give a positive nonsymmetric denominator.

## Why the implicit work needs the momentum equation

The material impulse in midpoint variables is

\[
I_j=\frac4h\bar J_j^T\lambda_j,
\qquad W_j=\bar v^TI_j,
\qquad \bar v=\frac{v_++v_-}{2}.
\]

For the accepted state, including the gravity impulse separately,

\[
\boxed{
\Delta T-W_{\rm ext}-\sum_jW_j
=\frac4h\bar v^Tr_x.
}
\]

The material equation supplies the elastic counterpart:

\[
\Delta U+\sum_jW_j
=r_\lambda^T\alpha^{-1}\Delta C.
\]

Thus a small constitutive residual alone is insufficient. A nonzero \(r_x\) does not necessarily cause energy error: its projection onto \(\bar v\) can vanish, as in the special radial pendulum case. Momentum and trajectory errors still require separate checks.

Alternatively, recording every accepted solver impulse gives the exact telescoping identity

\[
\delta W_i=\frac{v_i^{\rm before}+v_i^{\rm after}}2\cdot\delta I_i
=\delta T_i.
\]

This must include predictors and all correction impulses. It is algorithmic-path accounting; individual source attribution depends on update order, and it does not automatically equal minus the analytical potential change. The first-order expression \(v^{\rm before}\cdot\delta I\) omits \(\|\delta I\|^2/(2m)\).

## Measured GS evidence

The following recorded runs use \(h=1/60\), 10 seconds, and **10 inner linear GS sweeps per nonlinear iteration**. The two-spring case has one moving particle coupled to two nonlinear anchored springs. Energy errors are absolute joules; the position residual is \(\|M^{-1}r_x\|\).

| Case and update | Outer limit | Max energy error | Max position residual | Max constitutive residual |
|---|---:|---:|---:|---:|
| Two springs, drop \(r_x\) | 5 or 80 | 0.136570 | 7.74674e-5 | 9.30e-15 |
| Two springs, retain \(r_x\), drop curvature | 5 | 1.83578e-9 | 1.09746e-12 | 5.10e-16 |
| Two springs, retain \(r_x\), drop curvature | 80 | 1.18744e-10 | 9.98703e-15 | 4.81e-16 |
| Pendulum, retain \(r_x\), drop curvature | 5 | 1.19016e-13 | 8.51152e-16 | 2.22e-16 |

The retained-residual two-spring run uses at most seven outer iterations when allowed 80. Dropping \(r_x\) stalls with a per-step kinetic/material-work defect as large as 0.0352600 J despite an almost-zero constitutive residual. These are measured outcomes at these settings, not a general finite-sweep guarantee. See `results/work_balance/midpoint-gs/summary.csv` for all limits and stepsizes.

The cheaper **one-inner-sweep** variant also converges in this case: retaining \(r_x\) gives maximum energy errors 3.89235e-9 J with outer limit 5 and 1.13527e-11 J with limit 20, using at most seven outer iterations. Dropping it leaves 0.131530 J at limits 20 and 80. One sparse inner sweep and residual evaluation cost \(O(E+N)\); backtracking can add residual evaluations. The CSV covers 216 combinations, including one and ten inner sweeps.

A separate single radial-spring block gives the same distinction: at \(h=1/60\), five curvature-dropped iterations retaining \(r_x\) give maximum relative energy error 1.44e-13; dropping it gives 6.55e-7 even after 100 iterations. This single-block result does not establish coupled GS convergence. A linear scalar spring is a special case with constant gradient: midpoint XPBD agrees with implicit midpoint and preserves its quadratic energy to roundoff, while retaining phase error.

For the backward-Euler scalar spring at \(h=1/240\) over 20 seconds, exact impulse-work accounting has a kinetic residual below 4.2e-16 while mechanical energy loses 99.9758%. Work accounting alone does not repair the trajectory.

### Coupled free-body momentum check

A separate 144-run study reuses the same GS core on three free particles connected by three nonlinear springs, without gravity or pins. The initial velocity field combines translation and rotation; a second initialization adds deformation. The springs are initially strained, so “rigid” in the CSV labels the initial velocity field, not rigid-body motion throughout the run.

For the first initialization, \(h=1/15\), 10 seconds, **one inner GS sweep and outer limit 20**, the measured maximum absolute errors are:

| Update | Energy, J | Linear momentum, kg m/s | Angular momentum, kg m²/s |
|---|---:|---:|---:|
| Drop \(r_x\) | 0.0268623 | 1.74391e-12 | 0.0148378 |
| Retain \(r_x\), drop curvature | 8.35509e-12 | 4.52239e-13 | 6.36948e-12 |

The dropped-residual run has a constitutive residual below 1e-14 but fails the full residual criterion on every step. The retained-residual run meets that criterion on every step. Pairwise corrections preserve linear momentum to floating-point accuracy in both runs; angular momentum distinguishes this coupled case from the special single-spring and pendulum geometry. The program checks per-step momentum-error bounds derived from \(r_x\); these are residual-dependent bounds, not unconditional conservation assertions. Raw data are in `results/work_balance/coupled-momentum/coupled_momentum_summary.csv`.

## Hard pendulum: energy and length are not endpoint tangency

For \(c(r)=\tfrac12(\|r\|^2-\ell^2)\), the discrete gradient is \(\bar r\). A midpoint radial reaction does zero work when both endpoint lengths agree. However,

\[
\boxed{r_+\cdot v_++r_-\cdot v_-=\frac{h}{2m}\Delta T.}
\]

The hard-row extension \(\bar C=0\) assumes a consistent starting geometry: it implies \(C(x_+)=-C(x_-)\), so a pre-existing constraint error can alternate rather than disappear. Enforcing \(C(x_+)=0\) is a separate endpoint-constraint choice. The recorded hard-pendulum runs start on the constraint manifold.

With changing speed under gravity, Cartesian midpoint drift and fixed length cannot generally preserve both endpoint tangencies. The GS pendulum run above reaches maximum radial velocity **4.24280e-4 m/s** despite its small energy error. Dropping \(r_x\) also happens to preserve energy in this geometry because the residual is work-orthogonal; it does not converge the full equations.

The closed-form Cartesian reference at the same step gives energy change -2.58e-13 J and the same radial-velocity defect. Projecting its velocity tangent removes approximately **1.07214e-4 J** over 10 seconds and breaks Cartesian midpoint drift. An angular-coordinate discrete-gradient reference preserves tangency and energy but uses different kinematics; it is not the proposed GS method. Large-step singularities and branch failures are recorded in `pendulum_singularities.csv`.

## Extra dualization: what has and has not been derived

The nonlinear midpoint residuals above are generally **not** the direct gradient of a scalar Lagrangian in \((q,\lambda)\): their mixed derivatives contain \(J_{\rm end}\) and \(\bar J\), which differ. They define a concrete GS residual solve; a variational derivation requires a separately specified augmented problem.

With independent impulses \(I_{ij}=h\bar J_{ji}^T\ell_j\) and implicit work \(W_{ij}=I_{ij}\cdot\bar v_i\),

\[
\Delta U_j+\sum_iW_{ij}
=\Delta C_j\left(\frac{\bar C_j}{\alpha_j}+\ell_j\right).
\]

This determines the material law only away from \(\Delta C_j=0\). Retain the material law at zero work rather than dividing by \(\Delta C_j\). Impulse balance already makes the kinetic work row an identity, so adding it independently creates redundancy.

Actual implicit-work derivatives include \(\nabla_IW=\bar v\), \(\nabla_vW=I/2\), and the mixed Hessian \(I_3/2\). Eliminating velocity makes work quadratic across impulses. `formalism_notes.md` gives an augmented scalar optimization using independent positions, velocities, impulses, force coefficients, and reservoirs, with extra multipliers for drift and impulse balance. It is an explicitly labeled formulation candidate, not a proved equivalent or a tested solver.

## Files and reproduction

| Source or note | Purpose | Outputs under `results/work_balance/` |
|---|---|---|
| `midpoint_gs.hpp` | Reusable curvature-dropped midpoint/DG GS core | — |
| `midpoint_gs_study.cpp` | Target: nonlinear linearization, curvature drop, dual GS | `midpoint-gs/summary.csv` |
| `coupled_momentum_study.cpp` | Coupled free-body energy, linear/angular momentum | `coupled-momentum/coupled_momentum_summary.csv`, `run.log` |
| `kinetic_work_study.cpp` | Work accounting on unchanged XPBD trajectories | `kinetic-work/summary.csv` |
| `spring_study.cpp`, `spring_notes.md` | Linear/radial references and single-block linearized variants | `spring/spring_*.csv`, `run.log` |
| `pendulum_study.cpp`, `pendulum_notes.md` | Cartesian/angle references, tangency, singularities | `pendulum/pendulum_*.csv`, `run.log` |
| `formalism_notes.md` | Derivations, residual identities, variational boundary | — |

The integrated targets are enabled by default through `XPBD_BUILD_WORK_STUDIES`. From the repository root:

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j4
ctest --test-dir build --output-on-failure -j4
python3 experiments/work_balance/run_studies.py
```

The runner defaults to `build/` and `results/work_balance/`; `--build-dir` selects another build directory. Executables are `xpbd-midpoint-gs-study`, `xpbd-coupled-momentum-study`, `xpbd-kinetic-work-study`, `xpbd-spring-study`, and `xpbd-pendulum-study`. The spring, pendulum, and coupled-momentum programs write their CSV files themselves; the other two emit CSV to standard output. Small residuals and conservation errors should be read alongside trajectory errors, tangent-velocity error, solver limits, and rejected/singular cases.

The integrated Release build passed all seven CTest targets and all five study runs. An independent CSV scan confirmed finite numeric values and the expected row counts, including 216 midpoint-GS and 144 coupled-momentum cases.

## Primary sources

- [Macklin, Mueller, Chentanez, XPBD](https://mmacklin.com/xpbd.pdf), equations 12–17: curvature approximation and separate omission of the position residual.
- [Gonzalez, Time Integration and Discrete Hamiltonian Systems, 1996](https://web.ma.utexas.edu/users/og/PUBLICATIONS/paper_DisHamSys.pdf), sections 3.2 and 3.6–3.7: discrete chain rules, symmetry compatibility, and radial secant-force scheme. The conservative reference construction is established methodology.
- [Gonzalez, Mechanical systems subject to holonomic constraints, 1999](https://web.ma.utexas.edu/users/og/PUBLICATIONS/paper_Constrnts.pdf), equations 21–24 and pages 169/172: constrained conservation and the explicit limitation concerning endpoint velocity tangency.
