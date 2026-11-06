# Global energy: six linearized update schedules

This C++17 study compares direct endpoint energy constraints on an anchored 2D
spring and a hard pendulum. It has no work reservoir or velocity projection.
All local solves are linearized; schedules differ in when positions change,
when gradients are reevaluated, and how much of a frozen dual system is solved.
This is a research comparison, not a validated conservative integrator.

See [RESULTS.md](RESULTS.md) for measured outcomes. The [Jacobi versus GS rerun](JACOBI_VS_GS.md) adds independent simultaneous corrections to the [earlier frozen-schedule graphs](FROZEN_SCHEDULES.md), including both symmetric GS orders and immediate/deferred position application. The basic counterexample is
the spring released from rest: the joint method can preserve energy and close
every stationary equation by leaving the spring motionless.

## 1. Objective and physical energy

With constant acceleration \(a\),

\[
y=x_n+h v_n+h^2a,\qquad v(x)=(x-x_n)/h,
\]
\[
F(x)=\tfrac12(x-y)^T(M/h^2)(x-y)+U(x),\qquad
H(x)=\tfrac12v(x)^TMv(x)+U(x)-a^TMx-E_\star.
\]

The target is **fixed initial physical energy**, \(E_\star=E_0\).
For the spring, \(C=\|x\|-L,\ U=C^2/(2\alpha)\).
For the pendulum, \(C=0\) is a hard length constraint, \(\alpha=0,\ U=0\);
gravity supplies its potential energy.

Hard energy uses \(\min F\) subject to \(H=0\). Compliant energy uses

\[
\min_x F(x)+H(x)^2/(2\beta),\qquad \beta>0.
\]

The pendulum length constraint remains hard. The spring has a compliant material
relation, not an additional hard \(C=0\) constraint. Larger \(\beta\) weakens
the energy penalty; it does not conserve energy.

## 2. Dualization and residuals

Dualize the objective's elastic potential, keeping analytical position potential
inside \(H\). Use scaled multipliers \(\lambda\) and \(\eta=h^2\nu\):

\[
\mathcal L=I-\frac{\lambda C}{h^2}-\frac{\alpha\lambda^2}{2h^4}
-\frac{\eta H}{h^2}-\frac{\beta\eta^2}{2h^4}.
\]

For \(\alpha=0\), the material term imposes the hard constraint.
Let \(J=\nabla C^T,\ g=\nabla H,\ \tilde\alpha=\alpha/h^2,\ \tilde\beta=\beta/h^2\).
The scaled mixed residuals are

\[
R_x=h^2r_x=M(x-y)-J^T\lambda-g\eta,\quad
r_C=C+\tilde\alpha\lambda,\quad r_E=H+\tilde\beta\eta.
\]

CSV rx values report the **physical** norm \(\|r_x\|=\|R_x\|/h^2\).
No off-manifold substitution of \(\lambda\) for \(\nabla U(x)\) is made.
The equivalent gravity-centered gradient

\[
g=(M/h^2)\bigl(x-(x_n+h^2a)\bigr)+\nabla U(x)
\]

avoids manufacturing a nonzero gradient at a resting pendulum predictor.

## 3. Linearization and the two approximations

The exact top-left Jacobian is
\(A=M-\lambda\nabla^2C-\eta\nabla^2H\).
Every experiment drops curvature by replacing \(A\) with \(M\).
The kinetic energy Hessian is \(M/h^2\), not zero; its contribution is dropped.
Eliminating the position increment gives

\[
S=\begin{bmatrix}
JM^{-1}J^T+\tilde\alpha & JM^{-1}g\\
g^TM^{-1}J^T & g^TM^{-1}g+\tilde\beta
\end{bmatrix},
\]
\[
S\begin{bmatrix}\Delta\lambda\\\Delta\eta\end{bmatrix}
=-\begin{bmatrix}r_C\\r_E\end{bmatrix}
+\begin{bmatrix}J\\g^T\end{bmatrix}M^{-1}R_x,\qquad
\Delta x=M^{-1}(J^T\Delta\lambda+g\Delta\eta-R_x).
\]

XPBD-style variants additionally omit both appearances of \(R_x\).
Joint and frozen-GS controls retain them. Dropping curvature and dropping the
stationarity residual are distinct approximations.

## 4. Six schedules

Multipliers start at zero each timestep and accumulate within the solve.
Every run uses fixed sweep counts without early convergence exits.

| Name | One outer sweep | Algorithmic gradient evaluations |
|---|---|---:|
| baseline | Material scalar update, move | 1 |
| joint | Joint material/energy block, move once | 1 |
| immediate | Material, move; reevaluate; energy, move | 2 |
| frozen | Freeze both rows; material then energy; move once | 1 |
| reverse | Energy, move; reevaluate; material, move | 2 |
| symmetric | Material, energy, material; move and reevaluate between | 3 |

For frozen GS, denote the diagonal entries of \(S\) by \(a,d\), its off-diagonal
entry by \(b\), and its RHS by \(f_C,f_E\). Start total increments at zero:

\[
\Delta\lambda\leftarrow(f_C-b\Delta\eta)/a,\qquad
\Delta\eta\leftarrow(f_E-b\Delta\lambda)/d.
\]

The pending material displacement changes the energy RHS. Omitting that cross
term gives a different iteration. Apply the summed correction after 1, 5, or 20
inner passes. In retained-\(r_x\) controls, apply \(-M^{-1}R_x\) **once**.
For positive-definite \(S\), inner GS converges to the joint block, with
contraction factor \(b^2/(ad)\).

Immediate and deferred application with identical frozen gradients and correct
RHS updates are algebraically equivalent. Here immediate reevaluates between
updates, so it differs from frozen. Only joint/frozen have retained-\(r_x\)
controls. Symmetric uses three scalar solves per sweep versus two for immediate;
CSV counts expose unequal budgets. No optimized runtime benchmark is claimed.

## 5. Numerical policies

All corrections are full and undamped. There is no line search, artificial
denominator compliance, velocity rescaling, or random perturbation.
A numerically zero hard energy gradient with inconsistent RHS is reported and
skipped. An inconsistent rank-deficient joint block is also skipped.
For a redundant energy row, choose zero energy increment and solve the material
row. Continuing after a skip may leave a timestep unresolved.

The singular_updates counter counts unavailable attempted rows/blocks;
incompatible_updates counts their inconsistent RHS. The separate
rank_deficient_blocks counter records joint/frozen evaluations with a singular
full block. Frozen scalar GS can still execute with nonzero diagonal entries
of a rank-one block; this does not restore the block's rank.
Rank tests use 64 machine epsilons times a gradient scale, without changing
compliance.

The joint Schur complement is evaluated as nonnegative terms.
In a regular hard pendulum block, the two hard rows directly determine both
position coordinates. Computing those first avoids cancellation between huge
multiplier forces and \(R_x\). Retaining \(r_x\) changes multipliers but not
these positions. Numerical guard exits and compatible rank-deficient blocks
can still behave differently.

A run aborts if position/speed exceeds \(10^6\) in the chosen units, energy
error exceeds \(10^{12}E_s\), either scaled multiplier exceeds \(10^{18}\),
or a correction becomes undefined/nonfinite. Nothing is clipped.
The last completed timestep and rejected attempt are recorded separately.
These permissive guards prevent overflow; they are not physical acceptance
thresholds. Finite completed runs can fail badly.

Diagnostics include physical energy error \(H\), material residual \(r_C\),
energy row \(r_E\), physical stationarity \(r_x\), displacement, path length,
kinetic energy, and position error against a continuous reference.
For hard pendulums, a circle-intersection check distinguishes infeasible
endpoint constraints from singular local linearizations. An unresolved
timestep can make a later fixed-target problem infeasible.

The rows_closed_steps check requires \(10^{-8}L\) material and \(10^{-8}E_s\)
energy-row tolerances; kkt_closed_steps also requires
\(\|r_x\|<10^{-8}E_s/L\). Baseline counts ignore the energy row.
These check equations, not motion. RMS errors for aborted runs cover only their
accepted prefixes. Initial trace residuals are placeholders before any solve.
A late_mean_kinetic value is meaningful only with nonzero
late_samples; a run stopped before 7.5 seconds has no late-time observation.

## 6. Matrix and reproduction

There are 1,020 ten-second runs:

- Spring \(m=L=1,\ k=100,\ x_0=(1.25,0)\), initially at rest, moving radially
  with velocity \((-0.7,0)\), or tangentially with velocity \((0,1)\).
- Pendulum \(m=L=1,\ g=9.81\), angle 0.65 radians from downward vertical,
  initially at rest or with tangential speed 1.2.
- \(h=1/240,1/60\), and 1, 5, 20 outer sweeps.
- Energy softness \(\beta/E_s=0,0.01,1\), with \(E_s\) initial energy above
  minimum potential. Thus \(\beta\) has energy units.
- Retained \(r_x\) controls for joint/frozen; 1, 5, 20 frozen inner passes.

The continuous reference uses RK4 with 32 substeps, checked against 64 and
against physical energy drift. Refinement is a numerical accuracy check.
Independent tests verify production XPBD and analytic baseline parity,
energy gradients, the frozen spring root, singular rows, soft energy without
conservation, one-pass GS versus Jacobi, and joint/frozen updates against an
uneliminated mixed-system solve.

From the repository root:

~~~sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
ctest --test-dir build --output-on-failure
./build/xpbd-energy-ordering-study results/energy_ordering
python3 experiments/energy_ordering/analysis.py results/energy_ordering
~~~

Use the analysis option --no-plots for CSV checks without Matplotlib.
Plot generation requires Matplotlib. The study executable has no additional
compiled dependency; independent tests use the existing Eigen dependency.
summary.csv covers every run; traces.csv contains the 60 Hz / five-sweep slice;
reference_checks.csv records reference refinement. The analysis script writes
plots and derived tables.

## Why convergence is insufficient at rest

Direct expansion gives

\[
F(x)=H(x)-v_n^TM(x-x_n)/h+\text{constant}.
\]

At \(v_n=0\), the hard-energy objective is constant on the feasible set.
The radial spring's frozen branch has
\(\Delta x=0,\lambda=-h^2kC_n,\eta=h^2\), despite a nonzero elastic force.
It is not caused by an inaccurate block solve or omission of \(r_x\).
Material-first reaches the backward-Euler minimum of endpoint energy, where
\(g=0\) but \(H<0\); a hard linearized energy row cannot correct it.

With compliance at rest,
\(F+H^2/(2\beta)=(H+\beta)^2/(2\beta)+\text{constant}\).
This objective can prefer an energy deficit. Compliance regularizes equations;
it does not establish conservation or correct acceleration.
