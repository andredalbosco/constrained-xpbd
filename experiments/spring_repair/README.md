# Spring repair experiments

The current backward-Euler objective with hard endpoint energy has a turning-point
problem. Midpoint time integration repairs the harmonic spring; an averaged
central force with exact work balance also conserves energy and angular momentum
for the tested tangential 2D spring. Both use scalar linearized multiplier
updates with curvature dropped. The core solver and previous studies remain available.

## Why the current objective stalls

For scalar extension q, stiffness k, and mass m, the current equations use

\[
v_+=(q_+-q_n)/h,\quad
F=\frac{m}{2h^2}(q_+-q_n-hv_n)^2+\frac{k}{2}q_+^2.
\]

On the hard energy constraint \(\tfrac12mv_+^2+\tfrac12kq_+^2=E_\star\),

\[
F=E_\star-mv_nv_++\tfrac12mv_n^2.
\]

If the old state also has energy E_star, the two possible new velocities obey

\[
(m+kh^2)v_+^2+2khq_n v_+-mv_n^2=0.
\]

For nonzero old velocity the roots have opposite signs. The global minimum
of F chooses the new velocity with the same sign as the old one. It cannot
choose the reversal required for an oscillation. At rest, freezing and moving
have the same objective value. Both isolated feasible roots can satisfy the
stationarity equations; selecting a different root requires an additional
branch-selection rule rather than a better minimization of F.

## Solver-only controls

`controls_study.cpp` tests 448 combinations of ordering, relaxation, retained
or dropped r_x, and ordinary or half-force initial guesses. A further 52 runs
test the most promising bounded-motion candidate over 10/100 seconds, timesteps
1/30 through 1/240, and 5/20/80 outer sweeps.

That candidate uses EL, energy-increment relaxation 0.5, a half-force initial
guess, and dropped r_x. At 60 Hz and 20 sweeps it has relative energy errors
3.65e-6/4.04e-6 and position RMS errors 0.04849/0.06287 m for rest/radial velocity.
It is not a reliable repair: smaller timesteps do not show consistent motion
convergence; 80 sweeps drive energy error to roundoff without improving motion;
100-second runs have RMS position errors around 0.23–0.28 m. Five sweeps often
fail. Retaining r_x makes the selected rest/radial cases hit guards after 7/35
steps. The half-force guess introduces a stationarity residual that remains
unbalanced when r_x is dropped.

## Midpoint with the ordinary XPBD constraint solve

Use midpoint position z as the unknown and minimize

\[
F_{\mathrm{mid}}(z)=
\frac12\left(z-x_n-\frac h2v_n\right)^T\frac{4M}{h^2}
\left(z-x_n-\frac h2v_n\right)+U(z).
\]

This is ordinary XPBD at effective timestep h/2. For each spring,

\[
\widetilde\alpha=4\alpha/h^2,\qquad
\Delta\lambda=-\frac{C(z)+\widetilde\alpha\lambda}
{JM^{-1}J^T+\widetilde\alpha},\qquad
\Delta z=M^{-1}J^T\Delta\lambda.
\]

Reconstruct endpoints as

\[
x_+=2z-x_n,\qquad v_+=\frac{4(z-x_n)}h-v_n.
\]

This deliberately changes both the integration stage and velocity update.
The chord velocity (x_plus-x_n)/h is now the average velocity, not the endpoint
velocity used in mechanical energy. It is not velocity rescaling or an
energy projection. The scalar GS update remains linearized with curvature
dropped; no energy multiplier is needed for the harmonic spring.

At a solved step,

\[
\Delta K=-\nabla U(z)^T(x_+-x_n).
\]

For quadratic U the right-hand side equals -DeltaU exactly. With
\(a=h^2k/m\), the scalar map is

\[
q_+=\frac{(1-a/4)q_n+hv_n}{1+a/4},\qquad
v_+=\frac{(1-a/4)v_n-h(k/m)q_n}{1+a/4}.
\]

This conserves harmonic energy and reverses at turning points. Its numerical
frequency is \(2\arctan(h\sqrt{k/m}/2)/h\), so it still has phase error.

For a general 2D radial spring, U(x)=k(|x|-L)^2/2 is not quadratic in Cartesian
coordinates. Midpoint central forces preserve angular momentum, but their
energy defect is

\[
\Delta E=kL(r_+-r_n)\left(\frac{r_++r_n}{2|z|}-1\right).
\]

This defect vanishes for the tested collinear motions. It need not vanish for
tangential motion. Imposing a hard rod constraint at z would not automatically
preserve its endpoint length, so this spring change is not a validated pendulum
or cloth implementation.

At 60 Hz over ten seconds, one scalar material update per step:

| Initial motion | Midpoint max relative energy error | Midpoint position RMS | Standard XPBD position RMS |
|---|---:|---:|---:|
| Rest | 2.08e-13 | 0.02366 m | 0.15969 m |
| Radial | 2.67e-13 | 0.02451 m | 0.16606 m |
| Tangential | 7.83e-4 | 0.02219 m | 0.87536 m |

At 240 Hz, tangential energy error is 4.90e-5 and position RMS is 0.001394 m.
Radial ten-second phase lag decreases from 0.23052 to 0.01446 radians, matching
the discrete frequency formula. Angular momentum error at 60 Hz is below 8e-14.
One, five, and twenty sweeps give identical trajectories for this single spring;
one update already solves its midpoint stationarity equation.

![Midpoint comparison](../../results/spring_repair/midpoint/plots/midpoint_spring_h60_sweeps1.png)

## Exact work balance for the 2D radial spring

Let \(r_n=|x_n|\), \(r_+=|x_+|\), and \(\bar r=(r_n+r_+)/2\).
Replace the midpoint force by

\[
\bar f=-k\left(1-\frac L{\bar r}\right)z.
\]

Its work is exactly the negative potential change, because

\[
\Delta U=k(\bar r-L)(r_+-r_n)
=k\left(1-\frac L{\bar r}\right)z^T\Delta x.
\]

Solving \(M(v_+-v_n)=h\bar f\) with the same midpoint kinematics gives
\(\Delta K=\bar f^T\Delta x=-\Delta U\). The force is parallel to z, so
\(\Delta L=h\,z\times\bar f=0\) about the anchor.

To retain a scalar multiplier solve, define

\[
n=z/|z|,\quad C_w=|z|\left(1-\frac L{\bar r}\right),\quad
R=M(z-y)-n\lambda,\quad r_c=C_w+\frac{4\alpha}{h^2}\lambda,
\quad y=x_n+\frac h2v_n.
\]

With \(n_+=x_+/|x_+|\), the exact row gradient is

\[
d=\nabla_z C_w
=n\left(1-\frac L{\bar r}\right)+\frac{|z|L}{\bar r^2}n_+.
\]

Drop the geometric term \(\lambda Dn\) in the position Jacobian. For scalar
mass m the update is

\[
b=-R/m,\qquad
\Delta\lambda=\frac{-r_c-d^Tb}{d^Tn/m+4\alpha/h^2},\qquad
\Delta z=b+n\Delta\lambda/m.
\]

Reevaluate these rows between scalar visits. There is no coupled energy
multiplier, energy reservoir, nonlinear root solve inside a visit, or velocity
rescaling. This is a new discrete force/work formulation. It is **not the
dualization of the old objective**: the top multiplier direction n and bottom
row gradient d generally differ, making the linear system nonsymmetric.
For collinear same-ray motion it reduces to ordinary midpoint XPBD.

Conservation is tied to residual convergence. At any finite iterate,

\[
\Delta E=\Delta x^T\left(4R/h^2+knr_c\right),\qquad
\Delta L=(4/h)\,z\times R.
\]

The prototype checks these identities, the potential-work identity, and the
row gradient independently. Undefined directions and nonpositive scalar
denominators are reported explicitly rather than regularized. This does not
prove global convergence for large steps or collapsed springs.

On the tangential spring, five scalar visits at 60 Hz over ten seconds give
relative energy error about 2.53e-13, angular error 1.49e-13, and position RMS
0.02343 m. At 240 Hz position RMS is about 0.001472 m. The 100-second 60 Hz run
stays bounded with relative energy error about 4.20e-13; its position RMS grows
to 0.203 m from phase drift. Conservation does not eliminate discretization error.
One visit has small nonzero energy error, so the five-visit result should not
be attributed to a single approximate correction.

These tests cover an anchored spring, whose anchor exerts external force.
Linear momentum conservation for free pairs, multiple interacting springs,
and hard pendulum/cloth constraints require separate validation.

![2D spring work-balance comparison](../../results/spring_repair/work_midpoint/plots/work_midpoint_tangent_h60.png)

[All three initial conditions](../../results/spring_repair/work_midpoint/plots/work_midpoint_spring_h60.png)

## Reproduction

~~~sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
./build/xpbd-spring-midpoint-study results/spring_repair/midpoint
./build/xpbd-spring-work-midpoint-study results/spring_repair/work_midpoint
./build/xpbd-spring-controls-study results/spring_repair/existing_controls.csv results/spring_repair/existing_controls_targeted.csv
python3 experiments/spring_repair/plot_midpoint.py results/spring_repair/midpoint
python3 experiments/spring_repair/plot_work_midpoint.py results/spring_repair/midpoint results/spring_repair/work_midpoint
ctest --test-dir build --output-on-failure
~~~

All 15 repository CTests pass. The work-midpoint probe also passes ASan and
UBSan (LeakSanitizer disabled for this environment).

Plots require Matplotlib. The midpoint study validates analytical baseline and
midpoint steps, exact harmonic phase advance, conservation, stationarity, and
RK4 reference refinement. Its study also runs as a CTest. The controls check
that relaxation one reproduces the existing GS equations. Repository-release
and staging builds produce identical CSV files.

[Midpoint summaries](../../results/spring_repair/midpoint/summary.csv) ·
[Midpoint traces](../../results/spring_repair/midpoint/traces.csv) ·
[Numerical checks](../../results/spring_repair/midpoint/checks.csv) ·
[Solver controls](../../results/spring_repair/existing_controls.csv) ·
[Long-run/refinement controls](../../results/spring_repair/existing_controls_targeted.csv)

[Work-midpoint summaries](../../results/spring_repair/work_midpoint/summary.csv) ·
[Work-midpoint traces](../../results/spring_repair/work_midpoint/traces.csv) ·
[Work-midpoint checks](../../results/spring_repair/work_midpoint/checks.csv)
