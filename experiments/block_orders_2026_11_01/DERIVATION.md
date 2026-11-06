# Two-order block Gauss–Seidel experiment

This experiment implements the position-only energy-relaxed formulation agreed
in the conversation. The two methods differ only in the order of their material
and energy blocks. The primary policy recomputes analytical residuals and
gradients after **every block**, with accumulated multipliers inside the timestep.
The notation follows
[notation.md](/home/blacknoir/Projects/constrained-xpbd/notation.md).

This is a CPD-inspired energy-relaxation experiment. Original CPD uses PD local
projections and a coupled global constrained solve, including momentum rows. Here
elasticity is substituted into material multipliers after a Gauss–Newton
approximation and the resulting equations are visited in overlapping blocks.
Energy closure, iteration convergence, and accurate motion are separate results.

## Objective and energy targets

The constrained objective is

\[
\begin{aligned}
\min_{x,\rho}\quad&
\frac{\|x-x_n-hv_n\|_M^2}{2h^2}
+U(x)+\frac{\epsilon}{2}\rho^2,
\\
\text{subject to}\quad&
\frac{\|x-x_n\|_M^2}{2h^2}+U(x)
=(1-\rho)E_\star+\rho K_\star.
\end{aligned}
\]

There is no independent endpoint velocity, velocity compliance, or velocity
coupling multiplier. After accepting a step, return

\[
v_{n+1}=\frac{x_{n+1}-x_n}{h}.
\]

The incoming velocity \(v_n\) is the previously accepted velocity. It is fixed
during the timestep and need not be inferred from a second old position.

The target is recomputed from the accepted incoming state on every physical
step and then held fixed during all solver visits:

\[
E_\star=\frac12v_n^TMv_n+U(x_n).
\]

The total potential contains both spring energy and conservative gravity. Gravity
is not also added to the predictor or as work in the energy target. This target
choice allows cumulative physical energy drift when relaxation or finite
iterations change a step's output energy; report both incoming-step energy
closure and drift from the run-initial energy.

The reference \(K_\star\) is the minimum potential compatible with the anchored
scene, in the **same additive potential gauge** as \(U\). This is a deliberate
anchored-system adaptation, not original CPD's prescribed-momentum construction.
For one planar particle, gravity \(g\), and mass \(m\):

| Scene | Potential \(U(x)\) | Reference \(K_\star\) |
|---|---|---|
| Radial spring without gravity | \(\tfrac{k}{2}(\|x\|-\ell)^2\) | \(0\) |
| Radial spring with gravity | \(\tfrac{k}{2}(\|x\|-\ell)^2-mg^Tx\) | \(-m\|g\|\ell-(m\|g\|)^2/(2k)\) |
| Hard rod with gravity | \(-mg^Tx\), with \(\|x\|=\ell\) | \(-m\|g\|\ell\) |

The multiplier \(\eta\) is dimensionless. Both \(\epsilon>0\) and the effective
energy compliance \((E_\star-K_\star)^2/\epsilon\) have energy units. The parameter
\(\epsilon\) is fixed throughout each run. Because the energy target changes
between steps, the effective compliance need not remain fixed.

The relaxation variable \(\rho\) is unbounded. There is no implicit
\([0,1]\) restriction and no clipping of \(\rho\), \(\eta\), or \(\lambda\).

## Stationarity before spring substitution

Introduce only the energy multiplier initially:

\[
\begin{aligned}
\mathcal L={}&
\frac{\|x-x_n-hv_n\|_M^2}{2h^2}+U(x)
+\frac\epsilon2\rho^2
\\
&+\eta\left[
\frac{\|x-x_n\|_M^2}{2h^2}+U(x)
-(1-\rho)E_\star-\rho K_\star
\right].
\end{aligned}
\]

The stationarity residuals are

\[
\begin{aligned}
r_x={}&\frac{M(x-x_n-hv_n)}{h^2}+\nabla U(x)
+\eta\left[\frac{M(x-x_n)}{h^2}+\nabla U(x)\right],
\\
r_\rho={}&\epsilon\rho+\eta(E_\star-K_\star),
\\
r_\eta={}&\frac{\|x-x_n\|_M^2}{2h^2}+U(x)
-(1-\rho)E_\star-\rho K_\star.
\end{aligned}
\]

The exact position Hessian before substitution is

\[
(1+\eta)\left[\frac{M}{h^2}+\nabla^2U(x)\right].
\]

For a compliant spring, write \(U=U_{\mathrm{ext}}+C^2/(2\alpha)\), with
\(\alpha=1/k>0\). Drop energy-multiplier curvature and the spring's geometric
curvature:

\[
\eta\left[\frac{M}{h^2}+\nabla^2U(x)\right]\longmapsto0,
\qquad
\nabla^2\frac{C(x)^2}{2\alpha}
\approx\frac1\alpha\nabla C(x)\nabla C(x)^T.
\]

Constant gravity has zero Hessian. No energy residual, physical energy gradient,
or primal stationarity residual is discarded.

## Material substitution and approximate system

Use the chosen unscaled XPBD sign convention:

\[
\lambda=-\frac{C(x)}\alpha,
\qquad r_\lambda=C(x)+\alpha\lambda.
\]

After the approximation, substitute the complete linearized coefficient:

\[
-\frac{C(x)+\nabla C(x)^T\Delta x}{\alpha}
=\lambda+\Delta\lambda.
\]

This gives the constitutive row

\[
\nabla C(x)^T\Delta x+\alpha\Delta\lambda=-r_\lambda.
\]

The resulting system for one material constraint is

\[
\boxed{
\begin{bmatrix}
\dfrac M{h^2}&-\nabla C(x)&0&
\dfrac{M(x-x_n)}{h^2}+\nabla U(x)
\\[5pt]
\nabla C(x)^T&\alpha&0&0
\\[5pt]
0&0&\epsilon&E_\star-K_\star
\\[5pt]
\left[\dfrac{M(x-x_n)}{h^2}+\nabla U(x)\right]^T&0&E_\star-K_\star&0
\end{bmatrix}
\begin{bmatrix}\Delta x\\\Delta\lambda\\\Delta\rho\\\Delta\eta\end{bmatrix}
=-\begin{bmatrix}r_x\\r_\lambda\\r_\rho\\r_\eta\end{bmatrix}.}
\]

After substitution, evaluate

\[
\begin{aligned}
r_x={}&\frac{M(x-x_n-hv_n)}{h^2}
+\nabla U_{\mathrm{ext}}(x)-\lambda\nabla C(x)
+\eta\left[\frac{M(x-x_n)}{h^2}+\nabla U(x)\right],
\\
r_\lambda={}&C(x)+\alpha\lambda,
\\
r_\rho={}&\epsilon\rho+\eta(E_\star-K_\star),
\\
r_\eta={}&\frac{\|x-x_n\|_M^2}{2h^2}+U(x)
-(1-\rho)E_\star-\rho K_\star.
\end{aligned}
\]

Elasticity is represented by \(-\lambda\nabla C\) only in the objective part of
\(r_x\). The energy terms retain the analytical \(U\) and \(\nabla U\). For
multiple compliant constraints, replace that objective term by
\(-\sum_j\lambda_j\nabla C_j\) and give each constraint its constitutive row.

The displayed system has a nonsymmetric sign convention because its material
row uses \(r_\lambda=C+\alpha\lambda\); multiplying that row by \(-1\) gives the
equivalent symmetric saddle matrix. This is a sign convention, not a change in
the solve.

The full derivative of these substituted residuals would also contain
\(-\lambda\nabla^2C\) and
\(\eta[M/h^2+\nabla^2U]\) in the position block. Both are omitted here. This is
an approximate Newton system rather than the exact nonlinear Newton system.

For a hard rod, set \(\alpha=0\), use
\(C(x)=\|x\|-\ell\), and include **no** elastic energy in \(U\). The constitutive
row is then the linearized hard geometry constraint. This is a hard-constraint
extension, not evaluation of \(C^2/(2\alpha)\) at zero compliance.

## Material block: solve position and material multiplier together

Hold \(\eta,\rho\) fixed. At the current state, evaluate the full \(r_x\) and
the current \(\nabla C\), then solve the position and constitutive rows:

\[
\boxed{
\Delta\lambda=
\frac{-r_\lambda+h^2\nabla C(x)^TM^{-1}r_x}
{\alpha+h^2\nabla C(x)^TM^{-1}\nabla C(x)}.
}
\]

Recover the position correction:

\[
\boxed{
\Delta x=-h^2M^{-1}r_x
+h^2M^{-1}\nabla C(x)\Delta\lambda.
}
\]

Immediately apply \(x\leftarrow x+\Delta x\) and
\(\lambda\leftarrow\lambda+\Delta\lambda\). Neither is a freshly assigned
multiplier value. If the block is used inside a multiple-constraint system,
\(r_x\) contains **all** current material multipliers and the energy multiplier.

## Energy block: solve position, relaxation, and energy multiplier together

Hold material multipliers fixed. Solve the position, relaxation, and energy rows
jointly by elimination:

\[
\boxed{
\Delta\eta=
\frac{
r_\eta
-h^2\left[\dfrac{M(x-x_n)}{h^2}+\nabla U(x)\right]^TM^{-1}r_x
-\dfrac{E_\star-K_\star}{\epsilon}r_\rho
}{
h^2\left[\dfrac{M(x-x_n)}{h^2}+\nabla U(x)\right]^T
M^{-1}\left[\dfrac{M(x-x_n)}{h^2}+\nabla U(x)\right]
+\dfrac{(E_\star-K_\star)^2}{\epsilon}
}.
}
\]

Recover

\[
\boxed{
\begin{aligned}
\Delta x={}&-h^2M^{-1}r_x
-h^2M^{-1}\left[\frac{M(x-x_n)}{h^2}+\nabla U(x)\right]\Delta\eta,
\\
\Delta\rho={}&-\frac{r_\rho+(E_\star-K_\star)\Delta\eta}{\epsilon}.
\end{aligned}}
\]

Immediately apply \(x\leftarrow x+\Delta x\),
\(\eta\leftarrow\eta+\Delta\eta\), and
\(\rho\leftarrow\rho+\Delta\rho\). Eliminating variables to calculate this
block does not turn it into three sequential scalar visits: the energy
multiplier equation already includes both primal responses.

With exact arithmetic and fixed targets, this block leaves \(r_\rho=0\), since
that residual is affine in \(\rho,\eta\). The nonlinear position and energy
residuals generally do not become exactly zero in one approximate update.

The energy denominator can vanish if the analytical endpoint energy gradient
and \(E_\star-K_\star\) both vanish. Positive \(\epsilon\) alone is not a
guarantee against this degeneracy.

## Ordering and refresh policy

| Schedule | One sweep |
|---|---|
| EL | Joint energy block, immediate update and refresh, then joint material block, immediate update and refresh |
| LE | Joint material block, immediate update and refresh, then joint energy block, immediate update and refresh |

These are overlapping blocks because both contain position. With one frozen
linearization and correctly updated linear residuals, the same visits implement
Gauss–Seidel on the reduced multiplier equations together with restoration of
the position equation. Under the primary policy, analytical gradients and all
residuals are refreshed after each primal push, so the experiment is nonlinear
overlapping block Gauss–Seidel.

Every physical timestep cold-initializes

\[
x=x_n+hv_n,\qquad\lambda=0,\qquad\eta=0,\qquad\rho=0.
\]

Multipliers accumulate over the visits and sweeps in that timestep. They are
not carried to the next physical step in this experiment. After each complete
sweep, measure all residuals from the same refreshed state. Report results for
the stated finite sweep budget, or for a stopping rule that tests both
stationarity and the material/energy rows. Testing only the energy row cannot
establish convergence of the coupled formulation.

Full, undamped increments define the primary method. There is no line search,
clamp, unreported substep, alternate branch, or fallback integrator. A failed
denominator/rank check, nonfinite state, or declared numerical failure stops the
rollout. Preserve and report the accepted prefix and the first failed step.
An iteration-budget exhaustion must be distinguished from a numerically failed
visit and from a converged step.

The radial distance gradient is undefined at zero spring/rod radius. Encountering
that state is a declared geometry failure, rather than a reason to introduce an
unreported fallback direction.

## Frozen-system ordering diagnostic

For one material constraint, elimination of position and relaxation gives a
two-multiplier system. After choosing the symmetric row signs, its diagonal
entries are

\[
\alpha+h^2\nabla C(x)^TM^{-1}\nabla C(x)
\]

and

\[
h^2\left[\frac{M(x-x_n)}{h^2}+\nabla U(x)\right]^T
M^{-1}\left[\frac{M(x-x_n)}{h^2}+\nabla U(x)\right]
+\frac{(E_\star-K_\star)^2}{\epsilon}.
\]

Its off-diagonal entries are both

\[
-h^2\nabla C(x)^TM^{-1}
\left[\frac{M(x-x_n)}{h^2}+\nabla U(x)\right].
\]

Assuming positive diagonals, the scalar GS contraction factor over a frozen
full sweep is the squared off-diagonal entry divided by the product of the
diagonals. Cauchy–Schwarz bounds it by one. Positive compliance makes the bound
strict whenever it removes the relevant collinear-gradient degeneracy.
Therefore a nonsingular positive definite frozen reduced system cannot acquire
an EL-versus-LE divergence solely from reversing those two scalar multiplier
visits. Truncated solves can differ, and nonlinear refreshes can change the
system after every visit. This check distinguishes those effects.

An independent numerical algebra check compared the two block formulas above
against direct dense solves for 250 randomized positive-mass systems; the
largest absolute difference was \(2.28\times10^{-13}\). This checks the block
elimination signs and formulas only.

## What to expect at release and near the hard-energy limit

For an unforced spring released from rest, the cold initial iterate is \(x=x_n\)
with zero multipliers. Thus \(r_x=r_\eta=r_\rho=0\), but generally
\(r_\lambda=C(x_n)\ne0\). Its first **energy** visit is exactly a no-op. The
subsequent material visit moves the particle. LE first makes the material
correction and then evaluates energy at that new position. A finite-sweep
difference is therefore expected even before any long-time instability appears.

More generally, at release from rest the position-dependent implicit objective
equals the reconstructed physical endpoint energy:

\[
\frac{\|x-x_n\|_M^2}{2h^2}+U(x).
\]

On the relaxed energy constraint, the entire objective is

\[
E_\star-\rho(E_\star-K_\star)+\frac\epsilon2\rho^2.
\]

For a compliant spring without extra hard geometry constraints, if the preferred
energy level is attainable and the endpoint energy gradient is nonzero,
stationarity selects

\[
\rho=\frac{E_\star-K_\star}{\epsilon},\qquad
\eta=-1,\qquad
\frac{\|x-x_n\|_M^2}{2h^2}+U(x)
=E_\star-\frac{(E_\star-K_\star)^2}{\epsilon}.
\]

Positions on that same energy level remain indistinguishable to the objective.
Finite relaxation can therefore remove a singular energy correction while
leaving a motion/branch-selection ambiguity. The target level also need not be
attainable; no clipping to the physical energy floor is part of this method.

For a compliant spring with finite \(\epsilon\) and nonzero
\(E_\star-K_\star\), the unchanged release state is **not** a full stationary
point when the total physical potential gradient is nonzero: position
stationarity requires \(\eta=-1\), relaxation stationarity
then requires nonzero \(\rho\), and energy at the unchanged state requires
\(\rho(E_\star-K_\star)=0\). In the hard-energy limit, the unchanged state is a
stationary feasible branch; very strong penalties can approach it and produce
nearly frozen motion.

A hard rod has the additional reaction direction; at a rod equilibrium, its
reaction can balance a nonzero radial potential gradient. The preceding
nonstationarity statement is not a reason to move such an equilibrium.

For a scalar spring of stiffness \(k\), mass \(m\), initial extension \(a\), and
zero initial velocity, the joint **first frozen linear solve** predicts

\[
\Delta x=
-\frac{ka}
{m/h^2+k+(ka)^2\epsilon/(E_\star-K_\star)^2},
\]

when \(E_\star\ne K_\star\). This explicitly tends to zero as
\(\epsilon\to\infty\). It is a diagnostic for the joint linear reference, not
a claim that either finite nonlinear EL/LE sweep equals that solve.

For the hard pendulum, only endpoint length is constrained here. Reconstructed
velocity has a radial defect even if both old and new positions lie exactly on
the rod circle:

\[
x^T\frac{x-x_n}{h}=\frac{\|x-x_n\|^2}{2h}
\quad\text{when }\|x\|=\|x_n\|=\ell.
\]

Report this defect. Adding endpoint tangency while keeping the hard reconstructed
velocity relation and both endpoint lengths would force zero displacement, so
no hidden endpoint-tangency projection is included.

## Required interpretation of results

Compare the two orders using the same state, timestep, sweep budget, compliance,
energy targets, and failure rules. Record motion error against a continuous
reference, physical energy drift, relaxed energy closure, material error,
stationarity residual, multiplier magnitudes, \(\rho\), correction sizes, and
accepted duration. For the rod, also record radial-velocity defect. An accurate
energy number with frozen, phase-shifted, or discontinuous motion is not success.

Frozen-linear-system agreement validates the algebra and signs, not nonlinear
stability. Convergence of both orders to a common residual root validates that
comparison only at the tested state. An apparent bounded rollout with a fixed
sweep budget may still have large stationarity or constraint defects.
