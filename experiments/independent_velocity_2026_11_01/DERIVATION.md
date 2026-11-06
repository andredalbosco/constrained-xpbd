# Independent endpoint velocity with refreshed block Gauss–Seidel

This experiment keeps the original implicit position objective and adds a
compliant relation between endpoint velocity and position displacement. Endpoint
velocity is an independent returned state. The primary method uses positive
primal velocity compliance, with \(\beta_v\in\{0.1,1,10\}\), energy relaxation,
and the same two block orders EL and LE. A coupled approximate-system solve is
the JOINT reference. The position-only control eliminates velocity separately;
it is not evaluation of the primal penalty at zero compliance.

## Objective, targets, and accepted state

For \(\beta_v>0\) and \(\epsilon>0\), the objective is

\[
\begin{aligned}
\min_{x,v,\rho}\quad&
\frac{\|x-x_n-hv_n\|_M^2}{2h^2}+U(x)
+\frac1{2\beta_v}\left\|v-\frac{x-x_n}{h}\right\|_M^2
+\frac\epsilon2\rho^2,
\\
\text{subject to}\quad&
\frac12v^TMv+U(x)=(1-\rho)E_\star+\rho K_\star.
\end{aligned}
\]

The energy target is incoming total physical energy,

\[
E_\star=\frac12v_n^TMv_n+U(x_n),
\]

and is fixed throughout the timestep solve. The reference \(K_\star\) is the
scene's minimum potential in the same energy gauge. This is an anchored-scene
adaptation of CPD's energy relaxation, not its original momentum-based reference.
The parameter \(\epsilon\) remains fixed throughout a run, typically selected
as a ratio times initial energy above the potential floor. The primal variable
\(\rho\) is unbounded and is not clipped.

For a radial anchored spring,

\[
C(x)=\|x\|-\ell,\quad \alpha=1/k>0,\qquad
U(x)=\frac{C(x)^2}{2\alpha}+U_{\mathrm{ext}}(x),
\qquad U_{\mathrm{ext}}(x)=-mg^Tx.
\]

The potential floor is

\[
K_\star=-m\|g\|\ell-\frac{(m\|g\|)^2}{2k}.
\]

For a hard rod, \(C(x)=\|x\|-\ell\), \(\alpha=0\), and
\(U(x)=-mg^Tx\), without a spring-energy division by zero. The constrained
potential floor is \(-m\|g\|\ell\). Gravity appears in the potential, and is
not also added to the inertial predictor or energy target as external work.

After acceptance, return both computed \(x\) and computed \(v\). **Do not
overwrite \(v\) with \((x-x_n)/h\).** That overwrite changes both the energy
and the method. The mismatch between those two velocities is an output to
measure, rather than a defect silently projected away.

Without the energy constraint, stationarity in \(v\) sets
\(v=(x-x_n)/h\); eliminating it recovers the original implicit position
objective exactly.

## Stationarity before material substitution

Introduce the energy multiplier \(\eta\). The Lagrangian is

\[
\begin{aligned}
\mathcal L={}&
\frac{\|x-x_n-hv_n\|_M^2}{2h^2}+U(x)
+\frac1{2\beta_v}\left\|v-\frac{x-x_n}{h}\right\|_M^2
+\frac\epsilon2\rho^2
\\
&+\eta\left[\frac12v^TMv+U(x)
-(1-\rho)E_\star-\rho K_\star\right].
\end{aligned}
\]

The residuals before material substitution are

\[
\begin{aligned}
r_x={}&\frac{M(x-x_n-hv_n)}{h^2}+\nabla U(x)
-\frac{M}{\beta_vh}\left(v-\frac{x-x_n}{h}\right)
+\eta\nabla U(x),
\\
r_v={}&\frac{M}{\beta_v}\left(v-\frac{x-x_n}{h}\right)+\eta Mv,
\\
r_\rho={}&\epsilon\rho+\eta(E_\star-K_\star),
\\
r_\eta={}&\frac12v^TMv+U(x)-(1-\rho)E_\star-\rho K_\star.
\end{aligned}
\]

The exact position–velocity Hessian at this stage is

\[
\begin{bmatrix}
\dfrac{1+1/\beta_v}{h^2}M+(1+\eta)\nabla^2U(x)
&-\dfrac{M}{\beta_vh}
\\
-\dfrac{M}{\beta_vh}&\dfrac{M}{\beta_v}+\eta M
\end{bmatrix}.
\]

The \(\rho\rho\) entry is \(\epsilon\), with no primal cross term involving
\(\rho\). The energy equality linearization is

\[
\nabla U(x)^T\Delta x+(Mv)^T\Delta v
+(E_\star-K_\star)\Delta\rho=-r_\eta.
\]

## Curvature approximation and material substitution

Drop energy-multiplier curvature \(\eta\nabla^2U\) and \(\eta M\). Approximate
the objective's spring Hessian by

\[
\nabla^2\frac{C(x)^2}{2\alpha}
\approx\frac1\alpha\nabla C(x)\nabla C(x)^T.
\]

The omitted objective spring term is
\(C(x)\nabla^2C(x)/\alpha\). The gravity Hessian is zero. Keep **all** velocity
coupling curvature, including its position–velocity cross blocks.

Use the chosen unscaled XPBD material convention:

\[
\lambda=-\frac{C(x)}\alpha,\qquad r_\lambda=C(x)+\alpha\lambda.
\]

Substitute the complete linearized objective spring coefficient:

\[
-\frac{C(x)+\nabla C(x)^T\Delta x}{\alpha}
=\lambda+\Delta\lambda,
\qquad
\nabla C(x)^T\Delta x+\alpha\Delta\lambda=-r_\lambda.
\]

After this substitution, evaluate

\[
\boxed{
\begin{aligned}
r_x={}&\frac{M(x-x_n-hv_n)}{h^2}
+\nabla U_{\mathrm{ext}}(x)-\lambda\nabla C(x)
-\frac{M}{\beta_vh}\left(v-\frac{x-x_n}{h}\right)
+\eta\nabla U(x),
\\
r_v={}&\frac{M}{\beta_v}\left(v-\frac{x-x_n}{h}\right)+\eta Mv,
\\
r_\lambda={}&C(x)+\alpha\lambda,
\\
r_\rho={}&\epsilon\rho+\eta(E_\star-K_\star),
\\
r_\eta={}&\frac12v^TMv+U(x)-(1-\rho)E_\star-\rho K_\star.
\end{aligned}}
\]

For gravity, \(\nabla U_{\mathrm{ext}}=-mg\) remains outside the energy
multiplier term. Elasticity in the energy row and in \(\eta\nabla U\) is
evaluated analytically from position, rather than replaced by the independently
carried material multiplier. A hard rod uses the same constitutive row at
\(\alpha=0\) with no elastic potential.

The approximate system is

\[
\boxed{
\begin{bmatrix}
\dfrac{1+1/\beta_v}{h^2}M&-\dfrac{M}{\beta_vh}&-\nabla C(x)&0&\nabla U(x)
\\[5pt]
-\dfrac{M}{\beta_vh}&\dfrac{M}{\beta_v}&0&0&Mv
\\[5pt]
\nabla C(x)^T&0&\alpha&0&0
\\[5pt]
0&0&0&\epsilon&E_\star-K_\star
\\[5pt]
\nabla U(x)^T&(Mv)^T&0&E_\star-K_\star&0
\end{bmatrix}
\begin{bmatrix}\Delta x\\\Delta v\\\Delta\lambda\\\Delta\rho\\\Delta\eta\end{bmatrix}
=-\begin{bmatrix}r_x\\r_v\\r_\lambda\\r_\rho\\r_\eta\end{bmatrix}.}
\]

The material row can be multiplied by \(-1\) to obtain the equivalent symmetric
saddle sign convention. Both \(r_x\) and \(r_v\) are retained in every block.
After substitution the omitted material curvature is
\(-\lambda\nabla^2C\); the omitted energy curvature remains
\(\eta\nabla^2U\) and \(\eta M\). These omissions do not replace any
analytical residual or gradient.

## Exact inverse of the retained primal block

The retained position–velocity block and its inverse are

\[
A_0=
\begin{bmatrix}
\dfrac{1+1/\beta_v}{h^2}M&-\dfrac{M}{\beta_vh}
\\
-\dfrac{M}{\beta_vh}&\dfrac{M}{\beta_v}
\end{bmatrix},
\qquad
\boxed{
A_0^{-1}=
\begin{bmatrix}
h^2M^{-1}&hM^{-1}
\\
hM^{-1}&(1+\beta_v)M^{-1}
\end{bmatrix}.}
\]

For positive mass and \(\beta_v>0\), the block is positive definite because its
quadratic form is

\[
\left\|\frac{\Delta x}{h}\right\|_M^2
+\frac1{\beta_v}\left\|\Delta v-\frac{\Delta x}{h}\right\|_M^2.
\]

This positivity concerns the approximation, not the exact multiplier-dependent
Lagrangian Hessian or the stability of a nonlinear rollout.

## Material visit: coupled position, velocity, and material multiplier

Hold \(\eta,\rho\) fixed. The material increment is

\[
\boxed{
\Delta\lambda=
\frac{
-r_\lambda+\nabla C(x)^T\left[h^2M^{-1}r_x+hM^{-1}r_v\right]
}{\alpha+h^2\nabla C(x)^TM^{-1}\nabla C(x)}.}
\]

Recover the coupled primal increments:

\[
\boxed{
\begin{aligned}
\Delta x={}&-h^2M^{-1}r_x-hM^{-1}r_v
+h^2M^{-1}\nabla C(x)\Delta\lambda,
\\
\Delta v={}&-hM^{-1}r_x-(1+\beta_v)M^{-1}r_v
+hM^{-1}\nabla C(x)\Delta\lambda.
\end{aligned}}
\]

The material Schur denominator is independent of \(\beta_v\), but the residuals
and velocity response are not. In particular, a material multiplier change
immediately changes both position and independent velocity.

## Energy visit: coupled position, velocity, relaxation, and energy multiplier

Hold \(\lambda\) fixed. Eliminating the primal increments gives

\[
\boxed{
\Delta\eta=
\frac{
\begin{aligned}
&r_\eta
-\nabla U(x)^T\left[h^2M^{-1}r_x+hM^{-1}r_v\right]
\\[-2pt]
&\quad-(Mv)^T\left[hM^{-1}r_x+(1+\beta_v)M^{-1}r_v\right]
-\dfrac{E_\star-K_\star}{\epsilon}r_\rho
\end{aligned}
}{
h^2\left[\nabla U(x)+\dfrac{Mv}{h}\right]^TM^{-1}
\left[\nabla U(x)+\dfrac{Mv}{h}\right]
+\beta_vv^TMv+\dfrac{(E_\star-K_\star)^2}{\epsilon}
}.}
\]

The sum-of-squares form of the denominator is useful computationally; expanding
its cross terms and subtracting near-equal values would introduce unnecessary
roundoff cancellation. Recover

\[
\boxed{
\begin{aligned}
\Delta x={}&-h^2M^{-1}r_x-hM^{-1}r_v
-\left[h^2M^{-1}\nabla U(x)+hv\right]\Delta\eta,
\\
\Delta v={}&-hM^{-1}r_x-(1+\beta_v)M^{-1}r_v
-\left[hM^{-1}\nabla U(x)+(1+\beta_v)v\right]\Delta\eta,
\\
\Delta\rho={}&-\frac{r_\rho+(E_\star-K_\star)\Delta\eta}{\epsilon}.
\end{aligned}}
\]

This is a joint energy block even though elimination produces one scalar solve.
The equation for \(\Delta\eta\) already includes the position, velocity, and
relaxation responses.

The extra positive term \(\beta_vv^TMv\) prevents a vanished combined position
energy gradient from being the sole source of a singular denominator when
\(v\ne0\). It provides no such protection at rest. If velocity, potential
gradient, and the target difference all vanish, the row can still be rank
deficient. A compatible zero row can be skipped consistently; an incompatible
zero row must be reported as a failure.

## Iteration contract

At each physical timestep initialize

\[
x=x_n+hv_n,\qquad v=v_n,\qquad\lambda=\eta=\rho=0.
\]

Hold \(x_n,v_n,E_\star,K_\star,h,\epsilon,\beta_v\) fixed during that solve.
The two primary schedules are:

| Schedule | One sweep |
|---|---|
| EL | Solve the complete energy block, immediately update \(x,v,\rho,\eta\), refresh all analytical residuals and gradients, then solve the complete material block and immediately update \(x,v,\lambda\) |
| LE | Solve the complete material block, immediately update \(x,v,\lambda\), refresh, then solve the complete energy block and immediately update \(x,v,\rho,\eta\) |

Refresh all residuals and gradients after the final block too. Accumulate the
multiplier increments within a timestep. These are overlapping nonlinear block
Gauss–Seidel visits, since both blocks contain position and velocity. The JOINT
reference solves the displayed full approximate system and then refreshes.

Use the stated finite sweep budget. A converged flag tests both primal
stationarity residuals, the material row, relaxation stationarity, and energy
closure; a finite output with an exhausted budget is not automatically
converged. There is no hidden velocity overwrite, energy rescaling, damping,
line search, clamping, substepping, or fallback branch. Undefined distance
normals and nonfinite/runaway states are explicit failures with an accepted
prefix, rather than reasons to repair the method silently.

The primary pendulum comparison imposes length only. It does **not** include
endpoint tangency. Report both velocity mismatch and radial velocity, together
with actual physical energy, accepted position motion, stalls, and jumps.
These matter even when energy closure is excellent.

## Relation to the separately eliminated position-only control

The primal objective is undefined at \(\beta_v=0\); the matrix entries also
diverge. Do not pass zero as its parameter. The formal inverse limit is

\[
\lim_{\beta_v\to0^+}A_0^{-1}
=\begin{bmatrix}h^2M^{-1}&hM^{-1}\\hM^{-1}&M^{-1}\end{bmatrix}.
\]

It is rank deficient in the independent position–velocity coordinates and
generates \(\Delta v=\Delta x/h\) for compatible hard-coupled increments.
The combined stationarity residual is

\[
r_x+\frac{r_v}{h}
=\frac{M(x-x_n-hv_n)}{h^2}
+\nabla U_{\mathrm{ext}}-\lambda\nabla C
+\eta\left[\nabla U+\frac{Mv}{h}\right].
\]

Substituting the hard relation \(v=(x-x_n)/h\) recovers the position-only
residual and energy gradient. Likewise, the energy denominator becomes its
position-only expression when the \(\beta_vv^TMv\) term disappears. This
explains the relation between formulations without claiming that a finite
independent-variable algorithm is defined by setting its primal penalty to zero.

## Remaining nonlinear and kinematic failure mechanisms

The positive approximate matrix is not a guarantee of benign motion. The
physical energy gradient changes after both position and velocity updates, and
the omitted \(\eta\nabla^2U\) and \(\eta M\) can become large. A global energy
constraint also leaves multiple position/velocity distributions with the same
energy, so the solver may choose an unintuitive branch.

At a stationary point with no optional velocity constraint, the velocity row
implies the exact relation

\[
\boxed{
\frac{x-x_n}{h}=(1+\beta_v\eta)v.}
\]

Thus an energy multiplier can materially alter how returned velocity relates
to visible displacement. Near \(1+\beta_v\eta=0\), small displacement can
coexist with large independent velocity; when that factor is negative, their
directions are opposite. Monitor this factor and the actual mismatch. Dividing
by it to eliminate velocity would reintroduce a singular expression and is not
part of the algorithm.

The exact velocity Hessian is \((1/\beta_v+\eta)M\), so precisely that same
factor can indicate a loss of positive curvature omitted from the approximate
solve. This is a mechanism to investigate, not a new output-rejection rule.
The primary success criterion remains approximately conserved physical energy
and bounded, intuitive motion.

For an unforced spring released from rest, the cold first energy visit remains
a no-op: both primal residuals and the energy residual are zero while the
material residual is not. The first material visit creates both position and
velocity changes. Independent velocity therefore does not alone eliminate
every release degeneracy or the hard-energy stationary freeze branch.

## Optional endpoint-tangency extension, disabled in the primary comparison

For a rod anchored at the origin one can add \(x^Tv=0\), with the separately
named multiplier \(\lambda_t\) and Lagrangian term \(+\lambda_t x^Tv\).
Its stationarity contributions are \(+\lambda_t v\) in \(r_x\) and
\(+\lambda_t x\) in \(r_v\). If its multiplier curvature is also dropped,
its linearized row is

\[
v^T\Delta x+x^T\Delta v=-x^Tv.
\]

A standalone tangency block, holding other multipliers fixed, would use

\[
\Delta\lambda_t=
\frac{
x^Tv-v^T\left[h^2M^{-1}r_x+hM^{-1}r_v\right]
-x^T\left[hM^{-1}r_x+(1+\beta_v)M^{-1}r_v\right]
}{
h^2v^TM^{-1}v+2hv^TM^{-1}x+(1+\beta_v)x^TM^{-1}x
}.
\]

Its primal corrections are the ordinary residual responses above minus

\[
\left[h^2M^{-1}v+hM^{-1}x\right]\Delta\lambda_t
\quad\text{and}\quad
\left[hM^{-1}v+(1+\beta_v)M^{-1}x\right]\Delta\lambda_t
\]

in position and velocity, respectively. The denominator is positive for
\(\beta_v>0\) and nonzero rod position. Its schedule or binding with geometry
and energy must be specified as a **different experiment**. Adding it silently
to the primary method would confound the velocity-disentangling comparison.

### Optional joint rod-length and tangency block

If tangency is enabled, a stronger local choice is to bind it to the geometry
visit. Solve \((\Delta x,\Delta v,\Delta\lambda,\Delta\lambda_t)\) together,
holding \(\eta,\rho\) fixed. Here \(\lambda\) is the length multiplier with
the negative material sign; \(\lambda_t\) has the positive tangency sign.
Evaluate the full residuals including

\[
r_x\mathrel{+}=\lambda_t v,\qquad
r_v\mathrel{+}=\lambda_t x.
\]

Drop the multiplier-dependent tangency cross curvature \(\lambda_t\mathrm{Id}\)
alongside the previously dropped material and energy curvature. The joint local
system is

\[
\boxed{
\begin{bmatrix}
\dfrac{1+1/\beta_v}{h^2}M&-\dfrac{M}{\beta_vh}&-\nabla C(x)&v
\\[4pt]
-\dfrac{M}{\beta_vh}&\dfrac{M}{\beta_v}&0&x
\\[4pt]
\nabla C(x)^T&0&\alpha&0
\\[4pt]
v^T&x^T&0&0
\end{bmatrix}
\begin{bmatrix}\Delta x\\\Delta v\\\Delta\lambda\\\Delta\lambda_t\end{bmatrix}
=-\begin{bmatrix}r_x\\r_v\\r_\lambda\\x^Tv\end{bmatrix}.}
\]

The rod sets \(\alpha=0\). Keeping \(\alpha\) visible in this display also
shows the constitutive sign. After eliminating \(\Delta x,\Delta v\) and
flipping the tangency equation's Schur row, the symmetric two-multiplier system
is

\[
\boxed{
\begin{bmatrix}
\alpha+h^2\nabla C(x)^TM^{-1}\nabla C(x)
&-\nabla C(x)^T\left[h^2M^{-1}v+hM^{-1}x\right]
\\[5pt]
-\nabla C(x)^T\left[h^2M^{-1}v+hM^{-1}x\right]
&h^2v^TM^{-1}v+2hv^TM^{-1}x+(1+\beta_v)x^TM^{-1}x
\end{bmatrix}
\begin{bmatrix}\Delta\lambda\\\Delta\lambda_t\end{bmatrix}
=
\begin{bmatrix}
-r_\lambda+\nabla C(x)^T\left[h^2M^{-1}r_x+hM^{-1}r_v\right]
\\[5pt]
x^Tv-v^T\left[h^2M^{-1}r_x+hM^{-1}r_v\right]
-x^T\left[hM^{-1}r_x+(1+\beta_v)M^{-1}r_v\right]
\end{bmatrix}.}
\]

The primal corrections are

\[
\boxed{
\begin{aligned}
\Delta x={}&-h^2M^{-1}r_x-hM^{-1}r_v
+h^2M^{-1}\nabla C(x)\Delta\lambda
-\left[h^2M^{-1}v+hM^{-1}x\right]\Delta\lambda_t,
\\
\Delta v={}&-hM^{-1}r_x-(1+\beta_v)M^{-1}r_v
+hM^{-1}\nabla C(x)\Delta\lambda
-\left[hM^{-1}v+(1+\beta_v)M^{-1}x\right]\Delta\lambda_t.
\end{aligned}}
\]

The symmetric Schur matrix is positive definite for positive mass,
\(\beta_v>0\), and nonzero rod position. Its two constraint-gradient columns
are independent: the length column has zero velocity part, while tangency has
the nonzero velocity part \(x\). For a single isotropic mass, \(\alpha=0\),
and the distance constraint, its determinant simplifies to

\[
\frac{h^2}{m^2}\left\{
h^2\left[\|v\|^2-\bigl(\nabla C(x)^Tv\bigr)^2\right]
+\beta_v\|x\|^2\right\}>0.
\]

Immediately accumulate both multiplier increments and push both primal
variables. Refresh geometry, tangency, energy, and both stationarity residuals
after the visit. The energy block then holds both \(\lambda\) and
\(\lambda_t\) fixed and uses the same energy formulas above with the augmented
full \(r_x,r_v\). It may disturb geometry and tangency, which the next local
block addresses. A finite local update satisfies the two **linearized** rows;
the nonlinear rows still require iteration.

The optional orders are therefore energy then joint length/tangency, or joint
length/tangency then energy. This is SHAKE/RATTLE-inspired constraint treatment
inside the current objective, rather than the original SHAKE or RATTLE time
integrator. It remains disabled in the primary comparison.

## Independent algebra and implementation audit

The stated inverse and block formulas were checked against direct dense block
solves for 375 randomized positive-mass systems, using \(\beta_v=0.1,1,10\),
both hard and compliant material rows, and seed 773641. These are algebraic
checks, not physical trajectory experiments. Maximum errors were:

| Check | Maximum error |
|---|---:|
| \(A_0A_0^{-1}-\mathrm{Id}\), absolute | 6.37e-13 |
| Material-block solution, relative to direct solve | 2.00e-14 |
| Energy-block solution, relative to direct solve | 4.11e-15 |
| Positive-sum energy denominator, relative | 4.32e-16 |

Read-only inspection of `independent_solver.py` confirms these signs, the full
coupling inverse in both blocks, analytical energy gradients, the sum-of-squares
denominator, cold initialization, accumulated updates, and preservation of
independent returned velocity. This checks implementation consistency with the
derivation, rather than establishing stability or visual quality.
