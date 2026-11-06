# Original CPD equations and declared anchored restrictions

This derivation follows the primary CPD paper's Eqs. (7), (11)–(21) and
Algorithm 1. The current experiment distinguishes the original free-body core,
central-anchor restrictions, finite-stiffness pendulums, and an exact hard-rod
extension. Source boundaries and locally selected settings are recorded in
[SOURCE_FIDELITY.md](SOURCE_FIDELITY.md).

There are no independently optimized endpoint velocities and no scalar material
multiplier Gauss–Seidel visits in the original core. The multiplier eta retains
its plus-sign energy convention. The paper's relaxation variable is called rho
here, reserving alpha for material compliance in the wider project notation.

## PD objective and local/global representation

For free particles of positive masses, use the source predictor

\[
y_i=x_{n,i}+hv_{n,i}+\frac{h^2}{m_i}f_{\mathrm{ext},i}.
\]

Let \(B_jx=x_{a_j}-x_{b_j}\) denote a spring difference; it is an incidence
operator, not a renamed gradient. With positive spring weights \(k_j\),

\[
U_s(x)=\sum_j\frac{k_j}{2}\left(\|B_jx\|-\ell_j\right)^2.
\]

The objective, before energy/momentum constraints, is

\[
\frac1{2h^2}\|x-y\|_M^2+U_s(x).
\]

At each current iterate, the PD local projection is

\[
p_j=\ell_j\frac{B_jx}{\|B_jx\|}.
\]

Its frozen quadratic spring model is
\(k_j\|B_jx-p_j\|^2/2\), giving the current gradient and constant PD matrix

\[
\nabla_x\left[\frac1{2h^2}\|x-y\|_M^2+U_s(x)\right]
=\frac{M(x-y)}{h^2}+\sum_j k_j B_j^T(B_jx-p_j),
\]

\[
\boxed{A_{\mathrm{PD}}=\frac{M}{h^2}+\sum_jk_jB_j^TB_j.}
\]

This uses the full PD element identity metric. It does not replace a radial
spring Hessian by a rank-one Gauss–Newton matrix. At a freshly projected
current iterate, the displayed spring gradient equals the analytical gradient
of \(U_s\). Refresh the local projections before evaluating the next iterate's
gradient or physical energy.

For a two-particle spring in the plane, the matrix is explicitly

\[
A_{\mathrm{PD}}=
\begin{bmatrix}
(m_a/h^2+k)\mathrm{Id}_2&-k\mathrm{Id}_2\\
-k\mathrm{Id}_2&(m_b/h^2+k)\mathrm{Id}_2
\end{bmatrix}.
\]

## Source momentum and energy targets

Choose a fixed origin o during the timestep. For the free-body comparison,
take o to be the incoming center of mass and transform all angular quantities
consistently. Incoming targets use actual incoming velocities, not an endpoint
guess or a reconstructed velocity from some other trajectory.

\[
\begin{aligned}
P_\star={}&\sum_i m_i v_{n,i}+h\sum_i f_{\mathrm{ext},i},\\
L_{\mathrm{ang},\star}={}&
\sum_i(x_{n,i}-o)\times m_i v_{n,i}
+h\sum_i(x_{n,i}-o)\times f_{\mathrm{ext},i},\\
E_\star={}&\frac12\sum_i m_i\|v_{n,i}\|^2+U_s(x_n)
+h\sum_i f_{\mathrm{ext},i}^Tv_{n,i}.
\end{aligned}
\]

The angular equation uses the declared correction of the printed Eq. (12)
torque order: position cross force, consistent with position cross momentum.
Gravity is an external force in this source convention. Therefore \(E_\star\)
is a target for internal kinetic-plus-elastic energy, **not** kinetic plus
elastic plus gravitational potential.

In planar motion, angular momentum and inertia are scalars. With incoming COM
origin,

\[
I_n=\sum_i m_i\|x_{n,i}-o\|^2,\qquad
\boxed{K_\star=\frac{\|P_\star\|^2}{2\sum_i m_i}
+\frac{L_{\mathrm{ang},\star}^2}{2I_n}.}
\]

For the 3D generalization,

\[
I_n=\sum_i m_i\left[
\|x_{n,i}-o\|^2\mathrm{Id}_3-(x_{n,i}-o)(x_{n,i}-o)^T\right],
\]

and the second term is \(L_{\mathrm{ang},\star}^TI_n^+
L_{\mathrm{ang},\star}/2\). A pseudoinverse does not make an out-of-range
angular target feasible; report rank and range compatibility. For the selected
planar finite-separation scenes, the rotational inertia is a positive scalar,
so a scalar division suffices. This inertia calculation is distinct from the
rank treatment of the conservation Schur system.

The momentum-derived kinetic lower bound alone does not prove that the total
elastic-energy target is feasible. Finite epsilon also does not guarantee rho
is zero whenever strict conservation is feasible.

The planar inertia is a mass/geometric quantity, not an entry of the PD
stiffness matrix. At the incoming COM origin, the momentum-gradient Gram matrix
in the inverse-mass metric has the blocks

\[
\nabla P^TM^{-1}\nabla P=\frac{\sum_i m_i}{h^2}\mathrm{Id}_2,
\qquad
\nabla P^TM^{-1}\nabla L_{\mathrm{ang}}=0,
\qquad
\nabla L_{\mathrm{ang}}^TM^{-1}\nabla L_{\mathrm{ang}}=\frac{I_n}{h^2}.
\]

These identities explain the translation/rotation split in the compatible
kinetic reference. They do not assert that the corresponding blocks of the
PD-metric conservation Schur matrix are equal to these mass-metric blocks.

## Original constrained objective

With backward velocity reconstruction, define momentum constraints explicitly:

\[
\begin{aligned}
P(x)&=\frac1h\sum_i m_i(x_i-x_{n,i}),\\
L_{\mathrm{ang}}(x)&=\frac1h\sum_i(x_i-o)\times m_i(x_i-x_{n,i}).
\end{aligned}
\]

The original CPD problem is

\[
\boxed{
\begin{aligned}
\min_{x,\rho}\quad&
\frac1{2h^2}\|x-y\|_M^2+U_s(x)+\frac\epsilon2\rho^2,\\
\text{subject to}\quad& P(x)=P_\star,\\
&L_{\mathrm{ang}}(x)=L_{\mathrm{ang},\star},\\
&\frac1{2h^2}\|x-x_n\|_M^2+U_s(x)
=(1-\rho)E_\star+\rho K_\star.
\end{aligned}}
\]

The physical energy in this equality is always the analytical current spring
energy. Using a previous frozen projection inside its evaluation after moving
position would define a surrogate-energy variant.

Stack \(q=(x,\rho)\), and collect the left sides minus targets as
\(\mathbf c(q)\). In the plane it has four entries: two linear momentum,
one angular momentum, and the displayed relaxed energy row. Gradients are
columns in \(\nabla\mathbf c\), as in the project convention.

## Analytical constraint gradients

For Cartesian unit vector \(\mathbf u_a\), a linear-momentum row has

\[
\nabla_{x_i}\left[\mathbf u_a^TP(x)\right]
=\frac{m_i}{h}\mathbf u_a,
\qquad \nabla_\rho P=0.
\]

The angular constraint is affine in endpoint position despite its initial
appearance. The identity

\[
(x_i-o)\times(x_i-x_{n,i})
=(x_{n,i}-o)\times(x_i-x_{n,i})
\]

shows that its derivative is fixed within a timestep. In the plane,

\[
\boxed{
\nabla_{x_i}L_{\mathrm{ang}}(x)
=\frac{m_i}{h}
\begin{bmatrix}-(x_{n,i}-o)_y\\(x_{n,i}-o)_x\end{bmatrix},
\qquad\nabla_\rho L_{\mathrm{ang}}=0.}
\]

For a 3D angular component selected by \(\mathbf u_a\), the corresponding
gradient is \((m_i/h)\mathbf u_a\times(x_{n,i}-o)\). The origin o is held fixed
when differentiating, rather than recomputed from the trial endpoint.

For the relaxed energy row,

\[
\boxed{
\nabla_{x_i}\left[
\frac{\|x-x_n\|_M^2}{2h^2}+U_s(x)
-(1-\rho)E_\star-\rho K_\star\right]
=\frac{m_i(x_i-x_{n,i})}{h^2}+\nabla_{x_i}U_s(x),}
\]

\[
\boxed{\nabla_\rho\left[
\frac{\|x-x_n\|_M^2}{2h^2}+U_s(x)
-(1-\rho)E_\star-\rho K_\star\right]=E_\star-K_\star.}
\]

These gradients are recomputed at the current trial state. Only the linear
and angular momentum gradients are reusable across the appropriate iterations.

## Exact Eq. (21) algorithm contract

Use the constant approximate primal Hessian

\[
A=\begin{bmatrix}A_{\mathrm{PD}}&0\\0&\epsilon\end{bmatrix}.
\]

This drops nonlinear objective-projection curvature and all conservation
multiplier curvature. It keeps the **objective gradient**, not the previous
iterate's stationarity residual. At every global iteration compute

\[
\nabla_q\left[
\frac{\|x-y\|_M^2}{2h^2}+U_s(x)+\frac\epsilon2\rho^2\right]
=\begin{bmatrix}
M(x-y)/h^2+\nabla U_s(x)\\\epsilon\rho
\end{bmatrix}.
\]

Let the **fresh conservation multiplier** stack be
\(\boldsymbol\eta_{\mathrm{cons}}^{\mathrm{new}}
=(\eta_P^{\mathrm{new}},\eta_L^{\mathrm{new}},\eta^{\mathrm{new}})\), where
\(\eta_P\) is the two-component linear multiplier, \(\eta_L\) the planar
angular multiplier, and eta the energy multiplier. This does not reuse the
material lambda or the work multipliers mu and nu. The source KKT subproblem is

\[
\boxed{
\begin{bmatrix}
A&\nabla\mathbf c(q)\\
\nabla\mathbf c(q)^T&0
\end{bmatrix}
\begin{bmatrix}\Delta q\\\boldsymbol\eta_{\mathrm{cons}}^{\mathrm{new}}\end{bmatrix}
=-\begin{bmatrix}
\nabla_q\left[\dfrac{\|x-y\|_M^2}{2h^2}+U_s(x)+\dfrac\epsilon2\rho^2\right]\\
\mathbf c(q)
\end{bmatrix}.}
\]

Equivalently, form

\[
S=\nabla\mathbf c(q)^TA^{-1}\nabla\mathbf c(q),
\]

and solve

\[
\boxed{
S\boldsymbol\eta_{\mathrm{cons}}^{\mathrm{new}}
=\mathbf c(q)-\nabla\mathbf c(q)^TA^{-1}
\nabla_q\left[\frac{\|x-y\|_M^2}{2h^2}+U_s(x)+\frac\epsilon2\rho^2\right].}
\]

Then apply the full coupled update

\[
\boxed{
q\leftarrow q-A^{-1}\left\{
\nabla_q\left[\frac{\|x-y\|_M^2}{2h^2}+U_s(x)+\frac\epsilon2\rho^2\right]
+\nabla\mathbf c(q)\boldsymbol\eta_{\mathrm{cons}}^{\mathrm{new}}
\right\}.}
\]

Positions and rho are accumulated through this primal update. Conservation
multiplier values are freshly assigned from the new subproblem; they are not
incremented by the displayed freshly solved value. Every conservation row is
coupled in the global solve. There are no EL/LE visits in this core method.

Initialize \((x,\rho)=(y,0)\), fix the incoming targets during the entire
timestep, refresh local projections and analytical conservation gradients after
each update, and eventually return \(v_{n+1}=(x-x_n)/h\). The experiment's
post-update, fixed-16, and strict stationarity policies are declared separately
in SOURCE_FIDELITY.md. No damping, energy rescale, hard velocity projection,
conservation-row deletion, or hidden regularization changes Eq. (21).

The implementation normalizes Schur rows/columns symmetrically and attempts a
direct solve of all active rows. An eigenvalue rank estimate at relative 1e-12
is diagnostic only, not an automatic rejection threshold. Actual singular
factorization, a nonfinite solution, or an explicit geometry/runaway guard is
reported as failure. No least-squares, pseudoinverse, row deletion, or diagonal
repair changes the subproblem.

For endpoint diagnostics, use the most recently solved multipliers with the
fresh endpoint gradients:

\[
r_x=\frac{M(x-y)}{h^2}+\nabla U_s(x)
+\nabla_x\mathbf c(q)\boldsymbol\eta_{\mathrm{cons}}^{\mathrm{new}},
\qquad
r_\rho=\epsilon\rho+(E_\star-K_\star)\eta.
\]

A small conservation residual can coexist with nonzero endpoint stationarity;
the published constraint-only stopping condition does not rule this out.

## Mechanical-energy drift decomposition for prescribed constant loads

For constant external forces with potential
\(-\sum_i f_{\mathrm{ext},i}^Tx_i\), actual mechanical energy is

\[
\frac1{2h^2}\|x-x_n\|_M^2+U_s(x)
-\sum_i f_{\mathrm{ext},i}^Tx_i.
\]

Define the named energy-row residual without hiding its expression:

\[
r_\eta=
\frac1{2h^2}\|x-x_n\|_M^2+U_s(x)
-(1-\rho)E_\star-\rho K_\star.
\]

Substituting the source target gives the exact bookkeeping identity

\[
\boxed{
\begin{aligned}
&\left[\frac1{2h^2}\|x-x_n\|_M^2+U_s(x)
-\sum_i f_{\mathrm{ext},i}^Tx_i\right]
\\
&\quad-\left[\frac12\|v_n\|_M^2+U_s(x_n)
-\sum_i f_{\mathrm{ext},i}^Tx_{n,i}\right]
\\
&=r_\eta-\rho(E_\star-K_\star)
+\sum_i f_{\mathrm{ext},i}^T\left[hv_{n,i}-(x_i-x_{n,i})\right].
\end{aligned}}
\]

The three terms distinguish solve residual, energy relaxation, and old-velocity
power versus actual displacement work. Exact source-row closure removes only
the first. For configuration-dependent attachment/collision forces, a constant
force times displacement is not automatically their exact potential work; the
present identity is used only for the selected prescribed constant-load cases.

## Free-pair relative fixtures

For two masses \(m_a=m_b=2\), initialize

\[
x_a=+r/2,\quad x_b=-r/2,\qquad
v_a=+u/2,\quad v_b=-u/2.
\]

The COM is at rest and reduced mass is one. Relative kinetic energy is
\(\|u\|^2/2\), angular momentum is \(r\times u\), and spring energy is
\(k(\|r\|-\ell)^2/2\): exactly the corresponding mass-one anchored quantities.
The original P/L/E constraints remain valid because both particles are free.

Known loads \(f_a=+g\), \(f_b=-g\) produce relative acceleration g, zero net
force, external torque \(r\times g\), and old-velocity power \(g^Tu\).
Their potential is \(-g^T(x_a-x_b)\). This is an opposite-load scene, not
uniform gravity on the pair. It avoids assuming an unknown support impulse
while retaining the earlier anchored gravity-relative problem.

The full free-pair Eq. (21) problem restricted by its zero-COM momentum target
reduces algebraically to the mass-one central-anchor angular/energy problem
below, with identical relative positions, velocities, epsilon, and iteration
policy. Their relative trajectories should therefore agree to numerical solve
accuracy. This provides an implementation parity check without deleting the
free-body linear-momentum rows.

## Central-anchor restriction and hard-rod extension

For one particle attached centrally to an origin support, omit total linear
momentum conservation. Retain angular momentum about that anchor and use

\[
L_{\mathrm{ang},\star}=x_n\times mv_n+h\,x_n\times mg.
\]

The central support has zero torque about that origin. Thus no support-force
estimate is required for this particular angular target. This does not extend
to arbitrary off-origin supports without transforming the torques correctly.

The angular-only kinetic reference is

\[
\boxed{K_\star=\frac{L_{\mathrm{ang},\star}^2}{2m\|x_n\|^2}.}
\]

It is the minimum kinetic energy for that discrete angular target, using
\(L_{\mathrm{ang}}=x_n\times mv\) under backward kinematics. Do not retain a
translational kinetic term from a linear-momentum target that was disabled.
For zero gravity this restriction conserves the physically correct anchored
angular momentum; it is more informative than the historical energy-only pin
restriction.

The finite-stiffness pendulum cases use the same radial spring with k=10,000.
They have compliant length and remain ordinary elastic PD models with the
declared anchored angular restriction. For the separately labeled exact rod,
remove the elastic spring energy, set the position metric to \(M/h^2\), and add
the global hard distance row \(\|x\|-\ell=0\), with gradient \(x/\|x\|\).
Solve it jointly with angular and relaxed energy rows. This is an extension of
the CPD solve, not a source-authored exact rod method.

The hard-row multiplier is named \(\lambda_c\), with plus-sign contribution
\(+\lambda_c(\|x\|-\ell)\). Its primal column is therefore
\(+\nabla C\), and its constraint row is \(\nabla C^T\Delta x=-C\). The
anchored hard-rod fresh stack is \((\eta_L,\lambda_c,\eta)\); no material
compliance multiplier is introduced. The full free-body core continues to use
\((\eta_P,\eta_L,\eta)\).

Nonzero outgoing angular momentum induced by gravity can prevent the trivial
zero-displacement rod release that occurs in an energy-only source-power
restriction. It does not make mechanical energy exact: the source old-power
target and backward velocity are unchanged. Exact endpoint tangency is not
added, because length at both endpoints plus chord velocity would then force
zero displacement.

## Release degeneracy remains a testable obstruction

When all individual incoming velocities and external forces are zero, the
position objective equals the constrained internal endpoint energy:

\[
\frac{\|x-x_n\|_M^2}{2h^2}+U_s(x).
\]

Strict energy conservation therefore makes the objective constant over the
feasible P/L/E manifold. Zero displacement is feasible when P and L targets
are zero. Finite energy relaxation permits a nearby nonunique energy level,
but large epsilon can select a nearly stationary branch. An initial
constraint-only precheck can return the unchanged state without an update;
the primary post-update interpretation is intentionally reported separately.

Restoring momentum rows directly addresses the energy-only method's unforced
angular-momentum gain. It does not remove this release obstruction by algebra,
nor guarantee that nonlinear turning branches or radial-mode behavior become
intuitive. The authorized experiments must assess those motions directly,
using physical energy, support/external-work accounting, actual multiplier
effects, stalls, and discrete velocity changes.

For a regular release branch with zero momentum targets and K_star=0,
stationarity of the relaxed formulation gives

\[
\eta=-1,\qquad\rho=\frac{E_\star}{\epsilon},\qquad
\frac{\|x-x_n\|_M^2}{2h^2}+U_s(x)
=E_\star-\frac{E_\star^2}{\epsilon}.
\]

The finite-energy relaxation selects a nearby energy level, while positions on
that level remain indistinguishable to the release objective. The full P/L
rows do not add a direction that makes that objective prefer a broad physical
elastic oscillation over a nearly stationary branch.

## Extension to incoming pure rigid motion

With zero external force, suppose incoming velocities are one rigid planar
field,

\[
v_{n,i}=t+\omega\begin{bmatrix}
-(x_{n,i}-o)_y\\(x_{n,i}-o)_x
\end{bmatrix}.
\]

Direct expansion of the implicit objective gives the exact identity

\[
\boxed{
\begin{aligned}
&\frac{\|x-x_n-hv_n\|_M^2}{2h^2}+U_s(x)
\\
&=\frac{\|x-x_n\|_M^2}{2h^2}+U_s(x)
-t^TP(x)-\omega L_{\mathrm{ang}}(x)+\frac12\|v_n\|_M^2.
\end{aligned}}
\]

The constraints fix every position-dependent term on the right for a given
rho. Including its penalty leaves only

\[
(1-\rho)E_\star+\rho K_\star-t^TP_\star
-\omega L_{\mathrm{ang},\star}+\frac12\|v_n\|_M^2
+\frac\epsilon2\rho^2.
\]

This broadens the release obstruction to initially rigid-rotating stretched
states: the objective cannot rank different feasible radial/rotating endpoints
with identical P/L/energy. It does not prove a unique frozen-radius result or
universal mesh freezing. Nonrigid incoming velocity, external-force terms, and
the subsequent discrete chord state can change the argument. The full proof,
anchored specialization, multiplier signs, and those limits are in
[OBJECTIVE_DEGENERACY.md](OBJECTIVE_DEGENERACY.md).
