# Why restoring angular momentum can leave radial motion undetermined

For zero external force, an exact identity explains an additional limitation of
the original constrained implicit objective. If the incoming velocity is a pure
rigid translation/rotation field, the objective cannot distinguish endpoints
having the same prescribed linear momentum, angular momentum, and energy.
Release from rest is one special case. An initially stretched rotating spring
is another.

This is a statement about the objective's ability to rank feasible endpoints,
not a proof that every mesh, initial velocity, or solver trajectory must freeze.
Feasible endpoints may form discrete branches or continuous families; the same
identity applies to either kind of feasible set.

## Exact planar identity

Let o be the fixed incoming center of mass. Suppose each incoming velocity is

\[
v_{n,i}=t+\omega
\begin{bmatrix}-(x_{n,i}-o)_y\\(x_{n,i}-o)_x\end{bmatrix},
\]

where t is one translation velocity and omega is one angular velocity shared
by all particles. There is no incoming nonrigid velocity component. Positions
may nevertheless be deformed and contain substantial elastic energy.

Expanding the objective with zero external force gives

\[
\begin{aligned}
&\frac1{2h^2}\|x-x_n-hv_n\|_M^2+U_s(x)
\\
&=\frac1{2h^2}\|x-x_n\|_M^2+U_s(x)
-\frac1h\sum_i m_i v_{n,i}^T(x_i-x_{n,i})
+\frac12\|v_n\|_M^2.
\end{aligned}
\]

The incoming rigid field makes the cross term

\[
\begin{aligned}
\frac1h\sum_i m_i v_{n,i}^T(x_i-x_{n,i})
&=t^TP(x)
+\frac\omega h\sum_i m_i(x_{n,i}-o)\times(x_i-x_{n,i})
\\
&=t^TP(x)+\omega L_{\mathrm{ang}}(x).
\end{aligned}
\]

The final equality is exact under backward kinematics; the angular constraint
is affine in endpoint position with its lever arm fixed at the incoming state.
Consequently,

\[
\boxed{
\begin{aligned}
&\frac1{2h^2}\|x-x_n-hv_n\|_M^2+U_s(x)
\\
&=\left[\frac1{2h^2}\|x-x_n\|_M^2+U_s(x)\right]
-t^TP(x)-\omega L_{\mathrm{ang}}(x)
+\frac12\|v_n\|_M^2.
\end{aligned}}
\]

On strict P/L/energy constraints, every term on the right is fixed. The
position objective is therefore constant on the complete feasible set. It
does not prefer one radial branch, elastic motion, or distribution of endpoint
positions over another with those same conserved quantities.

## CPD relaxation leaves only scalar rho dependence

With the source energy relaxation and its penalty, substitution of all
conservation rows leaves the objective equal to

\[
\boxed{
(1-\rho)E_\star+\rho K_\star
-t^TP_\star-\omega L_{\mathrm{ang},\star}
+\frac12\|v_n\|_M^2+\frac\epsilon2\rho^2.}
\]

There is still no position dependence among endpoints satisfying those rows
for a given rho. If the preferred energy level is attainable, minimization in
rho gives

\[
\rho=\frac{E_\star-K_\star}{\epsilon},\qquad
\frac1{2h^2}\|x-x_n\|_M^2+U_s(x)
=E_\star-\frac{(E_\star-K_\star)^2}{\epsilon}.
\]

For a pure rigid incoming field with zero external force and COM origin,
K_star equals the incoming kinetic energy. The difference E_star-K_star is
then incoming elastic energy. Energy relaxation picks a nearby scalar energy
level; it does not supply a positional preference between feasible motions on
that level.

At regular points, one compatible stationarity-multiplier choice is

\[
\eta_P=t,\qquad\eta_L=\omega,\qquad\eta=-1,
\]

together with the stated rho. These signs follow the plus-sign conservation
Lagrangian. They are fresh conservation multiplier values, not material lambda
increments. Attainability must still be checked: when E_star=K_star, the
separate feasibility diagnostic shows that a rest-length spring with pure
tangential velocity can have no feasible endpoint at all.

## Anchored central spring

For a single central anchor, omit translation conservation and use angular
momentum about the anchor. An initially tangential incoming velocity has

\[
v_n=\omega\begin{bmatrix}-(x_n)_y\\(x_n)_x\end{bmatrix}.
\]

The same identity simplifies to

\[
\boxed{
\frac1{2h^2}\|x-x_n-hv_n\|_M^2+U_s(x)
=\frac1{2h^2}\|x-x_n\|_M^2+U_s(x)
-\omega L_{\mathrm{ang}}(x)+\frac12\|v_n\|_M^2.}
\]

Thus initial tangential motion at a stretched radius does not resolve the
objective degeneracy merely because angular momentum is now constrained. The
new angular row prevents the energy-only method's unforced angular-momentum
gain, but it can leave the selection between radial/rotating branches to the
algorithm and its metric. The observed nearly constant-radius rotating pilot
is consistent with this limitation; it is not implied uniquely by the identity.

## Limits of the explanation

- The identity is exact for the specified incoming rigid velocity field and
  zero external force. It is not applied here to gravity/opposite-load scenes.
- With incoming nonrigid velocities, the expansion contains additional terms
  that the global P/L constraints do not fix. The objective can then distinguish
  endpoints through those terms.
- The identity applies to the true source objective. PD's frozen local
  quadratic model and constant metric still select a particular numerical
  update, which can break ties algorithmically. A tie-breaking metric does not
  restore a physical preference absent from the constrained objective.
- A constant-radius polygonal position trajectory has a chord velocity with
  a small radial component. After the first timestep, that velocity need not
  remain a pure tangential rigid field. The identity therefore cannot be
  promoted to a proof of every subsequent frozen-radius step.
- Conserved global quantities and smooth motion remain useful properties under
  the user's current criterion. They do not guarantee that stored radial
  deformation evolves intuitively; that is assessed from the actual trajectories.

No new rollout was run for this derivation. It complements the source-fidelity
audit and standalone feasibility fixtures rather than changing the solver or
reclassifying trajectories by phase/frequency accuracy.
