# Gravity release requires order-one energy relaxation

For the compliant pendulum released from rest at its spring rest length, the
source targets require **rho strictly greater than one**, rather than a tiny
energy adjustment. This is an analytic diagnostic of the initial target and
penalty scales. No new solver call or rollout was made for this note.

The selected anchored scene is an explicitly declared restriction of original
CPD: omit linear momentum, retain angular momentum about the central anchor,
and use the source external-power energy target. The same obstruction to small
rho also has a native free-pair counterpart below.

## Initial targets

The actual release parameters are

\[
m=\ell=1,\quad k=10000,\quad g=(0,-9.81),\quad\theta_n=0.65,
\]

\[
x_n=(0.6051864057360395,-0.7960837985490559),\qquad v_n=0.
\]

Incoming spring energy and kinetic energy are zero. The source energy excludes
gravitational potential, and its old-velocity power increment is also zero:

\[
E_\star=0.
\]

Gravity nevertheless supplies a nonzero angular target:

\[
\tau=x_n\times mg=-5.936878640270549,\qquad
L_{\mathrm{ang},\star}=h\tau,\qquad I_n=m\ell^2=1,
\]

\[
K_\star=\frac{(h\tau)^2}{2I_n}>0.
\]

| Rate | Angular target | K_star |
|---|---:|---:|
| 30 Hz | -0.1978959547 | 0.01958140444 J |
| 60 Hz | -0.09894797734 | 0.004895351110 J |
| 120 Hz | -0.04947398867 | 0.001223837777 J |

## Why rho must exceed one

Under chord kinematics the angular row is exactly affine:

\[
L_{\mathrm{ang}}=x_n\times m\frac{x-x_n}{h}.
\]

It fixes the tangential component of chord velocity, measured against the
incoming lever arm, to \(L_{\mathrm{ang},\star}/(m\ell)\). Its kinetic energy
is therefore at least K_star. Nonnegative spring storage gives

\[
\frac{m}{2h^2}\|x-x_n\|^2
+\frac{k}{2}(\|x\|-\ell)^2\ge K_\star.
\]

The energy row is now

\[
\boxed{
\frac{m}{2h^2}\|x-x_n\|^2
+\frac{k}{2}(\|x\|-\ell)^2=\rho K_\star.}
\]

Consequently rho must be at least one. Equality would require simultaneously
minimum kinetic energy (zero incoming-radial velocity component) and zero
spring extension. But the resulting tangential chord ends at radius

\[
\sqrt{\ell^2+h^2\left[L_{\mathrm{ang},\star}/(m\ell)\right]^2}>\ell.
\]

The spring then stores positive energy. Every exactly feasible endpoint
therefore has **rho>1** for the finite positive spring stiffness and nonzero
gravity torque used here. This is not infeasibility of the full unbounded-rho
problem; rho is allowed to exceed one. It is a conflict with treating rho as
an intrinsically tiny relaxation or clamping it to [0,1].

## Consequence for the penalty and multiplier scale

Relaxation stationarity is

\[
\epsilon\rho+\eta(E_\star-K_\star)=0,
\qquad\boxed{\eta=\frac{\epsilon\rho}{K_\star}>\frac\epsilon{K_\star}.}
\]

The source penalty epsilon*rho²/2 is thus at least of order epsilon/2. The
campaign's local epsilon ratio 100 million corresponds to epsilon
200522974.1234 J for this scene. If the source energy/angular rows and rho
stationarity are met, the energy multiplier must exceed approximately
1.024e10, 4.096e10, or 1.638e11 at 30, 60, or 120 Hz, respectively. These are
analytic implications, not measured values from capped unconverged steps.

The Eq. (21) metric omits energy-multiplier curvature. A required large
multiplier can therefore amplify the importance of the dropped curvature and
make small energy-gradient directions or nonlinear changes numerically
difficult. This explains a target/penalty tension worth inspecting in the
capped pendulum steps; it does not prove that every such step must fail.
Required order-one rho is consistent with the paper's intended role for
energy relaxation when external-force momentum targets demand kinetic energy
absent from the old-power target.

## Free-body counterpart and scope

The same initial relative position and zero relative velocity can be assigned
to the two free mass-two particles, with known opposite forces +(0,-9.81) N
and -(0,-9.81) N. Then
P_star=0, incoming COM inertia is one, external torque is the same tau, and
E_star=0. The full original planar P/L/energy system has the same relative
rho>1 requirement without an anchor. That is an analytic counterpart, not an
additional run. The campaign's existing opposite-load fixture starts stretched,
so it does not instantiate this particular zero-elastic-energy initial state.

This note concerns source internal energy and its old-power target. It does
not claim that actual mechanical energy including gravity should rise from
zero, and it does not reinterpret small physical-energy error as successful
source constraint convergence. Those quantities remain separately reported.
