# Constrained PD: exact energy, surrogate energy, and release degeneracy

The primary control should constrain **physical** energy. Frozen-projection energy is a useful second control, but exact solution of that surrogate does not mean exact physical energy at a finite iteration count. Neither choice repairs the release degeneracy of the full implicit objective with endpoint velocity reconstruction.

PD introduces local geometric projections and a quadratic global solve; that is the mechanism borrowed here, not a new conservation guarantee. [Bouaziz et al., 2014, author page](https://tiantianliu.cn/papers/bouaziz14projective/bouaziz14projective.html). The published *Constrained Projective Dynamics* is a particular energy-momentum method with an energy relaxation variable, distinct from the strict controls below. [Kee et al., 2021](https://doi.org/10.1145/3450626.3459878). The prior [work-balance audit](/home/blacknoir/Projects/dynamics-research/archive/docs/work_balance/report.md) already identifies strict release degeneracy; the equations below independently specialize it to the requested two scenes.

## Objective and the exact obstruction

Let `d=x-x_n`, `A=m/h²`, and initially reconstruct `v=d/h`. Include gravity in `U(x)` rather than also shifting the predictor:

\[
F(x)=\frac A2\|d-hv_n\|^2+U(x),\qquad
H(x)=\frac A2\|d\|^2+U(x),\qquad
H_n=\frac m2\|v_n\|^2+U(x_n).
\]

Constrain `H(x)=H_n`, and for the pendulum also `C(x)=||x||²-l²=0`. Then

\[
F|_{H=H_n}=H_n-\frac m h v_n\cdot d+\frac m2\|v_n\|^2.
\]

At release `v_n=0`, **the objective is constant on the entire feasible set**. The old configuration is feasible, but need not be its only point. A frozen result can therefore be a fully converged constrained solution. A better nonlinear solver or PD preconditioner does not remove this degeneracy.

Using only part of the same objective does not automatically help. For

\[
F_\gamma=\frac A2\|d-hv_n\|^2+\gamma U(x),
\]

the restricted objective is `γH_n+(1-γ)A||d||²/2-(m/h)v_n·d+constant`. At rest, `γ<1` selects `d=0`, `γ=1` is flat, and `γ>1` favors large displacement. Inertia-only `γ=0` strengthens freezing. If gravity is instead absorbed into `y=x_n+hv_n+h²g`, its linear potential is already present in the inertial square up to a constant; do not count it twice or call its removal merely a solver change.

## Planar radial spring: an exact feasible-set oracle

For `U=k(||x||-l)²/2`, write `x=r n(θ)`, `r_n=||x_n||`, and `n(θ)·x_n=r_n cosφ`. Define

\[
B(θ)=A r_n\cos φ+kl,\quad D=A r_n^2+kl^2-2H_n.
\]

Every physical-energy point has

\[
r_\pm(θ)=\frac{B(θ)\pm\sqrt{B(θ)^2-(A+k)D}}{A+k},
\]

retaining only real, nonnegative radii. This gives a one-dimensional constrained global-step oracle: optimize the proposed objective over both branches and their discriminant-zero endpoints. A dense angular scan plus refinement is numerical global-search evidence, not a formal certificate.

On the original radial ray at rest, with `q_n=r_n-l`, the two roots are

\[
d_r=0,\qquad d_r=-\frac{2h^2(k/m)q_n}{1+h^2k/m}.
\]

Reject a root if its radius is negative; it no longer lies on that ray. The moving root gives `v=d_r/h≈-2h(k/m)q_n`, twice the physical leading release velocity, and displacement four times the physical leading displacement. Thus choosing the moving exact-energy branch alone does not yield consistent release mechanics under this velocity definition.

## Rigid pendulum: circle/line intersection

Let `U=-m g·x`, with gravity acceleration vector `g`, and assume `||x_n||=l`. On the rod circle, physical energy reduces to

\[
b\cdot x=c,\qquad b=x_n+h^2g,\qquad c=l^2-h^2H_n/m.
\]

Feasibility requires `|c|≤l||b||`. For `b≠0`, the at-most-two points are

\[
x_\pm=\frac c{\|b\|^2}b
\pm\sqrt{l^2-\frac{c^2}{\|b\|^2}}\frac{Rb}{\|b\|},
\]

where `R` rotates a planar vector by 90 degrees. Handle `b=0` separately: `c=0` leaves the whole circle feasible; otherwise there is no solution. At rest, one solution is `x_n`; the other is `2(c/||b||²)b-x_n`, with displacement `≈2h²g_t`. This again gives twice the correct leading release velocity.

Do not add endpoint tangency while retaining `v=(x-x_n)/h`:

\[
x\cdot v=\frac{\|x-x_n\|^2}{2h}
\]

for feasible old/new rod endpoints. Exact endpoint tangency then forces no motion. With nonzero old kinetic energy, that frozen reconstructed velocity also contradicts energy conservation. This is geometric incompatibility, not a PD convergence defect.

## Concrete constrained local/global controls

For the spring, update `p^j=l x^j/||x^j||` when `x^j≠0`. The global surrogate is

\[
Q_j(x)=\frac A2\|x-x_n-hv_n\|^2+\frac k2\|x-p^j\|^2,
\qquad \nabla^2Q_j=(A+k)I.
\]

**Physical-energy control:** minimize `Q_j(x)` subject to the true `H(x)=H_n`, using the radial parameterization above. The surrogate majorizes `F` and touches it at the local projection, so exact alternating minimizations decrease `F` on the unchanged feasible set. At rest `F` is flat; this is not a strict descent or motion guarantee. The constrained global step is nonlinear despite its quadratic objective.

**Frozen-energy control:** impose

\[
H_j(x)=\frac A2\|x-x_n\|^2+\frac k2\|x-p^j\|^2=H_n.
\]

This is a circle after completing the square, but

\[
H_j(x)-H(x)=kl\|x\|\,[1-\cos\angle(x,p^j)]\ge0.
\]

A point satisfying the surrogate can lose physical energy. Refreshing `p` exposes this mismatch. Record both residuals. At rest `Q_j=H_j`, so the global objective is flat on its surrogate circle too. A nearest-old-state tie rule is an explicit branch rule that can freeze.

For a hard rod, enforce its circle exactly and compare the physical-energy intersections; substituting a large spring penalty changes the problem. A quadratic proximal term centered on the current feasible iterate can freeze immediately. A fixed mechanical target is a different objective, and even the implicit target is equidistant from the two release roots in these symmetric examples. A PD matrix used only to precondition a nonlinear constrained solve changes neither its exact feasible set nor stationary points.

## Nondegenerate controls, explicitly changing the contract

Keep these separate from the endpoint-velocity experiment:

- **Midpoint kinematics:** `v=2d/h-v_n`. A scalar spring release then has moving energy root `d=-2h²(k/m)q_n/(4+h²k/m)`, with the correct leading displacement/velocity. The repaired quarter objective selects this mechanical root; energy alone still admits a frozen root at rest.
- **Midpoint discrete-gradient mechanics:** solve `d=h(v+v_n)/2`, `m(v-v_n)=-h bar∇U`, where for the radial spring `bar∇U=(k/2)(r+r_n-2l)(x+x_n)/(r+r_n)`. Its exact chain rule closes energy when momentum and kinematics also converge. Treat vanishing denominator/antipodal geometry explicitly.
- **RATTLE rod control:** independent endpoint velocity permits both rod position and endpoint tangency. Standard RATTLE is a constrained mechanical baseline, not an exact-energy claim; any subsequent energy projection must state its own objective and constraints.

Tests must separate physical-energy residual, surrogate residual, objective stationarity, endpoint geometry/tangency, branch selection, and refined-reference motion. Include rest, near-rest, compression, turning points, large angular steps, empty intersections, and repeated steps. A release that stays still with zero energy error is a failed motion test, even when it solves the chosen constrained objective exactly.
