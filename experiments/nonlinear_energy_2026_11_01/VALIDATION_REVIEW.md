# Independent reference validation

`references.py` supplies diagnostic integrators and continuous references. The discrete-gradient controls change the endpoint kinematics; neither is a solver for the original endpoint-velocity EL objective.

## Discrete-gradient equations

For the anchored planar spring, let `r_i=|x_i|`, `U_s=.5 k(r-l)^2`, and let `g` be an acceleration vector. The symmetric discrete gradient is

\[
D_s(x_0,x_1)=\frac{k}{2}(r_0+r_1-2l)\frac{x_0+x_1}{r_0+r_1}.
\]

Its dot product with `x1-x0` equals `.5 k(r0+r1-2l)(r1-r0)=U_s(x1)-U_s(x0)`. At equal nonzero endpoints it reduces to the ordinary spring gradient. The factor one-half is necessary. When `l=0`, the nonsingular expression is `k(x0+x1)/2`; a positive-rest-length spring has no unique continuous force at the origin.

The step solves

\[
x_1-x_0=h(v_0+v_1)/2,\qquad m(v_1-v_0)=-hD_s+hm g.
\]

Thus kinetic change, spring storage change, and exact constant-gravity work sum to zero on a converged root. `v1` is the endpoint velocity reconstructed from midpoint kinematics, **not** `(x1-x0)/h`.

For the rigid control, `x=l(sin(theta),-cos(theta))` and `v=l omega(cos(theta),sin(theta))`. Its independent angular velocity obeys

\[
\delta=h(\omega_0+\omega_1)/2,\quad
ml^2(\omega_1-\omega_0)=-hD_\theta U,
\]
\[
D_\theta U=-ml\,\operatorname{sinc}(\delta/2)
 [g_x\cos(\theta_0+\delta/2)+g_y\sin(\theta_0+\delta/2)].
\]

This preserves endpoint length, endpoint tangency, and mechanical energy on a converged root. It satisfies angular midpoint kinematics; Cartesian midpoint kinematics generally differs. The tested `.65 rad, omega=1.2, h=.08` step has Cartesian midpoint displacement defect `7.59e-4 m`, intentionally exposing this distinction. Angular increments remain unwrapped, including rotations beyond `2 pi`.

## Solvers and failure semantics

The spring position equation uses SciPy `hybr`, followed on failure by least squares for the same residual, then explicitly reported `hybr` retries from the old position and free-flight predictor. Each attempt is identified in `method` and `reason`; `nfev` counts all residual calls. No clamping, substeps, energy projection, or velocity rescaling is applied. A verified residual establishes a root, not a unique branch or a minimum of another objective.

The pendulum scalar equation has a guaranteed bracket centered on `h omega0` with radius `h^2 |g|/(2l)`. Brent's method finds a root; a sufficient uniqueness condition is `h^2 |g|/(4l)<1`, since the derivative of the angular discrete gradient has magnitude at most `ml|g|/2`. Larger steps retain an explicit uncertified-branch flag. A forced one-iteration solve returns `converged=False` and residual `0.136`, rather than accepting a finite candidate.

Rigid initial length or tangency defects larger than relative tolerance `1e-12` are rejected. The reported energy change uses actual input/output Cartesian states, including any smaller accepted roundoff mismatch. Positive-rest-length spring force at an exact origin is rejected rather than clamped. Failed numerical steps retain their candidate, residual, and reason; callers must check `converged` before advancing.

## Executed checks

Run `/tmp/xpbd-nonlinear-venv/bin/python references.py`; the saved output is `reference_validation.json`. Environment: Python 3.13, NumPy 2.5.4, SciPy 1.18.1.

| Check | Spring | Rigid pendulum |
|---|---:|---:|
| Max step energy defect over 200 steps, h=.02 | 1.83e-14 J | 1.45e-14 J |
| Forward/backward one-step state discrepancy | 0 to measured precision | 6.47e-16 |
| DOP853 refinement discrepancy over 2 s, max position/velocity norm | 1.23e-11 | 6.47e-12 |
| Endpoint rod-length / tangency defect | Not applicable | 1.11e-16 / 1.11e-16 |

The spring discrete chain rule also passes 256 deterministic random endpoint pairs with maximum discrepancy `1.42e-14`. Initial radial rod velocity is rejected. A zero-gravity rotation of 10 radians preserves its unwrapped angular increment and velocity.

Continuous references use DOP853 with dense output at requested times. Default tolerances `rtol=2e-12, atol=2e-14` were checked against `2e-13,2e-15` and `max_step=.01`. Integrator failure returns only the valid requested prefix, never extrapolated samples. These refinement checks validate the measured two-second cases, not every future stiff or singular trajectory.

At a common 2 s endpoint, step refinement gives:

| h | Spring position error | Spring velocity error | Pendulum position error | Pendulum velocity error |
|---|---:|---:|---:|---:|
| .1 | .06127 | .41262 | .02339 | .07853 |
| .05 | .01786 | .09515 | .005808 | .01988 |
| .025 | .004626 | .02311 | .001449 | .004985 |

Errors decrease toward second order while energy remains approximately conserved. Conservation alone does not certify phase accuracy.

The initial 50-case stress probe is retained as `reference_stress_validation.csv`: compressed/stiff/near-rest springs and bottom/top/near-rest/near-separatrix/high-speed rods, with h=.02,.08,.2,.5,1. Its original acceleration-predictor search failed on the stiff spring at k=4000,h=.2, leaving residual 1.089 and a large energy defect. This was a root-search failure: independent old-position and free-flight starts found the same residual-verified root. The explicitly reported deterministic retries were added for that demonstrated failure; the old failed record is retained rather than relabeled as a mechanics success. Large-step existence and small residual do not certify selection of a physically accurate branch.
