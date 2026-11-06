# Scalar spring: exact energy visit followed by XPBD material visit

The scalar coordinate q is spring extension, with mass m, stiffness k, potential U(q)=k q²/2. The timestep starts at (q_n,v_n), its inertial predictor is y=q_n+h v_n, and the reconstructed endpoint velocity is (q-q_n)/h. We compare the fixed initial energy target with a target reset to the current physical energy at the beginning of every timestep. These are different methods.

## Objective and exact energy constraint

Before introducing multipliers, the constrained implicit objective is

F(q)=m(q-y)²/(2h²)+kq²/2, subject to E_n(q)=m(q-q_n)²/(2h²)+kq²/2=E_target.

Set σ=kh²/m. Completing the square gives

E_n(q)=A(q-c)²/2+E_min,

where A=m/h²+k, c=q_n/(1+σ), and E_min=kq_n²/[2(1+σ)]. The feasible positions are

q_±=c ± sqrt(2(E_target-E_min)/A).

We compute these completed-square roots directly, rather than the cancellation-prone general quadratic formula. If E_target<E_min beyond roundoff, no real root exists: the timestep is reported as infeasible and is not silently advanced. A negative difference within a relative roundoff tolerance is clamped to zero.

The exact E visit picks the root nearest the current primal iterate q. If both are equally near within roundoff, it picks the root nearest the original inertial predictor y; the final tie uses the lower root. The tie is a numerical policy, not a physical consequence of energy conservation. We count ambiguous visits and branch changes explicitly. A converged energy visit can choose the zero-displacement root at rest.

## XPBD L visit and ordering

After deriving the original objective, introduce the compliant material multiplier λ with α_tilde=1/(kh²). The material constitutive row is q+α_tilde λ=0. One L visit uses

Δλ=-(q+α_tilde λ)/(1/m+α_tilde),

λ←λ+Δλ, q←q+Δλ/m.

Both variables change immediately; λ is accumulated across all outer sweeps of a timestep and reset only between timesteps. An outer sweep is exactly E, then L. The scalar spring row is affine, so its refreshed and frozen forms coincide. We save the post-E and post-L states: satisfying energy at E does not imply it survives L.

For the linear E comparator, H=E_n(q)-E_target and ∇E_n=m(q-q_n)/h²+kq. The visit is

Δη=-H / [(∇E_n)²/m], q←q+(∇E_n)Δη/m, η←η+Δη.

It is skipped only for a gradient indistinguishable from zero at floating-point precision; the event is counted. The comparator explicitly drops stationarity and multiplier curvature, as the historical raw EL algorithm does.

The exact E visit directly changes q to a nonlinear energy root. We do not pretend that it supplies a unique accumulated energy multiplier η consistent with the full KKT equations. In one dimension a nonzero-gradient displacement could be encoded as mΔq/∇E at that visit, but this becomes undefined at a zero gradient, and reusing it does not restore stationarity after the gradient changes. The exact-E/L method is an alternating primal projection and XPBD material update, not an exact constrained optimizer.

## Exact constrained-objective control

We also solve the original constrained objective exactly in this scalar example: evaluate F at both feasible roots and choose the smaller value, with equal objective values resolved by proximity to y, then the lower root. This control has no material iterations.

The identity

F(q)-E_n(q)=-(m v_n/h)(q-q_n)+m v_n²/2

shows that the objective is constant on the entire feasible set when v_n=0. At release from rest, q=q_n is feasible, so the deterministic closest-to-predictor rule selects no movement. Energy conservation by itself does not resolve the physical branch.

## Why more exact-E/L sweeps can freeze release

At release q_n=a,v_n=0 with E_target=ka²/2, the two energy roots are a and a(1-σ)/(1+σ). The first E visit chooses a. The first L visit gives a/(1+σ), exactly the midpoint of the two energy roots. The next nearest-root choice is a genuine tie; our tie policy selects a again. If later E visits keep that root r, the material multiplier recurrence is

λ_next=[σ/(1+σ)](λ-mr),

whose fixed point is λ=-mσr=-r/α_tilde. Therefore q after L converges to r=a. More sweeps can converge rapidly to the stationary, no-motion branch. We record this behavior instead of treating energy closure as successful dynamics.

## Experiment and measurements

Cases: m=1,k=100, q_0=.25, v_0=0 (release) or .8 (mixed). Timesteps 1/60,1/120,1/240,1/480 seconds; 1,4,16 E/L sweeps; five simulated seconds. The exact objective control is run once per case/timestep/target.

Compare against the analytic harmonic oscillator q(t)=q_0 cos(ωt)+(v_0/ω)sin(ωt), v(t)=-ωq_0 sin(ωt)+v_0 cos(ωt), ω=sqrt(k/m). Record normalized energy errors before and after every row, endpoint position/velocity RMS error, path length, largest one-row displacement, no-root events, ambiguous root choices, near-zero-gradient skips, and accepted physical time. Nonfinite or runaway states end a run; shorter prefixes are not counted as successful five-second results.
