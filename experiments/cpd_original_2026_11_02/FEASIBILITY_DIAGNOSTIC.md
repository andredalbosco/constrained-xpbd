# A standalone infeasible energy/angular target

The CPD energy-relaxation variable does not guarantee feasibility for every
spring state. This diagnostic gives a closed-form obstruction and tests the
unchanged, validated Eq. (21) implementation on it. It is outside the primary
13-scene campaign and is not a reproduction of an authored paper scene.

The source-based implementation is distinguished from the authors' binary/code;
the result concerns the equations and declared numerical choices used here.

## Exact obstruction

Consider two free particles of mass 2 each, connected by one rest-length-one
spring. At the beginning of the step,

\[
x_{n,a}=(1/2,0),\quad x_{n,b}=(-1/2,0),\qquad
v_{n,a}=(0,u/2),\quad v_{n,b}=(0,-u/2),\quad u\ne0.
\]

The particles have no external forces and reduced mass one. The linear target
is zero, the angular target is u, and the incoming inertia about COM is one.
Therefore

\[
E_\star=\frac{u^2}{2}=K_\star,
\qquad E_\star-K_\star=0.
\]

The interpolation is constant for every rho. No epsilon value can relax the
energy equality in this state.

Zero total momentum fixes COM displacement to zero. The angular constraint,
under backward velocity, fixes the tangential component of relative chord
velocity to u. Write its unconstrained radial component as w, so the relative
endpoint position is

\[
(1+hw,hu).
\]

The relative kinetic energy is \((w^2+u^2)/2\), and the internal energy exceeds
the required target by

\[
\boxed{
\frac{w^2}{2}
+\frac{k}{2}\left[\sqrt{(1+hw)^2+(hu)^2}-1\right]^2.}
\]

Both terms are nonnegative. Equality would require w=0 and zero spring
extension simultaneously. But w=0 gives endpoint length
\(\sqrt{1+(hu)^2}>1\). Thus the gap is strictly positive for every real w
when h>0, k>0, and u is nonzero. **No exact endpoint satisfies all three
P/L/energy targets, even with unbounded rho.**

This is an elastic/backward-kinematics incompatibility above the compatible
kinetic lower bound. It is separate from release-from-rest freezing and from
an arbitrary rank-rejection threshold.

## Certified scalar minimum

The scalar objective above is coercive. Every global minimizer must have gap
no larger than its value at w=0. Its kinetic term consequently confines all
global minimizers to

\[
|w|\le\sqrt{2\left\{
\frac{k}{2}\left[\sqrt{1+(hu)^2}-1\right]^2\right\}}.
\]

On this compact interval, the diagnostic bounds the exact second derivative
away from zero. The full scalar objective has second derivative

\[
1+kh^2\left[
1-\frac{(hu)^2}{\left((1+hw)^2+(hu)^2\right)^{3/2}}
\right].
\]

The lower bounds in the two fixtures are approximately 2.10617 and 1.02777.
Thus there is one minimum on an interval containing every global minimizer.
A bracketed derivative root, an independent bounded minimization, and separate
bounded searches on both interval halves agree on the minimum value. This is
not a claim of global optimality based on one unconstrained optimizer start.

| Fixture | h | u | k | Global minimum energy gap | Relation to raw paper tolerance |
|---|---:|---:|---:|---:|---|
| Strong | 1/30 | 2 | 1000 | 0.0011697260 J | Above 1e-4 |
| Small | 1/60 | 0.8 | 100 | 0.0000003843521 J | Below 1e-4 |

The minimizing radial chord velocities are -0.0350527853 and -0.00014412562,
respectively. The proof of strictly positive gap applies independently of the
specific values or the numerical minimization.

## Observed Eq. (21) step results

The production solver was validated before these calls. Both fixtures use the
same direct scaled Schur solve, constant PD metric, fresh multipliers, full
updates, source energy target, and no row deletion or fallback. The policies
are the declared post-update conservation tolerance, fixed 16, and strict
stationarity capped at 200.

Epsilon ratios 10,000 and 100 million produce identical state and residual
results in each fixture/policy, as expected when the energy row has zero rho
derivative. Rho stays zero. The table applies to both ratios:

| Fixture | Policy | Iterations | Accepted finite step? | Converged? | Maximum raw constraint residual |
|---|---|---:|---|---|---:|
| Strong | Paper post-update | 200 | Yes | No, iteration cap | 0.0011986494 |
| Strong | Fixed 16 | 16 | Yes | No | 0.0080221009 |
| Strong | Stationarity | 200 | Yes | No, iteration cap | 0.0011986494 |
| Small | Paper post-update | 1 | Yes | Yes under 1e-4 criterion | 0.0000036546283 |
| Small | Fixed 16 | 16 | Yes | No under strict diagnostic | 0.0000014180599 |
| Small | Stationarity | 200 | Yes | No, iteration cap | 0.00011125395 |

No actual singular solve or explicit runaway guard occurs in these twelve
one-step tests. The minimum estimated normalized Schur rank is four, the full
number of rows. The strongly incompatible fixture nevertheless cannot close
its energy row; at 200 iterations its P/L residuals are exactly zero in the
saved state while its energy residual remains 0.0011986494 J. Its maximum
normalized Schur condition is about 2.53e8, and the largest frozen subproblem
linear residual is approximately 2.55e-11. Accurately solving the successive
linear subproblems does not create a feasible nonlinear endpoint.

For the small fixture, the paper-style rule legitimately stops within its
declared approximate tolerance even though exact mathematical feasibility is
absent. This is not automatically a practical animation failure at that error
scale. It shows why a numerical convergence flag must not be read as proof of
exact feasible targets. Continuing strict iterations does not monotonically
improve the result: its 200-iteration residual is larger than after one or
sixteen updates.

## Reproduction and limits

Run `python3 feasibility_diagnostic.py` in the experiment environment. It writes
`feasibility_diagnostic_results.json`, containing the proof bounds, scalar
minimization cross-checks, complete target values, both epsilon ratios, all
policy outputs, and every serialized update.

This experiment establishes one limitation of the source interpolation, not a
general statement that CPD cannot produce useful energy/momentum-conserving
motion. It does not add a new parameter grid, alter the primary scene results,
or use phase/frequency error as a criterion. The key distinction is between a
kinetic lower bound and feasibility of the full elastic endpoint constraints.
