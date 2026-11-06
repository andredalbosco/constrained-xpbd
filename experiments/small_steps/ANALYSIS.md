# Small-step EL analysis

This study uses the original endpoint velocity, one frozen E-then-L sweep per substep, fresh gradients, and reset multipliers. The 240/480/960 Hz rates mean 4/8/16 substeps per 60 Hz frame. This changes the integration schedule; it does not perform more iterations of one large-step constrained problem.

## What resets imply

At the predictor \(y=x_n+h v_n+h^2a\), \(\lambda=\eta=0\), so \(R=0\). Write \(u=\nabla U_{elastic}=k C n\). Then the primal objective gradient is \(q=u\), whereas the energy gradient is
\[
g=m v_n/h+u.
\]
Gravity is already included in the predictor. The actual energy residual includes gravitational potential. The final material correction changes the energy after E has been visited, so the endpoint energy is not generally the value imposed by the energy row.

The displacement after the two scalar solves is
\[
\Delta x=-\frac{n C}{m A}+\frac{g-nB/A}{m}\Delta\eta,
\quad A=1/m+\tilde\alpha,\quad B=n^Tg/m.
\]
The first term is ordinary material XPBD, and the second is the energy increment after the material row reacts.

## Raw and effective directional gates are equivalent here

The raw test uses \(q^Tg\Delta\eta/m\). The effective test uses \(q^T(g-nB/A)\Delta\eta/m\). Because \(q=k C n\),
\[
q^T(g-nB/A)\Delta\eta/m
=\frac{\tilde\alpha}{A}\,q^Tg\Delta\eta/m.
\]
For a compliant spring the factor is positive. Both gates therefore make the same sign decision in exact arithmetic. This is specific to one radial spring, one predictor sweep, and reset multipliers; it is not a claim about arbitrary multiple-constraint iterations.

All matching raw/effective variants in the audited sweep produced identical physical metrics and acceptance counters. The extra effective test therefore adds no independent stabilization in this experiment.

## Compliance versus relaxation

At the predictor, soft energy replaces the hard energy multiplier increment by
\[
\Delta\eta_\beta=s_\beta\Delta\eta_0,
\quad
s_\beta=\frac{h^2\|g\|^2/m}{h^2\|g\|^2/m+\beta}
=\frac{m\|v_n+h u/m\|^2}{m\|v_n+h u/m\|^2+\beta}.
\]
Thus, for this single sweep, compliance acts exactly like a state-dependent relaxation of the hard correction. This equivalence concerns this update, not the full converged objective. With physical \(\beta=b(E_0-V_{min})\) fixed across step sizes, the reduction remains finite at nonzero velocity and tends to zero as \(O(h^2)\) at rest when \(u\ne0\). It suppresses the small-gradient singularity much more strongly than a constant relaxation factor.

A constant factor merely rescales the hard correction. It cannot guarantee safety: an \(O(1)\) erroneous correction remains \(O(1)\) for any fixed positive factor as \(h\to0\). Compliance allows a nonzero physical energy error; a small soft-row residual is not proof of mechanical energy conservation.

## What the magnitude cap does

The bound is proportional to \(S=h\|v_n\|+h^2\|f_n\|/m\). Away from turning points it scales as \(O(h)\); near rest it scales as \(O(h^2)\), the appropriate size for a force-driven position change. It therefore limits singular corrections that a directional test can accept. The bound applies to the energy contribution, not the entire position correction.

The directional test only checks a first derivative. Even a downhill correction can increase the actual objective because higher-order terms can outweigh the negative first-order change. Smaller steps do not uniformly cure that problem: the energy gradient can become correspondingly smaller near a cancellation.

## How to judge the results

Small steps can reduce ordinary discretization error while the solver continues to use an approximate, unconverged energy correction. More substeps are not equivalent to converging the constrained objective within one step. Energy, trajectory, reversals, and torque balance must be considered together. A path-length ratio near one does not exclude phase error, and a finite completed run does not exclude a very large energy excursion.

In the current 465-run sweep, all cases complete but no tested family gives accurate energy and motion across all five cases. At 480 Hz, the baseline has substantial energy loss but smaller position error than the conservation variants on the moderate resting swing and coupled case. Directional-only and fixed-relaxation variants still have large vertical spikes; these are not cured monotonically by increasing the substep count.

Some useful comparisons at 480 Hz (one sweep per substep):

| Variant | Vertical position RMS | Vertical max energy error | Coupled position RMS | Coupled max energy error |
|---|---:|---:|---:|---:|
| Baseline | 0.0593 | 87.54% | 0.1736 | 60.90% |
| Directional + cap 0.001 | 0.0536 | 6.46% | 0.6020 | 4.06% |
| Directional + energy beta ratio 10 | 0.1346 | 5.92% | 0.6450 | 3.94% |
| Directional + energy beta ratio 100 | 0.0457 | 11.36% | 0.5288 | 6.71% |
| Directional + energy beta ratio 1000 | 0.0409 | 65.02% | 0.2434 | 32.27% |

Percentages use initial energy above the minimum potential. Strong softening recovers more ordinary XPBD-like motion at the expense of energy control. At 960 Hz, beta ratio 100 improves to vertical RMS 0.0323 and coupled RMS 0.2909, with 3.78% and 2.37% energy error, but the coupled trajectory still misses one reference reversal. These are tradeoffs, not a completed general conservation method.

## Independent audit

The audited 465 rows are complete. For every matching raw/effective pair, physical metrics and counters are identical. Every accepted gated correction has nonpositive tested objective work; every cap respects its bound within floating-point tolerance. Attempt counters obey accepted + directional skips + singular = attempts. Scalar solve counts include the candidate energy division even when that candidate is rejected, so scalar solves = completed material solves + nonsingular energy attempts.

The 135 selected reruns match the sweep summary field-for-field. The 216,045 saved rows reconstruct endpoint velocity exactly from consecutive positions; reconstructed energy agrees within 1.80e-15 relative error. At 240 Hz, old one-sweep baseline/directional runs have identical trajectory-only metrics and acceptance counters. Reference-derived RMS/path metrics differ at roundoff scale because this harness uses 16 RK4 subdivisions, checked against 32, whereas the older harness used 32 checked against 64.
