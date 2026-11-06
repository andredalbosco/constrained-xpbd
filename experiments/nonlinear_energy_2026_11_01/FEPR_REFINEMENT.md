# FEPR refinement and branch diagnostics

The implemented planar FEPR fast projection passes an independent source-and-numerics audit, but accurate energy does **not** make it a reliable motion solver in these single-particle tests. Refinement is nonmonotone, and the stiff orbit contains large position-branch jumps. These findings concern the declared implicit-Euler base, auxiliary metric, and fast projection; they are not a general impossibility result for energy projection.

## Sources, scope, and validation

The algorithm is Dinev et al.'s [FEPR](https://users.cs.utah.edu/~ladislav/dinev18FEPR/dinev18FEPR.pdf), Equations 8 and 12–14, as specified in [FEPR_DERIVATION.md](FEPR_DERIVATION.md). It resolves energy and affine momentum constraints using a current-iterate projection. Convergence does not certify the fixed-center closest-point KKT conditions. The spring base is the exact radial implicit-Euler minimum. The rod version adds simultaneous geometry and tangency constraints to a RATTLE base and is explicitly an extension of the paper.

The [independent audit](fepr_independent_validation/SUMMARY.md) passed 540/540 checks. Maximum Jacobian finite-difference error was 2.32e-10 relative; two successive fast updates agreed with independently assembled increment KKT solves within 1.14e-12. Checks cover unrestricted auxiliary values, current-iterate versus fixed-center updates, independent endpoint velocity, target semantics, caps and regularization. The implementation self-test also passed eight fixed-energy one-second trajectories, with energy error below 9.70e-11. These establish implementation consistency, not physical trajectory accuracy.

## Scalar release: tighter constraint tolerance does not remove phase error

The planar motion stays on one axis, so the reference is analytic: m=1, k=100, l=1, x(0)=(1.25,0), v(0)=0, g=0. With q=x_1-l, q(t)=0.25 cos(10t), v_1(t)=-2.5 sin(10t). Each trajectory lasts five seconds. Phase is unwrapped atan2(-v_1/10,q), compared with 10t. All 30 targeted runs completed; no Schur regularization occurred.

Default epsilon=.001, fixed initial-energy target, raw constraint tolerance 1e-10:

| Steps/s | Position RMS | Final phase error (rad) | Maximum absolute energy error | First-step position error | First-step velocity error |
|---:|---:|---:|---:|---:|---:|
| 30 | .219487 | -2.68881 | 5.66e-11 | 1.55e-3 | 4.34e-2 |
| 60 | .035862 | -.337190 | 6.31e-12 | 1.14e-5 | 6.80e-4 |
| 120 | .051206 | -.436814 | 3.96e-11 | 7.40e-6 | 8.88e-4 |
| 240 | .048396 | -.431259 | 9.74e-11 | 5.86e-7 | 1.41e-4 |
| 480 | .042219 | -.400677 | 9.22e-11 | 3.86e-8 | 1.85e-5 |

Changing tolerance from 1e-7 to 1e-10 barely changes these phase/RMS values. Replacing the fixed target by the incoming step energy also barely changes them. Per-step targeting allows small energy errors to accumulate: the 60 Hz, 1e-7 run reaches 7.38e-6 absolute energy drift, versus 8.66e-8 maximum error with a fixed target. Neither explains the much larger trajectory discrepancy. Improving the first release step at high rates also does not guarantee small accumulated error around the whole orbit.

Full numbers, including auxiliary extrema and iteration counts: [results/fepr_refinement.csv](results/fepr_refinement.csv). The tested range does not establish an asymptotic convergence order or prove failure as h tends to zero.

## Why the auxiliary metric deserves attention

Let Delta v=v_n-v_* for this one-particle problem. The exact linear-momentum constraint is

\[
v-v_*=s\,\Delta v.
\]

Every fast update preserves this affine relation. Its increments therefore satisfy delta v=Delta v delta s. Restricting the local quadratic metric to that direction gives

\[
\frac12mh^2\|\delta v\|^2+\frac\epsilon2(\delta s)^2
=\frac12\left(mh^2+\frac\epsilon{\|\Delta v\|^2}\right)\|\delta v\|^2,
\qquad \Delta v\ne0.
\]

Thus the effective velocity-correction cost is not simply mh². For Delta v=O(h) and fixed epsilon, the auxiliary contribution grows as O(h^-2). At exactly Delta v=0 the momentum row fixes v=v_*; this is a restriction, not an infinite coefficient to evaluate numerically. The angular auxiliary adds another coupled restriction in planar orbit tests.

This identity is exact. The explanation that it steers energy repair toward position, especially near a vanishing potential-energy gradient, is a **mechanistic hypothesis** for the observed refinement behavior; these experiments do not isolate it as the sole cause.

Two already-run diagnostic epsilon controls show substantial sensitivity. At fixed energy and tolerance 1e-10:

| Steps/s | Default .001 position RMS | epsilon=1e-9 RMS | epsilon=.001h² RMS |
|---:|---:|---:|---:|
| 30 | .219487 | .202828 | .203107 |
| 60 | .035862 | .092173 | .091888 |
| 120 | .051206 | .039676 | .039276 |
| 240 | .048396 | .017245 | .016489 |
| 480 | .042219 | .007336 | .006384 |

These are sensitivity experiments, not a proposed replacement default. They worsen some lower-rate results, worsen the first release velocity error, and permit much larger auxiliary extrapolation (maximum |s| about 212 for epsilon=1e-9). No robustness claim follows from their finer-step improvement.

## Stiff planar orbit: actual branch jumps

The existing stress case has m=1, k=1000, l=1, x(0)=(1.1,0), v(0)=(0,.8), g=0. Its energy is 5.32 and exact angular momentum is .88. The following diagnostics read the existing [stress states](stress_results/fepr_states.csv) and [substep states](substep_results/fepr_states.csv); no new trajectory sweep was required.

| Internal steps/s | Output steps/s | Target | Position RMS | Maximum relative energy error | Final angular momentum | Largest recorded position step |
|---:|---:|---|---:|---:|---:|---:|
| 30 | 30 | fixed | 1.91371 | 6.19e-9 | 1.02400 | 2.00821 |
| 60 | 60 | fixed | .308534 | 1.81e-8 | 1.13774 | .07109 |
| 120 | 30 | per internal step | 1.36306 | 9.32e-8 | .84945 | 2.04731 |
| 240 | 60 | per internal step | .921400 | 2.63e-6 | 1.03330 | 1.87097 |
| 480 | 60 | per internal step | .116224 | 1.52e-6 | .92064 | .06294 |

The 240 Hz internal trajectory sampled at 30 Hz instead gives RMS .920442 and largest recorded step 1.82114; these are sampling differences for the same trajectory. The table is not a perfectly matched target/sampling refinement study. At 30 and 60 Hz, the recorded fixed/per-step controls have essentially identical motion, so changing the target alone does not fix those cases.

At 30 Hz, the **second** projected step moves 2.00821 metres at t=.066667. The position correction relative to the tentative state is 1.93126 metres. Recorded radii remain between .902978 and 1.100388 throughout the run. The energy bound implies |v|<=sqrt(2H/m)=3.2619, so physical travel over one 1/30-second interval cannot exceed .10873 metres. This jump cannot be interpreted as merely ordinary physical angular motion with a slightly wrong phase. Accurate endpoint energy does not constrain the path between disconnected selected positions.

Although auxiliary scalars were not saved in the stress traces, the nonzero momentum differences allow their endpoint values to be recovered algebraically from each saved old/new state and the deterministic base step:

\[
s=\frac{(v_n-v_*)\cdot(v-v_*)}{\|v_n-v_*\|^2},\qquad
t=\frac{L-L_*}{L_n-L_*}.
\]

The reconstructed linear-momentum relation has error below 3.4e-16. At the 30 Hz jump, s=188.992 and t=.058406. Across that trajectory, s ranges from -3.418 to 188.992, and t from .019873 to .405234. At 60 Hz, s ranges from -15.637 to 3.367 and t from .002480 to .383205. The fast method allows these extrapolations; it does not clip s or t or certify global proximity to the tentative state.

The original angular constraint is L=L_*+t(L_n-L_*), not L=L_n. Consequently energy accuracy can coexist with angular-momentum drift. The 30 Hz run spans L in [.748714,1.118812]; the 60 Hz run spans [.850124,1.213170]. This changes the radial/angular allocation of energy as well as the phase.

The separate [full fixed-center objective diagnostic](FEPR_FULL_ORACLE.md) examines this same second step in the per-step-target trace. Fast FEPR has original objective value 19.7244. Three full Eq1 local solves reach objective values .00304916, .00300783, and 1.59072, with constraint residuals below 4e-15 and whitened stationarity residuals below 1.7e-10. Their position jumps are .15549, .000163, and 1.86049 metres respectively. This directly confirms that the fast feasible iterate need not minimize the original distance. It also exposes different full-objective stationary branches: a full solve is not, by itself, a validated motion repair, and these local solves do not prove global optimality.

The substep traces retain only outer frames, so their per-internal-step s/t values cannot be recovered from those files. Their large recorded jumps still prove serious path inconsistency: at 120 internal Hz the 2.04731-metre jump appears at outer time 2.53333; at 240 internal Hz the 1.87097-metre jump appears at 1.81667. Recorded minimum radii are .900917 and .901712 respectively. Crossings of r=1 in the saved frames number 39, 46, 45, 46, and 49 for the five table rows; the substep counts are sampling-dependent and do not resolve all internal crossings. Unwrapped angular phase across large jumps is branch-ambiguous, so it is not used as a physical rotation-count metric here.

## Consequence for the comparison

Keep FEPR as a source-faithful, audited control with independent x and v, and report energy, trajectory error, angular momentum, branch jumps, and failures separately. Its completion and small energy residual in these cases do not establish reliable motion or monotone improvement from substepping. The hard-rod results must retain their extension label. A future investigation could test a higher-order tentative integrator, continuation/proximity safeguards, or a controlled auxiliary scaling, but each is an additional algorithmic choice requiring its own evidence; none has been silently added to these results.
