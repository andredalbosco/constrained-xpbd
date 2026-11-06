# Work identities and EL/LE audit

These are fixed-motion identity probes and short solver diagnostics, not evidence that an exact work identity gives accurate or stable dynamics. Only new files in this experiment directory were written; the existing constrained-xpbd checkout is used read-only.

## Reproduce

From this directory:

```bash
g++ -std=c++17 -O2 -Wall -Wextra -Wno-misleading-indentation \
  -I/home/blacknoir/Projects/constrained-xpbd/external/eigen \
  -I/home/blacknoir/Projects/constrained-xpbd/experiments \
  -I/home/blacknoir/Projects/constrained-xpbd/experiments/gradient_policies \
  -I/home/blacknoir/Projects/constrained-xpbd/experiments/energy_ordering \
  audits.cpp -o audits
./audits audits_results
```

The suppressed indentation warning comes from compact existing/common headers; no other compiler warnings were reported in the audit build. Output is deterministic with seed 42. The executable returns nonzero if any analytic identity's normalized error exceeds 1e-11.

## Fixed-motion audit

There are 1,000 probes for each of five environments: affine scalar spring, planar distance spring, rigid pendulum, elastic pendulum, and a two-mass serial spring chain. Spring stiffness is 100; lengths range from 0.6 to 1.4; planar directions can rotate by up to 2.8 radians. These deliberately broad displacements test identities, not time integration accuracy. Rod endpoints have unit length, old velocity is tangent, and old/new signed reaction strengths vary independently. Gravity work is separately checked with its constant force.

Ten work laws are evaluated: old-force spatial work, endpoint-force spatial work, average complete force, separately averaged strength and direction, discrete gradient with old/endpoint/averaged strength, entirely frozen start power, and old/endpoint-force estimates of kinetic work. Here the kinetic estimates mean f^T[(x-x_n)+h v_n]/2, using endpoint position-derived velocity. They are not asserted to be kinetic-exact.

`fixed_work_summary.csv` reports max absolute and normalized **storage** defects, with normalization 1+|delta U|+|W|. A nonzero entry for a kinetic estimate means it is not potential-exact; it is not a test of its intended kinetic approximation. `fixed_work_incidence.csv` retains each endpoint's work before summing over a spring. In the chain, an internal spring's work delivered to one endpoint can be nonzero even when its net storage change is zero.

Results:

- Averaged-strength discrete-gradient work closes every quadratic spring's storage balance to at most 1.34e-14 absolute in these probes, including incidence sums in the chain.
- For feasible rigid endpoints, separate averaging, all three discrete-gradient strength choices, and frozen tangent start power give zero rod work to floating-point precision (maximum 2.40e-15 for geometric work).
- Average complete force is not zero-work for a rigid pendulum when its reaction strength changes: the maximum absolute work here was 18.35. Old and endpoint spatial work also fail that geometric test.
- For the affine scalar spring, full-force averaging, separate averaging, and averaged-strength discrete work coincide. Endpoint spatial and endpoint discrete work coincide.
- Old-force spatial work still changes with the solved displacement: its change is f_n^T delta x. The entirely frozen start power h f_n^T v_n is independent of the new position. Both identities are checked separately.

The discrete-gradient potential identity assumes U=k C^2/2 and q_n=-k C_n, q=-k C with fixed k. The rod reaction is independent; it is not obtained from a zero rod elastic energy.

## Kinetic bookkeeping audit

`kinetic_impulse_probes.csv` contains 1,000 random two-impulse probes. With m(v-v_n)=I_a+I_b, using the common full-step average velocity gives exact kinetic work. After a deliberately unlogged velocity correction, the residual obeys

    delta K - sum_j I_j^T (v+v_n)/2
        = [(v+v_n)/2]^T [m(v-v_n)-sum_j I_j].

This is checked numerically. The maximum nonzero defect after the unlogged correction was 0.655 in these units. If the previously recorded work itself is left unchanged, the defect equals the correction's actual additional kinetic-energy change; this is a separate analytic check.

Sequential per-impulse work also telescopes exactly, but attribution depends on order. Applying A after B rather than before B changes A's work by I_a^T I_b/m. The common full-step midpoint allocation gives half of this cross term to each impulse. The largest order-dependent allocation change here was 6.98, despite identical final velocity and total work. Two-mass chain probes separately verify both bodies' kinetic balances with the shared internal impulse applied equally and oppositely.

## Scalar storage and release obstructions

For an affine scalar spring, prescribed old-force storage gives

    delta U + W_old = +k (C-C_n)^2/2,

and prescribed endpoint-force storage gives the negative of that expression. Therefore enforcing either storage balance exactly requires C=C_n. These are algebraic restrictions, not merely a solver failure.

`scalar_release_analytic.csv` checks the separate endpoint-kinematics obstruction from h=1/30 through 1/960. At rest, with s=k h^2/m, exact total energy has roots C=C_n and C=C_n(1-s)/(1+s). The nontrivial root approaches four times the physical position displacement and twice the physical velocity over the first release step. The backward-Euler trial C_n/(1+s) lies at zero energy gradient with a nonzero negative energy residual. Energy exactness therefore does not select the physical step.

## EL/LE implementation and diagnosis

`el.hpp` exposes the shared `el_step` interface. Dropped-r_x frozen/refreshed sweeps call the existing traditional implementations directly. Retained-r_x sweeps use the same existing row evaluation and mass-block Schur approximation with the residual correction included. Refreshed, raw, one-sweep retained-r_x results agree with the existing research-agenda implementation; frozen first sweeps agree with the existing dropped-r_x solver because their initial r_x is zero. Both comparisons are checked over all ten initial conditions and both schedules.

The 480 short compositions cover ten initial conditions across five environments, EL/LE, raw/bounded correction, retained/dropped r_x, per-sweep frozen/per-row refreshed linearization, and 1/4/16 sweeps. Each runs up to 24 steps at h=1/120 (0.2 seconds). The adapter's energy target is the **old energy of each step**, matching the new local balance experiments. Some historical studies instead target the run-initial energy; those trajectory results are not identical by definition.

Bounded correction retains the historical energy-repair criterion and 0.1 motion-scale cap, applied per sweep. It accepts a partial correction and does not exactly solve the constrained optimization. At multiple retained-r_x sweeps, the material baseline is recomputed from that sweep's current iterate and residual; this is an explicit extension of the one-sweep test.

`EL_scalar_release_visit_trace.csv` isolates the failure mechanism. For the raw one-sweep frozen EL scalar release at h=1/120:

1. The first energy visit sees zero residual and does nothing. The material visit takes the backward-Euler step and loses 0.02155 energy.
2. On the next predictor, the energy residual is -0.04266 while its gradient is only -0.17241. The raw energy row requests a position correction of -0.24741.
3. That step adds about 443.78 energy. The short trajectory subsequently crosses the scalar adapter's positive-anchor domain.

The energy correction is -H times the mass-inverse energy gradient divided by its squared mass-inverse norm. A small gradient with finite residual produces a large correction. Refreshed LE encounters the exact zero-gradient release state after its material visit and skips that singular energy row. Frozen raw LE instead cancels the release displacement exactly: its zero energy error accompanies a motionless spring. Increasing sweep count, retaining r_x, or bounding the correction does not uniformly resolve these distinct problems; the complete CSV retains failures and motion measurements.

The scalar EL adapter uses the existing distance-spring solver on its positive-anchor branch. Its shared step interface rejects a nonpositive predictor or completed-sweep position with an explicit domain failure. The short standalone diagnostic also stops on nonpositive final scalar position. It is not claimed to implement the affine scalar spring through an anchor crossing.

All **45 analytic/compatibility identity groups passed**. The largest normalized error was 2.26e-13, from cancellation when evaluating the theoretically zero release gradient at small h. The fixed-motion guarantees are strong identity checks; the 0.2-second EL trajectories are deliberately short diagnostics. Longer dynamics, reference errors, and composition rankings belong to the shared runner.
