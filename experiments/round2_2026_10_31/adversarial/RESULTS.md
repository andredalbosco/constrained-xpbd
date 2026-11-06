# Adversarial mechanics and pendulum compositions

1,110 trajectory configurations: 1049 complete, 61 failed. A completed trajectory may contain unconverged steps; these remain explicitly counted.

42 rows have uncertified references and are excluded from trajectory ranking. All pendulum references pass the two-resolution check.

## Pendulum just above the separatrix

Initial angle 0, angular speed 1.001 × 2√9.81, unit mass and length. Ten seconds at 30 Hz. This is deliberately sensitive to small energy changes.

| Method | Position RMS (m) | Max relative energy error | Max scaled tangency | Status |
|---|---:|---:|---:|---|
| rattle | 1.29 | 0.00577 | 2e-16 | complete |
| dg64 | 0.00567 | 8.37e-11 | 0.00297 | complete |
| dg_tangent | 0.422 | 0.00106 | 1e-16 | complete |
| dg_tangent_energy | 0.00576 | 1.45e-15 | 2e-16 | complete |
| sym_geom | 0.98 | 0.00546 | 1.33e-11 | complete |
| sym_energy | 0.00266 | 1.18e-10 | 4.28e-12 | complete |
| post_energy | 0.00283 | 6.53e-13 | 7.47e-14 | complete |
| angle_dg | 0.00862 | 4e-10 | 1e-16 | complete |
| rattle_accuracy | 0.000589 | 2.67e-06 | 1e-16 | complete |

RATTLE is a useful structural baseline, but at30 Hz the just-below-separatrix case changes from reference libration into numerical rotation (the candidate angle reaches22.37 rad while the reference stays between−3.053 and3.053 rad). Symmetric projection and the angular discrete-gradient control substantially reduce that specific error. This does not establish superiority at equal cost.

## Finite convergence and extreme steps

| Case / h | Method | Max relative energy error | Unconverged accepted steps | Rejected attempts |
|---|---|---:|---:|---:|
| radial_stiff / 0.033333 | dg64 | 6.61 | 1 | 0 |
| radial_stiff / 0.033333 | dg_adaptive | 6.68e-08 | 0 | 43 |
| elastic_stiff / 0.25 | dg64 | 2.32 | 3 | 0 |
| elastic_stiff / 0.25 | dg_adaptive | 5.04e-09 | 0 | 8 |
| rod_fast_spin / 0.25 | dg64 | 1.4e-10 | 0 | 0 |
| rod_fast_spin / 0.25 | dg_adaptive | 1.4e-10 | 0 | 0 |

Residual gating rescues energy when an incomplete nonlinear solve is the cause. It does not certify temporal accuracy: the fully converged fast-spin Cartesian DG trajectory still has 1.41 m position RMS at h=0.25. Increasing Newton iterations cannot repair the wrong temporal resolution.

## The inherited radial mode

For a rigid rod, exact old/new length and Cartesian midpoint kinematics imply $n^T(v_{n+1}+v_n)=O(h)$. If the old velocity contains a finite normal component, the small-step limit is a reflection $v_{n+1}\to v_n-2n(n^Tv_n)$. This preserves kinetic energy. Two half steps reflect twice, while one full step reflects once. The resulting velocity discrepancy tends to $2|n^Tv_n|$, rather than zero.

The probe starts from a genuine DG output, not only an artificial off-manifold input. After h=0.01 from pendulum release, radial velocity is 8.8082e-6. At the next trial h≈2.03e-6 it becomes −8.8098e-6; the scaled step-doubling discrepancy is 1.2459e-5 against tolerance 2.03e-10, while the solver residual is 3.03e-11.

The accuracy controller fails on 24/28 ten-second gravity-pendulum configurations; both free-spin cases at each step size work because they do not generate this radial defect. Two additional stiff-spring configurations hit the controller depth limit. Those latter failures are not attributed to the rigid hidden mode.

## Projection and reversibility

The symmetric projection experiment follows Hairer’s same-multiplier before/after construction around an implicit midpoint step of the acceleration-level pendulum ODE. Position, tangency, and optionally energy are enforced together in fixed dimensionless coordinates. Ordinary postprojection is a separate control. [Primary source](https://www.unige.ch/~hairer/preprints/symproj.pdf).

For release at 30 Hz, symmetric energy projection has local reversal discrepancy 1.37e-15; postprojection 2.11e-6; DG followed by tangent/energy repair 2.30e-4. RATTLE plus energy-only velocity rescaling encounters a backward step whose position requires negative kinetic energy. A small energy error or exact tangency does not by itself imply reversibility.

Angular DG uses the same potential in angle coordinates and midpoint angular velocity. Rod length and tangency hold by construction; the exact potential secant closes energy. It intentionally changes Cartesian midpoint kinematics. All 42 angle-control trajectories complete, including h=0.25 fast rotation; that case has 0.161 m RMS versus 1.41 m for Cartesian DG. This is a single-pendulum coordinate control, not a general chain solver.

## Numerical limits

- Relative energy uses initial excitation energy above the hanging minimum; it is not divided by signed total energy.
- Reversal compares one forward/backward step from each of the first 32 states. The CSV sentinel 1e9 means a backward solve failed, not a measured discrepancy of that magnitude.
- Midpoint and DG impulse defects are macro-step diagnostics; adaptive compositions obey those equations at each accepted substep, not necessarily at the macro-step.
- Accuracy subdivision accepts two half steps when their scaled discrepancy from one full step is below 1e-4 |h|. It has a depth cap14 and keeps the two-half-step state without extrapolation.
- The residual controller subdivides only on solver/geometry failure, to depth10. State-dependent subdivision is not claimed time reversible.
- Timings here were collected concurrently and include the first32 backward probes; use them only as diagnostic workload, not performance rankings.
- This experiment does not cover collisions, friction, damping, external time-dependent forcing, or rigid chains.

![New failure mechanisms and controls](figures/new_failures_and_controls.png)
