# Adversarial mechanics and pendulum compositions

This is a separate extension of `composition_2026_10_31`, with its model and mechanical kernels reused unchanged. `study.cpp` adds explicit convergence checks, subdivision, new projections, and a pendulum coordinate control. Read `RESULTS.md` for measured outcomes; `summary.csv` files retain every configured trajectory, including failed prefixes and finite unconverged steps.

## Mechanical questions

1. Does a finite Newton iterate actually satisfy the energy-conserving mechanical equations?
2. Does satisfying those equations suffice for accurate motion at a large step?
3. Can energy, rigid position, endpoint tangency, and time reversibility be composed?
4. Does adaptive step-doubling behave sensibly when Cartesian DG leaves an endpoint radial velocity?

The twenty cases include the previous ten initial conditions, pendulums just above/below the separatrix, fast rotation, a near-upright release, gravity-free rotation, a spring compressed to length0.02, high spring velocity, stiffness100,000, an elastic pendulum of stiffness10,000, and a two-spring chain with mass ratio100. Ten-second runs use h=1/120,1/30,1/10,1/4. Sixty-second runs use the seven rigid cases at1/120 and1/30. These are finite stress tests, not a stability theorem.

## New compositions

- **Cartesian DG:** midpoint position/velocity and momentum closure with the exact distance discrete gradient; fixed4/16/64 Newton budgets. A finite but unconverged state is deliberately retained and flagged.
- **Residual-gated subdivision:** reject a failed or insufficiently converged step, then take two half steps recursively. It does not estimate temporal truncation error.
- **Accuracy subdivision:** compare one step with two half steps; accept the latter when the scaled state discrepancy is below1e-4|h|. Maximum depth14. The accepted state is not Richardson-extrapolated. This exposes rather than suppresses inherited radial-mode failure.
- **Tangent postprojection:** project velocity into the endpoint rigid tangent space. Optional rescaling then enforces the previous step’s energy. A negative kinetic-energy target or absent velocity direction is an explicit failure. RATTLE plus energy rescaling is a separate control.
- **Implicit midpoint ODE and projection controls:** for the single rod, extend acceleration away from its constraint manifold as

  \[
  \dot x=v,\qquad
  \dot v=g-x\frac{\|v\|^2+x^Tg}{\|x\|^2}.
  \]

  Projection constraints are rod length, endpoint tangency, and optionally energy. Coordinates are fixed dimensionless \(y=(x/\ell,v/v_s)\), with \(v_s=\sqrt{E_s/m}\) from the initial excitation energy. This metric is fixed across forward/backward steps.

- **Symmetric projection:** following [Hairer’s primary paper](https://www.unige.ch/~hairer/preprints/symproj.pdf), use the same correction multiplier before and after the implicit midpoint step. Solve the integration equations and invariants jointly. Our experiment adds energy to the geometry/tangency manifold and compares ordinary postprojection under the same metric. The finite solver and chosen nonlinear branch remain part of the observed result.
- **Angular DG:** write the single pendulum as \(x=\ell(\sin\theta,-\cos\theta)\), \(v=\ell\omega(\cos\theta,\sin\theta)\). For \(d=\theta_{n+1}-\theta_n\),

  \[
  a_d=\frac{g_x\cos(\theta_n+d/2)+g_y\sin(\theta_n+d/2)}{\ell}
      \operatorname{sinc}(d/2),
  \quad d=h\omega_n+\frac{h^2}{2}a_d,
  \quad\omega_{n+1}=\omega_n+h a_d.
  \]

  This is the exact potential secant combined with midpoint angular velocity. It preserves position and tangency by construction, while changing the Cartesian displacement relation. It is deliberately a single-pendulum control, not a proposed general solution for arbitrary constraints.

The RATTLE control follows the position/velocity constraint distinction of [Andersen’s original paper](https://doi.org/10.1016/0021-9991(83)90014-1). Its old-normal position solve may have no real root when a tangential trial displacement exceeds the rod length; a finite iteration is not evidence of successful constraint enforcement.

## Validation and interpretation

`adversarial_study validate` checks quadratic-spring equivalence with implicit midpoint, analytic invariant derivatives, four projection variants, symmetric reversal, the missing-RATTLE-root case, adaptive recovery, and angular-DG energy/reversal. `radial_mode_probe` records the hidden-mode limit from both a true previous DG output and an injected normal velocity.

References use RK4 with two substep counts, choosing at least64 coarse substeps and more for large hω. Both trajectories must be finite. Certification requires maximum scaled position+velocity difference<1e-5 and relative reference energy drift<1e-6. All seven pendulum references pass. Seven ten-second spring/chain scenarios fail the refinement threshold; their42 method rows retain a visible `reference_certified=0` and are not used for trajectory ranking.

The energy scale is initial excitation energy above the hanging minimum, floored at1e-4. Root residuals for DG are normalized by the velocity scale. RATTLE length and tangency residuals are normalized separately by length and velocity. Projection equations are solved in dimensionless coordinates. Thresholds therefore describe each stated residual, not identical raw algebra across methods.

CSV fields explicitly retain: status/reason, completed/requested steps, reference certification, RMS trajectory errors, energy and geometry, macro-step midpoint/impulse compatibility, solver residual, local reversal, angular change, nonconverged accepted steps, solver evaluations, rejected attempts, and accepted substeps. Angular change is only an invariant diagnostic in gravity-free cases. A reversal value1e9 is the failed-backward-solve sentinel.

The coarse-step compatibility diagnostics generally differ from zero for a composition of accepted substeps; they do not measure a failure of those substeps. Concurrent wall time includes reversal probes while forward evaluation counts do not. Do not interpret these times as a solver speed ranking.

## Reproduce

From the round2 experiment directory, configure and build using the top-level CMakeLists. Then:

```bash
OMP_NUM_THREADS=8 ./build/adversarial_study adversarial/results
OMP_NUM_THREADS=4 ./build/adversarial_study adversarial/long_results long
OMP_NUM_THREADS=4 ./build/adversarial_study adversarial/angles_results angles
OMP_NUM_THREADS=4 ./build/adversarial_study adversarial/angles_long_results angles_long
./build/radial_mode_probe > adversarial/radial_mode.csv
MPLCONFIGDIR=/tmp/xpbd-round2-plots python3 adversarial/analyze.py
```

Analysis needs NumPy and Matplotlib. The local workstation run used `/home/blacknoir/Projects/mandelbrot-diffusion/.venv/bin/python` for those packages. C++ uses Eigen from the constrained-XPBD checkout. No GPU is needed for these small mechanical problems.
