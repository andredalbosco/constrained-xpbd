# Measured results

The complete run matrix and qualifications are in [README.md](README.md). The animation viewer is [viewer.html](viewer.html).

Executed 20,877 configured runs; 20,804 completed and 73 ended early. A completed run can still be unconverged or nearly frozen.

## A consistent comparison

Nominal parameters, 120 Hz, 10 seconds. XPBD uses four edge sweeps. The midpoint discrete-gradient solver uses at most sixteen Newton iterations with line search and matching momentum closure. These are prototype comparisons at different computational costs.

| Case | XPBD position RMS (m) | Midpoint DG position RMS (m) | Midpoint DG maximum relative energy error |
|---|---:|---:|---:|
| scalar_release | 0.141847 | 0.00594121 | 3.63e-11 |
| radial_orbit | 0.15904 | 0.00222297 | 2.92e-11 |
| rod_release | 0.119724 | 0.000517995 | 1.03e-12 |
| elastic_mixed | 0.0941984 | 0.00159904 | 6.42e-10 |
| chain_mixed | 0.158491 | 0.00420535 | 5.82e-15 |

RATTLE gives 0.0003659 m position RMS and 0.0001877 maximum relative energy error on the same rigid release, while enforcing endpoint tangency. Midpoint DG has smaller energy error but does not enforce that velocity condition.

All 45 selected midpoint-DG stress runs completed, with maximum relative energy error 4.56e-8; all ten 60-second runs completed, with maximum 6.60e-9. This is evidence for the tested parameter range and nonlinear solver, not a general convergence guarantee. Some long coupled references fail refinement certification; their trajectory RMS is excluded from rankings.

## What the experiments rule out

- Energy alone is not a motion criterion. Many independent-velocity and fixed-target projection runs retain energy while shifting orbit or phase substantially.
- Solver residuals alone are not a minimum certificate. The dense scalar search found a nonminimal stationary solution with an objective gap of 0.02563 despite tiny KKT residuals.
- The proposed soft velocity objective has its own bias. At scalar release, beta=1 gives a limiting initial acceleration factor of 1.5 for the endpoint penalty and 1.2 for the midpoint penalty. See [the derivation and measurements](analysis/objective_bias.md).
- Reordering alone does not fix EL. [The audit](audits_notes.md) reproduces a nearly vanishing energy-gradient denominator producing a large correction, and a different ordering conserving energy by freezing.

## Numerical techniques for later

[The primary-source research notes](paper_ideas.md) separate tested algorithms from untested transfers: Muon-inspired bounded spectral balancing; PPO/TRPO-inspired limits on local-model movement; Large Steps-inspired graph preconditioning; explicit-work augmented Lagrangian splitting; feasibility/stationarity filters. They target convergence and correction quality after the mechanical formulation is chosen.

## Validation

Three numerical validation executables pass. Derivative, work identity, impulse closure, incidence sharing, rod constraints, temporal convergence, and objective-landscape checks are recorded. The viewer was checked in Chrome: playback, environment switching, RATTLE selection, scrubbing, long-run switching, and failed-prefix overlays work; the final page reports no console errors. Screenshots are in output/playwright/.
