# Numerical findings

Observed **20,877 / 20,877** planned runs. **20,804** completed their full requested duration; **73** did not. Completion is not convergence: **2,674** independent-velocity runs pass the stated balance, stationarity and requested geometric thresholds. **647** completed trajectories have less than 10% of the reference path length.

All selections below require full completion and a certified reference. Each comparison holds case, timestep, duration and stress factors fixed. The balanced score is the maximum of position RMS / 0.01 m, relative energy error / 0.01, and absolute log(path ratio) / 0.1. These are disclosed illustrative scales, not a visual-quality assessment.

## Nominal 120 Hz, 10 second comparisons

| Case | Family | Configuration | Work / velocity target | Position RMS (m) | Max relative energy error | Path ratio | Balance residual | Stationarity residual |
|---|---|---|---|---:|---:|---:|---:|---:|
| scalar_release | xpbd | core_000568 | none / hard_endpoint | 0.141847 | 0.999753 | 0.236351 | 0.00689655 | 1.32394e-14 |
| scalar_release | independent | traces_000013 | start / soft_endpoint | 0.129538 | 0.13574 | 0.946303 | 0.00459135 | 0.00461104 |
| scalar_release | force_momentum | core_000588 | force_average / hard_midpoint | 0.00594121 | 3.62532e-11 | 0.999515 | 4.10978e-12 | 1.82491e-11 |
| radial_orbit | xpbd | core_001218 | none / hard_endpoint | 0.15904 | 0.580084 | 0.808524 | 0.00444995 | 1.33352e-14 |
| radial_orbit | independent | core_001169 | kinetic_start / hard_endpoint | 0.280159 | 5.66119e-08 | 0.891807 | 5.9866e-11 | 1.966e-11 |
| radial_orbit | force_momentum | core_001246 | separate_average / hard_midpoint | 0.00222275 | 6.41184e-06 | 0.999871 | 6.69797e-07 | 1.02445e-10 |
| rod_release | xpbd | core_002347 | none / hard_endpoint | 0.119724 | 0.610929 | 0.795178 | 0.00114873 | 1.31191e-14 |
| rod_release | independent | core_002067 | dg_start / soft_endpoint | 0.0143994 | 0.00730418 | 0.995713 | 0.000443313 | 0.000499323 |
| rod_release | force_momentum | core_002366 | force_average / hard_midpoint | 0.000517995 | 1.04606e-12 | 0.999924 | 3.15237e-14 | 1.07276e-13 |
| rod_release | rattle | core_002403 | none / verlet | 0.000365872 | 0.000187678 | 1.00005 | 4.97461e-06 | 1.69664e-12 |
| elastic_mixed | xpbd | core_002988 | none / hard_endpoint | 0.0941984 | 0.639146 | 0.733835 | 0.00209385 | 1.3669e-14 |
| elastic_mixed | independent | traces_000233 | kinetic_end / soft_endpoint | 0.0537666 | 0.00220506 | 0.972661 | 0.000680921 | 0.0005535 |
| elastic_mixed | force_momentum | core_003008 | force_average / hard_midpoint | 0.00156881 | 0.000213681 | 0.999763 | 1.28952e-05 | 1.89871e-10 |
| chain_mixed | xpbd | core_003641 | none / hard_endpoint | 0.15813 | 0.816207 | 0.605538 | 0.00372874 | 5.29696e-12 |
| chain_mixed | independent | core_003449 | force_average / soft_endpoint | 0.127713 | 0.00374999 | 0.9916 | 0.000777859 | 0.00102274 |
| chain_mixed | force_momentum | core_003669 | separate_average / hard_midpoint | 0.00420078 | 0.000105184 | 0.999589 | 6.60502e-06 | 2.28541e-10 |

Residual columns have method-specific meanings and normalizations; compare stationarity within a family. The independent solver's threshold-pass flag is provided separately in the CSVs. Its small energy error may coexist with a poor trajectory or nonzero optimization residual. An exact work identity does not imply the numerical solve reached its constraints.

## Solver evidence

The dedicated solver tier uses its actual duration from the manifest (0.1 seconds here); it diagnoses short-horizon iteration behavior and is not a long-horizon stability test. The convergence tier uses 64 iterations and its own one-second horizon. The long tier separately tests 60 seconds. Curvature, scheduling and iteration comparisons must retain those distinctions.

True edge blocks include both endpoint kinetic rows and keep shared vertex multipliers. In the fixed DG-average, penalized, soft-endpoint slice shown in solver_residuals.png, joint and overlapping edge blocks reduce balance residuals much more than storage/kinetic split passes on the two-spring chain. More iterations do not consistently repair the split schedule at the budgets tested.

## Scalar objective bias, independent of solve convergence

For release from rest with positive spring extension C_n, write omega²=k/m and v=−h omega² C_n s. Physical acceleration has s→1. Exact energy conservation gives x−x_n=−h² omega² C_n s²/2+O(h⁴). Substituting into the displayed independent-velocity objective gives the following leading nonzero stationary equations for the moving branch:

- Soft endpoint: (beta+1)s²−3s+2(1−beta)=0; its moving minimum has s=[3+sqrt(1+8 beta²)]/[2(1+beta)]. At beta=1, s=1.5.
- Soft midpoint, using the implemented unit-v residual v+v_n−2(x−x_n)/h: (beta+4)s²−6s+2(1−beta)=0; the moving minimum has s=[3+sqrt(1+6 beta+2 beta²)]/(beta+4). At beta=1, s=1.2; at beta=0.1, s≈1.042144.

The landscape experiment confirms this distinction: at 120→240 Hz, scalar-release velocity/reference ratios approach 1.49397→1.49848 for endpoint beta=1 and 1.19740→1.19935 for midpoint beta=1, despite approximately solved objectives and energy balances. Midpoint beta=0.1 gives 1.04118→1.04190. Therefore this particular finite-compliance objective carries a release-acceleration bias, even when solved accurately. These statements apply to the chosen global-energy objective and moving branch; they do not establish a universal limitation of every independent-velocity formulation. The frozen state can also be a local minimum for beta<1, so initialization and branch selection remain relevant.

## Interpretation limits

- Momentum-closed midpoint spring methods and RATTLE are reference candidates supported by matched trajectory evidence, not automatic replacements for the current solver.
- The work_kinematics_quality figure fixes four iterations, joint scheduling, beta=gamma=1, dropped multiplier curvature, and rigid tangency when relevant. It illustrates composition effects at that finite budget; it is not a ranking of converged optimization objectives.
- Candidate evidence is screened against the best certified XPBD position RMS in the same physical cohort. The candidate gate requires both position and velocity accuracy to stay within 5%, path agreement to remain comparable, and at least 20% improvement in either position or energy. Residual passage is reported independently.
- Family medians summarize heterogeneous tested configurations within one physical cohort. They do not establish a universal property of every member of a family.
- Failed prefixes remain in failures.csv and never enter rankings. Pending jobs are identified separately in coverage.csv.
- Timings are compute-only and affected by parallel contention; use a serial benchmark before making performance claims.
- No visual-quality claim is inferred solely from RMS, energy, path length, or residual statistics.
