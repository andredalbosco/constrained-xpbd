# Independent validation pass

Validated spring_study.cpp, pendulum_study.cpp, and the current midpoint_gs_study.cpp / midpoint_gs.hpp. Existing study outputs were preserved. Sanitizer artifacts were isolated in /tmp/xpbd-work-research.X6BItq/validation_spring, validation_pendulum, and validation_midpoint_gs.

## Sanitizers and optimized builds

All three full studies passed builds/runs using:

    c++ -std=c++17 -O1 -g -fsanitize=address,undefined -fno-omit-frame-pointer ...
    ASAN_OPTIONS=detect_leaks=0 ./study_binary

The spring and pendulum studies also passed their optimized -O2 builds and built-in numerical assertions. The pendulum build used -Wall -Wextra -pedantic after the finite-value guard below was added.

AddressSanitizer and UndefinedBehaviorSanitizer reported no findings. Default LeakSanitizer could not run because the execution environment uses ptrace; it reported that environmental limitation at process exit. Leak checking was therefore disabled explicitly. This is not evidence of a program leak, and no leak-validation claim is made.

## CSV integrity

Independent Python parsing verified consistent row widths and finite numeric fields:

| Output | Rows | Columns |
|---|---:|---:|
| spring_linear.csv | 24 | 20 |
| spring_radial.csv | 12 | 14 |
| spring_degeneracy.csv | 4 | 9 |
| spring_linearized.csv | 16 | 17 |
| pendulum_summary.csv | 48 | 22 |
| pendulum_trajectories.csv | 25,710 | 16 |
| pendulum_singularities.csv | 14 | 10 |
| midpoint GS sanitizer output | 216 | 17 |

The only nonfinite CSV values are 12 explicitly written NaN placeholders in the two expected pendulum singularity-error rows. These are the exact critical-step cases: the continuum degeneracy with zero initial speed and the antipodal/no-finite-multiplier case with nonzero speed. All 48 ordinary pendulum summary rows have status ok.

Sanitizer and optimized numerical outputs agreed bit-for-bit for the linear, degeneracy, linearized-spring, pendulum, and 216-row GS tables. The nonlinear radial spring tables differed by at most 1.12e-16 in a numeric field.

## Independent numerical checks

- The backward-Euler scalar energy matched the analytical decay \(E_N=E_0(1+\omega^2h^2)^{-N}\) with a maximum relative discrepancy of 5.21e-13.
- All 48 pendulum summaries satisfied the independently checked accumulated identity: mechanical-energy change equals reported constraint-impulse work plus reported gravity-work mismatch. Maximum discrepancy: 2.67e-15.
- The 216-row GS study completed all assertions. Maximum residual-to-kinetic-work identity error: 3.43e-13.
- For retained position residual, h=1/60, outer limit 80, and one inner GS sweep, all timesteps converged in each of the spring, pendulum, and two-spring runs. Maximum absolute energy errors were respectively 6.50e-13, 1.28e-13, and 1.14e-11.
- The conservative Cartesian pendulum still had maximum endpoint radial velocity 4.24280e-4 at h=1/60. Its near-zero energy error is not evidence of endpoint tangent-velocity preservation.

## Concrete guard fixed

The pendulum summaries accumulated maxima using std::max and checked them with greater-than thresholds. Without an explicit finite-value check, NaNs could be silently omitted from those maxima. The current runs had no such values, but this was a real future-test weakness.

At the root agent's request, advance() now rejects nonfinite position/velocity components, squared norms, mechanical energy, and every numeric work/solver metric before returning a state. Optimized and sanitizer full reruns passed after this change.

The root GS header already checks finite row data, residual totals, and accepted position/velocity state. The integrated header also rejects nonfinite StepReport work/identity fields. All seven Release CTest targets and the five-program result runner passed after this guard was added.

## Interpretation limits

- Energy, momentum, tangent velocity, geometric constraint, and solver residuals are separate checks. No claim was inferred merely from a small bookkeeping identity.
- The midpoint pendulum's large-step singularity probes are intentional failure/branch diagnostics, not successful physical-accuracy tests.
- The angular-coordinate DG pendulum reconstructs tangent Cartesian velocity, but generally does not satisfy Cartesian midpoint drift. Its inferred Cartesian non-gravity impulse is diagnostic bookkeeping; it is not automatically a radial force quadrature.
- Sanitizers test memory and undefined behavior. They do not establish mathematical correctness, convergence for arbitrary systems, or stability outside the tested cases.

## Integrated regression run

All seven CTest targets passed: baseline/passive, active reservoirs, unchanged-XPBD work audit, midpoint GS, spring references, pendulum references, and coupled momentum. The runner generated 216 midpoint-GS cases and 144 free coupled-momentum cases. Independent parsing verified the expected row counts and finite numeric fields in all main summary tables. The coupled-momentum study reused the independently reviewed and sanitized midpoint solver; its full run was validated in Release, without a separate coupled-fixture sanitizer run.
