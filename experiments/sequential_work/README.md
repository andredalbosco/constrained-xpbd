# Local reservoir update ordering: measured results

This benchmark compares ordinary XPBD with three reservoir update schedules over 10 seconds. The sequential version changes each spring by material lambda, then potential-work nu, then its two particle-owned kinetic mu multipliers. It improves some strongly coupled cases but does not conserve mechanical energy, linear momentum, or angular momentum. Weaker feedback approaches the dissipative XPBD baseline.

The [trajectory behavior study](behavior/README.md) adds detailed motion histories, the kicked 2D grid, cloth, and isolated momentum checks across feedback strengths and timesteps.

## Reproduce

From the repository root:

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --target xpbd-ordering-study -j 4
mkdir -p results/sequential_work
(cd results/sequential_work && ../../build/xpbd-ordering-study)

# Run the production solver with the new order and independent compliances.
./build/constrained-xpbd --scene spring --mode active --order sequential --beta 1000 --gamma 1000
```

See [FORMULATION.md](FORMULATION.md) for the scalar updates, compliance units, and shared-particle multiplier details.

The output is `ordering_results.csv`. An optional positive duration argument changes the default 10 seconds. All measurements below use 10 seconds.

## Schedules and coverage

- `paired`: jointly update (lambda, nu) for every spring, then update every particle mu once per sweep.
- `paired_local`: jointly update (lambda, nu), then update both endpoint mu values immediately for each spring.
- `sequential`: update lambda, relinearize and update nu, then update both endpoint mu values immediately.
- A mu belongs to a particle and uses the frozen force sum from all incident springs. It is shared across spring visits. Local schedules therefore revisit a triangle particle twice per sweep; `paired_local` is the matched-schedule comparison for `sequential`. All schedules update isolated particles once per sweep.
- Positive beta and gamma are the existing kinetic and potential reservoir compliances. Larger values weaken feedback; zero is not tested or supported.

The 684 cases cover four scenes, h = 1/240, 1/60, 1/15, and 1, 5, 20 sweeps. At each combination there is one XPBD baseline and three active orders for (beta, gamma) = (1,1), (10,10), (1000,1000), (1,1000), (1000,1), (1e6,1e6). The free pair has unequal masses and translation plus rotation; the triangle has three springs sharing particle multipliers. Exact initial states are in the source.

All cases completed the requested duration with finite numeric CSV cells and exactly stationary pins. Status counts are 591 `completed` and 93 `large_accepted_excursion`. The latter only means an accepted position norm exceeded 10 or speed exceeded 100; it is a magnitude flag, not a claim of divergence. In a freely translating scene it also depends on the origin. `completed` does not imply row convergence or stationarity.

## Representative measurements

Here h = 1/60 and there are five sweeps per step. Energy error is max over the trajectory of |E(t) - E(0)|, in energy units; it includes gravitational potential for the pendulum. The energy zero is arbitrary, so absolute error is used below. The single-spring scenes produce bit-for-bit identical `paired` and `paired_local` CSV metrics for every setting; those duplicate rows are omitted.

| Scene | Order | Max energy error, beta=gamma=1 | Max energy error, beta=gamma=1000 | Max stationarity residual, beta=gamma=1 |
|---|---|---:|---:|---:|
| spring | xpbd | — | 3.125 | — |
| spring | paired | 25.0761 | 3.125 | 3179.65 |
| spring | sequential | 14.8686 | 3.125 | 181.023 |
| pendulum | xpbd | — | 1.67294 | — |
| pendulum | paired | 14.9731 | 1.23398 | 246422 |
| pendulum | sequential | 2.33749 | 1.23398 | 2.11956 |
| free_pair | xpbd | — | 3.10546 | — |
| free_pair | paired | 3.02478 | 2.58266 | 1.4418 |
| free_pair | sequential | 3.03901 | 2.58265 | 0.503471 |
| free_triangle | xpbd | — | 1.69238 | — |
| free_triangle | paired | 2.99815 | 1.136 | 0.405578 |
| free_triangle | paired_local | 2.97276 | 1.13592 | 0.635495 |
| free_triangle | sequential | 3.05477 | 1.13592 | 0.627313 |

The stationarity diagnostic is the Euclidean norm of the physical position KKT residual: `M(x-y)/h^2 - sum_j grad(C_j) lambda_j/h^2 - sum_j grad(d_j) qnu_j/h^2 - sum_i grad(b_i) qmu_i/h^2`, restricted to free coordinates. These methods drop that residual during updates. Closing the material and feedback rows therefore does not close stationarity: for the spring at 20 sweeps, beta=gamma=1, `paired` has maximum material/feedback residuals 4.44e-16 / 3.46e-14 but stationarity residual 3183; `sequential` has 1.11e-16 / 2.17e-14 but stationarity residual 181. More sweeps do not remove this discrepancy.

The strong pendulum case trades failures: at five sweeps, paired has maximum length error 0.0135078 and 31 rejected updates; sequential has length error 0.0343605 and no rejected updates. With beta=gamma=1000, both have length error at roundoff and very similar trajectories. At beta=1, gamma=1000, the spring ends with energy gain 3315.2 (paired) or 2661.7 (sequential): splitting the updates is not a general stability fix.

For the free systems only, linear momentum is an invariant of the intended isolated dynamics. At the same h and sweep count:

| Scene | Order | Max momentum change, beta=gamma=1 | Max momentum change, beta=gamma=1000 | Max angular momentum change, beta=gamma=1 |
|---|---|---:|---:|---:|
| free_pair | xpbd | — | 4.0589e-13 | — |
| free_pair | paired | 5.3393 | 1.52191 | 4.00848 |
| free_pair | sequential | 5.26242 | 1.52198 | 2.63529 |
| free_triangle | xpbd | — | 1.1736e-12 | — |
| free_triangle | paired | 2.05273 | 0.345111 | 3.52743 |
| free_triangle | paired_local | 2.15479 | 0.345167 | 3.56586 |
| free_triangle | sequential | 2.15277 | 0.345166 | 3.17071 |

Particle kinetic corrections need not sum to zero, so endpoint-local scheduling does not restore momentum conservation. At beta=gamma=1e6 the free pair has momentum error 0.0010648, approaching but still above the ordinary XPBD value 4.06e-13. Angular momentum is also not conserved; the CSV includes its actual error for every case. P and L changes in anchored/gravity scenes are diagnostics, not conservation failures of an isolated system.

## What was verified

- Mandatory checks compare 600 ordinary XPBD spring steps against the analytic backward-Euler recurrence, check the stationary pin, and check ordinary XPBD free-system linear momentum. These pass before the sweep runs.
- An independent ledger reconstruction uses the previous accepted multiplier divided by the previous timestep squared, together with the previous accepted spring direction. It separately checks every particle and spring reservoir update.
- Across all active rows, the maximum independent local ledger residual is 1.4552e-11 and maximum extended-energy drift is 2.4663e-10. This closes E + sum(R) + sum(S), including gravity potential. It does not close mechanical E: reservoirs are explicitly assigned from accepted residuals even at finite iterations.
- The source checks finite state and diagnostics before accepting a step, reports exceptions/partial durations, and checks output-stream errors. All 684 cases finished without exceptions. Reservoirs and previous force history are carried between steps.
- ASan+UBSan passed the mandatory checks and all 684 settings over a 0.2-second smoke sweep with `ASAN_OPTIONS=detect_leaks=0`. LeakSanitizer is disabled because this execution environment uses ptrace. The full 10-second sweep was run with the optimized ordinary binary; it was not repeated under sanitizers.

## Integration checks

All 11 CTest targets pass, including seven independent sequential-work test suites. Those suites check the literal one-sweep order, shared particle multipliers, fixed and isolated particles, separate compliance responses, the weak-feedback limit, and accepted-state local accounts. The sequential tests also pass ASan+UBSan with leak detection disabled.

The original binary and updated binary produced identical summary fields (excluding wall time) for all four base scenes in XPBD, passive, and default paired-active modes over 2400 steps. New sequential runs also completed those four scenes for ten seconds at h=1/240, 20 sweeps, beta=gamma=1000:

| Scene | Maximum mechanical energy error | Maximum extended-energy error |
|---|---:|---:|
| spring | 6.43759 | 2.31e-14 |
| pendulum | 0.202757 | 2.64e-14 |
| 2D spring grid | 0.0302241 | 3.96e-14 |
| cloth | 10.9410 | 2.59e-13 |

These are finite-run observations, not stability or conservation guarantees. The timestep and sweep count differ from the representative ordering tables above.
