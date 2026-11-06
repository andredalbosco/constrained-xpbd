# Matched row-evaluation ceilings

`fair_cost.cpp` holds the existing objective, work target, initialization, beta/gamma, and solver definitions fixed. It tests `chain_mixed`, `chain3`, and `chain8` at `h=0.02,0.08`, with three methods, six schedules, and row-evaluation ceilings `100,300,1000,3000`: **432 comparisons**. No engine changes were made for this experiment. All six dense reference roots converged, had full row rank, and passed positive reduced-curvature certification.

The operation-budget comparison changes the ranking seen at equal scalar-visit counts. A local stationarity gather evaluates multiple incident rows; treating one such visit as equal cost to one dropped-stationarity or frozen visit overstates its efficiency. At a ceiling of 1000 evaluations, frozen GS plus restoration has the smallest full KKT residual in **30 of 36** matching case/timestep/schedule groups; at 3000 it also wins **30 of 36**. These are operation-count comparisons, not wall-clock speed claims.

For `chain_mixed`, `h=0.08`, material→storage→kinetic, ceiling 1000:

| Method | Actual row evaluations | Scalar visits | Coordinate writes | Full KKT residual | Distance to oracle state |
|---|---:|---:|---:|---:|---:|
| Refreshed, dropped stationarity | 1000 | 1000 | 7250 | 2.90e-3 | 3.52e-3 |
| Refreshed, local stationarity | 996 | 146 | 1058 | 2.94e-5 | 1.19e-4 |
| Frozen GS, full restoration | 1000 | 1000 | 8875 | 4.02e-11 | 1.57e-10 |

The dropped method's energy error here is only `9.44e-12`; its residual is stationarity. The local method's remaining error is primarily storage feasibility; stationarity is `4.73e-9`. Frozen restoration uses more coordinate writes than either alternative and 125 full restoration passes. Its extra arithmetic is not charged as additional row-gradient evaluations, so this ceiling alone is not a complete work metric.

Across all 36 cases/orders at each ceiling:

| Ceiling | Method | Median actual visits | Median full KKT residual | Residual below 1e-6 |
|---:|---|---:|---:|---:|
| 1000 | Dropped stationarity | 1000 | 6.03e-3 | 0/36 |
| 1000 | Local stationarity | 137 | 4.18e-4 | 1/36 |
| 1000 | Frozen restoration | 1062.5 | 8.08e-5 | 12/36 |
| 3000 | Dropped stationarity | 3000 | 6.00e-3 | 0/36 |
| 3000 | Local stationarity | 399 | 1.15e-4 | 10/36 |
| 3000 | Frozen restoration | 3187.5 | 1.49e-5 | 12/36 |

There is no universal winner. At `chain8,h=0.08,ceiling=3000`, the material→storage→kinetic full residuals are `0.0161` dropped, `0.00635` local, and `0.193` frozen. Frozen GS can continue to make substantial nonlinear errors even with many cheap frozen visits. Local stationarity reaches the `coordinate_limit` guard for the same eight-mass state under `edge` and `edge_reverse` at ceiling 3000; those two failed prefixes are retained in the data, not replaced by their last successful prefix. No trajectory claim follows from these same-state steps.

## Budget certification and limitations

Deterministic doubling and binary search select the largest successful visit prefix within the ceiling, with a cap of twice the ceiling. Every candidate restarts from the same initial state. Whole frozen-sweep preparation is charged upfront. All **430 nonfailure rows** have a recorded next prefix whose evaluation count exceeds the ceiling; the visit cap was never active. Both guard-failure rows identify the first failing prefix, together with the largest successful prefix. All 432 selected rows obey their evaluation ceiling. Cached count monotonicity is checked through the first failure, and failures stop search expansion.

CSV fields also record coordinate writes, restoration passes, selected solve time, setup time, and total replay-search time. Replay search is diagnostic overhead, not a proposed runtime algorithm. Final global diagnostics are excluded from the engine's operation counters and selected solve time, as in the main study. The fixed start-only predictor has its own setup solve. Different row families and incident gathers need different arithmetic, and these timings were not isolated; no wall-time performance claim is made.

Data: `results/fair_cost.csv`, `results/fair_cost_oracles.csv`, and `results/fair_cost.log`. The oracle gap combines per-particle position distance and mass/energy-normalized velocity distance, matching the main runner. It is distance to the solved discrete objective, not error against continuous mechanics.

## Reproduction

From this study directory, with `XPBD_ROOT` set to the existing checkout:

```bash
clang++ -std=c++17 -O2 \
  -DDUAL_GS_HOUR_PROBLEM_HEADER='"problem.hpp"' \
  fair_cost.cpp -o /tmp/dual_gs_fair_cost \
  -I"$XPBD_ROOT/experiments/hour_2026_10_31/local_solver" \
  -I"$XPBD_ROOT/experiments/round2_2026_10_31/local_work" \
  -I"$XPBD_ROOT/experiments/composition_2026_10_31" \
  -I"$XPBD_ROOT/external/eigen"
/tmp/dual_gs_fair_cost results > results/fair_cost.log 2>&1
```

The executed command used `/home/blacknoir/Projects/constrained-xpbd` for `XPBD_ROOT`; the binary was built outside the shared CMake build. Successful completion produced all 432 expected rows and six oracle rows.
