# Independent review of ordering and stationarity

Reviewed the final numerical-semantics CSVs regenerated on 2026-10-31: [one-step](results/one_step.csv), [controls](results/controls.csv), [scaling](results/scaling.csv), [visit traces](results/visit_trace.csv), and [short trajectories](results/short_rollouts.csv). These contain 4,158, 1,296, 540, 702, and 814 rows respectively. No new solver experiments were run for this review. This review uses operation counters; recorded timings are descriptive, not a robust performance benchmark. See [the cost review](RESULTS_COST_REVIEW.md) for the subsequent comparison at matched row-evaluation budgets.

## What is compared

Let `q` be the number of distinct lifted rows. Budgets are **actual row visits**, `b*q`, not completed sweeps. For the two-spring and N3 scenes, `q=8` and `q=12`. An edge traversal contains 9 and 14 visits because shared kinetic multipliers are revisited. Therefore a budget of `q` visits can truncate an edge traversal; low-budget results include this coverage effect. `edge_reverse` reverses edge traversal while preserving the order within each edge.

The full residual `R` is the maximum of each row-family infinity norm and the Euclidean stationarity norm `||r_z||`. Energy is `|H-H_old|/E`. These are different diagnostics. The dense oracle converged for all 33 primary case/timestep pairs; its positive reduced curvature supports a local minimum, not global uniqueness. “Root gap” below means distance to that particular oracle solution, using the driver's mixed position/normalized-velocity metric.

## Ordering at matched visits

For `refreshed_local_rx`, material→storage→kinetic is the strongest tested default for these two chains at `h=.02` and `.08`. At **16q visits**, full residuals are:

| Scene / h | Material→storage→kinetic | Kinetic→storage→material | Storage→kinetic→material | Edge | Edge storage first | Reverse edge |
|---|---:|---:|---:|---:|---:|---:|
| Two-spring / .02 | **1.06e-6** | 1.30e-4 | 5.12e-5 | 7.81e-6 | 1.84e-4 | 4.89e-5 |
| Two-spring / .08 | **3.97e-5** | 2.66e-4 | 6.10e-5 | 3.70e-4 | 1.89e-3 | 5.10e-4 |
| N3 / .02 | **1.17e-6** | 1.08e-4 | 4.29e-5 | 1.67e-4 | 3.04e-4 | 1.78e-4 |
| N3 / .08 | **7.93e-5** | 4.84e-4 | 5.22e-4 | 8.07e-4 | 6.17e-3 | 3.76e-3 |

The same ordering wins within each refreshed method at 1q and 4q visits for these four scene/timestep pairs. It is **not universal**: at two-spring `h=.2`, 16q visits, edge gives `R=1.17e-4` versus `8.87e-4` for material→storage→kinetic. With 10× stiffness at `h=.08`, 8q visits, storage→kinetic→material gives `2.89e-4` versus `2.24e-3`. Frozen full correction often favors storage→kinetic→material at larger budgets: N3 `h=.08`, 16q gives `4.61e-4` versus `8.35e-4` for material→storage→kinetic.

Even within the preferred order, retaining correction does not win every small budget. N3 `h=.08`, 4q visits: dropped correction gives `R=.00793`, local correction `.00887`, frozen correction `.01447`. At 16q, the values become `.00573`, `7.93e-5`, and `8.35e-4` respectively.

## The dropped-correction plateau is real

For two-spring `h=.08`, material→storage→kinetic, **128q=1,024 visits**:

| Update | Largest row residual | `||r_z||` | Energy error | Root gap |
|---|---:|---:|---:|---:|
| Dropped correction | 2.46e-12 | **2.90e-3** | 5.95e-12 | 3.52e-3 |
| Local correction | 2.29e-12 | 5.67e-14 | 5.55e-12 | 7.32e-12 |
| Frozen full correction | 2.57e-11 | 2.85e-13 | 6.21e-11 | 9.91e-11 |

The dropped method satisfies the material and balance rows but remains away from stationarity. Its `||r_z||` barely changes from `.00293` at 4q to `.00290` at 128q. This is not a lack of energy projection accuracy. The independently tested identity is

\[
r_z^{new}=r_z^{old}-(G^{new}-G^{old})y^{new}
\]

for a dropped-correction visit, with constant base metric. Refreshing gradients accumulates stationarity error that row satisfaction alone does not remove.

N3 also has a second, distinct limitation. At `h=.08`, local correction gives storage residual `7.93e-5` at 16q, `4.01e-5` at 128q, and `1.36e-6` at 512q; corresponding stationarity is `1.06e-6`, `3.64e-7`, and `1.24e-8`. This is slow balance propagation, not the dropped method's stationarity plateau. The 512q result still misses the row tolerance `1e-8`.

Across the 198 primary runs at the largest budget, the declared row/stationarity convergence test passes 56 dropped, 167 local-correction, and 156 frozen-correction runs. The local method has six failed runs at that budget. These counts include easy scalar/release cases and must not be presented as mesh-scale success rates.

## Locality is verified; visit cost differs

The N3 edge trace revisits **the same row 4**, the middle particle's kinetic multiplier, at visits 8 and 14. Its stored value changes `.0741885 → .0442486`. At visit 8, the middle x coordinate changes `.236146 → .235052` and its velocity `−.095187 → −.129589`: a kinetic multiplier moves position through the coupled inverse metric. Local correction does not force global stationarity to decrease on every visit; that visit changes the full norm `.0790 → .0880`.

At N3 `h=.08`, 16q=192 visits, dropped and frozen methods evaluate 192 rows, while local correction evaluates **1,426** because it gathers incident rows for local stationarity. Coordinate writes are 1,472, 1,792, and 1,472 for dropped, frozen, and local respectively; frozen also performs 16 full restoration passes. Thus equal visits do not establish equal compute. The local method is structurally local, but its accuracy advantage here is only at matched visits; the cost review supplies a separate row-evaluation comparison.

## Failures and limits

- N3 `h=.2`: every local-correction order ultimately hits `coordinate_limit`; first failing visits are 66, 177, 75, 45, 449, and 52 in the table's order. At 16q, edge-storage-first survives with `R=.244`, which is not convergence. The 21 primary failure rows repeat these six failed prefixes at different requested budgets.
- N32 `h=.08`: every local-correction schedule fails by 4q. At 256q, four frozen schedules fail; its best survivor has `R=1.94`, energy error `10.6`. Dropped correction survives but its best `R=.423` is also unacceptable.
- At N32 `h=.02`, 256q, the best local residual is `1.84e-6`; dropped is `.0172` and frozen `2.66e-6`. This is improved accuracy, not demonstrated convergence at a small budget. These chain constructors change length and stiffness with N, so this is not an isolated algorithmic complexity experiment.
- Halving relaxation is not a general improvement: N3 `h=.08`, 512q, material→storage→kinetic local correction worsens from `1.36e-6` to `1.14e-5`. Damping was not tested on the failing N3 `h=.2` or N32 `h=.08` cases in the controls file.

Short trajectories also prevent equating a lower one-step KKT residual with better physical phase. N3 `h=.08`, 4q, material→storage→kinetic completes 1.04 seconds with position RMS `.00466` dropped and `.00497` local; the converged dense trajectory has `.00554`. Their maximum energy errors are `.0284`, `.0229`, and `1.82e-10`, respectively. The lower dropped-method position error can reflect cancellation with discretization error. Reverse-edge local correction at this budget fails after `.48` seconds; that prefix remains recorded. Selecting an order solely by position RMS would discard the objective and conservation distinctions.

The provisional **visit-only** candidate is material→storage→kinetic with local correction on moderate steps. This is not the compute-budget recommendation: defer that choice to [the cost review](RESULTS_COST_REVIEW.md), whose subsequent matched row-evaluation comparison favors frozen restoration in most tested combinations. The data also justify investigating slow balance propagation and a separately named stabilization for hard cases. They do not justify a universal ordering, unconditional stability, or a converged low-sweep mesh solver.
