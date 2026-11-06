# Sparse multiplier-coordinate GS implementation

`dual_gs.hpp` implements the lifted quarter objective described in `MATH_GS_REVIEW.md`. It is a numerical research engine for the existing serial-chain `compose::Model`, not an adapter for the application's general pinned graph. Compliant chains can have arbitrary particle count. The inherited old-reaction and work predictor currently support only one planar rigid rod; the constructor explicitly rejects larger or nonplanar rigid systems. Scalar springs retain the common model's affine constraint `C=x-l`.

## Fixed formulation and setup

The primal variables are `z=(q,u,w)`, with `x=x_old+xs*q`, `v=vs*u`, `W=E*w`, `mref=mean(mass)`, `vs=sqrt(E/mref)` and `xs=h*vs`. The supplied energy scale is floored at `1e-10`, matching the existing research setup. The independent midpoint velocity penalty has coefficient `1/beta`; the normalized work-target penalty has coefficient `gamma`. These are new objective parameters, unrelated to the application's reservoir parameters. Both must be strictly positive; the primary experiments use `beta=gamma=1`.

Each compliant material row is `C/xs+alpha*lambda=0`, with `alpha=4*E/(k*xs^2)`. The saddle uses minus multiplier signs and the material quadratic `-alpha*lambda^2/2`. Eliminating the material multiplier recovers quarter spring energy exactly. All kinetic, storage, length and tangency rows have zero diagonal compliance. Rigid edges use one length row instead of a compliant material row, without a duplicate geometric row. Rigid zero-storage rows retain their incidence works but have no position dependence.

The constructor obtains the old acceleration and frozen `implicit_midpoint_dg` incidence-work target from the prior `predictions.hpp`. This setup includes a dense old-state stiffness solve for compliant chains. The single-rod setup includes the old constraint reaction. The target is computed solely from the old state; it is a penalty target, not an identity for the resulting trajectory. This setup cost is outside `Result::milliseconds` and must be measured separately for an end-to-end cost comparison.

Each call to `solve` initializes all duals to zero and all primal variables at the exact minimum of the quadratic base objective: physical displacement rate `s=v_old+h*(a_old+gravity)/4`, independent velocity `v=2*s-v_old`, and work `W=W_hat`. The separate gravity contribution follows from the quarter gravity term left in the base objective.

## Row ordering and exact schedule manifest

The unique row IDs are stored as: all compliant material rows **or** all rigid length rows; all particle kinetic rows; all spring/rod storage rows; all gravity-storage rows when gravity is present; all rigid endpoint-tangency rows. A kinetic multiplier is owned by its particle and shared across every incident-edge visit. Repeated rows accumulate into the same dual entry.

Write `L_j` for material/length, `S_j` for spring/rod storage, `G_j` for gravity storage, `K_i` for particle kinetic and `T_j` for tangency. Missing gravity/tangency rows are omitted. Let `L`, `S`, `G`, `K`, `T` mean their complete forward-index lists. The pass schedules are:

| Name | Exact row list |
|---|---|
| `material_storage_kinetic` | `L, S, G, K, T` |
| `kinetic_storage_material` | `K, S, G, L, T` |
| `storage_kinetic_material` | `S, G, K, L, T` |
| `edge` | For `j=0..N-1`: `L_j, S_j, G_j, K_j, K_(j-1)` if `j>0`, then `T_j` |
| `edge_storage_first` | For `j=0..N-1`: `S_j, L_j, G_j, K_j, K_(j-1)` if `j>0`, then `T_j` |
| `edge_reverse` | For `j=N-1..0`: `L_j, S_j, G_j, K_j, K_(j-1)` if `j>0`, then `T_j` |

Thus reversal changes edge traversal and preserves within-edge order. Tangency is last in each pass schedule, and last within each rigid-edge bundle. Gravity is a per-edge bundle row in interleaved schedules. An interior vertex's kinetic row appears twice; the final vertex's appears once. For a two-mass compliant chain with gravity, `edge` is exactly `L0,S0,G0,K0,L1,S1,G1,K1,K0`; the corresponding reverse schedule is `L1,S1,G1,K1,K0,L0,S0,G0,K0`.

## Static locality audit

Every row has a fixed symbolic support. Particle support always includes its coupled position and independent-velocity block, even if a derivative currently vanishes. Material/length rows include their endpoints; kinetic rows include their particle and all incident works; compliant storage rows include their endpoints and their own incidence works; gravity rows include one particle and its gravity work; tangency includes its endpoints. Rigid zero-storage rows contain only incidence works and return before evaluating edge geometry.

The constructor precomputes each row's incident-row gather from these supports. `eval_row` creates only a sparse list of gradient entries and vectors of spatial dimension. `restricted_gradient`, `base_grad_local`, `apply_inv_local` and `stationarity_local` allocate vectors of support size. They do not create a dense global row or scan all rows. Stationarity gathers only the precomputed incident rows, skipping exactly zero dual coefficients. For the current bounded-degree serial graph, support and incident-gather sizes are bounded independently of chain length. Public primal and dual arrays are global storage; they are not copied inside a refreshed visit.

The inverse metric is exact for the base objective: each particle-coordinate pair uses `(mref/m_i)*[[1,2],[2,beta+4]]`, and each work entry uses `1/gamma`. Both directions of position/velocity support closure are checked. Multiplier-weighted row curvature is intentionally omitted by every primary method; objective and multiplier gradients remain the actual nonlinear gradients.

`refreshed_drop` evaluates the current scalar row, updates its accumulated multiplier, and immediately pushes its local primal support. `refreshed_local_rx` additionally gathers the current stationarity residual on that support, including all neighboring dual contributions, and applies the derived restricted KKT correction. It does not reset other duals or update unrelated primal coordinates. Relaxation scales the complete local primal and dual increment together.

`frozen_full_rx` evaluates all rows once at the beginning of each outer sweep, gathers full stationarity, and applies one restoration pass through the same particle/work inverse blocks. Subsequent scalar visits use frozen sparse gradients and current accumulated increments; they are the exact coordinate-GS updates for that frozen Schur system. A restoration is not repeated at every row. For `omega=0.5`, **only the scalar row pushes are relaxed**; the once-per-sweep restoration remains full strength, preserving the frozen stationarity relation. This differs explicitly from relaxed local retained-stationarity visits.

The only whole-system numerical operations are constructor setup and frozen-sweep preparation/restoration. Whole-system `base_gradient`, `stationarity`, `residuals`, `dense_metric` and `dense_gradients` helpers also exist for validation and final diagnostics. The refreshed visitor does not call these helpers. Initial safety is global; each refreshed/frozen row update subsequently checks only its changed support and active multiplier. A frozen restoration requires a global safety check because it writes all primal variables.

## Counters, stopping and timing boundaries

- `visits` counts scalar row updates that reached a primal/dual write. A failed row model before a write is not a visit; its evaluation is counted. A zero increment is still a completed visit.
- `row_evaluations` counts active row evaluations, incident-row evaluations used to gather local stationarity, and all frozen-row preparation evaluations. A row may be evaluated twice when it is both active and a nonzero-dual incident contributor. This is reported rather than deduplicated invisibly.
- `coordinate_writes` counts support size plus one dual entry for every completed visit, including zero values, and `n` primal entries for every full restoration. It is not a floating-point operation or byte-transfer estimate.
- `stationarity_passes` counts the full restoration passes used only by `frozen_full_rx`. It does not count local stationarity gathers.
- `completed_sweeps` counts complete prescribed row lists. `max_visits` may stop midway through a list; any frozen preparation performed for that partial sweep remains included in its counters.
- `milliseconds` starts before initialization and schedule construction and ends after numerical updates. It excludes constructor/predictor setup and final whole-system diagnostics. The runner must report setup-inclusive timing separately. Memory allocation and safety checks inside this interval are included.

There is no line search, energy repair, projection fallback, per-visit convergence stopping, or dual warm start. Finite inexact iterates advance for the requested fixed budget. Nonfinite iterates, nonfinite row models, incompatible zero-gradient rows, undefined edge directions, or normalized primal coordinates exceeding `1e8` stop with a reason. No finite iterate is clipped or replaced. The guard bounds `z`, not dual magnitude. A diagnostic failure/nonfinite iterate gives infinite reported residuals rather than misleading zeros.

Final diagnostics report maximum absolute row residuals by family, Euclidean full stationarity norm, and relative step energy change. `converged` requires all row-family maxima below `1e-8` and full stationarity below `1e-7`; energy error alone is not a convergence test. These diagnostic operations are excluded from numerical counters. The methods always perform the requested budget unless a recorded failure occurs, even if a smaller budget would already satisfy the final tolerance.

## Validation boundary

A separate smoke program compiled the engine and exercised all ten initial conditions from the five scene families with all three methods and two reverse-edge sweeps. This only established executable coverage and near-zero base gradient at initialization. Independent `validation.cpp` and `oracle_validation.cpp` provide the substantive algebraic tests. Trajectory quality, propagation speed, finite-budget stationarity, and exact-oracle agreement are measured by the runner and must not be inferred from this static audit or from the exactness of material lifting.
