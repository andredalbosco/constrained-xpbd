# Proposed integration into constrained-xpbd

This is a future implementation proposal. The main solver has not been changed. It is based on the current `cxpbd::Particle`, `Spring`, `Model`, `Settings`, `Diagnostics`, and `step(Model&, const Settings&)` in `/home/blacknoir/Projects/constrained-xpbd/xpbd.hpp`, the scene constructors in `scenes.hpp`, and `parse_mode`/`run` in `main.cpp`.

## Candidate compositions

| Candidate | Composition to implement | Evidence and remaining boundary |
|---|---|---|
| First compliant prototype | Repaired quarter objective; independent endpoint velocity; soft midpoint kinematics; frozen start-only implicit-midpoint-DG work target; exact local energy balances | Midpoint initialization converged every step in all nine small case/regime combinations. This is an objective-based method distinct from existing reservoir feedback. |
| Reduced compliant solver | Same objective, with work eliminated into its graph quadratic; one energy row per component; Gauss–Newton work curvature followed by exact Newton | Full/reduced agreement is measured on 40 five-second trajectories. GN-then-Newton improves the harder recorded controls. Sparse arbitrary-graph implementation and sustained mesh trajectories remain proposals. |
| Rigid extension | Quarter predictor including the old constraint reaction; exact endpoint lengths and endpoint velocity tangency; independent velocity | Double/triple rods establish the reaction term's importance. Mixed rigid/compliant 3D graphs are untested. |
| Mechanical reference/control | Midpoint discrete-gradient force and exact midpoint kinematics; existing preconditioning variants | The preconditioning study solves these mechanical equations, not the quarter/work objective. Keep it as a separately labeled family. |
| Local work research option | Shared vertex/storage multipliers; coupled endpoint-position/velocity blocks; separate dual updates; stable incremental acceptance differences | Long-chain local AL remains expensive or inexact at practical tested budgets. It is not the initial production candidate. |

See [local results](local_solver/RESULTS.md), [work reduction](work_reduction/README.md), [rigid results](rigid/README.md), and [conditioning equations](preconditioning/README.md). Combining individually measured ingredients does not validate the resulting combination.

## Minimum adapter and public interface

Add a separate experimental entry point, for example `step_constrained(Model&, const ConstrainedSettings&, ConstrainedState&)`, with its own diagnostics. Keep `step` and the existing `Mode::{XPBD,Passive,Active}` regression behavior intact; `main.cpp::run` can dispatch an explicitly selected new method. Its current `parse_mode`, mode enumeration, output naming, and CSV schema all need deliberate extensions.

The adapter must map free particles to contiguous **3D** position and velocity unknowns and map each `Spring::{a,b}` through that index. Use `k=1/compliance` only for positive compliance; zero compliance supplies a hard distance row. Preserve rest lengths, masses, gravity, and fixed endpoints. Do not reuse the study's serial-chain indexing (`j` and `j-1`), one-common-rigidity flag, or implicit origin anchor as a general graph representation. The existing `make_grid()` and `make_cloth()` already contain cycles and multiple pins; the cloth is spatial.

Snapshot old positions and velocities before building the predictor or frozen target. Assemble potential forces from `spring_energy`'s analytical energy and gravity. Solve in scratch state, then commit both solved `Particle::x` and `Particle::v`; the current `step` assignment `v=(x-old_x)/h` would destroy the independent-velocity method. Stationary pins remain fixed with zero velocity, as required by `validate`. Store per-step work and optional warm multipliers in the new state object, not in `Particle::reservoir` or `Spring::reservoir`. Work is a step transfer, not a carried reservoir.

The first adapter should support compliant graphs and retain a dense full-variable KKT oracle for small fixtures. Add exact work elimination against that oracle before optimizing its graph operator. Cache graph topology only while endpoints, pins, or constraint types remain unchanged.

## Parameters and component balances

Existing `Settings::kinetic_reservoir_compliance` and `elastic_reservoir_compliance`, exposed as `--beta` and `--gamma`, govern the reservoir equations; larger values weaken feedback. They are not the new work parameters. Introduce explicit names such as `velocity_coupling_compliance` for the coefficient in `||v+v_n-2Δx/h||²_M/(2 beta_v)` and `work_preference_weight` for `gamma_w ||W-W_hat||²/2` in the declared normalized coordinates. Larger `gamma_w` strengthens this preference. Keep numerical AL penalty, tolerances, and factorization options separate from physical weights.

For each free particle, use a kinetic-energy node. For each compliant spring, use a potential-storage node with an incidence work variable for each free endpoint. Gravity can have one storage leaf per free particle. With the convention `Δe=B W`, an edge contributes positively at its receiving kinetic node and negatively at its storage node. A fixed endpoint has zero work and no kinetic degree of freedom. Two otherwise independent bobs attached to the same stationary pin must not become one energy component merely because they share that pin. Springs between two pins carry no dynamic transfer. A hard rod may carry endpoint transfers with zero summed storage change; its general-graph implementation still needs validation.

Build connected components from the actual dynamic work graph. Elimination requires one retained energy constraint per component, including isolated dynamic components. A single global row permits measured spurious exchange between disconnected systems. Keep the quadratic work preference after elimination; dropping it changes the objective. Cyclic graphs also require selecting the penalty-minimizing circulation, whereas tree work is unique at a feasible state.

Declare physical energy/work scales explicitly. With a total excitation normalization, N independent identical copies reduce the effective per-copy work weight to `gamma_w/N`. Component scales restore independence between disconnected components, but heterogeneous local weights and near-zero excitation still require a documented dimensional convention and tests. Row preconditioning must not silently change those physical scales or the final raw-residual stopping criterion.

## Rigid reaction and numerical policy

For hard constraints, compute the old physical acceleration from the acceleration-level constrained equations. With constraint-gradient columns G, solve `(Gᵀ M⁻¹ G) lambda = -Gᵀ M⁻¹ f - dot(G)ᵀ v`, then use `a_n=M⁻¹(f+G lambda)` in the quarter predictor. Account for pinned degrees of freedom and redundant constraints. `Spring::previous_lambda / Model::previous_dt²` is the existing history estimate; it is not an established replacement for this reaction in the new formulation, especially at initialization or after timestep changes. Never identify energy multipliers with material reaction multipliers without deriving that relation.

Enforce endpoint tangency explicitly for the rigid candidate. The measured hidden normal-velocity mode survives exact position and total-energy constraints. A single pendulum cannot validate a multi-rod reaction implementation: the omitted-reaction failure appeared in double/triple rods.

Start with mechanically informed midpoint initialization. Preserve the exact objective and gradient when using approximate numerical curvature, then test final KKT residuals and exact reduced curvature. Use stable incremental objective/constraint differences for acceptance. Graph and local spectral preconditioners are promising transfers, but the existing measured conditioning kernels solve a different mechanical Jacobian; their application to an indefinite constrained KKT system needs its own derivation and timing. A dense elimination or smaller KKT matrix alone establishes no speedup.

## Acceptance gates before enabling a new default

1. **Adapter equivalence:** reproduce the five small environments in 3D, including planar embeddings, exact pins, arbitrary edge ordering, and compliance conversion. Check forces, gradients, Hessian products, and work incidence against independent finite differences/oracles.
2. **Same formulation:** compare full versus reduced work objectives and reconstructed balances/multipliers; test cycles, disconnected objects, shared stationary pins, and replicated copies. Report raw component energy, local balance, and stationarity separately from scaled solver residuals.
3. **Motion:** certify refined references; measure position/velocity error, timestep order, rigid tangency, phase, modal transfer, and component energy. Include scalar turning points, high rotation, stiffness/mass contrasts, and double/triple rods. Small energy error is not a motion acceptance criterion.
4. **Finite budgets and failures:** distinguish convergence, saddle detection, budget exhaustion, line-search failure, and geometric degeneracy. Preserve completed time and failed prefixes. Commit/reject/fallback behavior must be explicit; a fallback to XPBD changes the method and must be logged.
5. **Repository scenes and cost:** exercise `make_spring`, `make_pendulum`, `make_grid`, and `make_cloth`, plus genuinely nonplanar motion. Compare setup-inclusive CPU time and memory at matched residual and physical-error targets. The serial-chain O(N) direct control is not an arbitrary-mesh performance result.

Moving pins, collisions, external time-dependent forces, changing topology, damping, redundant rigid networks, and mixed rigid/compliant 3D graphs are outside the demonstrated mechanical suite. Their external work or dissipation must enter the formulation before conservation claims extend to them. Existing XPBD/passive/active outputs should remain regression controls throughout integration.
