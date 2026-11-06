# Focused multiplier Gauss–Seidel experiment

The deliverable is a derived, implemented, measured local multiplier iteration. A dense solve is a reference for the same equations. No main-repository solver changes are planned.

## 1. Derive before implementing

Write the unconstrained mechanical objective first, then add kinetic/storage balances, then introduce their multipliers. Lift the quarter spring energy into an equivalent compliant constitutive row only after deriving the objective. Specify the full stationary equations and Newton system, then state the curvature approximation and every local elimination.

Cover the five existing examples explicitly: affine scalar spring, planar radial spring, rigid pendulum, elastic pendulum, and two-compliant-spring pendulum. Add a three-moving-mass chain so that an edge block cannot equal the whole problem; use N8/N32 only to examine propagation and cost. Gravity has separate work/storage leaves. Rigid length and tangency are separate constraints, and rod storage is zero.

Declare signs, normalization, frozen start-work target, old physical acceleration, each row's symbolic support, shared multiplier ownership, and the exact read/write sequence. Treat degeneracy as a recorded outcome, never silently remove an inconvenient equation.

## 2. Numerical methods to compare

- Refreshed multiplier GS with the usual dropped stationarity correction: update one accumulated multiplier and immediately push only its supported primal variables.
- Refreshed multiplier GS with a locally retained stationarity correction: solve the selected row and the stationarity equations restricted to its structural support, holding other multipliers fixed.
- Frozen linearization GS with the full stationarity correction applied once, followed by local multiplier pushes. This supplies a clean Schur-GS consistency control; its stationarity pass is explicit and counted.
- Dense joint solution of the equivalent objective-only formulation, with the material multipliers reconstructed and the complete lifted residual independently checked. This is a reference, not a local-method result.

The inverse primal metric consists only of constant per-particle position/velocity blocks and independent work entries. No full primal factorization is permitted inside a local row visit.

## 3. Scheduling experiments

First compare material–storage–kinetic, kinetic–storage–material, and storage–kinetic–material passes. Geometry uses the material slot; tangency is an explicit final row group.

Then compare edge-interleaved material/storage/endpoint-kinetic visits, swapping material/storage order and reversing traversal. Shared kinetic multipliers are single stored values and may be revisited by incident edges. Symmetric forward/reverse traversal is a separate schedule with its extra visits charged.

Report actual row visits, coordinate writes, residual/gradient work, completed sweeps, and wall time. A repeated-row schedule does not receive the same cost as a once-per-row sweep merely because both are called one iteration.

## 4. Validation gates

Before trajectories: material-lift equivalence; analytic gradients/Hessians against independent differences/oracles; local row update against a small restricted KKT solve; one frozen sweep against an independently assembled Schur GS sweep; immediate pushes versus explicitly accumulated frozen residuals; shared-multiplier ownership and untouched-coordinate checks; deterministic order traces; rigid rank/zero-work cases.

The implementation must actually evaluate and update local quantities. Dense state snapshots, full residuals, and factorization belong to validation and diagnostics, not unreported per-visit work. Any safety cap, line search, or fallback must be a separately named treatment with its effect counted.

## 5. Measurements

Start with same-state one-step tests and sweep budgets 1/2/4/8/16/32. Include converged controls, but judge animation candidates at small budgets. Then run short complete trajectories and selected longer survivors. Record failures and failed prefixes explicitly.

Separate material residual, kinetic/storage balance, geometry, tangency, stationarity, energy, reference-certified position/velocity error, and order sensitivity. Compare reference roots only when the reference has converged. Record low-budget inexact motion as inexact, not as a solved constrained objective.

Use release, mixed motion, turning/near-rest, and moderate stiffness/mass contrast. Avoid another broad cross-product of unrelated integrators. The final report must say which ordering and multiplier update are promising, which fail, and whether any useful result survives a genuinely local small-sweep implementation.
