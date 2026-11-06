# What to carry forward from this hour

The experiments separate three causes of bad motion: a mechanically inconsistent objective, missing position/velocity geometry, and a numerical solver that does not reach the intended constrained minimum. Exact total energy can hide all three. The most useful new tests also ask whether independent systems stay independent when their objectives are combined.

This hour extends the [second-round findings](../round2_2026_10_31/RESULTS.md), where the original independent-velocity objective was shown to double the acceleration component perpendicular to velocity in a generic moving state. The corrected quarter objective, start-only predicted work, endpoint tangency and numerical conditioning were therefore treated as separate ingredients, then combined in targeted controls.

## The strongest findings

| Finding | Evidence | Consequence |
|---|---|---|
| The old constraint reaction belongs in the predictor | Double/triple pendulum order tests: coefficient 1 gives approximately second order; 0,0.5,1.5,2 give approximately first order | Use the derived reaction coefficient; tuning beta or examining total energy cannot replace it |
| Energy conservation can hide the wrong mode | Triple pendulum, highest initial mode, 60 s at 30 Hz: omitted reaction moves almost all its energy into other modes, with total-energy error about 7e-9 | Check modal amplitude, phase and trajectories alongside conservation |
| A global energy row can couple disconnected systems | Fully converged independent scalar/radial springs exchange up to 0.6% under the quarter objective and 8.4% under a Verlet projection target | Use one energy closure per connected physical/work component |
| Total-energy work normalization is not replication invariant | N independent identical copies change each copy's effective work weight from gamma to gamma/N | Give work preferences local physical scales, or compensate the chosen normalization explicitly |
| Work variables and local multipliers can be eliminated exactly | Weighted graph-QP checks, AD/analytic derivative agreement, 40 matching full/reduced trajectories | Retain the work preference as a graph quadratic; elimination is a concrete alternative to shared local dual iteration |
| Correct arithmetic removes some local-solver stalls | Exact incremental line-search differences improve long AL convergence from 13/28 to 26/28 | Keep the same objective, but avoid subtracting nearly equal total objective values |
| Better conditioning needs a nonlinear work policy | Adaptive forcing improves fully converged horizons from 59 to 72/252 at budget 12; direct chain solves remain faster | Measure factor construction, model refreshes and convergence together |

Percentages for the disconnected spring comparison refer to each component's own initial excitation energy. The modal-transfer diagnostic uses the fixed initial high-mode quadratic energy as denominator; the conserved energy is the nonlinear mechanical energy.

## 1. Rigid chains distinguish two objectives that a single pendulum did not

For linear gravitational potential, the quarter objective is exactly equivalent, up to a constant, to a midpoint-potential objective whose predictor includes `h² a_reaction,n / 4`. The plain midpoint-potential objective omits this contribution. Position and velocity remain independent, with the midpoint velocity penalty and endpoint geometry/tangency constraints.

| Position RMS over 1 s | h=1/30 | h=1/60 | h=1/120 | h=1/240 |
|---|---:|---:|---:|---:|
| Double, omitted old reaction |0.02951|0.01421|0.006972|0.003454|
| Double, quarter/full reaction |0.002138|0.000536|0.000135|0.0000338|
| Triple, omitted old reaction |0.02958|0.01371|0.006611|0.003249|
| Triple, quarter/full reaction |0.003523|0.000887|0.000224|0.0000564|

All these comparisons have converged solves and certified references. The order failure therefore survives tighter iteration. The 120-run reaction-weight grid also shows that overweighting the reaction can suppress a particular mode transfer while still giving first-order trajectory error.

![Reaction and modal failure](rigid/figures/reaction_and_mode_transfer.png)

The hidden Cartesian normal-velocity reflection also generalizes to multiple rods. Exact endpoint lengths and midpoint kinematics reflect an inherited normal velocity as the step size tends to zero. Its kinetic energy is preserved. One full step and two half steps then approach different velocity states, defeating the usual adaptive estimator. All 11 tested Cartesian-DG adaptive chain runs hit the depth limit; feasible-state methods complete 30/33 comparable runs, with the remaining heavy-tip cases investigated separately.

Canonical angle/momentum methods supply useful small-system controls: geometry and tangency hold by construction. A symmetry-aware exact-work correction can preserve both energy and free-chain angular momentum. Its near-degenerate implementation needs stable Hamiltonian differences; an algebraically equivalent difference evaluation rescues 40 targeted near-relative-equilibrium runs. This remains a reference construction for two/three-link planar chains, not a local solver for arbitrary XPBD graphs. Its tiny-denominator guard and broader failures are documented.

[Derivations, momentum/phase comparisons, validation and limits](rigid/README.md) · [Recorded-motion viewer](rigid/viewer.html) · [60-second comparison video](rigid/animation/modal_transfer_60s.mp4), normal speed with labeled 16× angle magnification.

## 2. Work reduction solves an algebraic issue, not every nonlinear issue

Start from a mechanical objective plus a quadratic preference for start-estimated work, constrained by local kinetic and storage balances. Writing these balances with an energy-graph incidence matrix permits exact elimination of all work variables. The reduced problem retains a graph-Laplacian quadratic penalty and one energy constraint per connected component. Kinetic and storage multipliers can be reconstructed afterward.

The graph construction passes 240 independent weighted QP checks, including cycles and disconnected components. Tree forward/transpose operators pass checks through 16,384 energy nodes. The mechanical full and reduced formulations agree across 40 five-second trajectories, and their analytic/AD derivatives agree to about 3e-15 in four cases. The current mechanical Newton prototype stores the graph operator densely; the fast tree operator has not yet been integrated into that solve.

Exact elimination changes the path followed by a nonlinear solver. At aligned multipliers, 16 independent first-Newton-step comparisons agree below 4e-13, but full Newton subsequently permits intermediate local imbalance while the reduced solver reconstructs work after every state change. Aligned dual initialization can even be worse than zero duals.

| Same final objective, 16 hard one-step configurations | Converged within 128 evaluations |
|---|---:|
| Full variables, zero initial duals |15/16|
| Full variables, graph-aligned initial duals |12/16|
| Reduced, exact work curvature |14/16|
| Reduced, Gauss-Newton work curvature |15/16|
| Reduced, Gauss-Newton followed by exact Newton |16/16|

The last method solves the hardest recorded case in 67 combined evaluations and has positive exact reduced curvature. It preserves the objective and its KKT stopping conditions. A subsequent one-second N8/N32 trajectory grid completes 8/8 horizons with GN-then-Newton, 7/8 with full Newton and 3/8 with exact reduced Newton; all 8 physical references pass refinement and energy checks. The hardest 30 Hz completed trajectory still has 0.0278 m position RMS and 1.15 m/s velocity RMS despite 1e-10 energy error. These are promising short-run numerical results, not a claim of long-duration robustness. Absolute row tolerances are reported in their actual coordinates; row scaling is not silently treated as an identical raw-energy threshold.

[Exact reduction, normalization and disconnected-system experiments](work_reduction/README.md).

## 3. Local shared-multiplier iteration remains a distinct option

The hour tests a fixed quarter objective with explicit work, using start-only midpoint-DG targets, joint Newton, block augmented Lagrangian updates, several initializations, and larger chains. Midpoint initialization fixes the tested small joint-solver failures across all 9 case/regime combinations. A continuation control nevertheless finds a stationary saddle with residuals near 1e-13 and reduced curvature −0.979; KKT residuals alone are insufficient.

The local AL solver has both a genuine slow mode and an avoidable arithmetic stall. Direct incremental objective/constraint differences improve the long-run success count from 13/28 to 26/28. Even after this repair, a 32-mass chain with coupled edge blocks requires roughly 3,457–4,801 evaluations in the reported cases. A single global primal block takes 41–58. A joint first step and warm multipliers do not make subsequent local steps converge within 2,048 evaluations.

The linearized iteration spectrum explains part of that gap: at N32, original split updates have spectral radius 0.962 per three-pass/dual cycle; larger edge blocks at a suitable penalty reduce it to 0.766. These are local asymptotic rates, not global convergence guarantees. Small one/two-mass examples can make an edge block equal the entire system, so they do not establish large-system locality.

[Local solver results](local_solver/RESULTS.md) · [Turning points and the iteration spectrum](local_solver/TURNING_AND_LOCALITY.md).

## 4. The useful paper transfers are numerical, with measured costs

The experiments distinguish numerical coordinates from changes to the physical objective. Graph metrics inspired by Large Steps and local spectral balancing inspired by Muon preserve the converged mechanical equations. The earlier capped-step tests were motivated by trust-region ideas; they did not implement PPO's clipped surrogate, and no PPO-specific or mechanical guarantee is claimed.

The strongest new transfer is adaptive inexact-Newton forcing: avoid solving an outdated linear model too accurately before recomputing the nonlinear state. It improves aggregate convergence, although some metric/budget combinations regress. Caching factors also helps some cases and hurts others.

For 128 masses, mass ratio 100, h=.04, ten steps:

| Method / budget | Isolated total CPU | Converged steps |
|---|---:|---:|
| Diagonal scaling, 48 |0.849ms|0/10|
| Local 2×2 right scaling, 48 |0.904ms|10/10|
| Curvature graph + adaptive forcing, 12 |1.06ms|10/10|
| Exact serial-chain block Newton |0.259–0.270ms|10/10|

The serial-chain direct control remains strongest here. General graphs require their own cost study. Fixed-old-state step continuation increases completion from 13/36 to 34/36 hard horizons, but four cases select different coarse-step roots. One completed case conserves energy to 8e-14 while retaining 2.80 rad modal phase error. Better nonlinear convergence does not fix insufficient temporal resolution.

[Measured costs and physical-reference separation](preconditioning/RESULTS.md) · [Equations and primary-paper sources](preconditioning/README.md).

## Five-second extension: a later failure boundary

The final [larger-chain follow-up](work_reduction/longer_rollouts/RESULTS.md) extends the same finite-budget solver from one to five seconds. GN-then-Newton completes 6/8 horizons, versus 8/8 at one second. The stiff N32, mass-ratio-10 cases stop at 1.5667 s (30 Hz) and 4.25 s (60 Hz), after exhausting the 128-evaluation budget inside the KKT line-search process. The inherited failure label is `kkt_line_search_stalled`; it does not establish that a larger budget would succeed. Every accepted first-second prefix matches the earlier saved run, including evaluation counts.

The two stiff N32 physical references also fail five-second refinement certification, even at 512/1024 RK4 substeps per frame. Their longer-run physical RMS values are therefore not accepted as certified motion errors. The six other reference configurations pass the original refinement checks. These failures strengthen the requirement for sustained trajectories and independently certified references before integration.

## Combinations worth implementing next

| Direction | Ingredients to combine | What must still be demonstrated |
|---|---|---|
| Reduced compliant-work solver | Quarter objective; independent velocity; start-only midpoint work; component closure; graph elimination; GN then Newton | Sparse tree/graph application, sustained larger-system trajectories, physical normalization across heterogeneous components |
| Local XPBD-like work solver | Coupled endpoint/velocity blocks; separate shared-dual pass; stable incremental differences; suitable block preconditioning | Equal-time performance on larger graphs, rather than iteration counts alone |
| Rigid constrained integrator | Position and endpoint tangency; old-reaction predictor; symmetric invariant treatment | Robust mixed rigid/compliant graphs and phase accuracy at practical step sizes |
| Geometric reference method | Canonical coordinates; symmetry-aware discrete work/AVF; stable energy differences | Broader mass contrasts and degeneracies; this is not automatically a local XPBD algorithm |

The [composition matrix](COMPOSITION_MATRIX.md) separates actual combinations from proposed transfers. The [integration proposal](INTEGRATION.md) maps these candidates onto the actual repository state, including independent velocity commits, pinned graphs and parameter names. The first implementation should expose these choices independently and retain the current reference controls. Energy error, geometry, tangency, stationarity, reduced curvature, phase/modal transfer, component-energy exchange and setup-inclusive cost should remain separate diagnostics. The experiments are saved additions; the main project solver has not been changed.

The root CMake validation currently passes all 11 checks. The original and this hour's viewers pass data/JavaScript checks, but interactive playback is unverified: the supported browser was unavailable and automatic approval review rejected the Chrome fallback. Static scientific figures were inspected directly. The native 60-second MP4 passed metadata and full-decode checks, and representative frames were visually inspected.
