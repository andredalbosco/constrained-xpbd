# Eliminate work variables exactly, then audit the numerical consequences

Quadratic work preferences can be retained while eliminating every explicit work variable and all local kinetic/storage multipliers. The remaining constraint is energy conservation for each connected component of the energy graph. This is useful for the shared-multiplier problem, but does not make the resulting optimization vertex-local or guarantee easier Newton convergence.

The small full/reduced comparisons agree, including 40 five-second trajectories. Larger stiff chains expose nonlinear convergence differences despite exact formulation equivalence. A separate normalization experiment finds that dividing a quadratic work penalty by the square of the entire system's energy makes unrelated copies alter each other's effective penalty.

## Start with the constrained objective

Keep the mechanical position/velocity objective `Phi(x,v)` and let `e(x,v)` collect kinetic energy per particle, spring energy per spring, and gravitational potential energy per particle. Each work variable is an edge from an energy store to a kinetic-energy node. The signed incidence matrix `B` has +1 at its kinetic destination and -1 at its storage source. Write `Delta e = e(x,v)-e(x_n,v_n)`.

For a start-only target `W_hat`, fixed positive-definite work metric `Q`, and positive weight `gamma`, the constrained problem is

\[
\min_{x,v,W}\quad \Phi(x,v)
+\frac{\gamma}{2}(W-\widehat W)^TQ(W-\widehat W),
\qquad \Delta e(x,v)=BW,
\]

together with any geometric or tangency constraints. `W_hat` is frozen during the step. If targets or metrics depend on the unknown endpoint, their derivatives must also be included; the implementation here does not claim to cover that case.

Let `Z` contain one component-indicator column per connected component, and define

\[
d=\Delta e-B\widehat W,\qquad L=BQ^{-1}B^T.
\]

Feasibility requires `Z^T Delta e = 0`. At a feasible `(x,v)`, the unique least-penalty work is

\[
W_*=\widehat W+Q^{-1}B^TL^+d.
\]

Substitution gives the exactly equivalent constrained objective

\[
\boxed{\min_{x,v}\quad\Phi(x,v)+\frac{\gamma}{2}d^TL^+d,
\qquad Z^T\Delta e(x,v)=0.}
\]

The quadratic penalty survives. Merely discarding it would instead produce the different free-work formulation. The matrix `L` is a graph Laplacian; the connection with minimum-energy network flow is standard, for example in [Stanford's Laplacian notes](https://ee263.stanford.edu/lectures/laplacian.pdf). Its use for these local mechanical energy balances is the construction tested here.

Only now introduce multipliers. Preserve the earlier negative-balance convention: collect the kinetic multipliers `mu` and storage multipliers `nu` into `y`, and use the Lagrangian term `-y^T(Delta e-BW)`. The eliminated multipliers recover as

\[
y=-\gamma L^+d+Z\eta.
\]

Thus the reduced solve has one energy multiplier per connected component; local work multipliers are reconstructed afterward. The numerical KKT helper adds the balance rows with a plus sign, so its stored dual is the negative of this `y`. The code uses normalized energies and works, with `Q=I`, and explicitly checks reconstructed full stationarity using its declared sign convention.

## Chains, cycles, and disconnected systems

For the anchored compliant chain with gravity, the energy graph is a tree. Its leaves include gravity stores, not extra mechanical degrees of freedom. On a tree, feasible work is unique; the target selects states through its penalty but cannot independently change the feasible work at a fixed state. A subtree sum computes the flow, and an accumulated potential computes its transpose, both in linear time. `tree_operator.py` verifies both operators, including the transpose identity, up to 16,384 energy nodes.

Cycles allow circulation in work space. The quadratic metric selects the circulation closest to the target. `graph_audit.py` compares the weighted formula to an independent full equality-constrained quadratic solve on trees, cycles, complete graphs and disconnected graphs. All 240 randomized checks pass; maximum balance and full-KKT discrepancies are below 1e-12. These are graph algebra tests, not physical grid trajectories.

Disconnected components need separate energy rows. The counterexample has total energy change zero but changes of +1 and -1 in two components. No internal work assignment can realize it. One global energy row would incorrectly permit exchange between isolated systems.

The mechanical prototypes store the inverse operator densely. They do not yet use the linear-time tree operator inside the Newton solve. Exact elimination introduces nonlocal mechanical coupling even when its operator can be applied cheaply.

## Mechanical equivalence evidence

The mechanical objective is the repaired quarter objective, plus the independent-velocity midpoint penalty at beta=1. Work targets use the start-only stiffness-aware midpoint prediction and discrete work along its predicted displacement. Tests cover scalar springs, orbiting radial springs, elastic pendulums, two-spring chains, and larger anchored compliant chains. Rigid chains are tested in the separate rigid study.

The first implementation differentiates the full coordinate map using automatic-differentiation jets. A second analytic implementation uses

\[
\nabla\Phi_{\rm reduced}=\nabla\Phi+(\nabla\Delta e)\gamma P^T(W-\widehat W),
\]

and its exact Hessian, where the columns of `nabla Delta e` are the energy gradients and `P=B^+` for the uniform metric. Analytic and AD reduced objectives, gradients and Hessians agree to about 3e-15 in four independent cases. Twelve finite-difference checks of the original coordinate map also pass.

| Off-manifold extension | Converged reduced cases / 84 | Both full and reduced converged | Largest position difference in those pairs |
|---|---:|---:|---:|
| Balanced pseudoinverse | 78 | 78 | 2.81e-9 m |
| Omit first balance row | 81 | 81 | 2.80e-9 m |
| Omit last balance row | 81 | 80 | 1.17e-8 m |

This table uses the 128-evaluation budget. The complete 1,008-row study also includes budgets 16 and 64 and both backward-Euler and midpoint initializations. Both-converged filtering is essential: a failed full solve is not an equivalence oracle. Velocity differences are at most 1.05e-6 m/s in those pairs, consistent with the finite KKT tolerances.

The three extensions agree on the energy-feasible set but distribute energy imbalance differently during Newton iteration. They are numerical choices, not different feasible objectives. The reduced KKT dimension falls from 79 to 33 for an eight-mass chain, and from 319 to 129 for the 32-mass chain used in the later analytic study.

All 40 five-second trajectory runs converge on every step. Full and reduced RMS position errors agree at the reported precision across both tested step sizes. Paired RK4 references differ by less than 6.6e-11 in the stated scaled state diagnostic. The larger one-step grid gives 23/24 converged full solves and 22/24 reduced solves at budget128. It therefore does not support a universal convergence advantage for elimination. Follow-up controls in `nonlinear_results` separate initial multipliers, residual-work curvature, and energy-row scaling; consult `summary.json` for their final counts.

Recorded timings in these reduction runs overlap other research jobs. The original AD implementation also does redundant work. No runtime speedup is claimed from its timings or smaller KKT dimension.

### Elimination changes the nonlinear path

At a common state, aligning the full local multipliers with the eliminated-work stationarity condition makes the first full and reduced Newton position steps agree. All16 independent Schur-complement checks pass, with relative discrepancies below4e-13. Subsequent iterates differ: the reduced method reconstructs work nonlinearly after each position/velocity change, whereas full Newton takes an affine work update and permits intermediate local balance defects.

Across16 harder N8/N32 configurations at budget128, the full zero-dual solve converges15, the dual-aligned full solve12, exact reduced Newton14, and reduced Gauss-Newton15. Scaling the one energy row does not change these success counts. Gauss-Newton here drops only residual-work curvature from the numerical Hessian; it retains the objective, gradient, constraints, and final KKT stopping test. Every converged result in this grid has positive exact reduced curvature (minimum0.0335 in the chosen normalized coordinates).

In the N32, stiffness-times10, mass-ratio10, h=1/60 case, exact reduced Newton remains inexact after128 evaluations, while reduced Gauss-Newton converges in10 and the full zero-dual solve in7. At h=1/30, Gauss-Newton reaches energy residual3.1e-12 but stationarity7.6e-7, so it still fails the1e-7 stationarity threshold. A separate `polish_results` follow-up switches from Gauss-Newton to exact Newton after at most64 evaluations; it converges16/16 cases, including the hardest case in67 combined evaluations. The objective and stopping conditions are unchanged. `nonlinear_audit_results` and `polish_audit_results` replay the controls with raw residuals and failure reasons added; they reproduce the same success counts.

The aligned dual is algebraically correct but not automatically a better nonlinear initializer. Likewise, a smaller system does not establish better conditioning or branch selection. Row-scaled variants retain the same numerical stopping tolerance after dividing the energy row by the square root of the number of energy nodes; this is a different tolerance on the raw component-total energy. Audited CSVs explicitly retain raw component-total and reconstructed full balance residuals. Their timing and success counts are not an equal-physical-tolerance performance claim.

### Sustained larger-chain check

`rollouts.cpp` tests24 one-second trajectories: N8/N32, stiffness multipliers1/10, mass ratio10, and h=1/30 or1/60, with full, exact-reduced and GN-then-Newton solves. It stops at the first unconverged step and records NaN errors when no step is accepted. GN-then-Newton completes8/8 horizons, full Newton7/8, and exact reduced Newton3/8. Fully completed full/reduced comparisons agree in motion at the finite-tolerance level. `reference_audit.cpp` independently records paired RK4 refinement and energy checks; all8 references pass, with worst combined state discrepancy5.4e-6 and reference energy error3.7e-9.

Completion still does not establish accurate motion at a coarse step. The hardest N32 stiff/mass-contrast run at30 Hz completes with energy error1.0e-10 but position RMS0.0278m and velocity RMS1.15m/s over one second. At60 Hz these become0.00916m and0.693m/s. This supports the reduced numerical method as a candidate while retaining timestep accuracy as a separate issue. It is not evidence of long-duration robustness or equal-time superiority.

## Normalization must survive independent copies

The tested normalized objective has the form

\[
\frac{\Phi}{E_*}+\frac{\gamma}{2E_*^2}\sum_j(W_j-\widehat W_j)^2.
\]

For N disconnected identical copies with total excitation scale `N E_*`, the symmetric per-copy problem has effective work weight `gamma/N`. Thus adding an unrelated copy changes the optimization problem. This is a physical weighting change, not harmless rescaling of equality rows.

`normalization.cpp` solves this exact symmetry-reduced problem for N=1,2,8,32, using two targets and three nonlinear compliant systems, at 30 Hz for five seconds. It is not a benchmark of an explicitly assembled N-copy solver. All 42 trajectories converge at every step and have certified references.

| Case and target | Position RMS, one copy | Position RMS, 32-copy total normalization |
|---|---:|---:|
| Radial spring, frozen power | 0.189253 m | 0.021013 m |
| Elastic pendulum, frozen power | 0.078838 m | 0.015355 m |
| Two-spring chain, frozen power | 0.129592 m | 0.034741 m |
| Two-spring chain, predicted midpoint DG | 0.031188 m | 0.031296 m |

The apparent improvement for frozen power comes from weakening that imperfect target; it is not a benefit of adding copies. The better predicted target is less sensitive but not invariant. Multiplying the total-normalized gamma by N restores identical trajectories in this experiment. Per-component or per-energy-term physical scales are the more general direction; their heterogeneous-system design remains open.

![Elimination equivalence and normalization dependence](figures/elimination_and_normalization.png)

## A global energy row can exchange energy between disconnected objects

This is a separate issue from penalty normalization. `disconnected.cpp` explicitly assembles two independent mechanical objectives with no shared forces or degrees of freedom. It then compares one energy row per component with a single summed energy row, keeping the objective and its scale identical. There are no work variables in this experiment: it isolates the effect of the energy constraint itself.

The36 five-second runs combine three pairs of compliant systems, three objectives, two step sizes and two choices of energy closure. Thirty-two converge at every step; four original-objective/component runs retain inexact-step flags. Every quarter and projection comparison is fully converged, and paired physical references agree below1.4e-9.

For an independent scalar spring and radial spring at30 Hz, the quarter objective plus one total-energy row changes the components' excitation energies by up to0.148% and0.622%, while total-energy error stays1.1e-7. With projection toward a Verlet target, the changes reach1.99% and8.38%, while total-energy error stays2.7e-8. Per-component energy rows remove this exchange to the accumulated nonlinear tolerance. The original biased objective gives even larger exchange, but is not needed to demonstrate the effect.

The single multiplier lets the correction borrow energy from one disconnected system to compensate another. Total conservation is therefore insufficient as a composability criterion even when the mechanical objective is separable and the nonlinear solve converges.

![Spurious exchange under a global energy row](figures/disconnected_energy_transfer.png)

## Reproduction and next composition

Build the CMake targets `work_reduction_study`, `work_normalization`, `work_reduction_larger`, and `work_reduction_nonlinear`. The first executable accepts `OUTPUT validate|single|trajectories|all`. The others accept an output directory; `work_reduction_larger OUTPUT validate` runs just its independent analytic-oracle checks.

Run `graph_audit.py`, `tree_operator.py`, and `analyze.py` with NumPy/Matplotlib available. Original results, failed solves and per-trajectory traces remain in their named result directories. `summary.json` is generated from saved CSVs rather than transcribed run logs.

The concrete integration candidates are: retain start-only work targets; use graph elimination or graph-derived shared-dual initialization; preserve component energy closure; use local physical normalization; and compare exact versus Gauss-Newton work curvature while keeping the same final objective. A sparse implementation and physical cyclic-graph trajectories are still required before claiming this solves the general vertex-iteration problem.
