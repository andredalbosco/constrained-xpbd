# Fixed protocol

The primary comparison holds the objective, start-only work predictor, normalization, beta=gamma=1, physical state and initialization fixed. The material lift is exact; the three numerical methods differ in gradient refresh and treatment of stationarity.

1. Validate the equations and local support before running trajectories.
2. Same-state steps: all ten release/mixed examples from the five scene families, h=0.02, 0.08 and 0.2 seconds. Add the three-mass chain. Compare the three methods and six schedules at budgets of 1, 2, 4, 8, 16, 32 and 128 times the number of distinct scalar rows. Budgets count actual visits, so an interleaved repeated-kinetic schedule may finish mid-sweep. Record restoration passes, incident gradient evaluations and coordinate writes separately.
3. Additional controls: selected coarse-step cases at relaxation 0.5; stiffness times ten; a tenfold mass contrast; near-rest states. These remain named controls rather than unreported numerical safeguards.
4. Advance fixed low-budget inexact steps for one second. Report all complete and failed trajectories; a step with small constraint error but large stationarity error is not a solved constrained objective. Select representative promising cases for five-second runs only after the one-step and short-trajectory results are available.
5. Examine N=3, 8 and 32 chains for locality, propagation and setup-inclusive cost. The frozen work predictor contains an old-state linear solve. Local row visits do not make that setup free or turn the whole current implementation into an optimal large-system solver.

The oracle uses the independently implemented unlifted objective and exact Lagrangian Hessian. Only converged oracle roots are used for root-distance comparisons. Reduced curvature is checked separately: a root is not automatically a local minimum. Physical reference trajectories use paired RK4 subdivisions; their discrepancy is recorded before reporting trajectory error.

No trajectory projection, energy repair, fallback integrator or silent row deletion is allowed. Nonfinite iterates, collapsed edges or large-state runaway stop the run with a recorded reason. Finite inexact steps otherwise advance and retain their measured residuals. Multipliers restart at zero for each time step in this study.
