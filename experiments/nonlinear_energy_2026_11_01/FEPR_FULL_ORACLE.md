# One-step full FEPR objective diagnostic

This solves the fixed-center Eq1 objective on the stiff spring second step at30Hz. It does not replace the source-faithful fast FEPR algorithm. The state is read from the saved per-step-target trajectory. Variables are whitened by the exact D metric, equality rows are nonsingularly scaled, and SLSQP uses analytic Jacobians. Three starts are tested; this is local minimum evidence, not a proof of the global minimum.

| Solver/seed | Success | Objective | Raw constraint L1 | Whitened stationarity | Position jump (m) |
|---|---|---:|---:|---:|---:|
| Fast FEPR | yes | 19.724425 | 4.595e-11 | not enforced | 2.00820537 |
| Full / tentative | True | 0.00304915968 | 3.908e-15 | 2.133e-13 | 0.15548596 |
| Full / old_feasible | True | 0.00300782567 | 3.824e-15 | 5.099e-11 | 0.000163462675 |
| Full / fast | True | 1.59072191 | 2.901e-15 | 1.612e-10 | 1.86048615 |

This diagnostic distinguishes reaching the energy/momentum manifold from minimizing distance to the tentative state. No full-objective long trajectory or scalable local implementation is claimed. The observed local solution can motivate a safeguarded/fixed-center fallback, but this single-step result is not a validated general repair.

The lowest objective among these starts is the nearly stationary position with a substantial nonzero returned velocity. The different solutions remain locally stationary. Therefore replacing fast projection by a full local SQP solve is not yet a demonstrated motion fix: branch choice and position/velocity consistency still require testing.
