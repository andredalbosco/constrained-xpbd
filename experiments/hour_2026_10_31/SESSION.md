# One-hour continuation

User request: continue exploring for one hour.

Started 2026-10-31 18:04:18 UTC. The requested hour ends at 19:04:18 UTC (20:04:18 Europe/Paris). Numerical exploration and final consolidation span that interval; build portability checks and result packaging conclude the session.

## Completed studies

- `local_solver`: corrected explicit-work objective, turning-point stalls, shared multipliers, larger chains, warm starts, stable incremental line-search differences and iteration spectra.
- `preconditioning`: local and graph coordinate metrics, setup-inclusive cost, adaptive inexact-Newton forcing, cached factors and fixed-old-state continuation.
- `rigid`: multi-link pendulums, old reaction in the predictor, tangency/energy/reversibility, modal transfer, canonical momentum controls and stable near-degenerate energy differences.
- `work_reduction`: exact work-variable elimination, component energy closure, tree operators, objective equivalence, nonlinear controls, physical rollouts, replication normalization and disconnected-system energy exchange.

The preceding round-two package is also consolidated and preserved. A final bounded five-second follow-up extends larger-chain GN-then-Newton coverage: 6/8 complete, with two later finite-budget failures and two uncertified stiff references. Exact scope and any interrupted controls are recorded in `work_reduction/longer_rollouts/RESULTS.md`. Failed solves, partial trajectories, physical-reference certification and independent review findings remain in the saved records. The hour CMake suite passed 11 checks and the round-two suite passed 4 checks in their staging locations. Both saved canonical packages were then freshly configured and built in `/tmp`; all 15 checks passed again. Receipts and logs are stored with the packages.

## Locations and presentation

Work was staged under `/home/blacknoir/Global/constrained-xpbd-round2-2026-10-31` and `/home/blacknoir/Global/constrained-xpbd-hour-2026-10-31`. Archive destinations are new directories under `/home/blacknoir/Projects/constrained-xpbd/experiments/`: `round2_2026_10_31` and `hour_2026_10_31`. Each archive includes a SHA-256 manifest. No main solver changes, commits or pushes are part of this session.

The browser connection was unavailable. Automatic approval review rejected the isolated Chrome fallback; interactive viewer controls remain unverified. A native pendulum comparison MP4 was exported and checked independently. It plays at normal speed with prominently labeled 16× rod-angle magnification.

Personal memories were not changed. Reports distinguish measured outcomes from proposed integrations.
