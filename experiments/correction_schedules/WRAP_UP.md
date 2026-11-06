# Wrap-up state

The user requested wrapping up before the planned hour elapsed. The research
began at 10:55:53 UTC on 2026-10-25; wrap-up was requested around 11:22 UTC.

Completed primary work: 1,809 frozen schedule/filter runs, 270 selected trace
reruns, 1,728 parameter runs, 162 hundred-second runs, 16 fine-step holdouts,
90 circular-speed/direction holdouts, and 648 assembly stress fixtures. There
are nine primary figures and two circular-speed figures, all visually checked.

The final refreshed-row guard and its tests are complete. It compares the
gradient magnitude with a floating-point position-representability bound,
without changing the solve denominator or adding energy compliance. Guarded
`refresh.csv` has 1,728 runs and three safety stops. The unprotected dataset is
preserved as `refresh_unprotected.csv` (71 stops). These two datasets must not
be combined. The guarded refresh comparison is a follow-up, not part of the
1,809-run main matrix.

At 480 Hz, refreshed LE with whole-proposal repair/cap gives released-spring
RMS 0.00500 m and max energy error 0.146%; the frozen counterpart gives 0.00959
m and 0.154%. However eccentric-orbit RMS remains 0.9195 m, so refresh alone
does not fix mode redistribution. Other refreshed cases are in `refresh.csv`;
full plotting and equal-cost comparison remain unfinished.

Core tests, refreshed guard tests, and the stable-assembly probe passed both
optimized and Address/UndefinedBehavior Sanitizer builds. The saved-data audit
checks primary frozen runs, not all follow-up datasets. Primary traces and
reference convergence were audited independently.

No chain implementation or simulation was started. No background research
agents remain running. Next work should focus on branch selection and energy
distribution, then a baseline-validated 2- or 3-mass spring pendulum. A chain
graph preconditioner remains a proposal in LARGE_STEPS.md.
