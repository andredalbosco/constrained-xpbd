# Five-second work-reduction follow-up

21/24 requested trajectories have recorded terminal outcomes. All runs stop at the first rejected nonlinear step. Physical reference certification covers the full five-second horizon, independently of candidate completion.

| Method | Complete five seconds | Completed with certified reference | Recorded terminal outcomes |
|---|---:|---:|---:|
| reduced_GN_polish | 6/8 | 6/8 | 8/8 |
| full | 5 complete; 1 interrupted; 2 not started | 5 certified completions | 5/8 |
| reduced_exact | 3/8 | 3/8 | 8/8 |

The original 128/256 reference pair certifies 6/8 configurations. Both N32 stiffness×10 references fail over five seconds; their physical trajectory RMS remains **uncertified**, regardless of a candidate's solver convergence or energy error.

| Failed method | N / stiffness | h | Accepted time (s) | Reason | Rejected primal | Rejected stationarity |
|---|---|---:|---:|---|---:|---:|
| reduced_GN_polish | 32 / 10 | 0.0333333 | 1.56667 | kkt_line_search_stalled (128-evaluation budget exhausted) | 0.00358 | 0.16 |
| reduced_GN_polish | 32 / 10 | 0.0166667 | 4.25 | kkt_line_search_stalled (128-evaluation budget exhausted) | 0.0647 | 0.417 |
| reduced_exact | 8 / 10 | 0.0333333 | 0.0666667 | kkt_line_search_stalled (128-evaluation budget exhausted) | 0.123 | 0.116 |
| reduced_exact | 32 / 1 | 0.0333333 | 0.3 | kkt_line_search_stalled | 0.16 | 0.316 |
| reduced_exact | 32 / 1 | 0.0166667 | 0.15 | kkt_line_search_stalled (128-evaluation budget exhausted) | 0.0698 | 0.202 |
| reduced_exact | 32 / 10 | 0.0333333 | 0 | kkt_line_search_stalled (128-evaluation budget exhausted) | 4.67 | 17.5 |
| reduced_exact | 32 / 10 | 0.0166667 | 0 | kkt_line_search_stalled (128-evaluation budget exhausted) | 0.221 | 0.307 |

Targeted 512/1024-substep reference follow-up (the candidate-error baseline is still the original256-substep reference):

| N / stiffness | h | Position gap | Scaled velocity gap | Energy drift | Certified |
|---|---:|---:|---:|---:|---|
| 32 / 10 | 0.0333333 | 0.0256 | 0.395 | 3.41e-11 | 0 |
| 32 / 10 | 0.0166667 | 0.00682 | 0.0842 | 1.33e-12 | 0 |

First-second regression: 21/21 saved-prefix checks pass; 19 contain accepted states. Maximum scaled numeric discrepancy is 0. Zero-frame comparisons only check matching immediate failure. Per-step evaluation counts are compared separately and retained.

These are bounded follow-up observations. They do not establish long-duration universal success, a general-mesh result, or a runtime ranking. Full residuals, failed prefixes, and raw accuracy columns remain in the CSVs. See [REPRODUCE.md](REPRODUCE.md) for exact commands, dependencies, and the unchanged formulation.

## Wall-clock cutoff

All 8 GN-polish and all 8 exact-reduced configurations reached terminal outcomes. Five full controls completed (all four N8 cases and N32 stiffness1 at30Hz). Full N32 stiffness1 at60Hz was interrupted by the wall-clock cutoff before any accepted trace row was flushed; both full N32 stiffness10 controls were not started. These three missing outcomes are **not solver failures**. The original interrupted file is retained, but no unsaved state or accuracy is inferred. See `run_status.json` for every requested configuration.

The first full N32 terminal row arrived between an interim analysis and interruption; the final report was regenerated after all processes stopped. All21 terminal-outcome first-second comparisons pass exactly, including evaluation counts; two have no accepted states because of immediate solver failure. No prefix comparison is possible for the interrupted empty trace.
