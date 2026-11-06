# Local sweep stability diagnostic

**The N3 failure at `h=.2` is not explained solely by a poor initial guess.** The dense-oracle KKT root is linearly unstable for all three tested refreshed-local-correction orders. Near that root, material→storage→kinetic and storage→kinetic→material have spectral radius `1.165262`; reverse edge has `6.150813`. This is numerical local-stability evidence about this root and these orders, not a claim that no alternative root or stabilized iteration exists.

The diagnostic was derived in [STABILITY_MATH.md](STABILITY_MATH.md), then implemented in [stability_probe.cpp](stability_probe.cpp) without changing the engine. [Raw results](results/stability.csv) contain 108 rows: nine scene/timestep pairs, two methods, three schedules, and two finite-difference scales. Frozen correction was outside this focused probe.

## Refreshed local correction

Spectral radius of **one complete sweep**, at the oracle root:

| Scene | h | Material→storage→kinetic | Storage→kinetic→material | Reverse edge |
|---|---:|---:|---:|---:|
| Two springs | .02 | .862425 | .862425 | .840255 |
| Two springs | .08 | .861732 | .861732 | .841377 |
| Two springs | .20 | .806578 | .806578 | .786636 |
| N3 | .02 | .993248 | .993248 | .991970 |
| N3 | .08 | .991259 | .991259 | .990514 |
| N3 | .20 | **1.165262** | **1.165262** | **6.150813** |
| N8 | .02 | .993785 | .993785 | .992428 |
| N8 | .08 | .985084 | .985084 | **3.021608** |
| N8 | .20 | Unreliable FD | Unreliable FD | Nonfinite perturbed map |

Values below one indicate local attraction under the smooth-map assumptions. Values close to one explain slow asymptotic progress: N3 at `.02` needs about 148 sweeps for an e-fold reduction along its slow mode. Values above one indicate an expanding mode. At N3 `.20`, the pass orders' dominant eigenvalues are a complex conjugate pair approximately `−.056111 ± 1.16391i`; reverse edge has a real leading eigenvalue `−6.150813`.

The two pass schedules have the same spectrum because they cyclically rotate the same group products at a fixed point. This does not imply identical nonlinear trajectories or errors at a truncated visit budget. Reverse edge is a different product and revisits shared kinetic rows. Its sweep has 9, 14, and 39 visits for N2, N3, and N8 versus 8, 12, and 32 for a pass schedule. This table is not a matched-cost ordering comparison.

N8 at `.08` shows why ordering can change stability: both pass schedules are locally attractive, while reverse edge is repelling. Even the attractive pass maps have matrix 2-norms `5.88` and `6.87`, so some perturbations grow transiently. Local attraction does not certify the standard initial guess lies in the convergence region.

## Dropped correction

All 27 dropped-correction maps have spectral radius numerically one, with exactly 13, 20, and 55 eigenvalues near `+1` for N2, N3, and N8. These counts equal the normalized primal dimensions. They agree with the derived feasible-manifold identity: once all material/balance rows vanish, the dropped map fixes the state regardless of stationarity. This explains a stationarity plateau; it is not evidence of explosive linear instability.

## Validation and limits

All nine dense oracles converged, and their smallest reduced Hessian eigenvalues were positive; the minimum was `.173965`. A constrained local minimum can therefore repel this particular iteration. Across the 51 accepted maps, lifted root residuals were at most `1.90e-11`, the unperturbed sweep moved by at most `9.94e-11`, and the two finite-difference Jacobians differed relatively by at most `1.52e-8`. Spectral radii differed by at most `1.05e-7`. The independently assembled analytic derivative products agreed with finite differences to relative error `1.53e-8` or better. Scalar derivative/unit-mode identity checks also passed.

The three local-correction maps for N8 at `.20` are deliberately rejected as finite-difference measurements. The analytic row products indicate enormous amplification, and even the small oracle residual or rounding error becomes macroscopic during a full sweep. Two unperturbed sweeps violate the fixed-point-defect gate; the reverse-edge perturbation becomes nonfinite. Their raw values and statuses remain in the CSV, but the table does not promote them into precise spectral estimates.

Reproduce with a standalone C++17 build using the same Eigen/hour/round2/composition include paths as the existing CMake configuration, then run `stability_probe <output.csv>`. The measured run used `/tmp/xpbd_stability_probe results/stability.csv`. This probe adds no line search, damping, clipping, or solver fallback.
