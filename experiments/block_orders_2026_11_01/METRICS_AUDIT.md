# Metrics and reference audit

The complete campaign contains 468 runs and 208,735 trajectory rows. Independent
recomputation of the summary values and CSV physical quantities passed, including
the 5,292 serialized first-step visits. No metric or serialization bug was found.
The audit did not modify the solver, runner, or saved result data.

## What the summary columns measure

| Column or group | Population and interpretation |
|---|---|
| `completed` | All requested timestep budgets were advanced without a guard/rank failure. This does not assert convergence or good motion. |
| `accepted_time` | Length of the accepted prefix; excludes a failing candidate. |
| `rms_position` | RMS against the continuous reference over accepted endpoint positions, excluding the initial row and any failing candidate. A failed run's RMS is a prefix measurement. |
| `path_length`, `max_angular_error` | Accepted prefix only. Angular error is measured for nongravity springs, where angular momentum is physically conserved. |
| `converged_fraction` | Fraction of accepted steps meeting the final residual tolerance within the fixed budget. EL/LE/JOINT test all stated equations; the XPBD control tests position/material equations only. |
| `max_rel_energy` | Maximum absolute physical-energy error relative to the initial energy, normalized by initial energy above the potential floor. Includes a retained failing candidate. |
| `final_rel_energy` | Energy error of the last stored candidate, which can have `accepted=0`. |
| Other maximum residual, geometry, tangency, rho, eta diagnostics | Include the retained failing candidate. |
| `epsilon` | Ratio times the initial energy-above-floor scale; fixed for the entire run. It is not rescaled each timestep. |

Each timestep's relaxed energy equation uses the incoming accepted energy as its
target, held fixed during that solve. Therefore satisfaction of the relaxed
equation, or a small rho, does not establish conservation relative to the start
of the rollout. The CSV energy columns measure that separate cumulative drift.

The fixed-budget method can accept a finite step with substantial residual or
physical error. The runaway guard is a diagnostic cutoff, not a stability test.

## Saved-data checks

- 463 runs complete their requested budgets; 5 terminate with `runaway_guard`.
- Every failed trajectory retains exactly one rejected last candidate, marked
  `accepted=0`, with an explicit reason and the correct accepted-prefix time.
- 297 of the completed runs have at least one accepted step that fails the
  requested convergence tolerance. Completion and convergence remain separate.
- Independently recomputed RMS, path, angular error, acceptance time, convergence
  fractions, physical energy, initial-relative energy, geometry, tangency, and
  maximum diagnostics match the stored summaries.
- Independent dense-block replay of all 5,292 stored first-step visits reproduces
  the serialized states to a maximum scaled discrepancy of 1.02e-11 and physical
  residuals to 2.72e-10. This includes EL, LE, JOINT, and XPBD controls.

The failed runs are all spring-release cases at epsilon ratio 1e8:

| Order | Hz | Sweeps | Accepted prefix |
|---|---:|---:|---:|
| EL | 60 | 16 | 0.400000 s |
| LE | 60 | 16 | 0.383333 s |
| EL | 120 | 1 | 0.025000 s |
| EL | 120 | 16 | 0.083333 s |
| LE | 120 | 16 | 0.016667 s |

Some accepted prefixes already contain large errors. A small RMS over a short
prefix, conversely, must not be presented as a successful five-second result.

## Continuous reference cross-check

The new `reference.py` was compared at runtime with the earlier independently
validated `references.py` from
`experiments/nonlinear_energy_2026_11_01/` in the canonical constrained-xpbd
repository.
All six scenes were checked over 0.5 seconds at both 60 and 120 Hz.

| Quantity | Maximum difference across all twelve comparisons |
|---|---:|
| Position | 8.54e-12 m |
| Velocity | 2.37e-11 m/s |

Both implementations integrate the same physical ODE with DOP853; this is a
runtime cross-check against the earlier verified implementation, not a different
continuous model or a proof of global accuracy. The campaign also saves its own
five-second tolerance-refinement and energy-drift checks in
`results/reference_audit.csv`.

## Reproduce

After the campaign finishes, run from the experiment directory using its
NumPy/SciPy environment:

```bash
OPENBLAS_NUM_THREADS=1 python3 metrics_audit.py
```

The script resolves saved results relative to itself and modifies no result
files. It looks for the earlier verified reference beside this experiment, then
at its canonical repository path, then in a sibling staging experiment. It fails
explicitly if that comparison source is unavailable. Its JSON output identifies
the source actually used and contains the maximum differences.
