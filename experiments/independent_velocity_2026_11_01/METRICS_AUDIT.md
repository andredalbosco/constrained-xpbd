# Independent metrics audit

All saved trajectory mechanical quantities and summary statistics are recomputed
from CSV positions and independent velocities. The audit never calls the
production energy or increment functions, and it does not execute a rollout.

## What was checked

- Actual energy uses one half m times returned velocity squared plus the radial
  spring and gravitational potentials. Chord energy instead uses position
  displacement divided by h. Both are checked against their explicitly labeled
  columns; nonzero disagreement is allowed by the compliant formulation.
- Coupling mismatch, its h-scaled displacement measure, implicit coupling
  multipliers, radial/tangential kinetic energy, angular momentum, geometry,
  velocity/chord tangency, endpoint motion, and changes of velocity.
- Every recorded position and velocity stationarity residual and material,
  relaxation, and energy residual, reconstructed from the incoming accepted state.
  Energy targets use incoming actual energy and remain fixed during each solve;
  epsilon uses the initial energy scale and remains fixed throughout a run.
- Summary energy extrema, final energies, accepted-prefix path and tail travel,
  longest stalls, endpoint jumps, velocity disagreement, geometry/tangency,
  angular drift, and convergence fractions. Completed runs use the requested fixed
  visit budget. Completion does not imply residual convergence or plausible motion.
- The explicit guard conditions are checked on accepted states. Rejection and
  failed-candidate summary semantics are audited in the source and reconstructed
  when failed rows exist. This campaign contains no rejected rows, so it supplies
  no runtime failing-candidate example.

## Population conventions

`max_rel_energy` includes a retained rejected candidate when present.
`max_accepted_rel_energy`, most other maxima, and motion statistics cover accepted
states. `final_rel_energy` and `final_chord_rel_energy` describe the last accepted
state; `failed_trial_rel_energy` separately describes a rejected candidate.
This differs from the earlier position-only runner's final-energy convention.

The tail excludes the row exactly one second before the final accepted time.
Stall speed and tail-motion thresholds are descriptive diagnostics, not physical
acceptance requirements. A nonzero compliant velocity mismatch is not itself a
failure, and chord energy is not substituted for actual mechanical energy.

Serialized first-step traces are audited by independent dense block replay when
present. The primary campaign did not request `--trace-first`, so its saved trace
dictionary is empty and this receipt reports zero trace visits. The separate
solver validation covers 630 dense-replayed refreshed visits; this saved-data
audit does not claim those as serialized campaign traces.

The receipt also counts accepted rows with small or negative `1+beta_v*eta`.
These are diagnostics of the velocity stationarity relation, not an additional
acceptance rule.

Run with the experiment's NumPy environment:

```bash
OPENBLAS_NUM_THREADS=1 python3 metrics_audit.py --write-report
```

Use `--results /path/to/saved/results` for another saved campaign. The summary
block-displacement and block-velocity maxima, rank-defect visit counts, and CPU
timings cannot be reconstructed from endpoint CSVs alone; they remain recorded
solver diagnostics rather than independently verified endpoint statistics.

## Numerical receipt

```json
{
  "abs_one_plus_beta_eta_below_0p05_rows": 19608,
  "checks": 9377856,
  "counts": {
    "completed": 648,
    "completed_with_unconverged_steps": 511,
    "failed": 0,
    "runs": 648
  },
  "largest_scaled_errors": {
    "completed fixed budget visits": 0.0,
    "energy normalization above initial potential floor": 0.0,
    "fixed dimensional epsilon": 0.0,
    "fixed initial energy": 0.0,
    "potential energy floor": 0.0,
    "reconstructed physical r_eta": 3.552713678800501e-15,
    "reconstructed physical r_lambda": 0.0,
    "reconstructed physical r_rho": 8.881784197001252e-16,
    "reconstructed physical r_v_norm": 0.0,
    "reconstructed physical r_x_norm": 0.0,
    "reconstructed physical residual": 1.7780915628762273e-17,
    "summary accepted_time": 0.0,
    "summary angular_conservation_expected": 0.0,
    "summary completed": 0.0,
    "summary converged_fraction": 0.0,
    "summary failed_trial_rel_energy": 0.0,
    "summary final_angular_momentum": 0.0,
    "summary final_chord_rel_energy": 0.0,
    "summary final_radial_kinetic": 0.0,
    "summary final_rel_energy": 0.0,
    "summary final_tangential_kinetic": 0.0,
    "summary initial_angular_momentum": 0.0,
    "summary longest_stall_duration": 0.0,
    "summary max_abs_eta": 0.0,
    "summary max_abs_rho": 0.0,
    "summary max_accepted_chord_rel_energy": 0.0,
    "summary max_accepted_chord_velocity_jump": 0.0,
    "summary max_accepted_displacement": 0.0,
    "summary max_accepted_rel_energy": 0.0,
    "summary max_accepted_velocity_jump": 0.0,
    "summary max_angular_error": 0.0,
    "summary max_chord_tangency_error": 0.0,
    "summary max_displacement_velocity_mismatch": 0.0,
    "summary max_geometry_error": 0.0,
    "summary max_radius": 0.0,
    "summary max_rel_energy": 0.0,
    "summary max_residual": 0.0,
    "summary max_tangency_error": 0.0,
    "summary max_velocity_chord_mismatch": 0.0,
    "summary min_radius": 0.0,
    "summary path_length": 3.4023471759083023e-16,
    "summary position_velocity_travel_ratio": 4.440892098500626e-16,
    "summary reversed_velocity_fraction": 0.0,
    "summary rms_velocity_chord_mismatch": 0.0,
    "summary tail_path_length": 3.495833255062134e-16,
    "summary tail_stalled": 0.0,
    "summary tail_velocity_travel": 5.900258821223772e-16,
    "summary velocity_travel": 3.327818153947864e-16,
    "trajectory angular_momentum": 0.0,
    "trajectory chord_energy": 2.8447315718998525e-16,
    "trajectory chord_kinetic": 0.0,
    "trajectory chord_rel_energy": 7.927686285214008e-16,
    "trajectory chord_speed": 0.0,
    "trajectory chord_tangency": 0.0,
    "trajectory chord_velocity_jump": 0.0,
    "trajectory chord_vx": 0.0,
    "trajectory chord_vy": 0.0,
    "trajectory displacement": 0.0,
    "trajectory displacement_velocity_mismatch": 0.0,
    "trajectory energy": 4.4284795610783046e-16,
    "trajectory geometry": 0.0,
    "trajectory kinetic": 0.0,
    "trajectory pi_x": 0.0,
    "trajectory pi_y": 0.0,
    "trajectory potential": 2.220446049250313e-16,
    "trajectory radial_kinetic": 2.1320289155595755e-16,
    "trajectory radius": 0.0,
    "trajectory rel_energy": 1.5855372570428017e-15,
    "trajectory speed": 0.0,
    "trajectory tangency": 0.0,
    "trajectory tangential_kinetic": 8.881784197001252e-16,
    "trajectory velocity_chord_dot": 0.0,
    "trajectory velocity_chord_mismatch": 0.0,
    "trajectory velocity_jump": 0.0
  },
  "min_abs_one_plus_beta_eta": 3.999072184512897e-09,
  "negative_one_plus_beta_eta_rows": 15091,
  "numerical_solver_executed": false,
  "serialized_trace_runs": 0,
  "serialized_trace_visits": 0,
  "source": "/home/blacknoir/Global/constrained-xpbd-independent-velocity-2026-11-01/results",
  "status": "passed",
  "trace_first_requested": false,
  "trajectory_rows": 292248
}
```
