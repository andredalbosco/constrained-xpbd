# Independent FEPR audit

540 checks; 0 failures.

Read the primary-paper text through Equations 1, 8, 12, 13, and 14. The implementation matches the paper fast projection in planar form: four energy/momentum rows, independent x and v, unrestricted auxiliary scalars, D=diag(mI,h²mI,epsilon,epsilon), and a current-iterate distance center. It does not implement fixed-center full SQP.

| Check | Maximum measured error |
|---|---:|
| independent_constraint_rows | 2.22044605e-16 |
| finite_difference_jacobian | 2.31200031e-10 |
| old_state_per_step_feasibility | 2.2887834e-16 |
| fixed_target_old_state_not_automatically_feasible | 8.78654212e-16 |
| fast_update_full_increment_KKT | 1.13991472e-12 |
| iteration_count_no_hidden_extra_solve | 0 |
| default_metric_override_epsilon | 0 |
| finite_iteration_limit_explicit | 0 |
| fast_distinct_from_fixed_center_SQP | 0 |
| spring_release_converges | 0 |
| auxiliary_s_not_clipped | 0 |
| anchored_linear_momentum_not_forced_old | 0 |
| velocity_independent_after_FEPR | 0 |
| paper_raw_L1_stopping | 0 |
| singular_explicit_failure | 0 |
| singular_regularization_count | 0 |
| singular_return_not_hidden_fallback | 0 |
| zero_iteration_cap_explicit | 0 |
| feasible_free_drift_unchanged | 2.77555756e-16 |
| feasible_free_drift_zero_projection | 0 |
| rattle_base_length | 2.22044605e-16 |
| rattle_base_tangency | 2.22044605e-16 |
| rattle_base_reversibility | 6.22852247e-14 |
| hard_rod_extension_reports_identity | 0 |
| hard_rod_extension_length | 2.49289478e-12 |
| hard_rod_extension_tangency | 9.84334836e-12 |
| hard_rod_extension_energy | 3.22994964e-11 |

The additional rod length and tangency rows, together with the RATTLE base step, are an explicitly named extension. Passing this audit validates the row formulas and numerical update, not a claim that the paper evaluated hard rods. Finite-difference tests include nonzero s/t, and the KKT comparison checks two sequential updates so fixed-center SQP cannot pass accidentally. Raw dimensional residual and determinant thresholds match the declared policy; they are not unit-invariant.
