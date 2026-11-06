# Independent planar validation

2934 checks; 0 failures.

The nonlinear E visit is a numerical global closest-point search on the exact physical energy contour, using the single-particle mass metric. It is not a root along a frozen energy-gradient line.

| Check | Maximum measured error |
|---|---:|
| radial_contour_true_energy | 5.4058843e-14 |
| grid32_vs128_objective | 5.5511151e-17 |
| closest_true_energy | 4.8102532e-14 |
| rod_joint_length | 4.4408921e-16 |
| rod_joint_energy | 4.0412118e-14 |
| rod_rest_tie | 0 |
| rod_joint_objective | 0 |
| pd_unconstrained_implicit_minimum | 4.460365e-16 |
| pd_true_global_objective | 0 |
| pd_true_feasible | 1.3899992e-13 |
| pd_surrogate_frozen_energy | 1.0746959e-13 |
| pd_surrogate_overestimates_true | 0 |
| pd_surrogate_global_objective | 0 |
| scalar_nonlinear_position | 7.1917852e-09 |
| scalar_nonlinear_velocity | 8.6301423e-07 |
| raw_planar_historical_x | 1.2267964e-13 |
| raw_planar_historical_v | 5.8889782e-11 |
| rod_predictor_zero_energy_gradient | 0 |
| raw_rod_historical_x | 2.4825342e-16 |
| raw_rod_historical_v | 1.1916164e-13 |

All randomized cases are within the explicitly supported contour component that does not enclose the origin. Grid32 versus grid128 compares attained objective values, not a general proof of global minimization. The scalar comparison permits representation and angular search error and uses only the first two timesteps.

The planar raw EL implementation refreshes the material row after the energy displacement. Its rod control therefore compiles the original traditional_cases/relinearized.hpp implementation, not the frozen traditional_cases/solver.hpp schedule. These differ for a changing rod normal. For the first two positive-branch scalar spring steps, the material row is affine and the frozen and refreshed forms agree.
