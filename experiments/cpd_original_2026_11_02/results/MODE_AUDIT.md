# Radius and angular-row audit

All values below reuse measured trajectories at60Hz and epsilon ratio1e8. Relative radius of each free pair is its spring separation; no continuous-reference score is used.

| Full free-body case | Policy | 4–5s radius range | Final radius | Final radial KE | Angular drift over5s | Max actual E error |
|---|---|---:|---:|---:|---:|---:|
| free_orbit | paper | 0.00030114m | 1.10151496m | 1.3154e-05J | 2.0317e-14 | 1.7558% |
| free_orbit | fixed | 1.8841e-08m | 1.09999868m | 1.1755e-05J | 9.6256e-14 | 0.00010821% |
| free_orbit | stationary | 1.8841e-08m | 1.09999868m | 1.1755e-05J | 1.5765e-14 | 0.00010821% |
| free_stiff | paper | 3.6139e-05m | 1.10018032m | 1.1919e-05J | 5.3513e-14 | 0.3375% |
| free_stiff | fixed | 2.8278e-08m | 1.09999977m | 1.1755e-05J | 5.4956e-14 | 0.00026422% |
| free_stiff | stationary | 2.8278e-08m | 1.09999977m | 1.1755e-05J | 2.6201e-14 | 0.00026422% |

The orbiting masses keep rotating, but fixed16/full-KKT policies leave the stretched spring radius nearly constant. Correct energy and angular momentum therefore coexist with frozen radial motion.

## Angular-row removal

Kstar remains unchanged in this ablation. This isolates removal of the angular constraint; it is not the full original CPD method.

| Case | L row | Initial angular momentum | Final angular momentum | Max actual E error | Final tangential KE |
|---|---|---:|---:|---:|---:|
| spring_orbit | 1 | 0.88 | 0.88 | 0.00010821% | 0.320001J |
| spring_orbit | 0 | 0.88 | 1.29079 | 5.6173e-08% | 0.804587J |
| free_orbit | 1 | 0.88 | 0.88 | 0.00010821% | 0.320001J |
| free_orbit | 0 | 0.88 | 1.29079 | 5.6175e-08% | 0.804587J |
| spring_stiff | 1 | 0.88 | 0.88 | 0.00026422% | 0.32J |
| spring_stiff | 0 | 0.88 | 3.27786 | 5.1574e-07% | 5.26161J |
| free_stiff | 1 | 0.88 | 0.88 | 0.00026422% | 0.32J |
| free_stiff | 0 | 0.88 | 3.27786 | 5.1574e-07% | 5.26161J |
| free_release | 1 | 0 | 0 | 0.0003% | 0J |
| free_release | 0 | 0 | 0 | 0.0003% | 0J |

Removing the angular row permits spin-up while maintaining energy. Keeping it suppresses that spin-up but does not restore spring extension oscillation.
