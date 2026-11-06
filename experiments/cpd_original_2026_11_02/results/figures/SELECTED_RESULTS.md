# CPD measured comparison

Selected 60 Hz, epsilon ratio 1e+08.

Actual energy is K+elastic−sum(f·x), while the algorithm constrains H=K+elastic using incoming external power. Policies and extensions remain explicit.

| Case | Family | Policy | Complete? | Max actual E error | Relative travel | Tail travel | Raw residual | KKT | Stationary steps | Failure |
|---|---|---|---|---:|---:|---:|---:|---:|---:|---|
| spring_release | anchored | paper | 1 | 0.0003% | 3.75e-07 | 7.5e-08 | 5.33e-15 | 1.28e-09 | 99.67% |  |
| spring_release | anchored | fixed | 1 | 0.0003% | 3.75e-07 | 7.5e-08 | 2.66e-15 | 8.53e-16 | 100.00% |  |
| spring_release | anchored | stationary | 1 | 0.0003% | 3.75e-07 | 7.5e-08 | 3.11e-15 | 9.95e-16 | 100.00% |  |
| free_release | free | paper | 1 | 0.0003% | 3.75e-07 | 5.3033e-08 | 1.33e-14 | 6.42e-10 | 100.00% |  |
| free_release | free | fixed | 1 | 0.0003% | 3.75e-07 | 5.3033e-08 | 2.66e-15 | 8.53e-16 | 100.00% |  |
| free_release | free | stationary | 1 | 0.0003% | 3.75e-07 | 5.3033e-08 | 1.33e-14 | 6.42e-10 | 100.00% |  |
| spring_orbit | anchored | paper | 1 | 1.7558% | 3.9973 | 0.79902 | 4.83e-05 | 0.000151 | 0.00% |  |
| spring_orbit | anchored | fixed | 1 | 0.00010821% | 4.0001 | 0.80002 | 4.99e-15 | 6.09e-15 | 100.00% |  |
| spring_orbit | anchored | stationary | 1 | 0.00010821% | 4.0001 | 0.80002 | 5.77e-15 | 3.03e-10 | 100.00% |  |
| free_orbit | free | paper | 1 | 1.7558% | 3.9973 | 0.56499 | 4.83e-05 | 7.57e-05 | 0.00% |  |
| free_orbit | free | fixed | 1 | 0.00010821% | 4.0001 | 0.5657 | 4.67e-15 | 5.7e-15 | 100.00% |  |
| free_orbit | free | stationary | 1 | 0.00010821% | 4.0001 | 0.5657 | 6.66e-15 | 1.52e-10 | 100.00% |  |
| spring_gravity | anchored | paper | 1 | 4.254% | 4.9998 | 0.9584 | 1.19e-06 | 1.24e-05 | 0.00% |  |
| spring_gravity | anchored | fixed | 1 | 4.2543% | 4.9998 | 0.9584 | 9.83e-15 | 3.01e-15 | 100.00% |  |
| spring_gravity | anchored | stationary | 1 | 4.2543% | 4.9998 | 0.9584 | 9.27e-15 | 1.12e-13 | 100.00% |  |
| free_opposite_load | free | paper | 1 | 4.254% | 4.9998 | 0.67769 | 1.19e-06 | 6.21e-06 | 0.00% |  |
| free_opposite_load | free | fixed | 1 | 4.2543% | 4.9998 | 0.67769 | 1.33e-14 | 3.39e-15 | 100.00% |  |
| free_opposite_load | free | stationary | 1 | 4.2543% | 4.9998 | 0.67769 | 1.33e-14 | 5.62e-14 | 100.00% |  |
| spring_stiff | anchored | paper | 1 | 0.3375% | 3.9997 | 0.7999 | 6.01e-05 | 0.000204 | 0.00% |  |
| spring_stiff | anchored | fixed | 1 | 0.00026422% | 4.0001 | 0.80001 | 3.12e-14 | 5.87e-15 | 100.00% |  |
| spring_stiff | anchored | stationary | 1 | 0.00026422% | 4.0001 | 0.80001 | 3.3e-14 | 1.05e-11 | 100.00% |  |
| free_stiff | free | paper | 1 | 0.3375% | 3.9997 | 0.56561 | 6.01e-05 | 0.000102 | 0.00% |  |
| free_stiff | free | fixed | 1 | 0.00026422% | 4.0001 | 0.5657 | 3.21e-14 | 6.04e-15 | 100.00% |  |
| free_stiff | free | stationary | 1 | 0.00026422% | 4.0001 | 0.5657 | 5.76e-14 | 5.27e-12 | 100.00% |  |
| free_triangle_spin | free | paper | 1 | 0.084915% | 8.9736 | 2.0254 | 9.92e-05 | 0.000155 | 0.00% |  |
| free_triangle_spin | free | fixed | 1 | 0.0078024% | 9.7559 | 2.0266 | 0.00028 | 0.000239 | 97.33% |  |
| free_triangle_spin | free | stationary | 1 | 2.6398e-05% | 9.6368 | 2.0215 | 6.05e-09 | 1.09e-06 | 99.00% |  |
| pendulum_release_compliant | anchored | paper | 1 | 623.6% | 6.6555 | 1.2891 | 12.4 | 694 | 0.00% |  |
| pendulum_release_compliant | anchored | fixed | 1 | 14.221% | 6.3593 | 1.254 | 0.171 | 1.95 | 98.67% |  |
| pendulum_release_compliant | anchored | stationary | 1 | 623.25% | 6.6555 | 1.2892 | 12.4 | 694 | 99.33% |  |
| pendulum_release_hard | hard_rod_extension | paper | 0 | 2.4751% | 0.1893 | 0.1893 | 5.13e-05 | 7.86e+05 | 0.00% | singular_conservation_schur |
| pendulum_release_hard | hard_rod_extension | fixed | 0 | 2.4675% | 0.1893 | 0.1893 | 0.000127 | 3.77e+05 | 0.00% | singular_conservation_schur |
| pendulum_release_hard | hard_rod_extension | stationary | 0 | 2.464% | 0.18931 | 0.18931 | 0.0001 | 1.77e+05 | 0.00% | singular_conservation_schur |
| pendulum_kick_compliant | anchored | paper | 1 | 6.4262% | 5.3379 | 1.0589 | 0.049 | 0.661 | 0.00% |  |
| pendulum_kick_compliant | anchored | fixed | 1 | 12.497% | 5.3029 | 1.0523 | 0.172 | 18.8 | 99.00% |  |
| pendulum_kick_compliant | anchored | stationary | 1 | 137.61% | 5.4005 | 1.0705 | 2.05 | 48.9 | 99.00% |  |
| pendulum_kick_hard | hard_rod_extension | paper | 0 | 4.9832% | 0.32778 | 0.32778 | 9.32e-05 | 3.61e+07 | 0.00% | singular_conservation_schur |
| pendulum_kick_hard | hard_rod_extension | fixed | 0 | 4.9689% | 0.30761 | 0.30761 | 0.000131 | 8.71e+05 | 0.00% | singular_conservation_schur |
| pendulum_kick_hard | hard_rod_extension | stationary | 0 | 4.976% | 0.3076 | 0.3076 | 7.49e-05 | 3.37e+05 | 0.00% | singular_conservation_schur |

A completed run has reached its requested horizon; that flag does not establish plausible motion or stationarity.
