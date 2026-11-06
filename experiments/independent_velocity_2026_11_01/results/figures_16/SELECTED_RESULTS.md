# Independent velocity: selected measured runs

βᵥ=1, 60 Hz, 16 sweeps, epsilon ratio=1e+08.

Energy percentages normalize by initial mechanical energy above minimum potential. Actual E uses returned independent velocity; chord E substitutes the backward position difference. Their difference is diagnostic, not an imposed equality. Incomplete runs are prefixes.

| Case | Method | Complete? | Max actual E error | Max chord E error | Max v/chord mismatch | Tail travel | Min 1+βη |
|---|---|---|---:|---:|---:|---:|---:|
| spring_release | independent EL | 1 | 1.8485e-05% | 2.976% | 0.2762 m/s | 1.574 m | 0.4183 |
| spring_release | independent LE | 1 | 2.1727e-05% | 2.9508% | 0.2477 m/s | 1.499 m | 0.50071 |
| spring_release | independent JOINT | 1 | 0.00024042% | 2.6271% | 0.1491 m/s | 7.415e-08 m | -0.064222 |
| spring_release | derived-v EL | 0 | 3.5552% | 3.5552% | 0 m/s | 0.4934 m | — |
| spring_release | derived-v LE | 0 | 0.21378% | 0.21378% | 0 m/s | 0.4867 m | — |
| spring_orbit | independent EL | 1 | 2.4855e-07% | 0.89291% | 0.00885 m/s | 1.27 m | 0.98936 |
| spring_orbit | independent LE | 1 | 2.4855e-07% | 0.89291% | 0.00885 m/s | 1.27 m | 0.98936 |
| spring_orbit | independent JOINT | 1 | 2.4855e-07% | 0.89291% | 0.00885 m/s | 1.27 m | 0.98936 |
| spring_orbit | derived-v EL | 1 | 4.9744e-07% | 4.9744e-07% | 0 m/s | 1.249 m | — |
| spring_orbit | derived-v LE | 1 | 4.9744e-07% | 4.9744e-07% | 0 m/s | 1.249 m | — |
| spring_stiff | independent EL | 1 | 9.3707e-07% | 15.887% | 0.4346 m/s | 3.241 m | 0.79898 |
| spring_stiff | independent LE | 1 | 9.9288e-07% | 15.887% | 0.4346 m/s | 3.241 m | 0.79898 |
| spring_stiff | independent JOINT | 1 | 9.9285e-07% | 15.887% | 0.4346 m/s | 3.241 m | 0.79898 |
| spring_stiff | derived-v EL | 1 | 0.00076434% | 0.00076434% | 0 m/s | 3.191 m | — |
| spring_stiff | derived-v LE | 1 | 1.8972e-06% | 1.8972e-06% | 0 m/s | 3.191 m | — |
| spring_gravity | independent EL | 1 | 3.0553e-06% | 0.81323% | 0.03199 m/s | 1.261 m | 0.71644 |
| spring_gravity | independent LE | 1 | 3.0103e-06% | 0.81323% | 0.03199 m/s | 1.261 m | 0.71644 |
| spring_gravity | independent JOINT | 1 | 3.0077e-06% | 0.81323% | 0.03199 m/s | 1.261 m | 0.71644 |
| spring_gravity | derived-v EL | 1 | 1.6266e-05% | 1.6266e-05% | 0 m/s | 1.228 m | — |
| spring_gravity | derived-v LE | 1 | 6.4107e-06% | 6.4107e-06% | 0 m/s | 1.228 m | — |
| pendulum_release | independent EL | 1 | 0.079306% | 0.33144% | 0.04446 m/s | 1.249 m | 0.5031 |
| pendulum_release | independent LE | 1 | 0.00070544% | 0.31105% | 0.06392 m/s | 1.269 m | 0.51435 |
| pendulum_release | independent JOINT | 1 | 0.00014189% | 0.24639% | 0.03554 m/s | 0.5904 m | -0.15444 |
| pendulum_release | derived-v EL | 1 | 0.0027194% | 0.0027194% | 0 m/s | 1.839e-06 m | — |
| pendulum_release | derived-v LE | 1 | 0.15425% | 0.15425% | 0 m/s | 1.179 m | — |
| pendulum_kick | independent EL | 1 | 0.0059484% | 0.27204% | 0.04757 m/s | 1.065 m | 0.57728 |
| pendulum_kick | independent LE | 1 | 2.6025e-06% | 0.26722% | 0.04668 m/s | 1.065 m | 0.57137 |
| pendulum_kick | independent JOINT | 1 | 0.00039494% | 0.24214% | 0.03024 m/s | 0.3353 m | -0.19175 |
| pendulum_kick | derived-v EL | 1 | 0.0063766% | 0.0063766% | 0 m/s | 4.143e-08 m | — |
| pendulum_kick | derived-v LE | 1 | 0.14873% | 0.14873% | 0 m/s | 1.013 m | — |

Derived-velocity baselines are reused archived trajectories; they were not rerun.
