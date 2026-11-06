# Independent velocity: selected measured runs

βᵥ=1, 60 Hz, 4 sweeps, epsilon ratio=1e+08.

Energy percentages normalize by initial mechanical energy above minimum potential. Actual E uses returned independent velocity; chord E substitutes the backward position difference. Their difference is diagnostic, not an imposed equality. Incomplete runs are prefixes.

| Case | Method | Complete? | Max actual E error | Max chord E error | Max v/chord mismatch | Tail travel | Min 1+βη |
|---|---|---|---:|---:|---:|---:|---:|
| spring_release | independent EL | 1 | 4.2754% | 7.0382% | 0.26 m/s | 1.543 m | 0.46486 |
| spring_release | independent LE | 1 | 0.092374% | 2.8447% | 0.2475 m/s | 1.5 m | 0.046288 |
| spring_release | independent JOINT | 1 | 0.00030007% | 0.00030026% | 0.0001096 m/s | 7.514e-08 m | -0.00052481 |
| spring_release | derived-v EL | 1 | 2.5724e+07% | 2.5724e+07% | 0 m/s | 51.2 m | — |
| spring_release | derived-v LE | 1 | 1.4241e+07% | 1.4241e+07% | 0 m/s | 279.4 m | — |
| spring_orbit | independent EL | 1 | 6.2457e-05% | 0.8929% | 0.00885 m/s | 1.27 m | 0.98936 |
| spring_orbit | independent LE | 1 | 2.4855e-07% | 0.89291% | 0.00885 m/s | 1.27 m | 0.98936 |
| spring_orbit | independent JOINT | 1 | 2.4776e-07% | 0.89291% | 0.00885 m/s | 1.27 m | 0.98936 |
| spring_orbit | derived-v EL | 1 | 4.2984e-05% | 4.2984e-05% | 0 m/s | 1.249 m | — |
| spring_orbit | derived-v LE | 1 | 4.9744e-07% | 4.9744e-07% | 0 m/s | 1.249 m | — |
| spring_stiff | independent EL | 1 | 1.0385% | 15.661% | 0.4524 m/s | 3.257 m | 0.79052 |
| spring_stiff | independent LE | 1 | 0.0063669% | 15.872% | 0.4341 m/s | 3.241 m | 0.80164 |
| spring_stiff | independent JOINT | 1 | 0.26563% | 15.231% | 0.4276 m/s | 3.245 m | 0.80066 |
| spring_stiff | derived-v EL | 1 | 1.8885% | 1.8885% | 0 m/s | 3.161 m | — |
| spring_stiff | derived-v LE | 1 | 0.034499% | 0.034499% | 0 m/s | 3.192 m | — |
| spring_gravity | independent EL | 1 | 0.1975% | 0.80449% | 0.03214 m/s | 1.262 m | 0.70244 |
| spring_gravity | independent LE | 1 | 0.00023139% | 0.81323% | 0.03199 m/s | 1.261 m | 0.72273 |
| spring_gravity | independent JOINT | 1 | 0.0029095% | 0.8132% | 0.03197 m/s | 1.261 m | 0.71281 |
| spring_gravity | derived-v EL | 1 | 8.0819% | 8.0819% | 0 m/s | 1.287 m | — |
| spring_gravity | derived-v LE | 1 | 0.0012157% | 0.0012157% | 0 m/s | 1.228 m | — |
| pendulum_release | independent EL | 1 | 1.4864% | 1.4535% | 0.1501 m/s | 1.251 m | -0.03004 |
| pendulum_release | independent LE | 1 | 0.016161% | 0.27278% | 0.04893 m/s | 1.259 m | 0.37624 |
| pendulum_release | independent JOINT | 1 | 0.0021458% | 0.24432% | 0.029 m/s | 0.1616 m | -0.31692 |
| pendulum_release | derived-v EL | 1 | 9.4274e-05% | 9.4274e-05% | 0 m/s | 6.239e-08 m | — |
| pendulum_release | derived-v LE | 1 | 0.22804% | 0.22804% | 0 m/s | 1.237 m | — |
| pendulum_kick | independent EL | 1 | 1.434% | 1.388% | 0.1079 m/s | 1.02 m | 0.30991 |
| pendulum_kick | independent LE | 1 | 0.016004% | 0.31133% | 0.03717 m/s | 1.061 m | 0.47417 |
| pendulum_kick | independent JOINT | 1 | 0.00023179% | 0.13306% | 0.01073 m/s | 1.729e-07 m | -0.0015983 |
| pendulum_kick | derived-v EL | 1 | 10.081% | 10.081% | 0 m/s | 0.01266 m | — |
| pendulum_kick | derived-v LE | 1 | 0.37599% | 0.37599% | 0 m/s | 1.041 m | — |

Derived-velocity baselines are reused archived trajectories; they were not rerun.
