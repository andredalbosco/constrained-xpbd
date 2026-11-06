# Scalar exact-energy experiment results

91 analytic validation checks passed. 112 runs, requested duration 5 s.

Every E visit is followed immediately by the unchanged compliant XPBD L visit, with accumulated material multiplier. The exact-objective control solves the original scalar constrained objective directly.

## 480 Hz results

| Case | Target | Method | Sweeps | Status/time | Max energy error | Position RMS | Path length | Ambiguous roots |
|---|---|---|---:|---|---:|---:|---:|---:|
| release | initial | energy_linear | 1 | runaway_or_nonfinite, 4.2917s | 1.29317e+11 | 351.427 | 14863.2 | 0 |
| release | initial | energy_linear | 4 | complete, 5s | 81.1519 | 0.110867 | 8.00669 | 0 |
| release | initial | energy_linear | 16 | runaway_or_nonfinite, 0.46458s | 6.74957e+06 | 2.17266 | 37.3718 | 0 |
| release | initial | energy_exact | 1 | complete, 5s | 0.0210412 | 0.295019 | 0.499837 | 0 |
| release | initial | energy_exact | 4 | complete, 5s | 7.07701e-14 | 0.306998 | 8.85403e-15 | 1 |
| release | initial | energy_exact | 16 | complete, 5s | 0 | 0.306998 | 0 | 2400 |
| release | initial | objective_exact | 1 | complete, 5s | 0 | 0.306998 | 0 | 2400 |
| release | previous_step | energy_linear | 1 | runaway_or_nonfinite, 0.61042s | 7.39369e+06 | 206.294 | 530.578 | 0 |
| release | previous_step | energy_linear | 4 | complete, 5s | 145.867 | 2.20667 | 3.27954 | 0 |
| release | previous_step | energy_linear | 16 | runaway_or_nonfinite, 0.44792s | 0.476105 | 0.0688299 | 0.553738 | 0 |
| release | previous_step | energy_exact | 1 | complete, 5s | 0.943476 | 0.239787 | 0.620553 | 0 |
| release | previous_step | energy_exact | 4 | complete, 5s | 1.69757e-10 | 0.306998 | 2.12197e-11 | 2400 |
| release | previous_step | energy_exact | 16 | complete, 5s | 0 | 0.306998 | 0 | 2400 |
| release | previous_step | objective_exact | 1 | complete, 5s | 0 | 0.306998 | 0 | 2400 |
| mixed | initial | energy_linear | 1 | complete, 5s | 0.52019 | 0.271243 | 6.84864 | 0 |
| mixed | initial | energy_linear | 4 | runaway_or_nonfinite, 2.2396s | 2.57551e+10 | 143.774 | 3386.85 | 0 |
| mixed | initial | energy_linear | 16 | complete, 5s | 0.000166327 | 0.0788877 | 8.22021 | 0 |
| mixed | initial | energy_exact | 1 | complete, 5s | 0.011795 | 0.322188 | 0.0123173 | 0 |
| mixed | initial | energy_exact | 4 | complete, 5s | 9.94528e-13 | 0.322327 | 0.0124881 | 0 |
| mixed | initial | energy_exact | 16 | complete, 5s | 6.05869e-15 | 0.322327 | 0.0124881 | 2380 |
| mixed | initial | objective_exact | 1 | complete, 5s | 6.05869e-15 | 0.322327 | 0.0124881 | 2379 |
| mixed | previous_step | energy_linear | 1 | runaway_or_nonfinite, 1.2146s | 7.5181e+08 | 2158.8 | 6229.22 | 0 |
| mixed | previous_step | energy_linear | 4 | complete, 5s | 3.02632 | 0.429763 | 0.276625 | 0 |
| mixed | previous_step | energy_linear | 16 | runaway_or_nonfinite, 0.2875s | 6.11775 | 0.501816 | 0.450295 | 0 |
| mixed | previous_step | energy_exact | 1 | complete, 5s | 0.800026 | 0.259614 | 0.143974 | 0 |
| mixed | previous_step | energy_exact | 4 | complete, 5s | 1.77841e-10 | 0.322327 | 0.0124881 | 2380 |
| mixed | previous_step | energy_exact | 16 | complete, 5s | 2.6684e-14 | 0.322327 | 0.0124881 | 2380 |
| mixed | previous_step | objective_exact | 1 | complete, 5s | 2.6684e-14 | 0.322327 | 0.0124881 | 2379 |

Energy errors are fractions of initial physical energy, not percentages. Failed runs report their accepted prefix and must not be ranked as complete trajectories.

`summary.csv` contains every run; `states.csv` contains every accepted timestep; `visits_first8.csv` contains every row pair for the first eight timesteps of each run, including energy immediately after E and after L. Exact projections do not have a claimed KKT-consistent accumulated energy multiplier.
