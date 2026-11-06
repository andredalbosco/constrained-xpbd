# Measured results

Five-second runs. Energy is normalized by initial energy above the minimum potential. Incomplete runs remain explicitly failed; their prefix errors are not comparable to full-horizon accuracy. `raw_EL` refreshes the material row after E. FEPR uses the paper fast update and a separate hard-rod extension for rigid cases.

## spring_release

| Method | Hz | Sweeps | Full 5 s? | Max relative energy error | RMS position (m) | Path length (m) | Failure |
|---|---:|---:|---|---:|---:|---:|---|
| xpbd | 60 | 1 | yes | 0.9997 | 0.14207 | 1.9006 |  |
| raw_EL | 60 | 1 | NO | 306 | 0.2891 | 4.958 | runaway_guard |
| nonlinear_EL | 60 | 1 | yes | 0.1761 | 0.28818 | 0.49011 |  |
| nonlinear_EL | 60 | 8 | yes | 1.711e-10 | 0.307 | 2.1383e-11 |  |
| pd_true | 60 | 8 | yes | 0 | 0.307 | 0 |  |
| pd_surrogate | 60 | 1 | yes | 0 | 0.307 | 0 |  |
| pd_inertia | 60 | 1 | yes | 0 | 0.307 | 0 |  |
| fepr_step_target | 60 | 1 | yes | 2.362e-06 | 0.035862 | 7.9453 |  |
| fepr_base | 60 | 1 | yes | 0.9997 | 0.14207 | 1.9006 |  |
| dg | 60 | 1 | yes | 9.521e-14 | 0.011868 | 7.9721 |  |
| cpd | 60 | 1 | yes | 9.375e-06 | 0.307 | 1.1719e-06 |  |
| xpbd | 240 | 1 | yes | 0.8753 | 0.074167 | 4.9584 |  |
| raw_EL | 240 | 1 | NO | 1292 | 0.26455 | 0.375 | runaway_guard |
| nonlinear_EL | 240 | 1 | yes | 0.04246 | 0.29476 | 0.49935 |  |
| nonlinear_EL | 240 | 8 | yes | 0 | 0.307 | 0 |  |
| pd_true | 240 | 8 | yes | 0 | 0.307 | 0 |  |
| pd_surrogate | 240 | 1 | yes | 0 | 0.307 | 0 |  |
| pd_inertia | 240 | 1 | yes | 0 | 0.307 | 0 |  |
| fepr_step_target | 240 | 1 | yes | 8.89e-07 | 0.048396 | 7.9408 |  |
| fepr_base | 240 | 1 | yes | 0.8753 | 0.074167 | 4.9584 |  |
| dg | 240 | 1 | yes | 5.809e-13 | 0.00074363 | 7.9902 |  |
| cpd | 240 | 1 | yes | 3.75e-05 | 0.307 | 4.6874e-06 |  |

## spring_orbit

| Method | Hz | Sweeps | Full 5 s? | Max relative energy error | RMS position (m) | Path length (m) | Failure |
|---|---:|---:|---|---:|---:|---:|---|
| xpbd | 60 | 1 | yes | 0.5799 | 0.093692 | 4.3964 |  |
| raw_EL | 60 | 1 | yes | 0.06494 | 0.86234 | 6.2459 |  |
| nonlinear_EL | 60 | 1 | yes | 0.06478 | 0.86221 | 6.2458 |  |
| nonlinear_EL | 60 | 8 | yes | 1.511e-10 | 0.85661 | 6.2477 |  |
| pd_true | 60 | 8 | yes | 1.679e-14 | 0.84197 | 6.2372 |  |
| pd_surrogate | 60 | 1 | yes | 4.821e-06 | 0.84131 | 6.2367 |  |
| pd_inertia | 60 | 1 | yes | 1.056e-14 | 1.4715 | 1.2575 |  |
| fepr_step_target | 60 | 1 | yes | 3.693e-06 | 0.083826 | 5.323 |  |
| fepr_base | 60 | 1 | yes | 0.5799 | 0.093692 | 4.3964 |  |
| dg | 60 | 1 | yes | 5.632e-14 | 0.0044297 | 5.4328 |  |
| cpd | 60 | 1 | yes | 3.655e-09 | 0.84197 | 6.2372 |  |
| xpbd | 240 | 1 | yes | 0.4745 | 0.038151 | 4.8309 |  |
| raw_EL | 240 | 1 | yes | 0.02231 | 0.56449 | 5.9931 |  |
| nonlinear_EL | 240 | 1 | yes | 0.02231 | 0.56448 | 5.9931 |  |
| nonlinear_EL | 240 | 8 | yes | 2.181e-11 | 0.53489 | 5.9678 |  |
| pd_true | 240 | 8 | yes | 9.505e-14 | 0.49855 | 5.9353 |  |
| pd_surrogate | 240 | 1 | yes | 1.263e-09 | 0.49845 | 5.9352 |  |
| pd_inertia | 240 | 1 | yes | 1.164e-14 | 1.5703 | 0.72883 |  |
| fepr_step_target | 240 | 1 | yes | 3.808e-06 | 0.042194 | 5.4169 |  |
| fepr_base | 240 | 1 | yes | 0.4745 | 0.038151 | 4.8309 |  |
| dg | 240 | 1 | yes | 1.068e-12 | 0.00027757 | 5.4361 |  |
| cpd | 240 | 1 | yes | 0.0001101 | 0.49856 | 5.9353 |  |

## spring_gravity

| Method | Hz | Sweeps | Full 5 s? | Max relative energy error | RMS position (m) | Path length (m) | Failure |
|---|---:|---:|---|---:|---:|---:|---|
| xpbd | 60 | 1 | yes | 0.7825 | 0.11844 | 4.2576 |  |
| raw_EL | 60 | 1 | yes | 0.1465 | 0.37845 | 6.3013 |  |
| nonlinear_EL | 60 | 1 | yes | 0.1475 | 0.34276 | 6.3485 |  |
| nonlinear_EL | 60 | 8 | yes | 8.621e-10 | 0.30452 | 6.4348 |  |
| pd_true | 60 | 8 | yes | 1.585e-14 | 0.19396 | 6.6919 |  |
| pd_surrogate | 60 | 1 | yes | 2.036e-05 | 0.19353 | 6.6936 |  |
| pd_inertia | 60 | 1 | yes | 7.133e-15 | 0.50436 | 0.53389 |  |
| fepr_step_target | 60 | 1 | yes | 2.2e-06 | 0.13119 | 6.9228 |  |
| fepr_base | 60 | 1 | yes | 0.7825 | 0.11844 | 4.2576 |  |
| dg | 60 | 1 | yes | 6.736e-14 | 0.0070416 | 7.1052 |  |
| cpd | 60 | 1 | NO | 0.1249 | 0.24848 | 4.7197 | maximum_iterations |
| xpbd | 240 | 1 | yes | 0.5244 | 0.052933 | 5.7375 |  |
| raw_EL | 240 | 1 | yes | 0.04638 | 0.1236 | 6.9203 |  |
| nonlinear_EL | 240 | 1 | yes | 0.04638 | 0.12354 | 6.9204 |  |
| nonlinear_EL | 240 | 8 | yes | 2.825e-11 | 0.13441 | 6.8954 |  |
| pd_true | 240 | 8 | yes | 6.736e-14 | 0.085379 | 7.028 |  |
| pd_surrogate | 240 | 1 | yes | 2.165e-08 | 0.085345 | 7.028 |  |
| pd_inertia | 240 | 1 | yes | 8.718e-15 | 0.57901 | 0.2475 |  |
| fepr_step_target | 240 | 1 | yes | 5.092e-07 | 0.029346 | 7.0783 |  |
| fepr_base | 240 | 1 | yes | 0.5244 | 0.052933 | 5.7375 |  |
| dg | 240 | 1 | yes | 4.319e-13 | 0.00044107 | 7.1151 |  |
| cpd | 240 | 1 | yes | 0.03565 | 0.11256 | 6.9384 |  |

## pendulum_release

| Method | Hz | Sweeps | Full 5 s? | Max relative energy error | RMS position (m) | Path length (m) | Failure |
|---|---:|---:|---|---:|---:|---:|---|
| xpbd | 60 | 1 | yes | 0.6106 | 0.10746 | 5.0914 |  |
| raw_EL | 60 | 1 | yes | 0.02164 | 0.040273 | 6.4052 |  |
| nonlinear_EL | 60 | 1 | yes | 0 | 0.7409 | 0 |  |
| nonlinear_EL | 60 | 8 | yes | 0 | 0.7409 | 0 |  |
| pd_true | 60 | 8 | yes | 0 | 0.7409 | 0 |  |
| fepr_step_target | 60 | 1 | yes | 1.296e-06 | 0.0031457 | 6.4335 |  |
| fepr_base | 60 | 1 | yes | 0.0007507 | 0.00074187 | 6.4333 |  |
| dg | 60 | 1 | yes | 1.909e-14 | 0.00085673 | 6.4316 |  |
| cpd | 60 | 1 | NO | 0 | nan | 0 | maximum_iterations |
| xpbd | 240 | 1 | yes | 0.2201 | 0.032192 | 6.0433 |  |
| raw_EL | 240 | 1 | yes | 0.0007967 | 0.015046 | 6.4177 |  |
| nonlinear_EL | 240 | 1 | yes | 0 | 0.73977 | 0 |  |
| nonlinear_EL | 240 | 8 | yes | 0 | 0.73977 | 0 |  |
| pd_true | 240 | 8 | yes | 0 | 0.73977 | 0 |  |
| fepr_step_target | 240 | 1 | yes | 1.738e-07 | 0.000673 | 6.4335 |  |
| fepr_base | 240 | 1 | yes | 4.692e-05 | 4.6322e-05 | 6.4334 |  |
| dg | 240 | 1 | yes | 9.283e-12 | 5.3523e-05 | 6.4333 |  |
| cpd | 240 | 1 | NO | 0.01282 | 0.04178 | 1.3 | maximum_iterations |

## pendulum_kick

| Method | Hz | Sweeps | Full 5 s? | Max relative energy error | RMS position (m) | Path length (m) | Failure |
|---|---:|---:|---|---:|---:|---:|---|
| xpbd | 60 | 1 | yes | 0.596 | 0.08842 | 4.3444 |  |
| raw_EL | 60 | 1 | yes | 0.008706 | 0.10057 | 5.2349 |  |
| nonlinear_EL | 60 | 1 | yes | 0.006412 | 0.25797 | 4.8319 |  |
| nonlinear_EL | 60 | 8 | yes | 0.0002173 | 0.47461 | 4.6217 |  |
| pd_true | 60 | 8 | yes | 3.566e-15 | 0.61409 | 0.15922 |  |
| fepr_step_target | 60 | 1 | yes | 1.983e-06 | 0.0017082 | 5.4339 |  |
| fepr_base | 60 | 1 | yes | 0.0004189 | 0.00049959 | 5.4333 |  |
| dg | 60 | 1 | yes | 1.545e-14 | 0.00076463 | 5.4297 |  |
| cpd | 60 | 1 | NO | 0.04373 | 0.0057039 | 0.14486 | maximum_iterations |
| xpbd | 240 | 1 | yes | 0.2103 | 0.026264 | 5.1236 |  |
| raw_EL | 240 | 1 | yes | 0.0005512 | 0.022951 | 5.386 |  |
| nonlinear_EL | 240 | 1 | yes | 0.0004016 | 0.069316 | 5.2821 |  |
| nonlinear_EL | 240 | 8 | yes | 1.362e-05 | 0.14394 | 5.0997 |  |
| pd_true | 240 | 8 | yes | 2.14e-14 | 0.61347 | 0.15922 |  |
| fepr_step_target | 240 | 1 | yes | 2.978e-07 | 0.00034729 | 5.4324 |  |
| fepr_base | 240 | 1 | yes | 2.617e-05 | 3.1103e-05 | 5.4321 |  |
| dg | 240 | 1 | yes | 1.606e-11 | 4.7674e-05 | 5.4319 |  |
| cpd | 240 | 1 | NO | 0.01236 | 0.0015422 | 0.15559 | maximum_iterations |

## spring_stiff

| Method | Hz | Sweeps | Full 5 s? | Max relative energy error | RMS position (m) | Path length (m) | Failure |
|---|---:|---:|---|---:|---:|---:|---|
| xpbd | 30 | 1 | yes | 0.9416 | 0.15887 | 4.2385 |  |
| nonlinear_EL | 30 | 1 | yes | 0.7109 | 1.47 | 15.538 |  |
| nonlinear_EL | 30 | 8 | yes | 6.395e-05 | 1.4722 | 16.19 |  |
| pd_true | 30 | 8 | yes | 4.007e-15 | 1.4734 | 16.209 |  |
| fepr_step_target | 30 | 1 | yes | 9.107e-09 | 1.9137 | 10.74 |  |
| fepr_base | 30 | 1 | yes | 0.9416 | 0.15887 | 4.2385 |  |
| dg | 30 | 1 | yes | 8.348e-15 | 0.099267 | 10.187 |  |
| fepr_sub4 | 30 | 1 | yes | 9.317e-08 | 1.3631 | 18.375 |  |
| fepr_sub8 | 30 | 1 | yes | 2.629e-06 | 0.92044 | 15.844 |  |
| xpbd | 60 | 1 | yes | 0.935 | 0.11099 | 4.4028 |  |
| nonlinear_EL | 60 | 1 | yes | 0.3015 | 1.472 | 16.036 |  |
| nonlinear_EL | 60 | 8 | yes | 4.758e-09 | 1.4713 | 16.181 |  |
| pd_true | 60 | 8 | yes | 7.847e-15 | 1.4709 | 16.173 |  |
| fepr_step_target | 60 | 1 | yes | 5.315e-07 | 0.30853 | 10.501 |  |
| fepr_base | 60 | 1 | yes | 0.935 | 0.11099 | 4.4028 |  |
| dg | 60 | 1 | yes | 3.64e-14 | 0.10476 | 10.929 |  |
| fepr_sub4 | 60 | 1 | yes | 2.629e-06 | 0.9214 | 16.441 |  |
| fepr_sub8 | 60 | 1 | yes | 1.516e-06 | 0.11622 | 10.875 |  |

## pendulum_large

| Method | Hz | Sweeps | Full 5 s? | Max relative energy error | RMS position (m) | Path length (m) | Failure |
|---|---:|---:|---|---:|---:|---:|---|
| xpbd | 30 | 1 | yes | 0.9134 | 1.3269 | 9.3786 |  |
| nonlinear_EL | 30 | 1 | yes | 0.07614 | 0.98984 | 14.163 |  |
| nonlinear_EL | 30 | 8 | yes | 1.434e-05 | 1.2384 | 13.568 |  |
| pd_true | 30 | 8 | yes | 5.173e-17 | 1.3566 | 0.019064 |  |
| fepr_step_target | 30 | 1 | yes | 8.686e-08 | 0.028027 | 14.64 |  |
| fepr_base | 30 | 1 | yes | 0.005147 | 0.018839 | 14.642 |  |
| dg | 30 | 1 | yes | 1.791e-13 | 0.0010323 | 14.63 |  |
| fepr_sub4 | 30 | 1 | yes | 1.259e-07 | 0.0041946 | 14.63 |  |
| fepr_sub8 | 30 | 1 | yes | 7.614e-07 | 0.0018968 | 14.63 |  |
| xpbd | 60 | 1 | yes | 0.7642 | 1.2686 | 11.682 |  |
| nonlinear_EL | 60 | 1 | yes | 0.02047 | 0.5562 | 14.519 |  |
| nonlinear_EL | 60 | 8 | yes | 3.558e-06 | 0.75962 | 14.439 |  |
| pd_true | 60 | 8 | yes | 2.587e-16 | 1.3554 | 0.019064 |  |
| fepr_step_target | 60 | 1 | yes | 1.94e-07 | 0.01016 | 14.644 |  |
| fepr_base | 60 | 1 | yes | 0.001282 | 0.0046799 | 14.645 |  |
| dg | 60 | 1 | yes | 5.173e-15 | 0.00025806 | 14.642 |  |
| fepr_sub4 | 60 | 1 | yes | 7.614e-07 | 0.0018968 | 14.642 |  |
| fepr_sub8 | 60 | 1 | yes | 4.4e-07 | 0.0010042 | 14.642 |  |

## Scope

The single radial spring contour search currently supports the energy component not enclosing the anchor. Unsupported topology, missing contour/intersection, and nonfinite states report failure. These small cases do not establish robustness for coupled networks, contacts, or indefinite multibody energy surfaces. The strict PD controls are not the published relaxed-energy CPD algorithm. FEPR constraint convergence is not fixed-objective KKT convergence. Discrete-gradient controls change kinematics and use nonlinear root solves. No cross-method runtime claim is made from these Python prototypes.
