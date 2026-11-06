# Free coupled-spring momentum experiment

This experiment reuses `midpoint_gs.hpp` directly. It does not implement a second solver. `coupled_momentum_study.cpp` tests two initial motions, three timestep sizes, three solver variants, four outer iteration limits, and one or ten inner Gauss-Seidel sweeps: 144 ten-second runs. Every run completed and the independent momentum/residual bounds passed.

Build and run:

```sh
g++ -std=c++17 -O2 -Wall -Wextra -pedantic \
  -I/tmp/constrained-xpbd.PipwSV \
  -I/home/blacknoir/Projects/constrained-xpbd/external/eigen \
  coupled_momentum_study.cpp -o coupled_momentum_study
./coupled_momentum_study > coupled_momentum_run.log
```

`coupled_momentum_summary.csv` contains all results. Result opening/writing and finite-state checks are enforced. No project files were modified.

## Model

Three free particles have masses 1, 2, and 1.5, with positions (-0.7,-0.3,0), (0.9,-0.1,0.1), and (0.1,1,0.2). Three nonlinear distance springs form a triangle, with compliances 0.01, 0.015, and 0.008 and rest lengths respectively 0.9, 1.05, and 0.95 times their initial lengths. There are no anchors or external forces.

All particles initially share center-of-mass translation (0.4,-0.2,0.3) and rigid angular velocity (0.3,-0.5,1.1). The second case additionally gives particles 0 and 1 deformation velocities (0.2,0.05,-0.1) and (-0.1,-0.025,0.05), whose mass-weighted sum is zero. The nonzero initial spring strains cause deformation in both cases.

## Discrete invariants and the position residual

Let q_i=(x_i,+ + x_i,-)/2 and use the solver's scaled position residual

    rx_i = m_i(q_i - q_predictor,i) - sum_j D_ij lambda_j.

For no external forces, midpoint velocity reconstruction gives

    ΔP = (4/h) sum_i rx_i,
    ΔL = (4/h) sum_i q_i cross rx_i.

The midpoint discrete spring gradients are parallel to the difference q_a-q_b, so their pair impulses have zero total force and zero torque about q. The remaining kinematic term in ΔL vanishes because Δx_i=h vbar_i. Thus closing rx establishes linear and angular momentum conservation for this model, independently of whether a curvature matrix was used in the solver.

Each Gauss-Seidel pair correction is mass-weighted and opposite, so dropping rx still preserves total P to floating-point accuracy in these examples. It does not generally preserve L: directions from successive linearizations need not be central with respect to one common accepted midpoint configuration.

The test independently checks bounds on each step's P and L defect using the reported ||M^-1 rx||, particle masses, and accepted midpoint positions. Maximum total P drift across all 144 runs was 5.43e-11.

## Representative results: h=1/15, one inner GS sweep

Rigid translation/rotation initial velocity, simulated for 10 seconds. Errors are maxima relative to initial physical E, P, and L.

| Method | Outer limit | Energy error | Angular-momentum error | Scaled rx | Material residual |
|---|---:|---:|---:|---:|---:|
| drop rx | 1 | 0.479736 | 1.63e-12 | 1.43522e-4 | 0.0122757 |
| drop rx | 5 | 0.0268623 | 0.0148378 | 1.17451e-4 | 9.99e-12 |
| drop rx | 20 | 0.0268623 | 0.0148378 | 1.17451e-4 | 9.99e-15 |
| drop rx | 80 | 0.0268623 | 0.0148378 | 1.17451e-4 | 9.99e-15 |
| retain rx | 5 | 6.18e-8 | 4.39e-8 | 7.10e-11 | 5.52e-11 |
| retain rx | 20 | 8.36e-12 | 6.37e-12 | 8.98e-15 | 9.89e-15 |
| retain rx | 80 | 8.36e-12 | 6.37e-12 | 8.98e-15 | 9.89e-15 |

The dropped-residual method stops improving once the material rows converge. More outer iterations cannot close an equation that has been omitted. One outer sweep happens to preserve L because all update directions are frozen at the common predictor; its large energy error demonstrates that this is not full convergence.

## Representative results: 80 outer, 10 inner GS sweeps

| h | Method | Max E error | Max P error | Max L error | Max scaled rx |
|---|---|---:|---:|---:|---:|
| 1/240 | drop rx | 2.26081e-4 | 3.22e-11 | 1.30e-10 | 2.16e-9 |
| 1/240 | retain rx | 3.84e-10 | 1.39e-11 | 4.47e-10 | 9.96e-15 |
| 1/60 | drop rx | 0.00351879 | 3.80e-12 | 2.87e-8 | 5.40e-7 |
| 1/60 | retain rx | 3.66e-11 | 3.34e-12 | 2.79e-11 | 9.97e-15 |
| 1/15 | drop rx | 0.0346570 | 1.11e-12 | 6.69e-5 | 1.16e-4 |
| 1/15 | retain rx | 3.20e-12 | 4.38e-13 | 2.49e-12 | 9.64e-15 |

The smaller retained-rx errors at larger h reflect the fixed absolute residual stopping tolerance and the 4/h conversion from position residual to impulse; they are not evidence that larger steps are more trajectory-accurate. No trajectory-reference comparison is made by this experiment.

The added-deformation case gives the same distinction: at h=1/15 with 80 outer and 10 inner sweeps, dropping rx yields E error 0.0317583 and L error 6.30e-5; retaining rx yields E error 3.66e-12 and L error 3.67e-12.

## Independent root-solver review

I inspected `midpoint_gs_study.cpp` before helper extraction and compiled an independent binary, `midpoint_gs_study_review`, writing `midpoint_gs_review.csv`. All 108 original runs passed the source's assertions. Root subsequently added the one-inner-sweep sweep; the coupled experiment above independently exercises that case too.

The following equations and implementation details checked correctly:

- q is the spatial midpoint, with predictor x_old+h v_old/2+h²g/4.
- The compliant secant gradient is (r_new+r_old)/(|r_new|+|r_old|), and its material-row derivative with respect to q is the endpoint unit vector.
- Curvature dropping replaces the position Jacobian with M; it need not drop the positional residual.
- The frozen multiplier system is solved by actual incremental Gauss-Seidel updates, without Hessians or a dense solve.
- The kinetic-work defect is (4/h) sum_i vbar_i·rx_i; the program checks this identity separately from total energy drift.
- `retain_rx_symmetric` also approximates the material-row endpoint derivative by the discrete gradient. That is an additional Jacobian approximation, not merely curvature dropping, but it retains the same nonlinear residual and can converge to the same equations.

No concrete equation or loop bug was found. Hard-distance claims require a consistent initial constraint: the averaged hard row enforces C_new=-C_old, so exact length preservation follows from C_old=0; finite-iteration violations must continue to be reported.
