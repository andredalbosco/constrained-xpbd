# Direct global energy: measured update-ordering results

Changing the order changes the failure mode; it does not by itself produce correct energy-conserving dynamics. Hard joint solves can freeze a stretched spring exactly. Immediate material-then-energy updates move that spring but can grow catastrophically. Positive energy softness removes all guards in this matrix and gives useful pendulum trajectories, while allowing energy error. Retaining the position residual enables full stationarity convergence in the moving-spring control; even that convergence selects physically poor motion.

## Experiment and interpretation

These runs use direct global endpoint mechanical energy, no reservoirs, and backward-Euler velocity v_new=(x_new-x_old)/h. The target is the initial mechanical energy for the entire run. The material and softened energy rows are

```text
rC = C(x) + alpha/h² * lambda
rE = E(x, (x-x_old)/h) - E_initial + beta/h² * eta
R  = M(x-predictor) - grad(C)*lambda - grad(E)*eta
```

Physical energy error is H=E-E_initial, not the softened-row residual rE. Reported position residual is |R|/h². Energy errors are divided by Escale: initial spring elastic plus kinetic energy, or pendulum energy above the lowest gravitational potential. All model masses and rest lengths are one. For the hard pendulum rC is the actual length error; for a compliant spring it is the material equation residual, not geometric stretch.

| Method | One outer sweep |
|---|---|
| baseline | Ordinary XPBD material update. |
| joint | Coupled material/energy block at the current linearization. |
| immediate | Material update, then recompute and update the energy row. |
| frozen | Freeze both rows and derivatives, perform 1, 5, or 20 scalar passes on that one linear system, then apply the accumulated position update. |
| reverse | Energy update, then recompute and update the material row. |
| symmetric | Full material, energy, material updates, recomputing each row. This label does not assert second-order accuracy. |

Joint and frozen separately retain or drop R. Curvature terms are dropped in every tested update. There is no line search or hidden denominator regularization. The guard terminates a run without clipping or accepting the failed trial.

The 1020 runs cover five initial states, h=1/240 and 1/60, 1/5/20 outer sweeps, and beta/Escale=0,0.01,1. Each requests 10 seconds. The final results contain 930 ordinary completions, 48 completions with incompatible local updates, and 42 safety-guard aborts. All 660 positive-softness runs complete; all 42 guards occur in the 330 hard-energy runs. All recorded numeric fields are finite. Completion and finiteness do not establish acceptable dynamics.

## Hard energy: beta=0

The following slice uses h=1/60 and five sweeps, dropping R. RMS position error compares with the physical reference over the accepted interval. Guarded rows are explicitly marked and are not ten-second results.

| Case | Method | Max absolute H / Escale | RMS position error | Max position residual | Accepted duration / behavior |
|---|---|---:|---:|---:|---|
| spring_rest | joint | 0 | 0.307 | 0 | 10s; exactly frozen |
| spring_rest | immediate | 2.0039e+10 | 850.14 | 5.7252e+10 | 9.45s; guard |
| spring_radial | joint | 7.2936e-12 | 0.31305 | 66.895 | 10s; path only 0.5096 |
| spring_radial | immediate | 2.33e+08 | 186.6 | 6.8268e+10 | 4.283s; guard |
| spring_tangent | joint | 2.9402e-15 | 1.5103 | 3.2432 | 10s; wrong trajectory |
| spring_tangent | immediate | 4.4103e-15 | 1.5075 | 0.41324 | 10s; wrong trajectory |
| pendulum_rest | joint | 0.59878 | 0.20211 | 1.0075e+05 | 10s; incompatible local updates |
| pendulum_rest | immediate | 0.00064189 | 0.04254 | 7.6742 | 10s; length drift |
| pendulum_tangent | joint | 49.003 | 0.45885 | 4.6986e+06 | 10s; large energy spikes |
| pendulum_tangent | immediate | 0.00075388 | 0.11918 | 13.106 | 10s; length drift |
| pendulum_tangent | symmetric | 0.0025685 | 0.12007 | 9.1369 | 10s; energy drift |

The resting spring has zero motion, zero energy error, and zero full stationarity residual under hard joint updates, with either R choice. Across the matrix, 48 completed hard joint/frozen runs exhibit this exact freeze. This is a formulation degeneracy: from rest without gravity, the backward-Euler objective is the endpoint mechanical energy. Constraining that same energy makes the objective constant on the admissible energy surface. Full convergence therefore cannot choose the intended physical advance.

The tangent spring shows a different failure. Joint and immediate methods preserve energy to roundoff but produce RMS position errors near 1.51. Angular momentum about the fixed spring anchor starts at 1.25 and ends near 2.78223; the physical central spring should preserve it. Ordinary XPBD instead ends at 1.04017. The direct energy correction transfers energy into rotation while radial oscillation decays. Thus even excellent scalar energy numbers are not evidence of correct motion.

Hard pendulum results trade constraint errors. Immediate updates have maximum length errors 8.3595e-4 from rest and 1.0640e-3 with tangent velocity. The symmetric material-energy-material order restores length to roundoff but has maximum normalized energy errors 0.0030746 and 0.0025685. The five-sweep tangent joint trajectory has 29 timesteps where the exact hard endpoint length/energy intersection is infeasible, and energy spikes up to 49.003 times Escale.

![Hard energy ordering](../../results/energy_ordering/analysis/hard_energy_ordering.png)

## Positive energy softness

For beta/Escale=0.01 at the same timestep and five sweeps, the pendulum runs improve substantially. These are completed ten-second trajectories, not partial prefixes.

| Case | Method | Max absolute H / Escale | RMS position error | Max length error | Max position residual |
|---|---|---:|---:|---:|---:|
| pendulum_rest | joint | 0.0066351 | 0.06232 | 2.2204e-16 | 8.1951 |
| pendulum_rest | immediate | 0.0055969 | 0.046964 | 2.7596e-05 | 1.0918 |
| pendulum_rest | symmetric | 0.0056887 | 0.046972 | 1.1102e-16 | 1.0874 |
| pendulum_tangent | joint | 0.0060314 | 0.036991 | 2.2204e-16 | 6.0695 |
| pendulum_tangent | immediate | 0.0054824 | 0.025161 | 8.6828e-06 | 1.2899 |
| pendulum_tangent | symmetric | 0.0055036 | 0.025166 | 1.1102e-16 | 1.2884 |

For comparison, ordinary XPBD has normalized energy errors 0.83630 / 0.84357 and RMS position errors 0.18971 / 0.24729 for the resting/tangent pendulum. Softness 1 gives roughly 2.5% peak energy error and RMS position error 0.092–0.104 for these direct-energy methods. This is a tradeoff between a soft target and the resulting trajectory, not exact conservation.

The spring remains problematic. With softness 0.01, the joint resting-spring trajectory travels only 0.0022515 total distance and has zero late kinetic energy, while allowing 1.3061% energy error. Immediate updates travel 15.4388 but have 3.7426% maximum energy error. For the resting-spring case at softness 1, the methods allow about 12–16% energy error. No one ordering or softness is uniformly best across these cases.

![Soft energy ordering](../../results/energy_ordering/analysis/soft_energy_ordering.png)

## Position residual and frozen inner passes

Relinearizing scalar constraints without retaining R does not generally converge the complete position equations. For the tangent spring, h=1/60, beta=0:

| Joint variant | Outer sweeps | Max absolute H / Escale | Max position residual | Full KKT closed steps / 600 |
|---|---:|---:|---:|---:|
| drop R | 5 | 2.9402e-15 | 3.2432 | 0/600 |
| retain R | 5 | 8.8358e-11 | 0.00011943 | 590/600 |
| drop R | 20 | 2.6952e-15 | 3.2432 | 0/600 |
| retain R | 20 | 3.4302e-15 | 3.9968e-13 | 600/600 |

At 20 sweeps, retaining R closes the full equations on all 600 steps, but RMS position error is still 1.50749. Solving the constrained objective accurately does not repair its physical trajectory.

For frozen drop-R tangent-spring blocks at five outer sweeps, increasing inner passes from 1 to 5 to 20 reduces maximum material-row error from 1.11e-9 to roundoff, but the full position residual stays near 3.24. Scalar solve counts rise from 6000 to 30000 to 120000 while derivative evaluations remain 3000. More inner passes approximate the joint local block; they do not restore the omitted position residual, and they are not equal-work comparisons.

The hard pendulum has an additional geometric restriction: when its two scalar constraint gradients are independent, they fully determine the two position increments. Retaining R then changes the dual update, not that position increment. In this implementation the five-sweep hard tangent joint trajectories agree for both R choices. At 20 sweeps both variants abort after only 15 accepted steps (0.25s). Their accepted-prefix energy error is near roundoff; reporting that as ten-second success would be incorrect.

For the soft tangent pendulum, retaining R does matter: at softness 0.01 and 20 sweeps it lowers maximum position residual from 6.06946 to 2.1043e-7 and RMS position error from 0.036991 to 0.020705, while normalized energy error changes from 0.0060314 to 0.0051371.

![Frozen block and position residual controls](../../results/energy_ordering/analysis/residual_controls.png)

## Rank, feasibility, and validation

- `rank_deficient_blocks` measures the complete joint/frozen two-row linear system. A frozen scalar pass can have positive scalar denominators while the complete block is rank deficient. Scalar pass availability therefore does not demonstrate a repaired rank deficiency.
- `singular_updates` has method-specific scope: an unavailable energy scalar update for immediate/reverse/symmetric, or a dependent/unavailable joint block. `incompatible_updates` counts skipped local equations. These counters should not be compared as if they measured identical operations.
- `infeasible_steps` is separate: for a hard pendulum it checks whether the endpoint circle for fixed length intersects the endpoint energy circle at all. There are 82 runs with at least one such step, versus 48 runs with a rank-deficient joint/frozen block. Local rank and global endpoint feasibility are distinct.
- The reference is RK4 with 32 substeps per simulation step, checked against 64. Maximum refinement position difference is 1.4723e-10; maximum reference energy error is 2.1288e-11. This supports the reported trajectory comparisons at their displayed precision.
- The analysis validates unique settings, finite CSV cells, accepted-step counts, sequential trace indices, final times, and trace-versus-summary energy extrema. It checks all 1020 summaries and 95 trace groups containing 56719 states. Forty runs have no late-interval kinetic samples; `late_mean_kinetic_available` marks them unavailable rather than treating stored zero as freezing.
- Shared validation passed all 12 repository CTests, including the new ten-check energy-ordering test executable. The new tests also passed ASan+UBSan with leak detection disabled for the execution environment. The production XPBD solver and four base scene definitions remain unchanged.
- These are ten-second experiments with explicit finite-iteration errors and safety limits. They do not prove long-time stability, and no score here combines energy accuracy, physical trajectory accuracy, and solver convergence into a single success label.

## Reproduce and inspect

From the repository root:

```sh
./build/xpbd-energy-ordering-study results/energy_ordering
MPLCONFIGDIR=/tmp python3 experiments/energy_ordering/analysis.py results/energy_ordering
```

The analysis uses the Python standard library for CSV validation and summaries; plots require matplotlib. `--no-plots` skips that dependency and `--output` overrides the default results-directory `analysis` subdirectory.

- [Raw summary](../../results/energy_ordering/summary.csv) and [accepted-state traces](../../results/energy_ordering/traces.csv).
- [Representative slice](../../results/energy_ordering/analysis/representative.csv), [position-residual controls](../../results/energy_ordering/analysis/controls.csv), and [spring angular momentum](../../results/energy_ordering/analysis/spring_angular_momentum.csv).
- [Machine-readable audit](../../results/energy_ordering/analysis/audit.json).
