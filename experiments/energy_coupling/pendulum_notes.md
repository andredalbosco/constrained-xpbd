# Direct material/energy multiplier coupling: hard 2D pendulum

Independent C++17 experiment; the project repository is untouched. Run in this directory:

```sh
g++ -std=c++17 -O2 -Wall -Wextra -pedantic study.cpp -o study
./study > run.log
```

Outputs: `summary.csv` (640 cases), `one_step_traces.csv` (initial-step traces plus failed-step traces for the 80-iteration runs), and `run.log`. Both CSVs identify the energy-target policy explicitly. Experimental singularities and unstable rows are recorded as data; they do not make the program fail its independent validation checks.

## What was solved

Mass m=1, radius L=1, anchor at the origin. Gravity is (0,-9.81), except for the free-rotation case. Initial states are theta=0.65 rad from downward vertical with speed 0 or 1.2, bottom equilibrium with speed 0, and gravity-free rotation at speed 2. Sweeps use h in {1/240,1/60,1/15,1/5}, iteration limits {1,5,20,80}, and up to 10 s.

The primary energy target is the **previous actual mechanical energy**, recomputed from the accepted state each timestep. Consequently small per-step energy errors may accumulate, and accepting an incorrect step changes the next target. A complete control sweep instead holds the **initial mechanical energy** fixed for the whole run. The results below use the previous-actual policy unless explicitly labeled otherwise.

For the fixed-initial control, H_fixed=H_previous+(E_old-E_initial); the gradient G is unchanged. Both policy names appear in the CSV column `energy_target`.

    xhat = x_old + h v_old + h² g,
    v_new = (x - x_old)/h,
    C(x) = |x| - 1,
    H(x) = |v_new|²/2 - g·x - E_old,
    J = x/|x|,
    G = (x-x_old)/h² - g,
    rx = x-xhat - J lambda - eta G.

Both joint and sequential methods drop rx from the update and drop curvature. The full rx is still measured using accumulated multipliers and accepted gradients. Small C and H do not establish that this mixed residual has been solved.

The joint solve uses the true tangent Schur complement: with t perpendicular to unit J, s=t·G,

    Schur = s²,
    delta_eta = (-H + (J·G) C)/s²,
    delta_lambda = -C - (J·G) delta_eta.

Its position correction is evaluated as

    delta_x = -C J + t(-H + (J·G) C)/s,

which avoids cancellation of enormous parallel components. No diagonal regularization, compliance, pseudoinverse, or near-parallel cutoff is added. Exactly zero Schur is reported as singular. Small nonzero Schur is allowed to expose the actual large update.

Sequential projection first applies the material row, relinearizes, and then applies delta_x=-H G/|G|². It updates both accumulated multipliers. It is ordinary alternating nonlinear row projection, not an exact joint solve.

Variants:

- `xpbd`: validated ordinary one-constraint XPBD baseline.
- `joint_raw`: starts the joint solve at the external predictor.
- `joint_seed`: explicitly performs one material projection before starting the joint solve. This is a separate initialization experiment.
- `sequential`: material then energy, relinearized at each row.
- `joint_seed_backtrack`: same seeded joint solve with an explicitly separate, optional 30-halving merit backtrack. This cannot fix a missing linearization direction; it rejects the unusable step.

## Mandatory checks performed before the sweep

The executable fails if any of these checks fail:

- The ordinary XPBD one-step projection and velocity reconstruction agree with independent formulas; gravity-free one-step kinetic energy matches its analytic value.
- Every baseline scene remains finite and exactly unit-length to 1e-13 over 10 s at every tested timestep.
- At a resting predictor, G=0 and H=-h²|g|²/2<0.
- After material projection of a resting predictor, G is parallel to J within rounding error and H remains negative for the 0.65 rad case.
- The raw resting joint solve reports its singularity, and the first resting sequential solve remains unconverged after 80 iterations.
- A shifted fixed energy target shifts H by the required constant without changing its gradient, and a negative energy-circle radius squared is correctly classified as infeasible.

To avoid manufacturing a spurious tiny derivative at rest through subtraction, the code evaluates the algebraically identical complete-square form:

    p = x_old+h²g,
    H = (|x-p|²/h² - |v_old|² - h²|g|²)/2,
    G = (x-p)/h².

For fixed-initial targeting the constant E_old-E_initial is added to H. The code independently reports the original physical-energy expression and their floating-point difference. This exposes the exact zero derivative instead of hiding it behind cancellation noise. Final multipliers and rx are also checked explicitly for finiteness before a candidate can be accepted.

## Failure mechanism 1: stationary predictor and parallel gradients

At the predictor,

    G(xhat)=v_old/h,
    H(xhat)=-h²|g|²/2.

Thus G is exactly zero only when v_old=0. The joint linearization has a nonzero energy residual and a zero energy row. Even bottom equilibrium fails from this raw initialization, although a material-first projection finds its exact stationary solution immediately.

For the resting 0.65 rad pendulum, material projection produces x=unit(xhat). Since the energy center is xhat, G=(x-xhat)/h² is parallel to J. The projected joint Schur is zero in exact arithmetic. Rounding produces tiny nonzero perpendicular components: the unregularized seeded joint solve then makes enormous corrections before reaching an exactly singular floating-point state. At h=1/60 its largest trial displacement is 1.94e37 m, with no accepted timestep.

Sequential projection avoids the joint division by the tiny tangent component. Nevertheless it does not solve the resting first step: material projection restores radius 1, and the following energy projection moves along the same radial ray and violates the length again. Repeating the sweep cannot resolve this exact invariant-ray obstruction. Carrying this unconverged state to the next timestep changes the geometry and target, which is why a long sequential simulation can continue; that is not first-step convergence.

## Failure mechanism 2: a weak energy gradient before angular ill-conditioning

Moving case theta=0.65, initial speed 1.2, h=1/60, 80 iterations:

- `joint_raw` accepts 15 steps, with max physical energy error 4.19e-11 J and max length error 1.13e-14 m.
- At the next predictor, the speed is nearly zero: |G|=2.84445e-4. The row-normalized Gram condition number is only 1.00377; the gradient directions are healthy.
- H is still -0.013366125 J, so the first joint Newton correction is **46.9902 m**. The neglected quadratic curvature of H is decisive.
- Following that correction, the gradients become nearly parallel; the normalized condition grows to about 1.22e33, and a trial correction reaches 1.31e78 m before a rank-one failure.
- All 16 input states were geometrically feasible according to the independent circle-intersection test. This failure is not attributed to an impossible energy target.

This separates a small energy-gradient magnitude from a small angle between rows. Normalizing the rows can conceal the former; it does not make the finite Newton step reliable. Near v_old=0, the required linearized energy correction scales like h³|g|²/(2|v_old|).

The backtracked seeded variant rejects large unusable updates instead of accepting them, but it still cannot finish the run: for this case it accepts 17 timesteps before rejecting a step. It is recorded separately rather than presented as a cure.

## Fixed-initial-energy control

For the same moving case (theta=0.65, speed 1.2, h=1/60, 80 iterations), **raw joint coupling still accepts exactly 15 steps and then fails** when the target is held at the initial energy. All input states remain circle-feasible. Its last runaway trial is about 9.22e61 m, versus 1.31e78 m with previous-actual targeting. Those enormous terminal magnitudes depend on floating-point details and should not be compared quantitatively; the reproducible result is the same failure timing and mechanism.

| Target policy | Method | Accepted steps | Max accepted energy error | Max accepted length error | Unconverged steps | Max rx |
|---|---|---:|---:|---:|---:|---:|
| Previous actual | raw joint | 15 | 4.19e-11 J | 1.13e-14 m | 0 | 33.9865 |
| Fixed initial | raw joint | 15 | 4.14e-11 J | 1.15e-14 m | 0 | 33.9865 |
| Previous actual | sequential | 600 | 8.79e-9 J | 8.61e-13 m | 0 | 0.0419079 |
| Fixed initial | sequential | 600 | 8.49e-11 J | 1.48e-8 m | 9 | 0.0521065 |

Fixed targeting prevents accumulation of energy error in this sequential run, but it does not close the mixed residual, and nine steps fail the stated constraint tolerances. The backtracked seeded joint control lasts 82 steps before rejecting a step with the fixed target, so target policy can affect its trajectory even though it does not eliminate failure. Resting raw-predictor singularities are identical under both policies because the initial target energies coincide.

## Joint blow-up versus bounded sequential behavior

Gravity, h=1/60, 10 s target. Energy errors are against the initial actual energy; constraints are accepted-state physical errors. `rx` is the full omitted mixed residual.

| Initial state | Method/iterations | Accepted steps | Max energy error | Max length error | Max rx | Result |
|---|---|---:|---:|---:|---:|---|
| Rest | raw joint/80 | 0 | n/a | n/a | n/a | Zero-gradient singularity |
| Rest | seeded joint/80 | 0 | n/a | n/a | n/a | Runaway trial, then singular |
| Rest | sequential/20 | 600 | 0.00293865 J | 6.25091e-4 m | 0.00590063 | Bounded; 39 unconverged steps |
| Rest | sequential/80 | 600 | 0.000705892 J | 6.25091e-4 m | 0.0317842 | Bounded; 18 unconverged steps |
| Speed 1.2 | raw joint/80 | 15 | 4.19e-11 J | 1.13e-14 m | 33.9865 | Runaway trial on next step |
| Speed 1.2 | sequential/20 | 600 | 0.00390508 J | 8.91599e-4 m | 0.0118245 | Bounded; 43 unconverged steps |
| Speed 1.2 | sequential/80 | 600 | 8.79e-9 J | 8.61e-13 m | 0.0419079 | C/H converge; rx remains unclosed |

Therefore this experiment supports a concrete moving example where joint coupling fails dramatically while sequential projection remains bounded. It does **not** support saying sequential projection universally solves the intended mixed equations or is physically accurate.

## Feasibility and branch ambiguity

The two constraints are circle intersections:

    |x|=1,
    |x-p|=R,
    p=x_old+h²g,
    R²=h²|v_old|²+h⁴|g|²+2h²(E_target-E_old).

Negative R² is infeasible. Otherwise feasibility requires ||p|-1| <= R <= |p|+1. This is checked at every attempted timestep and logged as `infeasible_inputs`. For previous-actual targeting the last term vanishes, and in the free case the condition reduces to h|v_old|<=2. Low-iteration runaway states can make later targets infeasible; such cases are distinguished from the initially feasible moving failure above.

For a resting initial state, x=x_old always satisfies both constraints. The second intersection is the reflection of x_old across the direction p. At small h its displacement is approximately 2h² g_tangent, whereas the physical initially resting pendulum moves approximately h² g_tangent/2. Exact energy and length with the backward velocity reconstruction therefore do not select the accurate physical branch.

## Free rotation: constraints can converge without tangency or angular conservation

For g=0 and initial tangent speed 2, all three primary approaches can converge C/H at sufficiently many iterations. At h=1/60, raw joint keeps energy within 2.86e-11 J and length within 2.22e-16 m over 10 s, but endpoint radial velocity is 0.0333333 m/s and anchor angular momentum changes by 0.000277797.

This agrees with exact geometry: fixing speed via energy and using v_new=Delta x/h fixes chord length h|v|, not endpoint tangency. The first-step angular momentum is L_new=|v|sqrt(1-h²|v|²/4), which differs from the initial tangent value. Subsequent values can remain constant on the selected branch. Gravity-case angular momentum is not a conserved quantity and is not used as a pass/fail condition.

## Output interpretation

`accepted_steps` counts finite candidates actually carried to the next timestep. A maximum-iteration exit is accepted but separately counted in `unconverged_steps` if either C or H fails tolerance. Failed singular or rejected candidates do not contaminate the accepted-state maxima; their errors are in `max_attempted_*` and `last_candidate_*`. Thus zero accepted-state error with zero accepted steps is **not success**.

The classification separates exact singularity, runaway trial, accepted runaway, bounded nonconvergence, and bounded constraint convergence with rx still unclosed. The diagnostic runaway thresholds (radius above 10 m or speed above 100 m/s, or a trial update/radius above 10 m) only label output; they never clamp or alter the solver. Baseline XPBD is classified against the added H constraint too, so its `bounded_constraints_nonconverged` label means it does not enforce that extra energy constraint, not that its distance projection failed.

The row-normalized Gram condition number is (1+|cos(angle)|)²/sin²(angle). Undefined zero-gradient rows are assigned infinity. Near-parallel samples have |sin(angle)|<1e-6; this is a reporting threshold only. Minimum energy-gradient norm is reported separately.
