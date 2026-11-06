# Direct energy conservation: joint versus sequential multiplier updates

**The joint block can be much less robust, but coupling alone is not universally unstable.** The tests reproduce a moving pendulum whose undamped joint solve fails while sequential updates complete the run. A moving spring can also grow without bound at a small iteration budget. Controls distinguish energy-target drift, zero energy gradients, dependent rows, nonlinear overshoot, and convergence of only a subset of the equations.

These experiments use direct total-energy constraints, ordinary XPBD backward velocity reconstruction, linearization, curvature dropping, and omission of the position stationarity residual. They do not use the midpoint work-balance scheme from the other experiment directory. Ordinary XPBD baselines are checked before the experimental sweeps.

## Exact formulation being compared

With timestep endpoints \((x_-,v_-)\), constant gravity \(a\), and

\[
y=x_-+hv_-+h^2a,\qquad v_+(x)=(x-x_-)/h,
\]

use the implicit objective and direct energy constraint

\[
F(x)=\tfrac12\|x-y\|_{M/h^2}^2+U(x),\qquad
H(x)=\tfrac12\|v_+(x)\|_M^2+U(x)-(Ma)^Tx-E_\star=0.
\]

For the spring, \(C(x)=\|x\|-1\), \(U=C^2/(2\alpha)\). The pendulum has the hard constraint \(C=0\), \(\alpha=0\), and no elastic potential. A hard reaction multiplier is never substituted into elastic energy.

Dualize the objective's elastic energy only. Write \(\eta=h^2\nu\), \(g=\nabla H\), \(J=\nabla C^T\), and \(\widetilde\alpha=\alpha/h^2\). The mixed equations are

\[
r_x=M(x-y)-J^T\lambda-g\eta,\quad
r_C=C+\widetilde\alpha\lambda,\quad r_\nu=H.
\]

Dropping curvature and \(r_x\), the **joint** update is

\[
\begin{bmatrix}
JM^{-1}J^T+\widetilde\alpha & JM^{-1}g\\
g^TM^{-1}J^T & g^TM^{-1}g
\end{bmatrix}
\begin{bmatrix}\Delta\lambda\\\Delta\eta\end{bmatrix}
=-\begin{bmatrix}r_C\\H\end{bmatrix},\qquad
\Delta x=M^{-1}(J^T\Delta\lambda+g\Delta\eta).
\]

The **sequential** version updates the material row first, re-evaluates \(H\) and \(g\), then updates the scalar energy row. It relinearizes between rows; it is not merely another exact solver for the same frozen block. Iteration counts compare passes, not equal numbers of residual evaluations. The energy gradient always comes from analytical position energy, including away from the material constitutive relation.

Retaining the position residual alone does not remove the hard 2D pendulum failure: when its two constraint gradients are independent, the two equations for the position increment already determine that increment uniquely. The top stationarity equation changes multiplier increments only. At a rank-deficient point with an inconsistent right-hand side, it cannot make the linearized constraints feasible. For the spring one-pass comparison, both multipliers start at zero at the predictor, so the initial position residual is zero and retaining it would give the same first update. These observations do not establish equivalence for compliant multi-pass solves.

Both targets are tested: fixed initial energy and energy of the previous accepted state. They agree for exact conservation, but differ when finite iterations leave an energy error. Reusing the previous achieved energy can accumulate that error.

## Measured comparisons

All examples below run for 10 simulated seconds with \(h=1/60\), mass 1, and rest length 1. The spring starts at \(x=(1.25,0)\), \(v=(0,1)\), with \(k=1000\), so \(E_0=31.75\). The gravity pendulum starts at angle 0.65 radians with tangent speed 1.2 and gravity magnitude 9.81.

| Test | Joint update | Sequential update |
|---|---|---|
| Spring, 1 pass, previous-actual energy target | \(E(10)/E_0\approx1.624\times10^{19}\); radius reaches \(8.09\times10^8\) | \(E(10)/E_0\approx2.810\); remains bounded, energy row unresolved |
| Same spring, 1 pass, fixed initial target | Max relative energy error 0.1079; bounded | Max relative energy error 0.7653; bounded |
| Same spring, 5 passes, previous-actual target | Max relative energy error \(4.44\times10^{-14}\); rows converge | Max relative energy error \(5.41\times10^{-8}\); some material rows unresolved |
| Pendulum, 80 passes, previous-actual target | 15 accepted steps; next nonlinear solve runs away and fails | Completes 600 steps; max energy drift \(8.79\times10^{-9}\) J |
| Pendulum, 80 passes, fixed initial target | Also fails after 15 accepted steps | Completes 600 steps; max energy error \(8.49\times10^{-11}\) J, but 9 steps miss the row tolerance |

The spring's one-pass joint block has normalized condition number at most 2.66. Its runaway is therefore not evidence of an ill-conditioned linear block. The linearized energy constraint leaves nonlinear error; accepting that error and resetting the target compounds it. The implemented joint backtracking also accepts those updates because its combined residual merit decreases. Five passes remove the energy runaway in this example, but the mixed stationarity residual remains about 0.0261: energy/material row convergence is not optimization convergence.

In the moving pendulum's first failed previous-actual step, the energy gradient is only \(2.84\times10^{-4}\), although its angle-based normalized condition is about 1.004. Correcting a roughly \(-0.01337\) J energy residual with that tiny gradient proposes a 46.99 m displacement. Later trial iterates become nearly parallel and reach enormous magnitudes. All inputs up to this failure satisfy the exact circle-intersection feasibility test. The terminal magnitude is sensitive to roundoff; the useful finding is the failed nonlinear solve. Sequential updates avoid that failure here, but their stationarity residual reaches 0.0419. Fixed-target sequential updates retain an unresolved stationarity residual as well.

## Why rest and turning points matter

For a spring released from rest, \(x=y=x_-\) already has \(H=0\). The joint update can solve both rows with

\[
\Delta x=0,\qquad \lambda=-h^2 k C(x_-),\qquad \eta=h^2.
\]

The stretched spring freezes exactly. Its energy error is zero and its block is nonsingular. Conversely, the sequential material step reaches the backward-Euler energy minimum, where \(g=0\) while \(H<0\); the energy projection cannot repair this at first order. Unguarded floating-point division by a tiny residual gradient can create enormous corrections.

For the resting pendulum, \(y=x_-+h^2a\) is exactly the center of the energy circle: \(g(y)=0\), \(H(y)=-\tfrac12mh^2\|a\|^2<0\). The joint block is singular before its first correction. Projecting to the rod first makes the two gradients parallel, which is still singular. Sequential updates can remain trapped on the prediction ray with a nonzero row residual.

There is a deeper problem with the direct-energy-constrained backward-Euler objective. With gravity absorbed in the predictor,

\[
F(x)=H(x)-\frac1h v_-^TM(x-x_-)+\text{constant}.
\]

At \(v_-=0\), every energy-feasible position has the same objective value. The objective provides no preference between the stationary branch and a moving branch. This issue survives retaining \(r_x\); it is not caused solely by joint block inversion.

## What the block itself proves

Let \(A=JM^{-1}J^T\), \(b=JM^{-1}g\), \(d=g^TM^{-1}g\). Then

\[
\det S=(Ad-b^2)+\widetilde\alpha d\ge0.
\]

For positive compliance and nonzero \(g\), this block is positive definite. For a hard constraint, parallel gradients make it singular. A zero energy gradient also makes it singular. Schur complements in the implementations are evaluated from nonnegative orthogonal components to avoid interpreting determinant cancellation as a physical failure.

The frozen-block control shows that joint inversion solves a nearly dependent pair by taking a large correction, while truncated scalar GS can appear quiet simply because it leaves a large residual. For gradient angle \(10^{-4}\) and zero compliance, the joint displacement is about 3000; ten frozen GS passes give about 0.2 but leave a residual about 0.3. Coupling is not a generic linear instability of an otherwise well-posed positive-definite system.

## A different meaning of multiplier coupling

If one dualizes the combined \((1-\nu)U\) instead of only the objective's \(U\), the dual quadratic contains \(1/(1-\nu)\). The transformed multiplier is \(\kappa=(1-\nu)\lambda\), and its usual concave dual representation requires \(\nu<1\). At a rest-state degeneracy, \(\nu=1\) is exactly the stationary energy multiplier. That parameterization therefore becomes singular at a state already encountered in these tests. The objective-only dual used here avoids that division; it does not remove the underlying branch degeneracy. See [analytical_notes.md](analytical_notes.md) for the derivation.

## Reproduce and interpret the checks

From the repository root:

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j4
ctest --test-dir build --output-on-failure -j4
python3 experiments/energy_coupling/run_studies.py
```

The new targets are `xpbd-energy-block-study`, `xpbd-energy-spring-study`, and `xpbd-energy-pendulum-study`. Disable them with `XPBD_BUILD_ENERGY_STUDIES=OFF` if needed. Results go under `results/energy_coupling/{block,spring,pendulum}/`.

The studies check baseline dynamics and analytical failure identities, then record experimental successes and failures. A passing study means those checks and output generation succeeded; it does not mean every integration method was stable. The spring sweep includes 864 cases, both energy-target policies, 1/5/20/80 passes, two stiffnesses, and three initial motions. The pendulum compares raw joint, material-seeded joint, sequential, backtracking, and baseline updates with both energy targets. Diagnostic infinities and failed-trial values remain visible.

Source and detailed audit: [block_study.cpp](block_study.cpp), [spring_study.cpp](spring_study.cpp), [pendulum_study.cpp](pendulum_study.cpp), [spring_notes.md](spring_notes.md), [pendulum_notes.md](pendulum_notes.md), [analytical_notes.md](analytical_notes.md).

Validation of the integrated sources: all 10 CTest targets passed, including the existing baseline, reservoir, and work-balance checks. The new runner completed the 20 frozen-block, 864 spring, and 640 pendulum cases. Independent CSV checks verified row counts and schemas; successful states have finite physical diagnostics. ASan/UBSan runs of all three new programs reported no findings, with leak checking disabled because of the execution environment. Expected solver failures and singular-condition infinities remain in the data.

The baseline approximations follow [XPBD equations 12–18](https://mmacklin.com/xpbd.pdf). The additional direct-energy constraints, comparisons, and failure analysis above are the experiments in this directory.
