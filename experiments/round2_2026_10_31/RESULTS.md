# What the second research round changed

The strongest new result is a formulation error, not an iteration-order effect: the original independent-velocity objective gives the wrong acceleration in a generic moving state. Two repaired position objectives remove that error in the tested cases. Separately, the rigid Cartesian midpoint method has a hidden radial velocity mode that exact energy cannot detect. These two findings explain why additional solver iterations or tighter energy tolerances can fail to improve motion.

The experiments remain isolated from the project solver. [README](README.md) describes reproduction; [the recorded-motion app](viewer.html) compares334 method entries across16 scenarios. Native C++ executables perform the simulations; the HTML/Canvas app only replays recorded frames. The project itself is headless, with no existing native graphical interface. The viewer has passed data/schema and JavaScript syntax checks; browser playback verification is pending because the available in-app browser could not connect and automatic approval review rejected launching the fallback Chrome session.

## Evidence added

| Branch | New measurements | Boundary |
|---|---|---|
| Objective revision |22 formulations;350 principal trajectories,60 near-rest trajectories;660 local consistency cases;200 dense scalar searches;600 angle/speed probes;60 rod probes |13 long large-chain references fail refinement;9 principal runs contain unconverged steps |
| Adversarial mechanics |1,110 trajectories, including10/60-second runs,20 initial conditions, up to h=0.25 |61 failed prefixes;42 rows have uncertified spring/chain references; finite completion is separate from convergence |
| Numerical scaling |720 finite-budget configurations, chains up to32 masses;420 isolated timing repeats |8/36 dense reference scenarios stop early; errors compare the same discrete mechanics |
| Start-only work and local solves |4,320 work-estimate rows;1,824 one-step solver rows;320 frozen-target,88 local-AL,24 crossed-objective trajectories |Prediction error, KKT convergence and physical motion are reported independently |

All four optimized CTest suites pass. The cross-objective study additionally checks six gradient/Hessian/unchanged-constraint combinations, and scalar stalled states are replayed independently. Counts represent configured experiments and comparisons, including controls; they are not counts of distinct physical formulations.

## 1. Repair the mechanical target before tuning compliance

For a moving state, the old midpoint-target objective predicts

\[
a_{\mathrm{pred}}=2a_n-
\frac{v_n^TMa_n}{v_n^TMv_n}v_n,
\qquad a_n=-M^{-1}\nabla U(x_n).
\]

The force component perpendicular to velocity is doubled. Shrinking beta does not remove it. This limit assumes nonzero velocity; the turning-point limit is analyzed separately.

Two repaired position objectives are

\[
\Phi_Q(x)=\frac{1}{2h^2}
\|x-x_n-hv_n-\tfrac14h^2a_n\|_M^2+\tfrac14U(x),
\]

and

\[
\Phi_M(x)=\frac{1}{2h^2}\|x-x_n-hv_n\|_M^2+
U\!\left(\frac{x+x_n}{2}\right).
\]

They retain the independent-velocity midpoint penalty and hard endpoint energy. For the quadratic scalar spring their mechanical state is already a feasible unconstrained minimum; energy needs no correction. Moving-state probes also recover the correct local acceleration and approximately third-order local velocity error.

At120 Hz over10 seconds, position RMS in metres:

| Case | Original midpoint with beta(h) | Quarter objective | Midpoint potential | DG residual + energy |
|---|---:|---:|---:|---:|
| Scalar release |0.1696|0.00594|0.00594|0.00594|
| Orbiting spring |1.201|0.00221|0.00218|0.00222|
| Rigid pendulum release |0.0368|0.000518|0.000827|0.000518|
| Elastic pendulum |0.5646|0.00160|0.00156|0.00160|
| Two-spring chain |0.5515|0.00424|0.00385|0.00421|

These are complete trajectories with certified references. The old formulation also has tiny energy error, so conservation alone would have hidden the problem. The detailed derivation and angular-momentum comparison are in [objectives/README.md](objectives/README.md).

## 2. Endpoint tangency controls a hidden state mode

With exact old/new rod length and exact Cartesian midpoint kinematics, an inherited normal velocity reverses sign as h tends to zero. Its speed and kinetic energy can remain unchanged. One full step reflects once, two half steps reflect twice, so the usual adaptive error estimator can stop decreasing.

The probe uses an actual previous DG output. After a0.01-second release step, normal velocity is8.8082e-6; at the next h≈2.03e-6 it becomes−8.8098e-6. The scaled step-doubling discrepancy is1.2459e-5 despite a3.03e-11 mechanical residual. Endpoint tangency is therefore an independent state condition, not something to infer from energy or rod position.

The alternatives tested here are RATTLE; tangent and energy postcorrections; joint symmetric projection; and a single-pendulum angular DG control. Symmetric projection follows the same-multiplier before/after construction of [Hairer](https://www.unige.ch/~hairer/preprints/symproj.pdf); [RATTLE](https://doi.org/10.1016/0021-9991(83)90014-1) separately controls position and velocity constraints.

Near the separatrix at30 Hz, RATTLE changes the below-threshold reference swing into numerical rotation. Symmetric energy projection and angular DG preserve the qualitative motion in that test. At extreme h, exact conservation still permits large phase error. The angular control is specific to a single pendulum; its generalization is a subsequent experiment.

See [adversarial/RESULTS.md](adversarial/RESULTS.md) and the figure below.

![Adversarial controls](adversarial/figures/new_failures_and_controls.png)

## 3. Start-only work and the repaired objective can compose

Eight frozen work targets were tested, including full Taylor power, predicted-midpoint work, predicted discrete work, and stiffness-aware linear predictions. All are available from old state data. A target computed from a predicted endpoint remains an estimate for the optimized endpoint.

The stiffness-aware midpoint predictor solves

\[
(M+\tfrac14h^2B_n)\,\delta\hat x
=hMv_n+\tfrac12h^2f_n,
\quad B_n=\sum_j k_j\nabla C_j\nabla C_j^T,
\]

then evaluates DG work along that prediction. It greatly improves nominal work estimates and avoids the explicit predictor's tested stiff-step failure. The estimate requires a linear solve; its cost is not zero.

The direct composition test is more informative than predictor error alone. Fully converged one-second joint trajectories at60 Hz give:

| Case | Original objective + frozen power | Quarter objective + frozen power | Quarter objective + predicted midpoint DG |
|---|---:|---:|---:|
| Orbiting spring RMS(m) |0.08586|0.003991|0.0008308|
| Two-spring chain RMS(m) |0.2163|0.008144|0.001568|

The scalar cross-test initially stalls near turning points. Raising its budget from128 to4,096 evaluations produces the same stalls. Reinitializing the same constrained objective from the start-only midpoint predictor solves all four recorded failures in1–3 evaluations, with positive reduced curvature and agreement with a dense feasible-circle search. This is a solver initialization problem in these particular states, not a necessary rank loss of the constraints. [Detailed cross-test](local_work/CROSS_OBJECTIVE.md).

## 4. Shared kinetic multipliers can be updated locally, with costs and branch risks

The explicit-work AL prototype has one kinetic multiplier per vertex and one storage multiplier per spring. Spring blocks update both endpoint positions and their work values; vertex blocks update velocity and incident work. Shared multipliers receive a separate dual update. This is a concrete nonlinear block AL method, not a claim of convex ADMM convergence.

It converges in many one-step settings, but requires far more evaluations than the joint solve on these small systems. Some block orders converge to a different strict local minimum; small residuals do not identify the desired branch. Short rollouts retain inexact-step counts. [Local formulation and results](local_work/RESULTS.md).

## 5. Paper-inspired numerical transfers have mixed outcomes

[Large Steps](https://bnicolet.com/publications/Nicolet2021Large.html) motivates changing numerical coordinates while keeping the physical objective. [Muon](https://kellerjordan.github.io/posts/muon/) motivates spectral inspection; the tested bounded row balancing is our numerical analogue, not the optimizer itself. [TRPO](https://proceedings.mlr.press/v37/schulman15.pdf) motivates a metric-aware locality test; capped Newton does not inherit a mechanical guarantee.

On the32-mass, mass-ratio100, h=0.04 case under budget48:

| Metric | Median CPU for10 steps | Max modal solve error |
|---|---:|---:|
| Identity |0.576 ms|2.53e-4|
| Diagonal |0.573 ms|1.97e-5|
| Graph |1.23 ms|9.94e-6|
| Dense SVD rows |83.52 ms|2.73e-9|

This is setup-inclusive serial timing in a dense prototype. Diagonal scaling improves this finite-budget case cheaply; dense spectral balancing is a useful diagnostic but a poor runtime choice here. The trust-radius variant has3 paired cases with more than10× energy improvement and13 with more than10× regression. It is not ready to replace line search as a default. [Full results and limitations](preconditioning/RESULTS.md).

## What to pursue next

The promising ingredients are the repaired objective, start-only stiffness-aware work, endpoint tangency, and low-cost numerical scaling. The next tests should address their shared failure modes: initialization and nonlinear branches at turns, large-chain shared multipliers, energy-graph normalization, exact elimination of work variables, and multi-link geometry-aware integration. The one-hour continuation pursues those directions separately; none is silently counted as validated by this round.
