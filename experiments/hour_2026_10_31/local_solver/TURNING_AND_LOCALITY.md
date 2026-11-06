# Two distinct numerical failures

The scalar turning stall and the slow local-chain solve have different causes. Both occur while trying to solve a well-defined constrained objective, and neither is explained by a violated energy identity alone.

## Why the scalar initialization fails

For the harmonic scalar, write `q=x-l` and `omega²=k/m`. The two initial displacements are

\[
\delta x_{BE}=\frac{h v_n-h^2\omega^2q_n}{1+h^2\omega^2},
\qquad
\delta x_{mid}=\frac{h v_n-\tfrac12h^2\omega^2q_n}{1+\tfrac14h^2\omega^2}.
\]

Both initializers then reconstruct `v=2 delta x/h-v_n`, because this experiment uses soft midpoint kinematics. The BE position objective and that velocity reconstruction disagree most visibly near a turn.

One replayed failing state has `q_n=-0.24724823`, `v_n=-0.36990755`, `omega=10`, and `h=1/60`. The BE initialization gives positive displacement `+0.00068388` and velocity `+0.45197359`; the midpoint initialization gives negative displacement `-0.00271229` and velocity `+0.04443307`. The particle turns during the step, but its net displacement still points toward the turning point. The BE initializer starts on the wrong side of that net motion and adds energy0.01684 before the constrained solve.

For `q_n<0, v_n<0`, opposite displacement signs occur in the interval

\[
\tfrac12h\omega^2|q_n|<|v_n|<h\omega^2|q_n|.
\]

At the replayed state this is `0.20604 < 0.36991 < 0.41208`. This gives a concrete trigger rather than attributing every failure vaguely to stiffness.

The midpoint predictor with its own predicted DG target satisfies the scalar energy balances, midpoint velocity relation, and zero work-correction penalty. Its quarter-objective gradient is zero. It is already the relevant KKT solution. BE-initialized KKT instead approaches a nearly singular Newton system away from feasibility and stops in its line search. Increasing its cap128→4096 does not change the result. The balance Jacobian remains full row rank; the failure is not an unavoidable rank loss of the energy constraints.

The hour study tests that initializer across the scalar, planar spring, fixed rigid pendulum, elastic pendulum, and two-mass elastic chain. It removes every step failure in the tested nominal/k×10 one-second runs. This is evidence for a compatible initialization, not proof that arbitrary nonlinear states will converge.

## Why local AL needs many visits

At a converged solution `(z*,y*)`, let `J` be the balance Jacobian, `H` the Lagrangian Hessian, and `A=H+rho J^T J` the augmented-Lagrangian Hessian. An exact local block Newton update has the first-order error map

\[
\delta z\leftarrow\delta z-E_B A_{BB}^{-1}E_B^T(A\delta z+J^T\delta y).
\]

Composing three passes produces matrices `P,Q` with `delta z_new=P delta z+Q delta y`. The shared dual update gives

\[
\begin{bmatrix}\delta z_{new}\\\delta y_{new}\end{bmatrix}
=\begin{bmatrix}P&Q\\\rho JP&I+\rho JQ\end{bmatrix}
\begin{bmatrix}\delta z\\\delta y\end{bmatrix}.
\]

`spectrum.cpp` builds this iteration matrix from the exact same local blocks and Hessians used by the solver. All sampled block Hessians need no eigenvalue shift. Finite differences of the nonlinear full-step cycle agree with the derived matrix at relative error about1e-7 for N2/8/32. This checks the local linearization; it does not assert that an arbitrary nonlinear trajectory accepts full steps or excites the slowest mode.

For the nominal N32 chain, the original split blocks with rho10 have spectral radius0.9620 per outer cycle. Six decades of worst-mode reduction would need about357 cycles, or68,515 local block visits. Edge-full blocks with rho100 lower the radius to0.7656 and the corresponding estimate to4,966 visits. Raising rho to1000 worsens this local block iteration. A single global primal block with rho1000 has radius0.2078 and an estimate of26 global visits; each of those visits solves a much larger system.

The iteration spectrum shows why512 evaluations, warm duals, or a larger penalty alone cannot ensure a resolved chain solve. It also distinguishes the original split's position/velocity coupling from the remaining shared graph modes. It provides no claim that local blocks are faster than the dense joint solve on these small systems.

## An additional roundoff plateau

Long AL runs can stop making useful progress even after the slow modes have mostly decayed. The Armijo test compares total objective values that include a large negative gravitational constant. Near convergence, the true local decrease is smaller than the precision of that subtraction, leading to repeated rejected steps or ineffective tiny steps.

`incremental.hpp` evaluates the exact change of the same objective and balance rows. For squared terms it uses `(a+b)(a-b)`; for distance changes it uses

\[
\|e+\delta e\|-\|e\|=
\frac{\delta e\cdot(2e+\delta e)}{\|e+\delta e\|+\|e\|}.
\]

The augmented difference is then

\[
\Delta\mathcal L_\rho=\Delta F+
\sum_r(y_r+\rho g_r)\Delta g_r+
\frac\rho2\sum_r(\Delta g_r)^2.
\]

The final objective, constraints, gradients, Hessians, and stopping tolerances are unchanged. The difference formula passes120 checks across all five small environments and an N8 mass-contrast chain. Baseline and stable-difference results are stored in separate CSVs. This addresses arithmetic cancellation in globalization; it does not remove the slow iteration modes or turn a finite local solve into an exact one.
