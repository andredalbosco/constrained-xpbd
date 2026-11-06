# FEPR implementation contract

This implementation follows Dinev et al.'s **fast projection**, Equations 8, 13, and 14, specialized to one planar particle. It uses energy, both planar linear-momentum components, planar angular momentum, and both auxiliary variables. It is not a scalar energy rescale or a full nonlinear solution of the original closest-point objective.

Primary sources: [author-hosted paper](https://users.cs.utah.edu/~ladislav/dinev18FEPR/dinev18FEPR.pdf), [author project page](https://www.tiantianliu.cn/papers/dinev18FEPR/dinev18FEPR.html), [DOI](https://doi.org/10.1145/3197517.3201277). The locally available [paper text](/home/blacknoir/Projects/dynamics-research/archive/resources/papers/FEPR%20-%20Fast%20Energy%20Projection%20for%20Real-Time%20Simulation%20of%20Deformable%20Objects.md) was read through Sections 3, 3.1, and 3.2, including the approximation and regularization details. The author pages and targeted code searches did not reveal an author-maintained source release; no code-parity claim is made. The existing [research C++ implementation](/home/blacknoir/Projects/dynamics-research/archive/src/solvers/fepr.cpp) is a separate local implementation, not identified as author code.

## Objective and actual fast update

Let `(x_*,v_*)` be the completed base step, held fixed when constructing targets. Define `q=(x,v,s,t)`, `q_*=(x_*,v_*,0,0)` and

\[
D=\operatorname{diag}(mI_2,h^2mI_2,\epsilon,\epsilon).
\]

The paper's original objective is `min ||q-q_*||²_D/2`, subject to

\[
H(x,v)=H_{target},\qquad
P(v)=P_*+s(P_n-P_*),\qquad
L(x,v)=L_*+t(L_n-L_*),
\]

where `P=mv`, `L=m(x_1v_2-x_2v_1)`, and `H=U(x)+m||v||²/2`. The auxiliary scalars are unrestricted; they permit extrapolation as well as interpolation. The spring is attached to the origin, so linear momentum is not an invariant. The formulas control its change using the tentative and old values; the implementation does not incorrectly impose `P=P_n`. Likewise `L=L_n` is not automatically forced when the base integrator changed angular momentum.

The paper replaces fixed-center SQP by a projection centered at **the current iterate**. With row-oriented Jacobian `J=∇cᵀ`, each iteration is

\[
S=JD^{-1}J^T,\quad S\lambda=c(q),\quad
q\leftarrow q-D^{-1}J^T\lambda.
\]

The target momenta remain frozen; only the local distance center moves. This distinction is essential: convergence means the constraints are resolved, not that the original fixed-center objective's KKT equations are solved. There is no accumulated multiplier or retained fixed-objective stationarity correction in this fast update.

The default epsilon is `1e-3*m*l²`, using the rest length as the characteristic size; passing `epsilon` overrides it. This states the size convention explicitly and implements the paper's suggested mass/length scaling. At `m=l=1` it equals the paper's experimental `.001`. Default stopping is the paper's raw `sum(abs(c))<1e-7`, with a finite iteration cap. The determinant threshold is implementation-defined as `1e-12` because the paper does not specify “close to zero.” A small/nonpositive determinant adds `1e-7 I` to the Schur matrix, as in the paper. These raw dimensional thresholds are exposed parameters and are not advertised as unit-invariant. All regularizations are counted; no line search, clipping, fallback integrator, or energy rescale is hidden in this path.

## Derivatives

For `r=||x||>0`, the spring has

\[
U=\frac k2(r-l)^2-mg\cdot x,\quad
\nabla U=k(r-l)x/r-mg.
\]

The energy row gradient is `(∇U,mv,0,0)`. Linear-momentum rows have velocity block `mI`, s-column `-(P_n-P_*)`, and zero other blocks. The angular row has x-gradient `m(v_2,-v_1)`, velocity-gradient `m(-x_2,x_1)`, and t-derivative `-(L_n-L_*)`. Nonzero-rest-length spring gradients at `r=0` are undefined; the implementation records that case rather than inventing a direction.

## Base steps and the hard-rod extension

For the anchored radial spring, the tentative step is the exact global minimum of its implicit-Euler/PD objective. Set `y=x_n+hv_n+h²g`, `A=m/h²`. For nonzero `y`,

\[
x_* = \frac{A\|y\|+kl}{A+k}\frac{y}{\|y\|},\qquad
v_*=(x_*-x_n)/h.
\]

At `y=0`, the base objective has no preferred direction; the old radial direction is the explicit tie rule. FEPR subsequently treats x and v independently.

The paper handles attachments through spring potentials; it does **not** supply a hard rod position/tangency algorithm. `rigid=True` is therefore explicitly labeled **FEPR planar projection plus hard-rod extension**. It uses `U=-mg·x`, a standard RATTLE tentative step, and appends

\[
c_C=\tfrac12(x\cdot x-l^2),\qquad c_T=x\cdot v
\]

with gradients `(x,0,0,0)` and `(v,x,0,0)`. The same fast projection then has six rows and six unknowns. It may be rank deficient, branch sensitive, or fail within the cap. Keeping the two extra rows simultaneous avoids projecting geometry afterward and silently spoiling energy or momentum.

RATTLE first solves `x_*=x_n+hv_n+h²g/2+αx_n` on the rod circle, taking the root connected to `α=0`, then uses `v_half=(x_*-x_n)/h` and projects `v_half+hg/2` onto the tangent at `x_*`. A negative quadratic discriminant or incompatible initial geometry/tangency is an explicit baseline failure. No substepping is silently inserted.

## API and target semantics

`fepr_step(x,v,h,m,k,l,g,rigid=False, *, energy_target=None, ...)` returns a dictionary containing `x`, `v`, `converged`, raw `residual`, named residuals, `counts`, `reason`, base-state diagnostics, auxiliary scalars, and method label. Scalar `g` means downward magnitude; a 2-vector is used directly. `baseline_step` and `project_fepr` are independently callable for controls and tests.

`energy_target=None` uses the incoming step energy, matching the paper's conservative step. An explicit float imposes a fixed run target. The paper's old-state feasibility argument applies to the per-step target; it does not automatically apply to a different fixed target. The projection returns its endpoint velocity directly and never reconstructs it from the final position.

## Validation and observed limits

The implementation self-test passed its finite-difference derivatives (maximum relative error 5.60e-11), full increment-KKT comparison (1.45e-14), and eight fixed-target one-second trajectories (maximum absolute energy error 9.70e-11). The independent [source and numerical audit](fepr_independent_validation/SUMMARY.md) passed 540/540 checks, including two sequential projection updates, unrestricted auxiliaries, and target semantics.

The completed trajectory studies nevertheless show nonmonotone refinement and large position-branch jumps for a stiff spring. See [FEPR_REFINEMENT.md](FEPR_REFINEMENT.md) for quantitative evidence. Energy tolerance alone is not evidence of correct motion; the comparison also reports geometry/tangency, projection failures, phase and reference error, base-step costs, and the original energy target throughout a trajectory. No algorithmic changes were made in response to these findings.
