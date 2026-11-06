# Large Steps transfer: precondition the original frozen dual residual

For a small spring-pendulum chain, use its material Schur block as the first graph preconditioner. It is tridiagonal and can be solved in linear time. Keep the energy row and the physical residuals unchanged. This is a solver experiment, not an additional smoothing energy or a postprocessing pass over positions. The proposed scope is a few connected pendulum masses; no cloth, grid, or large-system simulation is needed.

## What the primary sources actually do

Nicolet, Jacobson, and Jakob replace the optimization metric using a shifted mesh Laplacian. Their differential coordinates satisfy `u = (I + λL)x`; differentiating through the inverse map gives two inverse applications for ordinary gradient descent in `u`, hence a position step proportional to `(I + λL)^−2 ∇Φ`. The objective itself remains unchanged. See equations 14 and 17–20 of the [paper](https://rgl.s3.eu-central-1.amazonaws.com/media/papers/Nicolet2021Large.pdf) and the [authors' project page](https://rgl.epfl.ch/publications/Nicolet2021Large).

The inspected reference implementation builds `I + λL` from a combinatorial Laplacian by default, offers a cotangent option, and also supports `(1−α)I + αL`. Its conversion code caches a Cholesky or CG solver for a matrix; the autodiff backward pass solves with that same symmetric matrix. Sources: [geometry.py](https://github.com/rgl-epfl/large-steps-pytorch/blob/master/largesteps/geometry.py), [parameterize.py](https://github.com/rgl-epfl/large-steps-pytorch/blob/master/largesteps/parameterize.py), and [solvers.py](https://github.com/rgl-epfl/large-steps-pytorch/blob/master/largesteps/solvers.py). `AdamUniform` uses one shared maximum second-moment scale rather than independent coordinate scales; see [optimize.py](https://github.com/rgl-epfl/large-steps-pytorch/blob/master/largesteps/optimize.py).

The transfer below is our derivation. Applying one inverse graph operator to a dual residual is a first-power preconditioner; it should not be presented as reproducing the paper's differential-coordinate optimizer.

## Keep this frozen system fixed

Let `J` contain the material constraint gradients as rows, let `g = ∇E`, and let `d = [Δλ; Δη]`. Use the existing scaled stationarity residual `r_x` and compliance convention:

\[
M\Delta x-J^T\Delta\lambda-g\Delta\eta=-r_x,
\qquad
J\Delta x+\widetilde\alpha\Delta\lambda=-r_C,
\qquad
g^T\Delta x+\widetilde\beta\Delta\eta=-r_E.
\]

Define

\[
A=\begin{bmatrix}J\\g^T\end{bmatrix},\quad
T=\operatorname{diag}(\widetilde\alpha,\widetilde\beta),\quad
b=-M^{-1}r_x,\quad
q=\begin{bmatrix}r_C\\r_E\end{bmatrix}+Ab.
\]

If the experiment drops `r_x`, set `b=0` in **every** compared schedule. Elimination gives

\[
S d=-q,\qquad S=AM^{-1}A^T+T,
\qquad \Delta x=b+M^{-1}A^Td.
\]

At a pending increment `d`, the original dual residual is `ρ = S d + q`. A preconditioned residual iteration is

\[
d\leftarrow d-\omega P^{-1}\rho.
\]

For finite invertible `P` and nonzero `ω`, its zero-update points are exactly `ρ=0`. This statement concerns the chosen frozen system: dropping `r_x` is still an approximation to the full constrained objective. A directional rejection rule can introduce additional zero-update points, so test and report that rule separately.

Do not replace `S` with `P` in residual evaluation, add a graph penalty to `r_C`, or smooth `x` after reconstruction. An approximate inverse changes the iterative path; those alternatives change the equations or the applied solution.

## First chain implementation: exact material-block preconditioning

Split the unchanged Schur matrix as

\[
S=\begin{bmatrix}H&B\\B^T&D\end{bmatrix},\quad
H=JM^{-1}J^T+\widetilde\alpha,\quad
B=JM^{-1}g,\quad D=g^TM^{-1}g+\widetilde\beta.
\]

For a chain, consistently orient edge `e=(i,i+1)` and put its unit direction `n_e` into `J`. Write `w_i=1/m_i` for a free vertex and zero for a fixed vertex. Then

\[
H_{ee}=w_i+w_{i+1}+\widetilde\alpha_e,\qquad
H_{e,e+1}=-w_{i+1}\,n_e^Tn_{e+1}.
\]

All other off-diagonal entries vanish. Positive material compliance makes `H` positive definite. Factor this tridiagonal matrix once per frozen linearization in `O(number of springs)` time and memory, and reuse that factor for the substep's residual applications.

Start with the block-diagonal preconditioner `P=diag(H,D)`. One application solves `H z_C=ρ_C` and `D z_E=ρ_E`, then returns `d←d−ω[z_C;z_E]`. Also compare the two block-triangular applications: material block then energy (LE), and energy then material block (EL). These are graph-wide versions of the same schedules, using the original cross term `B`; they do not form an inverse of the nearly singular fully coupled block.

Keep the existing compatibility checks for zero-gradient or dependent hard energy rows. In particular, the exact energy Schur complement

\[
\sigma=D-B^TH^{-1}B
\]

can be tiny. An improved material solve cannot make an incompatible energy row solvable or make the global energy constraint select the physically correct trajectory. The bordered tridiagonal exact solve is useful as an oracle, not evidence that coupling is safe at finite precision.

## Optional reusable topology-only approximation

For comparison with the paper's reusable graph operator, construct the scalar oriented incidence matrix `Q` for the same chain and eliminate fixed-vertex columns. Let `W=diag(w_i)` and

\[
H_0=QWQ^T+\widetilde\alpha.
\]

This replaces neighboring direction dot products by one. It is exact for a consistently oriented straight chain and remains positive definite for positive compliance. It has the correct mass and boundary terms; a plain line-graph Laplacian with endpoint degree one generally does not. Bent or folded chains can make this approximation poor, which the frozen residual test will expose.

For fixed masses, compliance, topology, and substep size, `H_0` can be factored once. It also extends to a sparse shared-vertex constraint graph. A conservative tuning axis is the SPD interpolation

\[
P_C(\theta)=(1-\theta)\operatorname{diag}(H)+\theta H_0,
\qquad 0\le\theta\le1.
\]

Use `diag(P_C(θ),D)` only as the inverse-application metric; retain the original `S` in `ρ`. For uniform coefficients, `H_0` has the shifted-graph form behind Large Steps. Do not automatically square its inverse: if `P_C` already approximates the material Schur block, one inverse is the natural linear-solver preconditioner.

## Small experiment that isolates the mechanism

1. Freeze identical states for straight, bent, and folded spring-pendulum chains of 3, 5, and 8 masses. Preserve the objective, endpoint energy, compliance, `r_x` policy, and one linearization per substep.
2. Compare ordinary scalar GS, diagonal/Jacobi, exact chain material-block preconditioning, and topology-only preconditioning at equal counted work. Use the same acceptance filter for every schedule. Run a no-filter frozen-system diagnostic separately so rejection is not mistaken for preconditioner convergence.
3. Report the **original** scaled dual residual before and after application; material and energy row residuals separately; the energy displacement before and after material response; predicted versus actual energy change; and the final physical state error. Use normalized residuals for comparisons across material and energy units.
4. Check that a converged solve matches an independent solve of the unchanged frozen system when it is compatible. Test near-zero `σ` explicitly and report incompatibility instead of silently adding compliance.
5. Only then run the existing 4/8/16-substep trajectories and long-duration holdouts. Lower frozen residual is evidence of a better linear solve, not a proof of energy-conserving physical motion.

Graph preconditioning retains the reconstruction `M⁻¹JᵀΔλ`; internal material corrections therefore still come from equal-and-opposite central spring directions at the frozen geometry. This preserves the corresponding force-nullspace identities. It does **not** by itself prove discrete angular-momentum conservation, and the energy multiplier's contribution must be audited separately.

Source inspection performed 2026-10-25. No reference package was installed, no rendering performance claim was tested, and no solver code was changed for this note.
