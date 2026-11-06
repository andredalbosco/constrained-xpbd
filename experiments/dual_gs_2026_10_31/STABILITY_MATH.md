# Local stability of a complete multiplier sweep

This diagnostic keeps the objective, constraints, multipliers, supports, and numerical visit formulas unchanged. It examines the existing ordered map near a dense-oracle KKT root. It does not introduce a solver or establish a global convergence guarantee.

## Fixed point and row derivatives

Write the normalized primal and dual variables as `w=(z,y)`, the constant base metric as `A=∇²f₀`, and row primitives as `b_a(z)`. Let

\[
g_a=\nabla b_a,\quad G=[g_1\ldots g_q],\quad
r_a=b_a+D_a y_a,\quad r_z=\nabla f_0-Gy,
\]
\[
H=\nabla^2 f_0-\sum_a y_a\nabla^2 b_a.
\]

Only material rows have positive `D_a`; all balance and geometry rows have zero compliance. The dense solution is first checked for `r_z=0` and every `r_a=0` in the lifted engine. At such a root every existing visit is the identity on that point.

For the structural support of row `a`, let `P_a` select its coordinates and define

\[
B_a=P_a^T(P_a A P_a^T)^{-1}P_a,\qquad
d_a=g_a^TB_ag_a+D_a>0.
\]

The support is closed under the per-particle position/velocity metric. Thus `B_a` is precisely the local inverse used by the engine, embedded in the full normalized primal space. With `e_a` the dual coordinate vector, define

\[
u_a=\begin{bmatrix}B_ag_a\\e_a\end{bmatrix},\quad
v_a^T=\begin{bmatrix}g_a^T&D_a e_a^T\end{bmatrix},\quad
R=\begin{bmatrix}H&-G\end{bmatrix}.
\]

The dropped-correction visit is

\[
\delta y_a=-r_a/d_a,\quad
z^+=z+B_ag_a\delta y_a,\quad y^+=y+e_a\delta y_a.
\]

Its derivative at the root is

\[
J_a^{drop}=I-u_av_a^T/d_a.
\]

The locally retained correction uses

\[
\delta y_a=(-r_a+g_a^TB_ar_z)/d_a,\quad
z^+=z-B_ar_z+B_ag_a\delta y_a.
\]

Its derivative at the root is

\[
J_a^{local}=I-
\begin{bmatrix}B_a\\0\end{bmatrix}R
+\frac{u_a}{d_a}\left(-v_a^T+g_a^TB_aR\right).
\]

Derivatives of `g_a` and `d_a` multiply zero update numerators at the root, so those terms vanish. Constraint curvature still appears through `H` in the derivative of stationarity; dropping curvature inside the local solve does not remove it from the nonlinear iteration's stability.

## Ordered sweep and spectral criterion

For the actual scheduled visit sequence `a₁,…,a_m`, including repeated visits to a shared kinetic multiplier,

\[
T=V_{a_m}\circ\cdots\circ V_{a_1},\qquad
J_T=J_{a_m}\cdots J_{a_1}.
\]

The order is therefore part of the iteration, even at a fixed root. A complete edge sweep has more visits than a once-per-row pass. This diagnostic compares **one complete traversal** of each order, not equal visit budgets. There is one useful exception to expecting different spectra: for these compliant-chain cases, material→storage→kinetic and storage→kinetic→material are cyclic rotations of the same group products. At the common fixed point their square matrix products have the same eigenvalues. Their transients, nonlinear paths, and states at a truncated visit budget can still differ.

For a continuously differentiable map at an exact fixed point, `ρ(J_T)<1` establishes local asymptotic attraction, while an eigenvalue with modulus greater than one establishes linear instability. A measured value very near one is inconclusive numerically; it can also imply extremely slow convergence. Spectral radius does not bound transient growth in a nonnormal map, and attraction near one root says nothing about whether the standard initial guess reaches its basin or whether another KKT root is attractive.

The dropped map has a special structural limitation. Any `(z,y)` satisfying all `b_a(z)+D_a y_a=0` is a fixed point, even if `r_z≠0`. For any tangent perturbation with

\[
\begin{bmatrix}G^T&D\end{bmatrix}\delta w=0,
\]

every `J_a^{drop}` and hence the complete sweep acts as the identity. Unit eigenvalues are therefore expected on this feasible manifold. They demonstrate that balance/material satisfaction alone cannot select a stationary point; they do not by themselves predict explosive divergence.

## Numerical measurement

The probe directly calls the unchanged `visit_refreshed` in schedule order, starting from arbitrary normalized `(z,y)`. For each coordinate it computes

\[
(J_{\epsilon})_{:j}=
\frac{T(w_*+\epsilon e_j)-T(w_*-\epsilon e_j)}{2\epsilon},
\qquad \epsilon\in\{10^{-5},10^{-6}\}.
\]

It reports the oracle's convergence, lifted residual before and after the unperturbed sweep, fixed-point defect, full-state dimension, largest eigenvalue modulus and real/imaginary parts, count of approximately unit eigenvalues, matrix 2-norm, and agreement between finite-difference scales. A finite-difference result is not accepted when the oracle fails its convergence gate, a perturbed map fails, the unperturbed sweep moves by more than `1e-7`, the relative Jacobian disagreement exceeds `1e-3`, or spectral radii differ by more than `1e-3*max(1,ρ)`. These are diagnostic reliability gates, not solver tolerances. The matrix 2-norm diagnoses possible transient amplification, not convergence.

For an independent cross-check, the probe also assembles the derived row-derivative product. It obtains the material-eliminated Lagrangian Hessian from the independent hour evaluator, then reconstructs the lifted Hessian using

\[
H_{lift}=H_{eliminated}-\sum_{a\in material}\frac{g_ag_a^T}{D_a}.
\]

This follows by differentiating the eliminated quarter spring energy: its extra Hessian relative to the lifted Lagrangian at `y_a=-b_a/D_a` is `g_a g_a^T/D_a`. The derivative product describes the zero-residual root linearization; at an approximate root with a highly amplifying sweep, it need not agree with finite differences of the actual floating-point sweep. Those cases are explicitly marked rather than silently treated as a precise numerical derivative.

Cases are the nominal two-spring scene and the existing `chain_case(3)`/`chain_case(8)` constructors, at `h=.02,.08,.2`. Orders are material→storage→kinetic, storage→kinetic→material, and reverse edge. The two methods are refreshed dropped correction and refreshed local correction, with relaxation one. This finite-difference study is numerical evidence; it is not an interval-certified proof.
