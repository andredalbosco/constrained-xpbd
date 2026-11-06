# Multiplier-coordinate GS for the repaired independent-position/velocity objective

Math review only; no implementation or numerical result is asserted here. The primary proposal is to lift the quarter spring energy into compliant material multipliers, then compare three explicit numerical updates of that same formulation: refreshed dropped stationarity, refreshed locally retained stationarity, and frozen GS with one complete stationarity-restoration pass. The objective-only form is an equivalence oracle.

## 1. Objective first, before introducing multipliers

Let the old state be `(x_n,v_n)`, the step be `h`, and `a_n` be its physical acceleration. For rigid systems, `a_n` includes the old constraint reaction. A compliant spring has extension `C_j(x)` and energy `U_j(x)=k_j C_j(x)^2/2`. Let `U_g` be gravitational potential and let `W_hat` be a work prediction computed entirely from the old state and frozen during this step.

Fix positive scales `E` and `m_ref`, and define

\[
v_s=\sqrt{E/m_{ref}},\qquad x_s=h v_s,\qquad
q=(x-x_n)/x_s,\quad u=v/v_s,\quad w=W/E.
\]

The primal variable is `z=(q,u,w)`. With `beta>0`, `gamma>0`, the objective is

\[
\Phi(z)=\frac{1}{2Eh^2}
 \|x-x_n-hv_n-\tfrac14h^2a_n\|_M^2
 +\frac{U_g(x)+\sum_j U_j(x)}{4E}
 +\frac{1}{2\beta E}\|v+v_n-2(x-x_n)/h\|_M^2
 +\frac\gamma2\sum_{ij}(w_{ij}-\widehat w_{ij})^2.
\]

For each free particle and each potential store, impose respectively

\[
k_i(z)=\frac{K_i(v_i)-K_i(v_i^n)}{E}-\sum_j w_{ij}=0,
\qquad
s_j(z)=\frac{U_j(x)-U_j(x_n)}{E}+\sum_i w_{ij}=0.
\]

Gravity may have its own storage leaf and transfer per particle. Work incidence signs are fixed by these equations. A fixed endpoint contributes zero work and is not an independent kinetic node. For a rigid edge, its storage energy is zero; its transfer row can still enforce zero net work. Position and endpoint-tangency constraints are additional hard rows where required.

Positions, velocities and actual transfers are independent. In particular, the kinetic rows have no position dependence, and storage rows have no velocity dependence. Neither the work target nor its predictor is differentiated with respect to the unknown endpoint.

This review holds the physical objective and its normalization fixed across numerical variants. Changing the energy scale without compensating the work weight changes physical weighting; it is not a solver option.

## 2. Lift the spring energy after formulating the objective

Write `a=1/4`, and for every compliant spring define

\[
c_j(z)=C_j(x)/x_s,\qquad
\alpha_j=\frac{E}{a k_j x_s^2}
=\frac{4m_{ref}}{k_jh^2}>0.
\]

For fixed `z`,

\[
\max_{\lambda_j\in\mathbb R}
 \left[-\lambda_j c_j(z)-\tfrac12\alpha_j\lambda_j^2\right]
=\frac{c_j(z)^2}{2\alpha_j}=\frac{U_j(x)}{4E},
\qquad \lambda_j^*=-c_j/\alpha_j.
\]

Let `f_0` be `Phi` with the quarter **spring** potential removed. It retains inertia, the velocity penalty, the work penalty, and quarter gravity. The lifted Lagrangian, with the requested negative-multiplier convention, is

\[
\mathcal L=f_0(z)
-\sum_j\lambda_jc_j(z)-\frac12\sum_j\alpha_j\lambda_j^2
-\sum_i\mu_i k_i(z)-\sum_j\nu_j s_j(z)
-\text{additional hard-row multiplier terms}.
\]

Do not also retain quarter spring energy in `f_0`: that would count it twice. Material multipliers `lambda`, kinetic multipliers `mu`, and storage multipliers `nu` have different roles even when their gradients involve the same edge.

All these multipliers are unrestricted real variables. There are no inequality constraints and **no complementarity conditions**. The negative `-alpha lambda^2/2` term is concave in the material multiplier; its negative saddle diagonal produces a **positive** compliance contribution in the Schur system. A hard rod is represented by a hard length row, not by evaluating a spring formula at infinite stiffness.

Collect the scalar functions `c_j`, `k_i`, `s_j`, and geometric rows into a vector `b(z)`, and collect their multipliers into `y`. Let `D` be diagonal, with `D_jj=alpha_j` for material rows and zero for hard rows. Define the matrix of gradient columns

\[
G(z)=[\nabla b_1(z)\;\cdots\;\nabla b_m(z)].
\]

Stationarity and row residuals are

\[
r_z=\nabla f_0-Gy=0,\qquad r=b(z)+Dy=0.
\]

Thus material stationarity gives `c_j+alpha_j lambda_j=0`; kinetic and storage rows remain hard. Eliminating material multipliers gives exactly the stationarity of the original objective. Work stationarity provides a useful sign check:

\[
\gamma(w_{ij}-\widehat w_{ij})+\mu_i-\nu_j=0,
\qquad w_{ij}=\widehat w_{ij}+(\nu_j-\mu_i)/\gamma.
\]

## 3. Full Newton system and its local metric approximation

At the current iterate, the exact primal Hessian is

\[
H=\nabla^2 f_0-\sum_a y_a\nabla^2 b_a.
\]

The full Newton equations are

\[
\begin{bmatrix}H&-G\\G^T&D\end{bmatrix}
\begin{bmatrix}\Delta z\\\Delta y\end{bmatrix}
=\begin{bmatrix}-r_z\\-r\end{bmatrix}.
\]

Multiplying the second block row by minus one gives the symmetric saddle form

\[
\begin{bmatrix}H&-G\\-G^T&-D\end{bmatrix}
\begin{bmatrix}\Delta z\\\Delta y\end{bmatrix}
=\begin{bmatrix}-r_z\\r\end{bmatrix}.
\]

For constant gravity, `f_0` is quadratic and its Hessian `A=nabla^2 f_0` is constant. In the `(q_i,u_i)` coordinates its particle blocks and inverses are

\[
A_i=\frac{m_i}{m_{ref}}
\begin{bmatrix}1+4/\beta&-2/\beta\\-2/\beta&1/\beta\end{bmatrix}\otimes I,
\qquad
A_i^{-1}=\frac{m_{ref}}{m_i}
\begin{bmatrix}1&2\\2&\beta+4\end{bmatrix}\otimes I.
\]

Each work diagonal is `A_ww=gamma`, with inverse `1/gamma`. The particle determinant before the mass factor is `1/beta>0`; therefore `A` is SPD for positive masses, beta and gamma. Its inverse consists only of per-particle position/velocity blocks and scalar work entries. A kinetic-row push can consequently move both position and velocity; a storage/material push can also move both.

Replacing `H` by `A` explicitly drops **all multiplier-weighted curvature**, including material, kinetic, storage and geometric curvature. It changes the numerical direction, not the objective or final KKT conditions. It is not exact Newton for nonlinear rows.

Eliminate the primal increment from the approximate system:

\[
\Delta z=-A^{-1}r_z+A^{-1}G\Delta y,
\qquad
S\Delta y=t,
\]
\[
\boxed{S=G^TA^{-1}G+D,\qquad t=-r+G^TA^{-1}r_z.}
\]

The objective-only oracle has the same derivation with `nabla Phi` replacing `nabla f_0`, only hard multipliers retained, and `D=0`. Its Hessian includes quarter spring curvature. The lifted and eliminated objectives agree, but their finite approximate Newton/GS paths need not agree.

## 4. Exact GS of one frozen approximate Newton system

Freeze `z_0,y_0,G_0,r_z0,r_0,A,D` for one outer sweep. Let `d` denote the **total** multiplier increment relative to `y_0`, initially zero. For a visit to row `a`, solve

\[
d_a^{new}=\frac{t_a-\sum_{b\ne a}S_{ab}d_b}{S_{aa}},
\qquad \delta=d_a^{new}-d_a.
\]

Use the latest available `d_b`: that is the Gauss–Seidel dependency. Repeated row visits replace the total coordinate by its newly solved value; they do not add that total again.

An implementation need not assemble `S`. First perform the stationarity restoration **exactly once**:

\[
z\leftarrow z_0-A^{-1}r_{z0},\qquad y\leftarrow y_0.
\]

At each frozen row visit, evaluate its frozen linear residual

\[
\ell_a=b_a(z_0)+\nabla b_a(z_0)^T(z-z_0)+D_{aa}y_a.
\]

Then apply

\[
\delta=-\ell_a/(\nabla b_a(z_0)^TA^{-1}\nabla b_a(z_0)+D_{aa}),
\quad y_a\leftarrow y_a+\delta,
\quad z\leftarrow z+A^{-1}\nabla b_a(z_0)\delta.
\]

These immediate pushes are exactly equivalent to the frozen multiplier GS solve, including repeated rows and different orderings. Deferred reconstruction gives the same iterate in exact arithmetic. The invariance being preserved is the **frozen stationarity equation**

\[
r_{z0}+A(z-z_0)-G_0(y-y_0)=0.
\]

Do not subtract `A^-1 r_z0` again at every row. Do not substitute the nonlinear value `b_a(z)` for its frozen linear expression while calling the method frozen GS. Both changes define different algorithms. Recompute actual nonlinear residuals and gradients at the next outer sweep.

The once-per-sweep restoration can be carried out as an independent local pass over all particle/work metric blocks. It covers the system but requires neither a dense factorization nor a global primal row solve. The subsequent gradient pushes have local structural support.

If `S` is SPD, each exact coordinate update decreases the frozen quadratic `d^T S d/2-t^T d` by `S_aa delta^2/2`, and cyclic GS solves that frozen system in the convergence limit. This says nothing by itself about nonlinear convergence or physical phase. A small fixed number of visits leaves a measurable frozen row residual. Applying a final common damping factor to the whole primal/dual step leaves the corresponding fraction of the stationarity residual; it must be measured, not declared restored.

## 5. Refreshed dropped-stationarity XPBD-style visits

Recompute the visited row and its gradient at the latest state, then use

\[
\delta y_a=-\frac{r_a}{\nabla b_a^TA^{-1}\nabla b_a+D_{aa}},
\qquad \delta z=A^{-1}\nabla b_a\,\delta y_a.
\]

Accumulate `y_a` and immediately push `z`. This is the direct XPBD-style scalar update for the lifted equations. It is not GS of a single frozen linear system: each visit changes the nonlinear row model.

The clean zero-dual initialization is the exact unconstrained minimizer of `f_0`. In physical variables, writing `s=(x-x_n)/h`, it is

\[
s=v_n+\tfrac14h(a_n+g),\qquad
v=2s-v_n,\qquad W=\widehat W,\qquad y=0,
\]

where `g` is gravitational acceleration. This includes quarter gravity already present in `f_0`; adding the full gravity predictor again would be wrong. It starts with `r_z=0`, although material and hard balances generally are not satisfied.

Initialization does not keep true stationarity zero under nonlinear row visits. Since `f_0` has constant Hessian `A`, the dropped-stationarity push satisfies the exact identity

\[
\boxed{r_z^{new}=r_z^{old}-[G(z^{new})-G(z^{old})]y^{new}.}
\]

For affine rows, `G` is constant and zero stationarity is preserved. For nonlinear rows it drifts. Even an affine scalar spring retains nonlinear kinetic/storage rows. Furthermore, any feasible state with nonzero `r_z` is a fixed point of the dropped update: every row numerator is zero, so nothing restores stationarity. Its acceptance criteria must therefore report hard/compliant row residuals **and** `r_z` independently.

Warm multipliers generally destroy the zero-residual initialization. They require a corresponding primal initialization/restoration or an explicitly measured initial `r_z`; they cannot silently inherit the fresh-zero-dual argument.

## 6. Refreshed retention: global diagnostic versus genuinely local update

### Global-primal diagnostic

At each row, retaining the complete current stationarity residual gives

\[
\delta y_a=
\frac{-r_a+\nabla b_a^TA^{-1}r_z}
 {\nabla b_a^TA^{-1}\nabla b_a+D_{aa}},
\qquad
\delta z=A^{-1}\nabla b_a\delta y_a-A^{-1}r_z.
\]

This is a valid single-active-multiplier approximate KKT step. It changes every primal metric block having a stationarity residual, even when the visited constraint is local. It must be labeled a **global-primal diagnostic**, not a local XPBD visitor. Its restoration is recomputed at each visit; it is distinct from accidentally applying a frozen restoration repeatedly.

### Local retained-stationarity visitor

For each row choose a fixed **structural** primal support `S_a`, closed under the per-particle `(q_i,u_i)` metric coupling. Include the involved particle coordinates even if the instantaneous gradient happens to vanish. Let `A_S` be the restricted metric, and define

\[
p=A_S^{-1}r_z[S_a],\qquad
d=A_S^{-1}\nabla b_a[S_a].
\]

With every other primal coordinate and multiplier fixed, solve

\[
\begin{bmatrix}A_S&-\nabla b_a[S_a]\\
\nabla b_a[S_a]^T&D_{aa}\end{bmatrix}
\begin{bmatrix}\delta z[S_a]\\\delta y_a\end{bmatrix}
=\begin{bmatrix}-r_z[S_a]\\-r_a\end{bmatrix}.
\]

Thus

\[
\boxed{\delta y_a=
 \frac{-r_a+\nabla b_a[S_a]^Tp}
 {\nabla b_a[S_a]^Td+D_{aa}},\qquad
\delta z[S_a]=d\,\delta y_a-p.}
\]

This is genuinely local restricted KKT/quasi-Newton, not a global-primal row block. It approximately restores stationarity on the selected variables and the selected row. Changes to neighboring rows are expected and are handled by subsequent visits.

Suggested fixed supports:

- Material or length row: both free endpoints' full position/velocity blocks.
- Kinetic row at particle `i`: its full position/velocity block and **all** incident work variables, including gravity.
- Spring storage row: both free endpoint position/velocity blocks and that spring's incidence works.
- Rigid zero-storage row: its incidence works only; its identically zero potential has no position/velocity support.
- Gravity-storage row: its particle position/velocity block and gravity work.
- Endpoint tangency row: both participating position/velocity blocks.

To compute `r_z[S_a]`, gather contributions from **all** multipliers whose gradients touch those variables, not only the visited row. This includes adjacent springs and shared kinetic multipliers. A work coordinate needs its owning kinetic and storage multipliers. Such incident-row gathering is sparse; no dense global matrix is mathematically required.

If every row is visited, the supports cover every primal variable, and every individual primal/dual increment vanishes, the local equations imply all row residuals and all components of `r_z` vanish. Conversely a KKT point produces zero updates. This fixed-point consistency is not a convergence theorem; overlapping local saddle updates can stall, cycle or diverge. Net zero displacement over an entire cycle is insufficient: intermediate updates can cancel. Monitor actual full residuals and individual visit changes.

Disconnected/unconstrained primal blocks still need coverage, potentially with primal-only stationarity visits. A numerically zero gradient is not permission to remove its structural variables. For local damping, scale both `delta z[S_a]` and `delta y_a` together; recompute actual residuals afterward.

## 7. Shared multipliers and meaningful schedules

Maintain one `mu_i` per kinetic row, one `nu_j` per storage row, and one material `lambda_j` per spring. A shared vertex's multiplier is read and accumulated from the same memory at every visit. Repeated visitation is valid; duplicating the multiplier per incident spring changes the formulation unless a separate consensus system is introduced.

Use explicit row lists, for example:

1. **Separated:** all material rows, all storage rows, then each kinetic row once.
2. **Material/storage pairs plus vertex pass:** visit `(material_j,storage_j)` per edge, then each kinetic row once.
3. **Interleaved:** per edge visit its material/storage rows and the kinetic rows of both free endpoints; shared kinetic rows recur with vertex degree.
4. **Kinetic-first interleaving:** visit the incident kinetic rows before the material/storage pair.
5. Reverse or genuinely palindromic versions of a specified list.

For rigid systems, specify the placement of length and tangency visits explicitly. For frozen GS, every occurrence uses the same frozen model and the latest total multiplier increments. For refreshed methods, every occurrence reevaluates the row and gathers current shared stationarity terms. Do not compare a degree-weighted repeated-vertex list to a single vertex pass as if they had equal cost.

The material row, spring storage row, and particle kinetic row must not be collapsed into a single unspecified “constraint visit.” Count scalar row solves, particle/work updates, incident-gradient gathers and stationarity-restoration passes. Symmetric lists and repeated rows cost extra work. A fixed row-visit budget may end mid-pass; report that boundary rather than treating all sweep counts as equivalent.

## 8. Rank and singular-row boundaries

For any vector `t`,

\[
t^TSt=\|A^{-1/2}Gt\|^2+t^TDt.
\]

Consequently `S` is SPD precisely when `ker(G) intersect ker(D)={0}`. Positive material compliance regularizes its multiplier coordinates, but does not cure dependencies among hard rows. Examples include redundant rigid constraints or adding a total-energy row that is already the sum of all kinetic/storage rows.

With explicit work and finite positive gamma, kinetic/storage row gradients normally retain nonzero work entries even when velocity or force vanishes; this avoids some isolated scalar zero denominators. It does not ensure full rank: at an energy critical point, the sum of all energy-row gradients can vanish, leaving a combined null direction. Hard zero-storage rod rows can be dependent in special states. An undefined distance direction at coincident endpoints is a separate geometric failure.

If a hard row has zero gradient and nonzero residual, its linearized equation is incompatible: report it, rather than divide by a small number or silently declare it solved. If both are zero, the multiplier is locally undetermined. A zero multiplier increment and the available primal stationarity restoration solve that linearized row, but second-order nonlinear feasibility still needs checking. Identically zero rows should be removed structurally; state-dependent degeneracy should be diagnosed. A rank-aware dense oracle can distinguish these cases during validation.

`gamma=0` makes the work block noninvertible; free work needs its separate exact elimination/component-closure treatment. The finite positive-beta metric does not cover a hard kinematic constraint by simply setting beta to zero. Rigid and disconnected graphs require their own consistent row sets; one global energy closure cannot substitute for independent component balances.

## 9. Recommendation before coding

Use the lifted saddle objective as the main family, with the objective-only formulation as a small-system equivalence oracle. Implement three named visitors first:

- `refreshed_drop`: scalar XPBD-style row updates, initialized at the exact `f_0` minimum with zero multipliers; track stationarity drift.
- `refreshed_local_rx`: the structural-support local retained-stationarity update above, gathering all incident multiplier contributions.
- `frozen_full_rx`: freeze one approximate Newton model per outer sweep, perform one all-block stationarity pass, then exact GS of that frozen Schur model by local immediate pushes.

The refreshed global-primal retained-residual update is an optional diagnostic, explicitly outside the local family. Do not replace local visitors by dense full-primal KKT blocks. Preserve identical objective, scales, hard rows, initial state and row-list accounting across the three main choices. Use the full exact Hessian only in the oracle/curvature diagnosis unless an additional numerical variant is named.

Before trajectory comparisons, verify lifted/eliminated objective stationarity, signs of work/material relations, frozen GS versus an independently assembled Schur solve, immediate/deferred equivalence, repeated-row accumulation, support coverage at zero gradients, and the exact dropped-stationarity drift identity. Then compare explicit ordering lists at matched row/incident-update budgets and setup-inclusive wall time. Converged KKT, local minimum, conserved energy and faithful motion remain separate outcomes.
