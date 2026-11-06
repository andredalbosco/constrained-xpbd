# Independent rigid-rod audit for multiplier Gauss–Seidel

Scope: the **same** objective and five constraints in `hour/local_solver/problem.hpp`, for one anchored planar rod, positive gravity, beta=gamma=1 by default. No new physical method or code is proposed here. The implementation currently rejects multiple rigid rods; the conclusions below do not authorize removing that guard.

## 1. Feasibility and counting

Write the unknowns as `(x,v,W_r,W_g)` in 2D. There are six scalar unknowns and five scalar constraints:

| Row | Physical equation | Work signs |
|---|---|---|
| Length P | `||x|| - L = 0` | — |
| Tangency T | `n(x) dot v = 0` | — |
| Rod storage R | `W_r = 0` | Rod has no potential storage |
| Gravity storage G | `U_g(x)-U_g(x_n)+W_g = 0` | Positive delivered work lowers potential |
| Kinetic K | `K(v)-K(v_n)-W_r-W_g = 0` | Positive delivered work raises kinetic energy |

The three balance rows imply total energy conservation. They are **not generically redundant**: their sum is the nonzero function `Delta K + Delta U`, rather than an identity. Adding another total-energy row would be exactly redundant. Do not add a second copy of length when the rigid material multiplier uses alpha=0: use the existing length row/multiplier.

Eliminating W gives `W_r=0`, `W_g=-Delta U_g`, and geometry+tangency+energy. In angle coordinates, for gravity `(0,-g)`,

\[
x=L(\sin\theta,-\cos\theta),\quad
v=L\omega(\cos\theta,\sin\theta),\quad
\omega^2=\omega_n^2+\frac{2g}{L}(\cos\theta-\cos\theta_n).
\]

The admissible radicand must be nonnegative. These constraints leave a one-dimensional energy curve at regular states. They contain no h: the old feasible state itself demonstrates existence for every finite h. The objective selects a point and a branch on this set. Thus **the current explicit-W, soft-midpoint formulation is not overdetermined**.

The frozen work estimate is a penalty target, not an extra equality. Eliminating W from the same objective leaves the gravity term `gamma (Delta U_g + W_hat_g)^2/(2 E^2)` and a constant rod-work penalty. Keep this term when using elimination as an independent reduced control.

## 2. Endpoint velocity targets: soft versus hard

A soft prescribed velocity target does not change feasibility. If a vector `v_bar != 0` is imposed as a hard equality, tangency restricts x to the two points `x = +/- L R90(v_bar)/||v_bar||`; energy must also match at one of those points. An arbitrary prescribed vector therefore usually has no feasible endpoint. If `v_bar=0`, feasibility instead requires `U_g(x)=E_initial`; this is impossible for rotational energies above the maximum potential. A prescribed position on the circle similarly allows only the two tangent velocity signs with the energy-prescribed speed, or no velocity if its potential exceeds available energy.

Hard endpoint reconstruction is particularly restrictive:

- `v=(x-x_n)/h` plus endpoint tangency and equal lengths gives `x dot x_n=L^2`, hence **x=x_n and v=0**.
- For `v=2(x-x_n)/h-v_n`, write `delta=theta-theta_n` and assume the old velocity is tangent. Tangency gives `2(1-cos delta)/h = omega_n sin delta`. Besides the zero-chord branch, the nonzero branch has `tan(delta/2)=h omega_n/2` and `omega_next=omega_n`. Energy then also requires zero gravitational work, which is generically incompatible with that angle change.
- The zero-chord midpoint branch remains `x=x_n, v=-v_n`, an unphysical velocity reversal that satisfies geometry and energy. At old rest, hard midpoint reconstruction freezes the position even away from vertical equilibrium.

These are reasons to preserve the **existing soft** midpoint penalty. They are not blockers for the present objective. In that penalty, the conditional target `2(x-x_n)/h-v_n` depends on the unknown x; it is not an additional prescribed endpoint velocity equation.

## 3. Local metric and row derivatives

Use exactly the evaluator coordinates `q=(x-x_n)/x_s`, `u=v/v_s`, `w=W/E`, with `x_s=h v_s`, `v_s=sqrt(E/m)` for a single mass. Put `n=x/||x||`, `P=I-n n^T`,

\[
a=\frac{x_s}{\|x\|}P u,\qquad
b=-\frac{m x_s}{E}g.
\]

For each spatial coordinate the rigid objective Hessian and inverse are

\[
H_{qu}=\begin{pmatrix}1+4/\beta&-2/\beta\\-2/\beta&1/\beta\end{pmatrix},
\qquad
B_{qu}=H_{qu}^{-1}=\begin{pmatrix}1&2\\2&\beta+4\end{pmatrix}.
\]

The work part is `H_ww=gamma I`, `B_ww=I/gamma`. Gravity is linear and contributes no objective Hessian. Consequently this B is the **exact constant objective inverse Hessian** for the rigid problem. It is not the inverse of the full Lagrangian Hessian, which additionally contains multiplier-weighted constraint curvature and can be indefinite. A GS method using B should be described as objective-metric/SQP-style iteration, not exact-curvature Newton.

In the coordinate order `(q,u,w_r,w_g)`, the dimensionless row gradients and their primal directions are:

| Row | `j_i = grad c_i` | `B j_i` |
|---|---|---|
| P | `(n,0,0,0)` | `(n,2n,0,0)` |
| T | `(a,n,0,0)` | `(a+2n,2a+(beta+4)n,0,0)` |
| R | `(0,0,1,0)` | `(0,0,1/gamma,0)` |
| G | `(b,0,0,1)` | `(b,2b,0,1/gamma)` |
| K | `(0,u,-1,-1)` | `(2u,(beta+4)u,-1/gamma,-1/gamma)` |

A geometry multiplier therefore generally pushes both position and velocity. A kinetic multiplier pushes position, velocity and both work variables. Updating only the variable suggested by a row's name changes this metric treatment. Omitting the position derivative `a` from tangency is also not its full linearization.

## 4. A minimum valid GS update, with explicit signs

Use the parent's convention

\[
\mathcal L=F-\lambda_P c_P-\tau c_T-\nu_r c_R-\nu_g c_G-\mu c_K,
\qquad r=\nabla F-J^T\ell.
\]

At an outer linearization, freeze J and B and initialize the primal increment `d=-B r`. Visit scalar rows in the selected order, always using the **latest shared d**:

\[
s_i=j_i^TBj_i,\qquad
\delta\ell_i=-\frac{c_i+j_i^Td}{s_i},\qquad
d\leftarrow d+Bj_i\delta\ell_i.
\]

Accumulate each multiplier increment, then commit/relinearize with a documented globalization rule. For a row block I replace the scalar division by

\[
(J_I B J_I^T)\delta\ell_I=-(c_I+J_I d),\qquad
d\leftarrow d+B J_I^T\delta\ell_I.
\]

This is actual multiplier GS on the frozen linearized Schur system, expressed through primal pushes. It is distinct from primal coordinate descent followed by a simultaneous AL multiplier update. The `-B r` predictor is essential for optimizing the same F: pure feasibility projections started arbitrarily are not sufficient. If using plus signs in the Lagrangian, negate all multipliers and use minus primal pushes consistently.

Freeze old state, old physical acceleration and work targets for the entire step. Freeze J/B for a controlled linear-GS pass; re-evaluation after every row is a different nonlinear-GS control and needs its stationarity treatment documented. Do not accumulate multipliers against changing normals and then declare success from constraints alone. Check both all constraint residuals and full final `grad F-J^T ell`.

## 5. Equilibrium rank defect: positive scalar denominators do not detect it

At `v=0` and `g` parallel to n, put `b=b_n n`. The gradients obey

\[
j_K+j_R+j_G-b_n j_P=0.
\]

The five-row Jacobian has rank four. Nonetheless every scalar denominator is positive:

\[
s_P=1,\quad s_T=\beta+4,\quad s_R=1/\gamma,\quad
s_G=\|b\|^2+1/\gamma,\quad
s_K=(\beta+4)\|u\|^2+2/\gamma.
\]

Thus a check of individual denominators cannot reveal this dependence. The geometry block is itself regular: at rest its Schur matrix is `[[1,2],[2,beta+4]]` with determinant beta. The degeneracy appears across geometry and energy rows. Near equilibrium expect a small full Schur eigenvalue, slow GS and potentially drifting nonunique multipliers despite healthy diagonals.

At the hanging zero-energy equilibrium, the feasible geometry/energy set is locally an isolated point, but its first-order linearization has extra null directions. At upright rest it has crossing separatrix branches. At a nonvertical turning point, v=0 alone does **not** destroy rank: the tangential gravitational gradient is nonzero. An infeasible vertical-rest iterate can have an inconsistent linearized energy equation; no amount of linear GS then repairs that linearization without a nonlinear displacement/reinitialization/globalization step.

Keep every physical residual. For a compatible singular control, a rank-revealing solve/minimum-norm multiplier increment or an explicitly recorded dual gauge can select among equivalent multipliers. Never silently delete the energy, rod-work or tangency equation. Report rank, compatibility residual, complete stationarity and multiplier norm. A separate total-energy row or duplicate rigid-length row would create artificial dependence at all states and must be excluded by construction.

## 6. Minimum ordering controls and checks

Use fixed groups: **geometry** `[P,T]`; **storage/material work** `[R,G]`; **kinetic** `[K]`. The rigid alpha=0 material multiplier is P, not another visit under another name. Compare all six group permutations with the same metric, start, multiplier persistence, globalization, number of visits and stopping criteria. Keep within-group order fixed initially. Then compare a two-row `[P,T]` block with its scalar counterpart; label that as a block-size change, not merely ordering. A reversed or symmetric pass doubles visits unless its budget is adjusted.

Before trajectories, check: (1) the five row values/derivatives and push signs against the evaluator; (2) one frozen GS pass against explicit Schur GS; (3) all scalar denominators versus full rank at bottom rest, top rest, nonvertical rest, small nonzero velocity and zero gravity; (4) complete KKT residuals after every converged solve; (5) feasibility against the angle energy curve and branch continuity. Without gravity omit `W_g` and G; at rest K and R then have dependent gradients directly. Do not manufacture a zero-gravity storage variable just to keep array sizes fixed.

**Blocker:** none for the current single-rod objective with soft midpoint kinematics. Generic hard endpoint-velocity/reconstruction equalities are incompatible additions, and exact-equilibrium multiplier uniqueness cannot be demanded. Multi-rod implementation remains outside this evaluator's supported scope.
