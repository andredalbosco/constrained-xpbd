# Five explicit-work objectives and their local multiplier support

This is a mathematical specification, not an implementation or a numerical result. It matches the final `fraction=1` objective and row normalization in [hour/local_solver/problem.hpp](/home/blacknoir/Projects/constrained-xpbd/experiments/hour_2026_10_31/local_solver/problem.hpp), with the frozen `implicit_midpoint_dg` target from [predictions.hpp](/home/blacknoir/Projects/constrained-xpbd/experiments/round2_2026_10_31/local_work/predictions.hpp). The five scenes are defined in [common.hpp](/home/blacknoir/Projects/constrained-xpbd/experiments/composition_2026_10_31/common.hpp). The three-mass benchmark below uses the existing hour `chain_case(3)` constructor; it is an extension of the original five-scene menu.

## 1. Common objective, variables, and signs

Use particle indices i=1,...,N, spatial dimension d, fixed anchor x_0=0, and fixed anchor velocity v_0=0. The old physical state is (x_n,v_n). During a timestep h, x, v and every incidence work W_{ji} are independent unknowns. W_{ji}>0 means energy delivered from store j to particle i. An internal spring has a separate work value for each moving endpoint; these two values are not assumed opposite. Their sum accounts for the spring's storage change.

Let M=diag(m_i I_d), K_i(v_i)=m_i||v_i||²/2, and U=Σ_j U_j. Gravity has its own store G_i(x_i)=−m_i gᵀx_i and work W_{gi}; it is not included in a spring's work row. A rod has a zero-energy store U_R=0, plus separate geometric constraints.

The fixed trajectory normalization is

\[
E=\max(10^{-8},H(x_{\rm init},v_{\rm init})-H_{\min}),\quad
v_s=\sqrt{E/\bar m},\quad x_s=h v_s,\quad \bar m=N^{-1}\sum_i m_i.
\]

Here “init” denotes the trajectory initial state, not the fixed anchor. The small-scene hour driver uses this floor; its positive-excitation chain drivers pass the unfloored positive excitation. E is held fixed across that trajectory, not recomputed from the candidate state. Changing this total scale changes the relative physical work preference; it is not merely numerical row scaling.

Define a_i=(x_i−x_{n,i})/h−v_{n,i}−h a_{n,i}/4 and b_i=v_i+v_{n,i}−2(x_i−x_{n,i})/h, where a_n is the **full old physical acceleration**, including old rod reactions when present. Both a_n and the work targets Ŵ are frozen during the step. With beta=gamma=1, minimize

\[
\boxed{\Phi(x,v,W)=\frac1{2E}\sum_i m_i(\|a_i\|^2+\|b_i\|^2)
+\frac{U(x)}{4E}+\frac1{2E^2}\sum_{(j,i)}(W_{ji}-\widehat W_{ji})^2.}
\tag{1}
\]

For each kinetic node and store, impose

\[
k_i=\frac{K_i(v_i)-K_i(v_{n,i})-\sum_{j\ni i}W_{ji}}E=0,
\qquad
s_j=\frac{U_j(x)-U_j(x_n)+\sum_{i\in j}W_{ji}}E=0.
\tag{2}
\]

Gravity stores are included in s with distinct multipliers ν_{gi}. Rod position C and tangency T are defined in §5. The requested sign convention is

\[
\boxed{\mathcal L=\Phi-\sum_i\mu_i k_i-\sum_j\nu_j s_j
-\sum_R\lambda_R C_R-\sum_R\tau_R T_R.}
\tag{3}
\]

The numerical helper in the earlier hour suite used plus-sign row multipliers: its stored dual is the negative of the multipliers in (3). Energy multipliers μ,ν and geometric multipliers λ,τ are distinct; none is identified with the material reaction without an additional derivation.

### Objective derivatives in physical coordinates

Let a,b stack the vectors above. Every unlisted Hessian block is zero:

\[
\begin{aligned}
\nabla_x\Phi&=\frac{M(a-2b)}{Eh}+\frac{\nabla U}{4E},&
\nabla_v\Phi&=\frac{Mb}E,&
\nabla_W\Phi&=\frac{W-\widehat W}{E^2},\\
\nabla^2_{xx}\Phi&=\frac{5M}{Eh^2}+\frac{\nabla^2 U}{4E},&
\nabla^2_{xv}\Phi&=-\frac{2M}{Eh},&
\nabla^2_{vv}\Phi&=\frac ME,&
\nabla^2_{WW}\Phi&=\frac I{E^2}.
\end{aligned}
\tag{4}
\]

The balance derivatives are

\[
\nabla_{v_i}k_i=\frac{m_i v_i}E,\quad
\nabla^2_{v_i v_i}k_i=\frac{m_i I_d}E,\quad
\partial_{W_{ji}}k_i=-\frac1E;
\qquad
\nabla_xs_j=\frac{\nabla U_j}E,\quad
\nabla^2_{xx}s_j=\frac{\nabla^2U_j}E,\quad
\partial_{W_{ji}}s_j=\frac1E.
\tag{5}
\]

All work-containing second derivatives of a balance row vanish. Gravity has ∇G_i=−m_i g and ∇²G_i=0. Frozen Ŵ contributes no state derivative.

Thus, before rod terms,

\[
\begin{aligned}
r_x&=\nabla_x\Phi-\frac1E\sum_j\nu_j\nabla U_j,\\
r_{v_i}&=\frac{m_i}E(b_i-\mu_i v_i),\\
r_{W_{ji}}&=\frac{W_{ji}-\widehat W_{ji}}{E^2}+\frac{\mu_i-\nu_j}E.
\end{aligned}
\tag{6}
\]

Add −λ∇C−τ∇T to the applicable state residuals. In particular,

\[
W_{ji}=\widehat W_{ji}+E(\nu_j-\mu_i)
\tag{7}
\]

at work stationarity. The exact Lagrangian Hessian is

\[
\nabla^2\mathcal L=\nabla^2\Phi-\sum_i\mu_i\nabla^2k_i
-\sum_j\nu_j\nabla^2s_j-\sum_R\lambda_R\nabla^2 C_R
-\sum_R\tau_R\nabla^2T_R.
\tag{8}
\]

For example its compliant xx block is 5M/(Eh²)+Σ_j(1/4−ν_j)∇²U_j/E, and its vv block is diag(m_i(1−μ_i)I_d/E). Neither is guaranteed positive definite. Positive definiteness of a chosen numerical metric must not be mistaken for positivity of this exact Hessian.

### Material lift: an actual compliant material visitor

For each compliant spring, write its physical extension as C_j(x) (affine for the scalar spring, distance minus rest length otherwise), and introduce a **material** multiplier λ_j^m. Keep λ_R^g for a hard rod's geometric multiplier, so neither is confused with storage multiplier ν_j. Define

\[
\overline C_j=C_j/x_s,\qquad
\alpha_j=\frac{4E}{\kappa_jx_s^2}=\frac{4\bar m}{h^2\kappa_j},\qquad
f_0=\Phi-\sum_{j\;\mathrm{compliant}}\frac{U_j}{4E}.
\tag{M1}
\]

Thus f_0 retains the inertial term, independent-v midpoint penalty, work penalties and quarter **gravity** potential, while omitting only compliant spring potentials. Use the lifted saddle Lagrangian

\[
\boxed{\mathcal L_{\rm lift}=f_0
-\sum_j\lambda_j^m\overline C_j-\frac12\sum_j\alpha_j(\lambda_j^m)^2
-\sum_i\mu_i k_i-\sum_j\nu_j s_j
-\sum_R\lambda_R^g C_R-\sum_R\tau_R T_R.}
\tag{M2}
\]

Its material stationarity equation is

\[
\overline C_j+\alpha_j\lambda_j^m=0,\qquad
\lambda_j^m=-\overline C_j/\alpha_j.
\tag{M3}
\]

Maximizing the concave material-dual term at fixed state gives \(\overline C_j^2/(2\alpha_j)=\kappa_j C_j^2/(8E)=U_j/(4E)\), so eliminating λ^m recovers (1) exactly. For a stretched spring λ^m is negative with this sign convention. Its positive endpoint direction is ∇C, so a negative increment pulls that endpoint inward. Storage balances still contain the **full** spring energy change ΔU_j/E; lifting the quarter objective potential does not replace that energy by a quarter-sized storage change.

In normalized z coordinates, ∇²f_0 is exactly the constant block metric A_0 consisting of (12) for each particle and identity on normalized W. Gravity is affine and contributes no Hessian. The full lifted primal Hessian is

\[
\nabla^2_{zz}\mathcal L_{\rm lift}=A_0
-\sum_j\lambda_j^m\nabla^2\overline C_j
-\sum_i\mu_i\nabla^2 k_i-\sum_j\nu_j\nabla^2 s_j
-\sum_R\lambda_R^g\nabla^2 C_R-\sum_R\tau_R\nabla^2 T_R.
\tag{M4}
\]

The spring energy's gradient-outer-product curvature is recovered through eliminating the compliant material dual. One must not compare only the primal block of (M4) with the unlifted Hessian and call the difference an omitted physical term. Approximating (M4) by A_0 is instead an explicitly declared numerical curvature approximation.

Let the columns of G be the gradients of (C̄_j,k_i,s_j,C_R,T_R), y=(λ^m,μ,ν,λ^g,τ), and A=diag(α_j,0,...,0). Define the full residuals

\[
r_z=\nabla f_0-Gy,\qquad
r_c=(\overline C_j+\alpha_j\lambda_j^m,\ k_i,\ s_j,\ C_R,\ T_R).
\]

For the exact lifted Hessian H, the Newton system is

\[
\begin{bmatrix}H&-G\\G^T&A\end{bmatrix}
\begin{bmatrix}\delta z\\\delta y\end{bmatrix}
=-\begin{bmatrix}r_z\\r_c\end{bmatrix}.
\tag{M5}
\]

With the SPD approximation H=A_0, its Schur operator is GᵀA_0^{-1}G+A. An individual material row has the familiar compliant denominator ∇C̄ᵀA_0^{-1}∇C̄+α; kinetic and storage rows have zero compliance. For one row at primal stationarity,

\[
\delta\lambda_j^m=
-\frac{\overline C_j+\alpha_j\lambda_j^m}
{\nabla\overline C_j^TA_0^{-1}\nabla\overline C_j+\alpha_j},
\qquad \delta z=A_0^{-1}\nabla\overline C_j\,\delta\lambda_j^m.
\tag{M6}
\]

Away from stationarity, the Schur right-hand side is −r_c+GᵀA_0^{-1}r_z. The retained primal residual and all shared multipliers therefore matter. Material λ → kinetic μ → storage ν, or another explicitly specified order, is now a sequence of actual compliant/hard row visitors on this lifted formulation. This is a mathematical candidate for dual GS; no nonlinear convergence or equivalence of arbitrary visit schedules is asserted.

For a material edge between particles a and b, normalized position derivatives are (−n_j,+n_j), with the anchored endpoint omitted. Its normalized position Hessian has signed blocks ±x_s(I−n_jn_jᵀ)/r_j. The affine scalar material gradient is 1 with zero Hessian. By (12), a material multiplier moves **both** endpoint x and v, but no W directly:

\[
\delta\xi_i=(\bar m/m_i)\nabla_{\xi_i}\overline C_j\,\delta\lambda_j^m,
\qquad \delta u_i=2\delta\xi_i.
\]

| Example | Additional material rows and support |
|---|---|
| Affine scalar spring | One C̄=(x−ℓ)/x_s, α=4m/(h²κ), gradient 1 in ξ; λ^m moves x,v |
| Planar radial spring | One C̄=(||x||−ℓ)/x_s; gradient n in ξ; curvature x_s(I−nnᵀ)/r; λ^m moves x,v |
| Rigid pendulum | No finite-compliance spring-material lift: U_R=0. Its hard C_R uses λ_R^g with α=0, plus the distinct tangency τ and zero-store ν_R |
| Elastic pendulum | One radial material row as above; gravity stays affine in f_0 and has its own storage ν_g |
| Two-spring pendulum | λ_1^m has support x_1; λ_2^m has support x_1,x_2. Their metric responses include the corresponding v blocks |
| N3 extension | λ_1^m touches particle 1; λ_2^m touches 1/2; λ_3^m touches 2/3. The middle-particle response overlaps, while μ_2 is still one shared kinetic multiplier |

The unconstrained minimizer of f_0 at zero duals is explicit. With s_i=(x_i−x_{n,i})/h,

\[
s_i^*=v_{n,i}+\frac h4(a_{n,i}+g),\qquad
v_i^*=2s_i^*-v_{n,i},\qquad W^*=\widehat W.
\tag{M7}
\]

Use g=0 when gravity is absent. This is the zero-dual base target, not generally a feasible state. It differs from the stiffness-aware midpoint initializer used in the hour experiments; choosing it as a GS initialization is a new numerical choice to be assessed, not an existing measured result.

### Variable and row counts before and after lifting

The primal coordinates (x,v,W) stay unchanged. Each compliant spring adds one material multiplier and one compliant material residual; it does not add another primal position or work value.

| Scene | Primal coordinates | Unlifted hard rows | Lifted material rows | Total dual/residual rows after lift | Full unlifted / lifted KKT dimensions |
|---|---:|---:|---:|---:|---:|
| Scalar spring | 3 | 2 | 1 | 3 | 5 / 6 |
| Planar radial spring | 5 | 2 | 1 | 3 | 7 / 8 |
| Rigid pendulum | 6 | 5 (3 balances, C, T) | 0 | 5 | 11 / 11 |
| Elastic pendulum | 6 | 3 | 1 | 4 | 9 / 10 |
| Two-spring pendulum | 13 | 6 | 2 | 8 | 19 / 21 |
| Three-spring chain | 20 | 9 | 3 | 12 | 29 / 32 |

### The frozen start-only midpoint-DG target

For compliant edges e_j=x_j−x_{j−1}, let c_j=||e_j||−ℓ_j, except that the scalar affine spring uses c=x−ℓ. Assemble the old-state stiffness surrogate

\[
B_n=\sum_j\kappa_j\nabla c_j(x_n)\nabla c_j(x_n)^T,
\quad
(M+h^2B_n/4)\delta\widehat x=hMv_n+h^2 f_n/2,
\quad \widehat x=x_n+\delta\widehat x,\quad
\widehat v=2\delta\widehat x/h-v_n.
\tag{9}
\]

Here f_n=−∇U(x_n). B_n is a positive semidefinite surrogate, not the complete spring Hessian. This predictor requires a linear solve, but uses no optimized endpoint information. For a radial edge set

\[
\overline n_j=\frac{e_{n,j}+\widehat e_j}{\|e_{n,j}\|+\|\widehat e_j\|},\quad
\widehat f_j=-\frac{\kappa_j}{2}
(\|e_{n,j}\|+\|\widehat e_j\|-2\ell_j)\overline n_j.
\tag{10}
\]

Then Ŵ_{j,j}=f̂_jᵀδx̂_j and Ŵ_{j,j−1}=−f̂_jᵀδx̂_{j−1} for a moving parent. The anchor supplies no work variable. Ŵ_{gi}=m_i gᵀδx̂_i. Equation (10) makes the sum of the predicted endpoint works equal −ΔU_j along the **predicted** displacement. It does not impose a DG identity at the optimized endpoint. The rod's distinct predictor is specified below.

## 2. Why a kinetic multiplier can move position

Use the hour coordinates

\[
z=(\xi,u,w),\quad \xi=(x-x_n)/x_s,\quad u=v/v_s,\quad w=W/E.
\tag{11}
\]

If D=diag(x_s I,v_s I,E I), every derivative above transforms as ∇_z f=Dᵀ∇_{x,v,W}f and ∇²_z f=Dᵀ(∇²_{x,v,W}f)D. This specifies all normalized Hessian blocks without changing the objective.

A natural **per-particle numerical metric**, retaining the inertial and velocity-coupling curvature while leaving potential/constraint curvature out of that metric, is

\[
A_i=\frac{m_i}{\bar m}
\begin{bmatrix}5I_d&-2I_d\\-2I_d&I_d\end{bmatrix},\qquad
A_i^{-1}=\frac{\bar m}{m_i}
\begin{bmatrix}I_d&2I_d\\2I_d&5I_d\end{bmatrix},\qquad A_w=I.
\tag{12}
\]

It is SPD. It is **not** the full objective Hessian in (4): a spring between moving particles also has off-particle xx curvature. A different local metric is possible, but must be declared.

For a fixed quadratic model with L=Φ−αc, changing one multiplier by δα gives the primal response δz=A^{-1}∇_z c δα when the other duals are held fixed and the primal model is re-solved. Consequently:

| Row / multiplier | Nonzero row derivatives | Response through (12), before any global correction |
|---|---|---|
| k_i / μ_i | v_i and **all** incident W_{ji} | ξ_i changes by 2(m̄/m_i)(m_i v_s v_i/E)δμ_i; u_i changes by 5(m̄/m_i)(m_i v_s v_i/E)δμ_i; every incident w decreases by δμ_i |
| s_j / ν_j, compliant spring | Positions of its moving endpoints and its endpoint W values | At each endpoint let q_i=(x_s/E)∇_{x_i}U_j: δξ_i=(m̄/m_i)q_iδν_j, δu_i=2(m̄/m_i)q_iδν_j; its endpoint w values increase by δν_j |
| s_{gi} / ν_{gi} | x_i and W_{gi} | Same storage formula with q_i=−x_s m_i g/E; both x_i and v_i move, and w_{gi} increases |
| Zero rod store / ν_R | W_R only | Only w_R increases; no direct state response |
| C / λ | Rod endpoint positions | Both position and velocity through the x/v cross block |
| T / τ | Rod endpoint positions and velocities | Both state blocks; no direct work coordinate |

The magnitude and final sign depend on the solved multiplier increment; the table is a derivative, not a claim that every correction goes in the positive listed direction.

Under (12), a position-only row has δu_i=2δξ_i and therefore preserves the current soft-midpoint residual b_i during its linearized push. A kinetic row generally changes that residual. Thus “μ is a velocity multiplier” does not imply that its update leaves position unchanged. Conversely a storage row has no direct velocity derivative, but its response generally moves velocity. The response support remains local under (12). Using the exact globally coupled inverse Hessian would generally spread that response beyond the row's endpoint support.

A scalar row correction at primal stationarity in this fixed model uses δα=−c/(∇cᵀA^{-1}∇c). Away from stationarity, the primal residual must also be included: δz=−A^{-1}r+A^{-1}Σ_a∇c_a δα_a and the linearized balances determine the dual increments. Reapplying the complete −A^{-1}r independently on every row would be a different iteration. This document does not select a GS scheduling rule.

## 3. Affine scalar spring

There is one scalar particle, one anchored spring, and no gravity. Write κ for stiffness to distinguish it from the kinetic row k:

\[
U(x)=\tfrac12\kappa(x-\ell)^2,\qquad
\nabla U=\kappa(x-\ell),\qquad \nabla^2U=\kappa,\qquad
a_n=-\kappa(x_n-\ell)/m.
\]

Unknowns are (x,v,W), objective (1), and the complete constraints are

\[
k=\frac{m(v^2-v_n^2)/2-W}E=0,\qquad
s=\frac{\kappa[(x-\ell)^2-(x_n-\ell)^2]/2+W}E=0.
\tag{13}
\]

In physical variable order (x,v,W),

\[
\nabla k=(0,mv,-1)^T/E,\quad \nabla s=(\kappa(x-\ell),0,1)^T/E,
\quad
\nabla^2k=\operatorname{diag}(0,m,0)/E,\quad
\nabla^2s=\operatorname{diag}(\kappa,0,0)/E.
\]

The objective Hessian is explicitly

\[
\nabla^2\Phi=
\begin{bmatrix}
5m/(Eh^2)+\kappa/(4E)&-2m/(Eh)&0\\
-2m/(Eh)&m/E&0\\0&0&1/E^2
\end{bmatrix}.
\]

The full Lagrangian is Φ−μk−νs. Its xx and vv entries are obtained by subtracting νκ/E and μm/E. Both rows share W; μ additionally responds through the local x/v metric, while ν additionally responds through the spring extension.

Here (9) is scalar: (m+h²κ/4)δx̂=hmv_n−h²κ(x_n−ℓ)/2, and Ŵ=−κ[(x̂−ℓ)+(x_n−ℓ)]δx̂/2. This equals −ΔU on the predicted chord. There is no absolute value of x: crossing the anchor does not turn this affine spring into a radial one.

Nominal scene: m=ℓ=1, κ=100; release x=1.25,v=0, or moving release v=0.8. H_min=0.

## 4. Planar radial spring

One planar particle is anchored by a compliant spring, with no gravity. Put r=||x||>0, n=x/r, P=I−nnᵀ, c=r−ℓ. Then

\[
U=\tfrac12\kappa c^2,\quad
\nabla U=\kappa c n,\quad
\nabla^2U=\kappa\left(nn^T+\frac c r P\right),\quad
a_n=-\kappa(r_n-\ell)n_n/m.
\tag{14}
\]

Unknowns are x∈R²,v∈R²,W∈R. Use (1) with this U, and

\[
k=[m(\|v\|^2-\|v_n\|^2)/2-W]/E=0,
\quad s=[\kappa(c^2-c_n^2)/2+W]/E=0,
\quad \mathcal L=\Phi-\mu k-\nu s.
\]

The row derivatives are (5) with ∇U and ∇²U from (14): k has only the vv block mI/E; s has only the xx block κ(nnᵀ+cP/r)/E. The spring Hessian can have negative tangential curvature in compression; do not replace it by κnnᵀ in an **exact** Hessian. That replacement is used only for the stated predictor/metric approximation.

Both rows touch one work coordinate. Through (12), μ moves both x and v along the current v direction; ν moves both along the current radial energy gradient. The potential and its classical derivatives are not differentiable at the collapsed edge r=0 for ℓ>0; that is a separate geometric failure.

Nominal orbit: m=ℓ=1, κ=100, x=(1.1,0),v=(0,0.8). Release uses x=(1.25,0),v=0. H_min=0. The target is (9)–(10), with B_n=κn_n n_nᵀ.

## 5. Rigid pendulum with gravity

One planar particle has a rod of length ℓ to the origin and g=(0,−g_0). Its stores are U_R=0 and G(x)=−mgᵀx. Unknowns are (x,v,W_R,W_g). The objective is (1) with U=G and both work penalties. The complete five constraints are

\[
\begin{aligned}
k&=[m(\|v\|^2-\|v_n\|^2)/2-W_R-W_g]/E=0,\\
s_R&=W_R/E=0,\qquad s_g=[-mg^T(x-x_n)+W_g]/E=0,\\
C&=(\|x\|-\ell)/x_s=0,\qquad T=n^Tv/v_s=0,\quad n=x/\|x\|.
\end{aligned}
\tag{15}
\]

Use \(\mathcal L=\Phi-\mu k-\nu_Rs_R-\nu_gs_g-\lambda C-\tau T\). No spring potential is assigned to the hard rod. Its explicit zero-storage equation sets its net work to zero; this is separate from enforcing endpoint length and tangency.

On an initially feasible state, the full old acceleration is

\[
a_n=g-\left(n_n^Tg+\frac{\|v_n\|^2}{\ell}\right)n_n.
\tag{16}
\]

It satisfies x_nᵀa_n+||v_n||²=0. For an off-tangent input, the hour helper replaces ||v_n||² by the squared tangential component when estimating the reaction; this is not a substitute for making the input physically feasible. Equation (16) is the physical acceleration on the intended manifold. Omitting the reaction changes the quarter predictor and the objective.

The balance derivatives follow (5): ∇s_R is only its +1/E work entry; ∇_xs_g=−mg/E with +1/E on W_g; both storage Hessians are zero. Put r=||x||, P=I−nnᵀ and b=nᵀv. The geometric derivatives in physical coordinates are

\[
\begin{aligned}
\nabla_xC&=n/x_s,&\nabla^2_{xx}C&=P/(x_s r),\\
\nabla_xT&=Pv/(v_s r),&\nabla_vT&=n/v_s,\\
\nabla^2_{xv}T&=P/(v_s r),&\nabla^2_{vv}T&=0,\\
\nabla^2_{xx}T&=\frac{-bI-vn^T-nv^T+3bnn^T}{v_s r^2}.
\end{aligned}
\tag{17}
\]

Apply (11) for normalized derivatives. In particular ∇_ξC=n, ∇²_ξξC=x_sP/r, ∇_uT=n and ∇²_ξuT=x_sP/r. These match the hour evaluator's distance and unit-normal constraints; replacing them by squared distance and xᵀv preserves the feasible set but changes row scales, multipliers, and off-manifold derivatives.

The start-only rod predictor uses θ=atan2(x_{n,x},−x_{n,y}), ω=(x_{n,x}v_{n,y}−x_{n,y}v_{n,x})/ℓ² and α=(g_x cosθ+g_y sinθ)/ℓ. Set θ̂=θ+hω+h²α/2, ω̂=ω+hα, then map to x̂=ℓ(sinθ̂,−cosθ̂), v̂=ℓω̂(cosθ̂,sinθ̂). The existing target averages old/predicted reaction strengths and multiplies the exact distance discrete direction. Its chord work Ŵ_R is zero in exact arithmetic because both predicted lengths equal ℓ; Ŵ_g=mgᵀ(x̂−x_n).

μ touches both work leaves and v; ν_R only the rod work; ν_g the gravity work and x; λ and τ have their separate state supports in (17). Despite ∇_xk=0, μ moves x via the metric. Because s_R fixes W_R=0 and Ŵ_R=0, (7) gives ν_R=μ at stationarity. This equality concerns energy bookkeeping multipliers, not the rod's physical reaction or λ.

Nominal state: m=ℓ=1, g_0=9.81, θ=0.65; release ω=0 or kick ω=1.2. H_min=−mg_0ℓ. Exact midpoint kinematics is not added: together with endpoint length/tangency it would impose extra restrictions on the endpoint speeds. Here that relation is penalized softly.

## 6. Elastic pendulum

Replace the rod by the radial spring (14), retaining gravity. Unknowns are (x,v,W_s,W_g), and

\[
U=\tfrac12\kappa(r-\ell)^2-mg^Tx,\qquad
\nabla U=\kappa(r-\ell)n-mg,\qquad
\nabla^2U=\kappa(nn^T+(r-\ell)P/r).
\]

Use (1) and the three constraints

\[
k=[\Delta K-W_s-W_g]/E=0,\quad
s_s=[\Delta U_s+W_s]/E=0,\quad
s_g=[\Delta G+W_g]/E=0,
\quad \mathcal L=\Phi-\mu k-\nu_s s_s-\nu_gs_g.
\tag{18}
\]

There are **no** hard length or tangency rows. Radial velocity is physical. The spring row has the radial gradient/Hessian (14); the gravity row has constant position gradient and zero Hessian; the kinetic row has its vv block and both −1/E work entries. Thus μ has one shared kinetic meaning across both incidences, while ν_s and ν_g are separate store multipliers.

The acceleration in (1) is a_n=g−κ(r_n−ℓ)n_n/m. The target uses (9), B_n=κn_n n_nᵀ, and separate spring/gravity work predictions. The nominal mixed scene has m=ℓ=1, κ=100, θ=0.45, r=1+9.81/100+0.05, v=−0.2 n+0.5 t with n=(sinθ,−cosθ), t=(cosθ,sinθ). The vertical release has θ=0, r=1+9.81/100+0.2, v=0. H_min=−mg_0ℓ−(mg_0)²/(2κ).

## 7. Two compliant springs, then a shared-vertex extension

For two moving planar particles, e_1=x_1 and e_2=x_2−x_1. Let U_j=κ_j(||e_j||−ℓ_j)²/2, G_i=−m_i gᵀx_i and U=U_1+U_2+G_1+G_2. The works are

\[
(W_{11},W_{21},W_{22},W_{g1},W_{g2}),
\]

where W_{21} is energy delivered by spring 2 to particle 1. The complete constraints are

\[
\begin{array}{ll}
k_1=(\Delta K_1-W_{11}-W_{21}-W_{g1})/E,&
k_2=(\Delta K_2-W_{22}-W_{g2})/E,\\
s_1=(\Delta U_1+W_{11})/E,&
s_2=(\Delta U_2+W_{21}+W_{22})/E,\\
s_{g1}=(\Delta G_1+W_{g1})/E,&
s_{g2}=(\Delta G_2+W_{g2})/E.
\end{array}
\tag{19}
\]

Each equals zero. Minimize (1) with these five W penalties, and use \(\mathcal L=\Phi-\mu_1k_1-\mu_2k_2-\nu_1s_1-\nu_2s_2-\nu_{g1}s_{g1}-\nu_{g2}s_{g2}\). There are 13 primal variables and six rows, without geometry/tangency constraints.

For each spring write

\[
\nabla_{e_j}U_j=\kappa_j(r_j-\ell_j)n_j,\qquad
\nabla^2_{e_j}U_j=\kappa_j[n_jn_j^T+(r_j-\ell_j)(I-n_jn_j^T)/r_j].
\]

Each derivative in the following blocks is with respect to its own edge vector e_j. Then

\[
\nabla_x U=\begin{bmatrix}\nabla U_1-\nabla U_2-m_1g\\\nabla U_2-m_2g\end{bmatrix},\qquad
\nabla^2_{xx}U=
\begin{bmatrix}\nabla^2U_1+\nabla^2U_2&-\nabla^2U_2\\-\nabla^2U_2&\nabla^2U_2\end{bmatrix}.
\tag{20}
\]

The spring-2 row's position gradient is (−∇U_2,+∇U_2)/E; its xx Hessian is [[∇²U_2,−∇²U_2],[−∇²U_2,∇²U_2]]/E. Spring 1 contributes ∇²U_1/E at x_1 only. Each kinetic row contributes m_i I/E at v_i only. Gravity rows are affine. These, together with (4), give the full objective and Lagrangian Hessians.

The old accelerations are a_{n,1}=g−∇U_1(x_n)/m_1+∇U_2(x_n)/m_1 and a_{n,2}=g−∇U_2(x_n)/m_2. The predicted spring force f̂_2 is applied positively to particle 2 and negatively to particle 1, giving separate Ŵ_{22} and Ŵ_{21}; their sum is −ΔU_2 at the predicted state.

| Multiplier | Direct row support | Response through the per-particle metric |
|---|---|---|
| μ_1 | v_1,W_{11},W_{21},W_{g1} | x_1,v_1 and those three W values |
| μ_2 | v_2,W_{22},W_{g2} | x_2,v_2 and those two W values |
| ν_1 | x_1,W_{11} | x_1,v_1,W_{11} |
| ν_2 | x_1,x_2,W_{21},W_{22} | Both endpoint positions and velocities, both spring-2 works |
| ν_{gi} | x_i,W_{gi} | x_i,v_i,W_{gi} |

There is one μ_1, shared when either adjacent spring is visited. Making a separate “μ for particle 1 on spring 1” and “μ for particle 1 on spring 2” adds variables and requires explicit consensus constraints; it is not the formulation above.

Nominal mixed scene has both m=ℓ=1, κ=100. Link angles are (0.45,−0.3), lengths (1+2g_0/100+0.05, 1+g_0/100+0.10), and relative velocities (−0.2 n_1+0.5 t_1, 0.3 n_2+0.8 t_2); positions and velocities are accumulated down the chain. The larger-release variant uses angles (0.65,−0.45), extra stretches (0.10,0.20), and zero velocity. The chain minimum is

\[
H_{\min}=-\sum_{j=1}^N\left(M_j g_0\ell_j+\frac{(M_jg_0)^2}{2\kappa_j}\right),\qquad M_j=\sum_{i=j}^Nm_i.
\]

### N=3 extension: a genuinely shared middle vertex

Add spring 3 between particles 2 and 3, with works W_{32},W_{33} and gravity leaf W_{g3}. The eight work variables are (W_{11},W_{21},W_{22},W_{32},W_{33},W_{g1},W_{g2},W_{g3}). Equations (19) extend to

\[
\begin{aligned}
k_1&=(\Delta K_1-W_{11}-W_{21}-W_{g1})/E,\\
k_2&=(\Delta K_2-W_{22}-W_{32}-W_{g2})/E,\\
k_3&=(\Delta K_3-W_{33}-W_{g3})/E,\\
s_1&=(\Delta U_1+W_{11})/E,\\
s_2&=(\Delta U_2+W_{21}+W_{22})/E,\\
s_3&=(\Delta U_3+W_{32}+W_{33})/E,
\qquad s_{gi}=(\Delta G_i+W_{gi})/E.
\end{aligned}
\tag{21}
\]

There are 20 primal variables and nine rows. Objective (1), Lagrangian (3), and derivative rules (4)–(5) apply without alteration. The potential Hessian is block tridiagonal, with diagonal blocks ∇²U_1+∇²U_2, ∇²U_2+∇²U_3, ∇²U_3 and off-diagonal blocks −∇²U_2,−∇²U_3. Accelerations follow the adjacent spring-force differences plus gravity.

The **single** μ_2 updates W_{22}, W_{32}, W_{g2} together and responds through both x_2 and v_2. ν_2 moves particles 1/2 and their spring-2 works; ν_3 moves particles 2/3 and their spring-3 works. Their visits overlap at particle 2 even though their W variables differ. Under (12), storage rows also couple through the shared position response; a shared-dual iteration must retain those effects.

For the actual N3 benchmark, use `hour_local::chain_case(3)` from the linked hour evaluator: nominal m_i=1, ℓ_i=1/N=1/3, κ_i=100N=300 and gravity (0,−9.81). For i=1,2,3 let M_i=Σ_{j=i}^3m_j, θ_i=0.35 sin(0.7i), n_i=(sinθ_i,−cosθ_i), t_i=(cosθ_i,sinθ_i), and

\[
r_i=\frac13+\frac{M_i9.81}{300}+\frac{0.03}{3},\quad
x_i=\sum_{j=1}^i r_j n_j,\qquad
v_i=\sum_{j=1}^i\frac{0.12}{3}t_j.
\]

This constructor also supports a stiffness multiplier and alternating mass contrast. These length/stiffness scalings differ deliberately from the original two-spring scene; N3 is a shared-vertex solver fixture, not a claim that only particle count changes. The mathematical specification itself supplies no new simulation result.

## 8. Redundancy, zero gradients, and limits of a local interpretation

Summing all rows k_i and s_j cancels every W and gives (H(x,v)−H(x_n,v_n))/E. Therefore **adding a separate total-energy row to the complete local balances is exactly redundant**. For disconnected dynamic energy graphs, the same statement holds per component; one global sum does not replace all component closures.

The existing local rows are not generically redundant merely because their W-only incidence columns have a null direction. Their x/v derivatives carry the component energy equation. Shifting every μ and ν by the same constant changes state stationarity by that constant times ∇H/E, so it is not a general multiplier gauge. Likewise, graph cycles create possible work circulations; the positive work penalty selects among them.

For the single anchored rod, s_R can be removed **only together with eliminating W_R=0** (and the corresponding constant penalty). Its geometric rows remain. Summing k,s_R,s_g with an additional energy row would again be redundant. For a rod between two moving endpoints, zero storage would impose W_left+W_right=0, not zero work separately at each endpoint; that mixed/general-graph case is not derived by the single-rod code.

At v_i=0, ∇_{v_i}k_i=0. The kinetic row still has nonzero work derivatives, so its full gradient is not zero. A μ-only response through (12) then changes W but cannot create velocity or position to first order. A nonzero spring/gravity storage gradient can still move x and, through the metric, v. Thus release at a turning point is not automatically a rank-deficient full constrained problem.

At a zero-energy equilibrium of the scalar/radial spring, both v=0 and ∇U=0. The combined energy gradient vanishes: the two balance gradients have opposite W entries and become dependent. Energy's second-order terms then matter. For a hanging or upright rigid rod at rest, ∇U is normal to the geometric manifold, so the energy row is dependent **after accounting for the position constraint**, even though its Cartesian gradient is nonzero. Elastic gravity equilibria similarly have cancellation between spring and gravity gradients.

A nonzero residual with zero gradient cannot be corrected by dividing by its first-order squared gradient. A zero residual at a stationary state may instead be physically correct. The objective, full KKT residual, row rank, second-order constrained curvature, and the old state's physical energy must distinguish these cases. Neither local dual updates nor preconditioning remove a true loss of first-order information. No claim of convergence or branch selection follows from the support formulas alone.
