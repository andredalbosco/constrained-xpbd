# Constrained XPBD notation

Working reference for the current formulation discussion, established 2026-11-01.
Use these definitions in subsequent derivations. If a convention changes, update
this document and state the conversion explicitly. Historical experiment files
retain their own conventions; the conversion table below identifies the main
differences. This document defines notation, not a newly selected or validated solver.

Prefer explicit expressions to auxiliary substitutions. Write the velocity
difference as \(v-(x-x_n)/h\), the energy residual as
\(\tfrac12v^TMv+U(x)-E_\star\), and the target difference as
\(E_\star-K_\star\). The optional notation \(E(x,v)\) denotes total mechanical
energy, never an energy residual. Expand the implicit objective when presenting
a new formulation.

## State, indices, and derivatives

| Symbol | Meaning / convention |
|---|---|
| \(n\), \(h>0\) | Timestep index and timestep duration. |
| \(i\), \(j\) | Particle index and material/geometry constraint index. |
| \(x_n,v_n\) | Accepted position and velocity at the start of the timestep. |
| \(x,v\) | Current trial endpoint variables; accepted values become \(x_{n+1},v_{n+1}\). |
| \(M\), \(m_i\) | Mass matrix on free coordinates and particle mass. Pinned coordinates are eliminated. |
| \((x-x_n)/h\) | Velocity reconstructed from position; write this expression directly. |
| \(v-(x-x_n)/h\) | Velocity-coupling residual; write this expression directly. |
| \(x^{(r)},v^{(r)},\lambda^{(r)}\) | Solver iterate \(r\) within one timestep. Use \(r\), not spring stiffness \(k\), as the iteration index. |
| \(\Delta x,\Delta v,\Delta\lambda,\Delta\eta,\Delta\pi\) | Solver increments. The timestep displacement is written \(x-x_n\), not \(\Delta x\). |
| \(\lambda^-,\lambda^+\) | When used for a local visit: multiplier immediately before/after that visit; \(\lambda^+=\lambda^-+\Delta\lambda\). Previous-timestep values carry a timestep subscript. |
| \(\|z\|_M^2=z^TMz\) | Squared mass norm; \(\|z\|_{M^{-1}}^2=z^TM^{-1}z\). |
| \(\nabla_x f,\nabla_v f\) | Column gradients. State the differentiation variables; write \(\nabla f\) when only \(x\) is involved. |
| \(\nabla^2_{xx}f,\nabla^2_{xv}f\) | Hessian blocks. Use nabla notation, rather than renaming derivatives \(g\) or \(J\). |
| \(\nabla\mathbf c=[\nabla c_1,\ldots,\nabla c_m]\) | For a stacked constraint vector \(\mathbf c\), gradients are **columns**; its linearization is \(\mathbf c+(\nabla\mathbf c)^T\Delta q\) after defining the primal stack \(q\). |

## Objective, energies, and constraints

| Symbol | Definition / meaning |
|---|---|
| \(C_j(x)\) | Spring extension or position-constraint residual. Standard distance row: \(C_j=\|x_a-x_b\|-\ell_j\). |
| \(\ell_j\), \(k_j\) | Rest length and spring stiffness. For a compliant distance spring, \(k_j=1/\alpha_j\). |
| \(u_j(x)=C_j(x)^2/(2\alpha_j)\) | Elastic energy of compliant spring \(j\), with \(\alpha_j>0\). |
| \(U_s(x)=\sum_j u_j(x)\) | Total elastic potential. A hard rod has a constraint and a reaction multiplier, not a spring energy with division by zero. |
| \(U_{\mathrm{ext}}(x)\) | Conservative external potential. For constant gravity \(g\), \(U_{\mathrm{ext}}=-\sum_i m_i g^Tx_i\). |
| \(U(x)=U_s(x)+U_{\mathrm{ext}}(x)\) | Total potential energy. |
| \(F(x)=\|x-x_n-hv_n\|_M^2/(2h^2)+U(x)\) | Original implicit-Euler objective, in energy units. Conservative external forces are represented in \(U\). |
| \(\Phi(x,v)=\|x-x_n-hv_n\|_M^2/(2h^2)+U(x)+\|v-(x-x_n)/h\|_M^2/(2\beta_v)\) | Independent-velocity primal objective, defined for \(\beta_v>0\). |
| \(K(v)=\tfrac12v^TMv\) | Total kinetic energy. Particle contribution: \(K_i(v_i)=\tfrac12m_i v_i^Tv_i\). |
| \(E(x,v)=\tfrac12v^TMv+U(x)\) | Optional notation for total mechanical energy. Prefer the full expression in formulation equations. |
| \(E_\star\) | Prescribed energy target, held fixed during a timestep solve. State whether it is run-initial energy, incoming-step energy, or a target including external work. |
| \(\tfrac12v^TMv+U(x)-E_\star\) | Physical energy residual before relaxation; write it explicitly. With derived velocity it is \(\|x-x_n\|_M^2/(2h^2)+U(x)-E_\star\). |
| \(c(x)\) | A separately declared hard geometry residual, when needed. Example: \(c=(x^Tx-\ell^2)/2\) for a rod anchored at the origin. This is not identically the distance residual \(C=\|x\|-\ell\). |
| \(t(x,v)\) | A separately declared velocity-tangency residual. For that anchored rod, \(t=x^Tv\). |
| \(\mathcal L\) | Lagrangian/saddle objective. Distinct from schedule letter \(L\). |

All derivatives of physical energy use the analytical potential, even when the
objective's spring energy has been dualized or a PD projection has been frozen.
In particular,

\[
\nabla_x E(x,v)=\nabla U(x),\qquad
\nabla_v E(x,v)=Mv,
\]
\[
\nabla_x\left[\frac{\|x-x_n\|_M^2}{2h^2}+U(x)-E_\star\right]
=\nabla U(x)+\frac{M(x-x_n)}{h^2}.
\]

The derivative changes because the variables differ, not because the energy
function changed.

## Compliance and CPD relaxation

| Symbol | Meaning / convention |
|---|---|
| \(\alpha_j\ge0\) | Material compliance; zero denotes a hard constraint in the dual formulation. |
| \(\widetilde\alpha_j=\alpha_j/h^2\) | Timestep-scaled material compliance used in XPBD. |
| \(\beta_v\ge0\) | Dimensionless velocity compliance in the dual formulation. The primal penalty requires \(\beta_v>0\). It is unrelated to legacy work-reservoir \(\beta\). |
| \(\beta_v M^{-1}\) | Compliance matrix for \(v-(x-x_n)/h\). |
| \(\rho\) | CPD energy-relaxation variable. It is a primal unknown, not an energy multiplier; no \([0,1]\) bound is implicit. |
| \(\epsilon>0\) | Weight of \(\epsilon\rho^2/2\), with energy units. |
| \(K_\star\) | Fixed reference energy in the CPD interpolation. Original CPD derives it from prescribed momenta; it is not \(K(v)\). An anchored/gravity adaptation must define its reference explicitly. |
| \(\alpha_E\ge0\) | Energy compliance, with energy units. For \(\alpha_E>0\), the direct penalty is \([\tfrac12v^TMv+U(x)-E_\star]^2/(2\alpha_E)\). |
| \(\alpha_E=(E_\star-K_\star)^2/\epsilon\) | The particular energy compliance obtained by eliminating CPD's \(\rho\) when \(E_\star\ne K_\star\). If the targets coincide, this relaxation leaves the energy constraint hard. |

The CPD energy row is
\[
\frac12v^TMv+U(x)=(1-\rho)E_\star+\rho K_\star.
\]

This is the notation for the energy-relaxation component. Full original CPD
also has momentum constraints and its own external-work target convention.

## Multipliers: fixed signs and scaling

All Lagrangians below use the unscaled objective in energy units. The current
material convention, explicitly chosen by the user, is
\(\lambda_j=-C_j(x)/\alpha_j\) for a compliant spring. This is an unscaled
multiplier with the XPBD sign: its objective-gradient contribution is
\(-\lambda_j\nabla C_j\). The timestep-scaled XPBD multiplier is listed in
the conversion table. Introduce this coefficient by substitution
after the stated curvature approximation; no maximization notation is needed.

| Symbol | Role | Lagrangian contribution / defining equation |
|---|---|---|
| \(\lambda_j\) | Unscaled material multiplier, with the XPBD sign | \(\lambda_j=-C_j(x)/\alpha_j\), or \(C_j(x)+\alpha_j\lambda_j=0\). Its objective-gradient contribution is \(-\lambda_j\nabla C_j\). |
| \(\eta\) | Unscaled energy multiplier | \(+\eta[\tfrac12v^TMv+U(x)-E_\star]\) for hard energy; replace the target by \((1-\rho)E_\star+\rho K_\star\) with CPD relaxation. |
| \(\pi\) | Unscaled vector velocity-coupling multiplier | \(+\pi^T[v-(x-x_n)/h]-\tfrac12\beta_v\pi^TM^{-1}\pi\); \(v-(x-x_n)/h-\beta_vM^{-1}\pi=0\). |
| \(\mu_i\) | Unscaled kinetic-work multiplier, one shared value per particle | \(-\mu_i b_i\), with \(b_i\) defined below. |
| \(\nu_j\) | Unscaled potential-work multiplier, one per spring | \(-\nu_j d_j\), with \(d_j\) defined below. |
| \(\lambda_c,\lambda_t\) | Multipliers for separately declared \(c,t\) rows, if used | Use \(+\lambda_c c+\lambda_t t\); these named multipliers are unscaled. Do not reuse \(\mu,\nu\) for them. |

For a distance spring, \(\lambda_j\) has force units and its physical spring
force is \(+\lambda_j\nabla C_j\). The vector \(\pi\) has momentum units.
The multipliers \(\eta,\mu,\nu\) are dimensionless with these residuals.

After the Gauss-Newton spring approximation, substitute the complete linearized
spring coefficient, including its position increment:
\[
-\frac{C_j(x)+\nabla C_j(x)^T\Delta x}{\alpha_j}
=\lambda_j+\Delta\lambda_j.
\]
Use the constitutive-residual convention
\[
r_\lambda:=C(x)+\alpha\lambda,
\qquad r_{\lambda,j}:=C_j(x)+\alpha_j\lambda_j.
\]
Thus
\[
\nabla C_j(x)^T\Delta x+\alpha_j\Delta\lambda_j=-r_{\lambda,j}.
\]
Use \(r_\lambda\) in the single-spring system's right-hand side, alongside
\(r_x,r_\rho,r_\eta\).
If the current multiplier already equals \(-C_j(x)/\alpha_j\), the right-hand
side is zero. For independently iterated multipliers, retain that residual.
The analytical energy constraint and its gradient are evaluated from position;
the spring-coefficient substitution in the objective does not replace them.

If the energy slack has already been eliminated, use
\[
+\eta\left[\frac12v^TMv+U(x)-E_\star\right]-\frac{\alpha_E}{2}\eta^2,
\qquad r_\eta=\frac12v^TMv+U(x)-E_\star-\alpha_E\eta.
\]
Use either this term or explicit \(\rho\), not both for the same relaxation.

## Three velocity formulations

Energy, work, geometry, or momentum constraints are added explicitly to the
chosen formulation. They are not implied by its name.

| Formulation | Variables for velocity treatment | Objective / saddle expression | Status of \(\beta_v\) |
|---|---|---|---|
| Independent velocity, dualized coupling | Primal \(x,v\); dual \(\pi\) | Stationarity of \(\|x-x_n-hv_n\|_M^2/(2h^2)+U(x)+\pi^T[v-(x-x_n)/h]-\tfrac12\beta_v\pi^TM^{-1}\pi\) | \(\beta_v\ge0\); zero imposes \(v=(x-x_n)/h\). |
| Independent velocity, primal coupling | Primal \(x,v\) | \(\min_{x,v} \|x-x_n-hv_n\|_M^2/(2h^2)+U(x)+\|v-(x-x_n)/h\|_M^2/(2\beta_v)\) | Requires \(\beta_v>0\). |
| Position only, derived velocity | Primal \(x\) | \(\min_x \|x-x_n-hv_n\|_M^2/(2h^2)+U(x)\); define \(v_{n+1}=(x_{n+1}-x_n)/h\) | No \(\beta_v\), \(v\) unknown, or \(\pi\). |

The first two describe the same optimization problem for positive \(\beta_v\).
Their chosen iterative solvers need not produce the same finite-iteration state.
The third is obtained by eliminating the hard relation from the first at zero
compliance. In the independent-velocity formulations, the returned \(v\) remains
an independent variable; overwriting it with \((x-x_n)/h\) changes the method.

## Residuals and solver terminology

| Symbol / term | Meaning / convention |
|---|---|
| \(r_x=\nabla_x\mathcal L\), \(r_v=\nabla_v\mathcal L\) | Unscaled primal stationarity residuals. They are not energy or geometry residuals. |
| \(r_{\lambda,j}=C_j+\alpha_j\lambda_j\) | Material constitutive row; equals \(-\partial\mathcal L/\partial\lambda_j\) in the reference expression below. |
| \(r_\pi=v-(x-x_n)/h-\beta_vM^{-1}\pi\) | Velocity-coupling row. |
| \(r_\eta\) | Energy row: \(\tfrac12v^TMv+U(x)-E_\star\) for hard energy; \(\tfrac12v^TMv+U(x)-(1-\rho)E_\star-\rho K_\star\) with explicit slack; \(\tfrac12v^TMv+U(x)-E_\star-\alpha_E\eta\) after slack elimination. State which representation is active. |
| \(r_\rho=\epsilon\rho+\eta(E_\star-K_\star)\) | Stationarity of explicit CPD slack. |
| \(\nabla^2_{qq}\mathcal L\) | Exact primal Lagrangian Hessian, with \(q\) the explicitly declared primal stack; distinct from the full primal-dual system matrix. |
| \(A\), \(A_0\) | Declared approximation to the primal Hessian, and its position-only counterpart. State variables, scaling, and retained curvature. |
| \(S\) | Declared multiplier Schur matrix, not a physical stiffness or an energy. |
| \(p_j\) | PD local projection variable for element \(j\); its iterate is \(p_j^{(r)}\). |
| \(E\), \(L\) in a schedule | Energy-multiplier visit and material/length-multiplier visit. \(EL\) means energy then material; \(LE\) means the reverse. |
| \(\pi\)-visit | Velocity-coupling multiplier visit, if it is a separate block. State its position in the schedule, e.g. \(E\to L\to\pi\). |
| \(L_{\mathrm{ang}}\) | Angular momentum. Never use bare schedule letter \(L\) for it in the same derivation. |
| Accumulated multiplier increment | \(\lambda^{(r+1)}=\lambda^{(r)}+\Delta\lambda\), and similarly for other accumulated multipliers. Each local update immediately applies its specified primal correction. |
| Fresh multiplier solve | A newly solved multiplier value, explicitly labeled as such; it is not automatically an increment to add. |
| Frozen / refreshed gradients | Held fixed during the stated solve/sweep, or recomputed after the specified primal update. Declare which. |

Dropping multiplier-dependent curvature, dropping \(r_x\) or \(r_v\), and choosing
Gauss-Seidel instead of a joint solve are different operations. A schedule name
alone does not specify these choices, initialization, or the treatment of \(\pi\).

For a sign check, the dualized-velocity formulation with dualized spring energy
and explicit CPD slack is

\[
\begin{aligned}
\mathcal L={}&\frac{\|x-x_n-hv_n\|_M^2}{2h^2}+U_{\mathrm{ext}}(x)
-\sum_j\left(\lambda_j C_j(x)
+\frac{\alpha_j}{2}\lambda_j^2\right)\\
&+\pi^T\left(v-\frac{x-x_n}{h}\right)-\frac{\beta_v}{2}\pi^TM^{-1}\pi\\
&+\eta\left[\frac12v^TMv+U(x)-(1-\rho)E_\star-\rho K_\star\right]
+\frac{\epsilon}{2}\rho^2.
\end{aligned}
\]

Here the objective's \(U_s\) has been replaced by its dual representation;
the analytical \(U_s\) remains inside the physical energy constraint. Consequently,

\[
r_x=\frac{M(x-x_n-hv_n)}{h^2}+\nabla U_{\mathrm{ext}}
-\sum_j\lambda_j\nabla C_j
-\frac{\pi}{h}+\eta\nabla U,
\qquad r_v=\pi+\eta Mv.
\]

This formula fixes signs and units; it does not select this combination over the
other formulations. Extra geometry or work rows add their stated contributions.

## Work-balance notation

| Symbol | Meaning / convention |
|---|---|
| \(\Delta K_i\) | Timestep kinetic-energy change: \(\tfrac12m_i(\|v_i\|^2-\|v_{n,i}\|^2)\). |
| \(\Delta U_j\) | Timestep spring-energy change: \(u_j(x)-u_j(x_n)\). |
| \(W_{ij}\) | Work done by spring/constraint \(j\) on particle \(i\), positive when it supplies kinetic energy. Its approximation must be specified separately. |
| \(W_i^{\mathrm{ext}}\) | External work in the particle kinetic balance, including gravity if handled as work there. |
| \(b_i=\Delta K_i-\sum_jW_{ij}-W_i^{\mathrm{ext}}\) | Kinetic-work balance residual. |
| \(d_j=\Delta U_j+\sum_iW_{ij}\) | Spring potential-work balance residual. |
| \(-\sum_i\mu_i b_i-\sum_j\nu_j d_j\) | Fixed sign convention for these work constraints. |

These definitions do not choose free, endpoint, kinetic-exact, potential-exact,
or discrete work. A work formula and any auxiliary variables must be stated.
When total energy includes gravitational potential, do not also add gravitational
work to its total-energy target. The kinetic balance may still include that work.

## Conversions from earlier messages and existing files

| Earlier notation / convention | Translation here |
|---|---|
| Standard XPBD material multiplier used earlier and in existing code | Write it as \(\lambda_{\mathrm{XPBD}}=h^2\lambda\) when comparing with the current convention. Its row is \(C+(\alpha/h^2)\lambda_{\mathrm{XPBD}}=0\). |
| Earlier convention \(\lambda_s=C/\alpha\) | \(\lambda_s=-\lambda\). The current convention is \(\lambda=-C/\alpha\). |
| Scaled minus-sign energy multiplier called \(\eta\) in older EL experiments | Call it \(\eta_{\mathrm{old}}\) when comparing: \(\eta_{\mathrm{old}}=-h^2\eta\). Their position push has a plus sign. |
| Older position residual scaled by \(h^2\) | Equals \(h^2r_x\) when variables, forces, and other terms match. Never silently exchange scaled and unscaled residuals. |
| Work multipliers stored as \(q_\mu,q_\nu\) in the root implementation | \(q_\mu=h^2\mu\), \(q_\nu=h^2\nu\). Their signs follow the negative work terms above. |
| Root README's \(K_i,U_j\) for energy **changes** | \(\Delta K_i,\Delta U_j\) here. Existing files are historical; this table makes the change explicit. |
| \(H\) used for energy in older derivations | Use \(E(x,v)\) or the full mechanical-energy expression here. Write residuals by subtracting the target explicitly. |
| Abbreviations for velocity residuals or differences between energy targets | Expand them as \(v-(x-x_n)/h\) and \(E_\star-K_\star\). |
| Predictor \(x_n+hv_n+h^2g\), with gravity omitted from objective potential | Equivalent up to a constant to the inertial predictor \(x_n+hv_n\) plus \(U_{\mathrm{ext}}=-\sum_i m_i g^Tx_i\), for constant gravity. Do not include gravity in both places. |
| Original CPD's energy excluding gravity, with gravity supplied through an external-power target | Label that energy \(E_{\mathrm{CPD}}\) and its target \(E_{\star,\mathrm{CPD}}\). It is a different energy-accounting convention, not automatically the total mechanical energy and target defined above. |
| CPD paper's relaxation variable \(\alpha\) | \(\rho\) here; reserve \(\alpha_j\) for material compliance. |
| Bare \(\beta\), \(\gamma\) in root work-reservoir experiments | Kinetic/potential reservoir compliances, both in energy units. They are not \(\beta_v\). |
| \(\mu,\nu\) used for rod geometry/tangency in an earlier proposal | Use named \(\lambda_c,\lambda_t\); reserve \(\mu,\nu\) for work balance. |
| \(g,J\) used as derivative abbreviations in older files | Write the derivatives explicitly with nabla notation. Reserve \(g\) for gravitational acceleration. |

The relevant existing conventions can be inspected in the
[root README](/home/blacknoir/Projects/constrained-xpbd/README.md),
[energy-coupling derivation](/home/blacknoir/Projects/constrained-xpbd/experiments/energy_coupling/README.md),
[nonlinear-energy derivation](/home/blacknoir/Projects/constrained-xpbd/experiments/nonlinear_energy_2026_11_01/DERIVATION.md),
and [CPD specialization](/home/blacknoir/Projects/constrained-xpbd/experiments/nonlinear_energy_2026_11_01/CPD_DERIVATION.md).
