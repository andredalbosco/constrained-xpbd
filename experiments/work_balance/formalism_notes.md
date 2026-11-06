# Work choices: identities, compatible integrators, and limits

Research notes accompanying the integrated work-balance experiments. All algebra below assumes constant positive particle masses and distinguishes physical timestep endpoints (-,+) from solver iterates (k,k+1). Pinned particles require external reaction accounting. Conservation identities do not imply accurate trajectories or convergence of the nonlinear solve.

## 1. Endpoint impulse work

Let `vbar=(v+ + v-)/2`, `Dv=v+ - v-`, and choose source impulses `I_a` (material, external, balance corrections, etc.). Define

    rho = M Dv - sum_a I_a
    W_a = vbar^T I_a.

Then the exact algebraic identity is

    DT = sum_a W_a + vbar^T rho.

Particlewise, `DT_i=sum_a W_ia+vbar_i dot rho_i`. Thus endpoint impulse work satisfies kinetic balance exactly only when the attributed impulses sum to the actual momentum change (or the residual happens to be orthogonal to vbar). This does not yet impose `sum_i W_ij=-DU_j`, momentum conservation, or angular momentum conservation.

For original XPBD kinematics, `v+=(x+-x-)/h`, predictor `xhat=x-+h v-+h^2 M^-1 fext`, and final-gradient material estimate `I_j=J_j(x+)^T lambda_j/h`,

    r_x = M(x+-xhat) - sum_j J_j(x+)^T lambda_j
    rho = r_x/h
    DT - Wext - sum_j W_j = vbar^T r_x/h.

If active balance corrections are present, either include their impulses as separate sources or include them in rho. Re-evaluating the final gradient times the accumulated multiplier is generally different from the sum of actual solver increments:

    J(x+)^T sum_k Dlambda_k != sum_k J(x^k)^T Dlambda_k.

Dropping r_x is precisely what prevents the final-gradient force estimate from necessarily reproducing the actual momentum change. Original XPBD explicitly makes this approximation; solving only constitutive residuals more accurately does not necessarily remove r_x for nonlinear constraints.

## 2. Exact kinetic accounting through sequential updates

For any accepted velocity update, define `DI_i=m_i(v_i^{k+1}-v_i^k)`. Then

    DW_i^k = (v_i^k+v_i^{k+1})/2 dot DI_i
           = v_i^k dot DI_i + ||DI_i||^2/(2 m_i)
           = T_i^{k+1}-T_i^k.

For a standard XPBD position update, `DI_i=m_i Dx_i/h=J_ji(x^k)^T Dlambda_j/h`. Include line-search scaling in the accepted impulse. Sum these works over *every* accepted update, including external predictor kicks and mu/nu corrections. The result telescopes to the total physical endpoint kinetic change, to arithmetic precision, independent of whether a nonlinear residual converges.

For a simultaneous paired update `DI=DI_material+DI_balance`, use the same before/after mean velocity for both source contributions; their sum is the exact update work. Treating the two parts as fictitious sequential subupdates changes attribution.

Attribution is order-dependent. For fixed impulses a then b on one mass,

    W_a = v0 dot a + ||a||^2/(2m)
    W_b = v0 dot b + a dot b/m + ||b||^2/(2m).

Reversing the order moves the cross term to the other source. Endpoint impulse work assigns half of that cross term to each source. The total agrees in either case; computed impulses can themselves depend on solver order.

This is an accounting identity along an algorithmic path, not generally a quadrature of the actual material force along the physical trajectory. In particular, the external predictor kick work is not automatically `fext dot (x+-x-)` after later position corrections. Discarding the quadratic impulse term gives only first-order correction work. Exact kinetic accounting by itself neither removes numerical material dissipation nor proves physical energy conservation.

## 3. Compatible physical-step work

Use midpoint drift and an impulse-consistent force:

    Dx = h vbar
    M Dv = h fbar.

Then `DT=fbar dot Dx` exactly. For `U=.5 C^T alpha^-1 C`, define

    Cbar = (C(x+)+C(x-))/2
    Jbar Dx = C(x+)-C(x-)
    D_U = Jbar^T alpha^-1 Cbar
    fbar_internal = -D_U.

The discrete chain rule gives `D_U dot Dx=DU`, hence `DT+DU=0` for an isolated solve. Constant gravity can be included through its linear potential or its exact work `f_g dot Dx`.

Linear momentum additionally requires `sum_i fbar_internal,i=0`. Angular momentum additionally requires `sum_i xbar_i cross fbar_internal,i=0`, where `xbar=(x++x-)/2`. A generic discrete gradient satisfies the energy chain rule but need not satisfy these symmetry conditions. They must be built into its construction.

Ordinary implicit midpoint uses `grad U(xbar)`, which obeys the exact chain rule for U quadratic in Cartesian positions. U quadratic in nonlinear constraint values C is not sufficient. Consequently midpoint alone is not a generic energy-conserving nonlinear XPBD replacement.

## 4. Midpoint position variable and compliance scaling

Let `y=(x++x-)/2`, `yhat=x-+(h/2)v-`, with a constant-force predictor contribution `(h^2/4)M^-1 fext` if desired. Then

    x+ = 2y-x-
    v+ = (4/h)(y-x-)-v-
    M(y-yhat) = (h^2/4) fbar_internal.

For ordinary midpoint forces, the incremental objective is

    .5 ||y-yhat||_M^2 + (h^2/4) U(y),

and the XPBD-scaled compliance is **4 alpha/h^2**. The `h^2/4` is the force/time factor, not the scaled compliance.

For the discrete-gradient force, one can instead introduce

    lambda_y = -(h^2/4) alpha^-1 Cbar
    M(y-yhat)-Jbar^T lambda_y=0
    Cbar+(4 alpha/h^2)lambda_y=0.

These are compatible mixed equations, but generally not the ordinary midpoint variational objective. With x- fixed, `d Cbar/dy=J(x+)`, whereas the force direction is Jbar. Therefore copying `Jbar M^-1 Jbar^T+4 alpha/h^2` as the exact Newton Schur block is incorrect in general. Differentiating the actual equations also produces derivatives of Jbar. An approximate local iteration can still be designed, but conservation rests on convergence of the actual equations, not on the name XPBD.

## 5. Radial spring example

Let `r=x_a-x_b`, `s=||r||`, `C=s-l0`, `k=1/alpha`. For `s++s->0`, choose

    gbar = (r++r-)/(s++s-)
    fbar_a = -k ((s++s-)/2-l0) gbar
    fbar_b = -fbar_a.

Then `gbar dot Dr=Ds`, so spring work is exactly `-DU`. Equal/opposite forces preserve linear momentum, and gbar is parallel to midpoint separation, so midpoint torque vanishes. This extends to radial potential phi(s) by replacing `k Cbar` with the secant `[phi(s+)-phi(s-)]/(s+-s-)`, using its derivative limit when radii coincide.

Both radii zero is a nonsmooth/undefined distance-gradient case. Antipodal endpoints give gbar=0; large timesteps can introduce branch/conditioning issues. Symmetry and conservation identities do not establish unique physical roots or accuracy.

## 6. Hard pendulum and velocity tangency

For a fixed anchor and `c(r)=.5(||r||^2-l^2)`, the discrete gradient is `rbar`. Use

    Dr=h vbar
    m Dv=h (m g+eta rbar)
    c(r+)=0, with c(r-)=0.

Reaction work is `eta rbar dot Dr=eta Dc=0`; with constant gravity, mechanical energy is conserved at converged equations. About the anchor, the reaction has zero midpoint torque, so angular impulse is that of gravity alone. Linear momentum is not conserved because of anchor reaction.

However endpoint tangency is not automatically preserved. Length and midpoint drift imply the exact identity

    r+ dot v+ + r- dot v- = (h/2) vbar dot Dv = h DT/(2m).

With initially tangent velocity and changing speed, the final velocity has a nonzero radial component, generally O(h^2). Energy plus endpoint length is therefore not the same as preserving the intrinsic tangent phase space. Adding an endpoint tangent projection changes energy unless the full scheme is redesigned.

If the compliant mixed row is instead extended directly to alpha=0, its equation is Cbar=0, implying Cnew=-Cold. This agrees with the endpoint constraint for a consistent initial state, but pre-existing position error can alternate. Enforcing Cnew=0 is a separate row choice and must be distinguished from simply setting compliance to zero.

## 7. Residual-to-invariant errors

For no hard constraint, define

    r_q = Dx-h vbar
    r_p = M Dv+h D_U-h fext
    delta_U = DU-D_U dot Dx.

Then

    DE-h vbar dot fext = vbar dot r_p + D_U dot r_q + delta_U.

For constant external force and external work defined as `fext dot Dx`, replace the right side by `vbar dot r_p+(D_U-fext) dot r_q+delta_U`.

With hard discrete constraints `Gbar Dx=Dc` and reaction `Gbar^T eta`, let `r_p=M Dv+h D_U-h fext-h Gbar^T eta`. Then

    DE-h vbar dot fext
      = eta dot Dc + vbar dot r_p
        +(D_U-Gbar^T eta) dot r_q + delta_U.

At exact drift and chain rule, energy error is `vbar dot r_p+eta dot Dc`. Large reaction multipliers can amplify small constraint residuals.

For total force fbar and `r_p=M Dv-h fbar`,

    DP = h sum_i fbar_i + sum_i r_p,i
    DL = h sum_i xbar_i cross fbar_i
         +sum_i xbar_i cross r_p,i
         +sum_i r_q,i cross (m_i vbar_i).

These identities identify the appropriate stopping diagnostics; small constitutive residual alone is insufficient.

## Primary-source audit

- [Macklin, Mueller, Chentanez: XPBD (2016)](https://mmacklin.com/xpbd.pdf), equations 12–17: geometric Hessian approximation and separate omission of the position/momentum residual; paper expressly distinguishes the resulting method from an exact solve of the implicit equations.
- [Gonzalez: Time Integration and Discrete Hamiltonian Systems (1996)](https://web.ma.utexas.edu/users/og/PUBLICATIONS/paper_DisHamSys.pdf), sections 3.2 and 3.6–3.7: discrete chain rule, symmetry-compatible derivatives, and central radial secant-force scheme in equation 3.43. The radial construction above is established energy-momentum methodology, not a new theorem.
- [Gonzalez: Mechanical systems subject to holonomic constraints (1999)](https://web.ma.utexas.edu/users/og/PUBLICATIONS/paper_Constrnts.pdf), equations 21–24: constrained discrete-gradient scheme and its conservation identity. Page 169 explicitly distinguishes preserved positions from potentially non-tangent endpoint velocities; page 172 exhibits velocity-constraint oscillations.
- [Simo, Tarnow, Wong (1992)](https://www.sciencedirect.com/science/article/pii/004578259290115Z): foundational energy-momentum conserving nonlinear dynamics algorithms. Publisher abstract verified; not relied upon for implementation details here.

## 8. GS target: actual residuals versus variational dualization

The DG midpoint equations can be solved with nonlinear GS and curvature dropping. Keeping the position residual does not require exact Newton. For one scalar material row, set

    r_x=M(q-yhat)-Jbar^T lambda
    r_lambda=Cbar+(4 alpha/h^2)lambda.

Using `D_q Cbar=J_end` and dropping derivatives of Jbar gives

    S=J_end M^-1 Jbar^T+4 alpha/h^2
    Dlambda=(-r_lambda+J_end M^-1 r_x)/S
    Dq=M^-1(Jbar^T Dlambda-r_x).

For many constraints this is a nonsymmetric reduced linearization; use block/coordinate GS with the actual couplings and evaluate the complete residual afterwards. Unlike ordinary Gram-form XPBD denominators, S need not be positive on arbitrary large-step configurations. Retaining r_x and dropping curvature are distinct choices. For the midpoint variable, accepted correction impulse is `DI=4M Dq/h`, not the original XPBD `M Dx/h`.

These residuals are generally not the direct KKT gradient of a scalar Lagrangian in `(q,lambda)`: cross derivatives would require `J_end=Jbar`. The endpoint DG force can have nonzero curl. This does not invalidate nonlinear GS for the equations, but one should not label the method a variational dual solve without a different, fully specified formulation.

### Implicit impulse work and degeneracy

Take independent source impulses `I_ij=h Jbar_ji^T ell_j` and `W_ij=I_ij dot vbar_i`. Midpoint drift implies

    sum_i W_ij=ell_j DC_j
    d_j=DU_j+sum_i W_ij=DC_j(Cbar_j/alpha_j+ell_j).

Thus the elastic work equation implies the material law only when `DC_j!=0`. At zero constraint change it loses the material law altogether; do not divide by DC. Retain a material equation and regard its work consequence as a diagnostic, or handle the redundant equations deliberately. With exact impulse balance, the kinetic work equation is already an identity; adding it as an independent constraint introduces redundancy (a kinetic reservoir instead remains unchanged).

Actual derivatives of implicit work are

    grad_I W=vbar,  grad_v W=I/2,  Hessian_(I,v) W=Identity/2.

If velocity is eliminated by impulse balance,

    W_ij=vold_i dot I_ij
         +I_ij dot (Iext_i+sum_k I_ik)/(2m_i),

which is quadratic in impulses and couples different material sources. Freezing its impulse or dropping its cross derivatives is another approximation. Exact kinetic/work curvature cancels when kinetic balance is an identity; retaining one half and discarding the other can create spurious stiffness.

### A genuine augmented scalar objective, not claimed equivalent

To illustrate what extra dualization actually means, an explicitly work-constrained optimization can use independent `(x,v,I,ell,S)` and a chosen positional XPBD objective `Phi(x)`:

    minimize Phi(x)+sum_j S_j^2/(2 gamma_j)
    subject to
      D=x-xold-h(v+vold)/2=0,
      P_i=m_i(v_i-vold_i)-Iext_i-sum_j I_ij=0,
      A_ij=I_ij-h Jbar_ji(xold,x)^T ell_j=0,
      d_j=U_j(x)-Uold_j+sum_i I_ij dot vbar_i+S_j-Sold_j=0.

Introduce respective multipliers `(a,b,eta,nu)`. The scalar Lagrangian is

    L=Phi+sum S^2/(2gamma)+a dot D+sum_i b_i dot P_i
      +sum_ij eta_ij dot A_ij+sum_j nu_j d_j.

In particular its actual stationarity includes

    grad_Iij L=-b_i+eta_ij+nu_j vbar_i,
    grad_vi L=-(h/2)a_i+m_i b_i+(1/2)sum_j nu_j I_ij,
    grad_Sj L=S_j/gamma_j+nu_j,
    grad_ellj L=-h sum_i Jbar_ji eta_ij,

plus the position equation, which includes derivatives of Jbar and of U. This can be linearized and solved using block GS with explicitly chosen Hessian approximations. It is a real optimization model, but **not proved equivalent to the simple DG residual solver or to ordinary XPBD**. Its force law, degeneracy at zero work, and any independent constitutive equation must be specified. Adding all constitutive and work rows without handling their dependencies can make KKT singular. These are formulation choices, not findings from the spring/pendulum experiments.
