# Nonlinear energy visits and constrained projective dynamics

These experiments answer whether accurately solving energy cures the original E-before-L failure. They preserve the original endpoint velocity `v=(x-x_n)/h` until a method explicitly changes it (FEPR and the discrete-gradient controls). All simulations use a single fixed anchor; no moving support, contact, damping, or external work injection is present. Gravity is a conservative potential.

## Original objective and energy

Set `y=x_n+h v_n+h²g`. For the radial spring,

\[
F(x)=\frac{m}{2h^2}\|x-y\|^2+\frac k2(\|x\|-\ell)^2,
\quad
H(x)=\frac{m}{2h^2}\|x-x_n\|^2+\frac k2(\|x\|-\ell)^2-mg\cdot x.
\]

The energy constraint is `e(x)=H(x)-E_star=0`. For the rod, omit the spring energy and impose `||x||=l`. We use run-initial physical energy as the planar target. The scalar study also tests incoming-step targets separately.

\[
\nabla e=\frac m{h^2}(x-x_n-h^2g)+k(r-\ell)\frac xr,
\qquad
\nabla^2e=\frac m{h^2}I+k\left[(1-\ell/r)I+(\ell/r)nn^T\right].
\]

The kinetic derivative always includes the current reconstructed velocity. At the free predictor alone it reduces to `m v_n/h`. Numerically, `x-(x_n+h²g)` preserves an exactly zero free-fall gradient; subtracting two separately rounded accelerations need not.

## The original linear E visit and the replacement

The linear E visit solves

\[
\Delta\eta=-\frac{e(x)}{\nabla e^TM^{-1}\nabla e},\qquad
x\leftarrow x+M^{-1}\nabla e\,\Delta\eta.
\]

It ignores the finite-step energy curvature. The scalar release's second timestep can therefore request an amplitude-sized correction even for a vanishing timestep.

The replacement solves a nonlinear projection:

\[
x_E=\arg\min_{H(z)=E_\star}\frac12\|z-x\|_M^2.
\]

For the affine scalar spring, enumerate the two exact quadratic roots. For the planar radial spring, enumerate the exact energy contour and minimize its distance to the current iterate. For the rod, energy alone is a circle, so its nearest point is analytic. Equal-distance branches use the stated old-state preference; a zero-radius or missing contour is handled explicitly. The spring search uses a 32-point periodic scan plus scalar refinement of every sampled local minimum; comparison to 128 points and dense feasible candidates checks its accuracy, without claiming a formal global certificate.

This is a full nonlinear E projection. It is **not** a full constrained XPBD KKT solve. Its displacement is not generally `M^-1 nabla e(final) eta` for the accumulated historical multiplier. We do not report such a fictitious multiplier.

Next perform the ordinary fresh material row, retaining its multiplier within the timestep:

\[
\tilde\alpha=1/(kh^2),\quad
\Delta\lambda=\frac{-C(x_E)-\tilde\alpha\lambda}{1/m+\tilde\alpha},\quad
x\leftarrow x_E+\frac{n(x_E)}m\Delta\lambda,\quad
\lambda\leftarrow\lambda+\Delta\lambda.
\]

For the rod use `alpha_tilde=0`. Every row immediately pushes the primal variable. All timestep multipliers start at zero. Repeat E then L for the declared sweep count, reconstruct velocity once at the final position, and record physical energy both after E and after L. Later material visits can reopen the energy residual.

The raw planar comparison also refreshes the material normal after E. This coincides with the historical frozen-normal EL for the scalar counterexample; it is a **fresh-row variant**, not bitwise reproduction of the historical nonlinear radial/rod scheme. The separate C++ historical control executes the original implementation unchanged.

## Constrained PD global steps

Borrow PD's local projection and quadratic position objective:

\[
p^j=\ell\frac{x^j}{\|x^j\|},\qquad
Q_j(x)=\frac m{2h^2}\|x-y\|^2+\frac k2\|x-p^j\|^2.
\]

The ordinary global minimizer is

\[
z_j=\frac{(m/h^2)y+kp^j}{m/h^2+k}.
\]

The primary constrained variant minimizes `Q_j` subject to **true** physical energy `H(x)=E_star`. For one mass its Hessian is a scalar times the identity, so it is the nearest point on the exact energy contour to `z_j`. Then refresh `p` and repeat. For a hard pendulum, enforce both rod and energy exactly by enumerating their two circle intersections and comparing the inertial objective; this is the exact constrained single-rod control, not a large-stiffness spring approximation.

The frozen-energy control substitutes `k||x-p^j||²/2` into energy as well. Both objective and energy then have quadratic Hessian `(m/h²+k)I`; completing the square makes the global solve a circle projection. Refreshing the local projection changes that energy, by

\[
H_j(x)-H(x)=k\ell\|x\|[1-\cos\angle(x,p^j)]\ge0.
\]

Thus a finite-pass exact surrogate solve can lose physical energy. It is deliberately named `pd_surrogate` and its physical residual remains visible.

The `pd_inertia` control omits the spring term from the objective but keeps it in the energy constraint. Gravity remains in `y`. This tests using only part of the PD objective; it is a changed model, not a preconditioner.

These strict PD-inspired controls are distinct from Kee et al.'s published CPD relaxation formulation; see the separate CPD derivation/implementation if included. Source for the local/global construction: [Bouaziz et al. (2014)](https://tiantianliu.cn/papers/bouaziz14projective/bouaziz14projective.pdf).

## Why exact energy can still freeze

An exact algebraic identity gives

\[
F(x)=H(x)-m v_n\cdot\frac{x-x_n}{h}+\text{constant}.
\]

From rest, `F` is constant on the energy surface. The old position is feasible and the objective supplies no preference for correct acceleration. For the scalar extension `q`, the two roots at rest are

\[
q_{n+1}=q_n,\qquad q_{n+1}=\frac{1-kh^2/m}{1+kh^2/m}q_n.
\]

The moving root's reconstructed velocity has twice the correct leading acceleration increment. Picking a different root therefore also changes the consistency problem, rather than fixing it through numerical accuracy. The corresponding rod intersections show the same rest ambiguity. Adding endpoint tangency while retaining endpoint reconstruction forces zero motion, since on the rod

\[
x_{n+1}\cdot v_{n+1}=\frac{\|x_{n+1}-x_n\|^2}{2h}.
\]

FEPR instead projects independent position and velocity variables. The discrete-gradient diagnostic instead changes the kinematic rule to a midpoint relation and enforces compatible momentum/work. Their favorable energy residuals must be assessed together with motion, constraint residuals, and convergence failures.
