# Paper-inspired directions for constrained XPBD

31 October 2026. Research notes for a later discussion. **The transfers proposed below are untested; they have not been added to the active benchmark.** The existing benchmark's DG, RATTLE, and energy-only FEPR-style controls are distinct from these proposed extensions. Sources were checked against author papers, author repositories, or original publication pages.

The useful separation is: **mechanical equations and objective; work law; variable/constraint scaling; linear solver; globalization; stopping criterion.** Improving the last four cannot make an incorrect mechanical objective correct. A frozen state can be feasible and stationary, so energy residual alone cannot certify useful motion.

## Suggested order

| Priority | Technique to explore | Failure it targets | First controlled experiment |
|---|---|---|---|
| 1 | Feasibility/optimality filter with a physical trust region | Oversized energy corrections; constraint reduction that worsens stationarity | Scalar release and elastic turning points; retain the same objective and work law |
| 2 | Mass-aware, graph-aware preconditioning | Slow propagation and ordering sensitivity in the spring chain | Same nonlinear problem, identity versus graph metric, equal evaluations |
| 3 | Augmented Lagrangian with explicit work variables | Local energy/work blocks competing or oscillating | Joint versus cyclic/symmetric updates with primal and dual residuals |
| 4 | Fixed-target projection with momentum structure | Energy-correct but mechanically implausible corrections | Extend the energy-only FEPR-style control using physically applicable momentum balances |
| 5 | Bounded spectral balancing inspired by Muon | Several nearly parallel constraint directions dominating others | Tiny systems with SVD diagnostics before any Newton–Schulz approximation |

PPO supplies an analogy for limiting exploitation of a stale local model; a physical trust region is the more direct first implementation. The DG control remains the test of whether consistent mechanics already solves the problem without an extra optimization technique.

## 1. Muon: separate update geometry from update scale

**Source mechanism.** Muon builds a momentum matrix and approximately replaces its singular values by order-one values using Newton–Schulz iterations. For ideal polar orthogonalization, if B=UΣVᵀ, the transformed update is UVᵀ. The practical quintic iteration is approximate; its official implementation explicitly does not claim an exact polar factor. [Original author explanation](https://kellerjordan.github.io/posts/muon/), [official implementation](https://github.com/KellerJordan/Muon/blob/master/muon.py).

**Proposed transfer.** Inspect the singular directions of the *nondimensionalized constraint derivative*, not an arbitrarily reshaped state update. For rows of G=∇cᵀ, try the bounded row transformation

\[
T=(GG^T+\varepsilon I)^{-1/2},\qquad
\widetilde c=T c,\quad \widetilde G=T G,
\]

holding T fixed during each local model and transforming residuals and multipliers consistently. Compare diagonal scaling, exact small-system SVD scaling, and no scaling before considering fast approximate matrix functions.

**What does not transfer.** x, v, and W are vectors with physical units, not neural-network weight matrices. Reshaping them changes the proposed polar update arbitrarily. Polar normalization would also erase physically meaningful force magnitudes. Whitening preserves rank: it cannot repair an exactly vanishing energy gradient. Across-timestep optimizer momentum would introduce additional temporal state and must not be confused with physical momentum.

**Composition:** numerical metric axis only; initially keep the objective, physical mass, work law, and hard constraints fixed.

## 2. PPO and TRPO: limit movement while a local model is trusted

**Source mechanisms.** PPO maximizes a clipped probability-ratio surrogate:

\[
L^{\rm clip}=\mathbb E\left[\min\left(r A,\operatorname{clip}(r,1-\epsilon,1+\epsilon)A\right)\right],
\qquad r=\pi_\theta/\pi_{\rm old}.
\]

The clipping removes an incentive for certain large policy changes; it is not a hard bound on every update. TRPO constrains policy KL divergence, using a local quadratic metric and a line search in its practical solver. [PPO, Eq. 7](https://arxiv.org/pdf/1707.06347), [TRPO, Sections 6–7](https://proceedings.mlr.press/v37/schulman15.pdf).

**Proposed transfer.** Freeze a startpoint work surrogate over an inner solve, then limit the correction while that surrogate is used:

\[
\frac{\|\Delta x/h\|_M^2+\|\Delta v\|_M^2}{E_{\rm scale}}
+\sum_{ij}(\Delta W_{ij}/s_{W,ij})^2\leq\Delta^2,
\]

where E_scale is a positive characteristic energy, s_W are positive work scales, and Δ is dimensionless. Compare actual versus predicted reduction in the chosen residual/merit model to shrink or expand Δ. Refresh the work model in an explicitly recorded outer iteration.

**What does not transfer.** Signed work can vanish or change sign, so W/W_old is not a valid analogue of a probability ratio. There is no policy distribution or KL metric here. RL improvement bounds do not prove stability or conservation of the mechanical integrator. Merely clipping the energy correction can leave the energy constraint unsatisfied; that residual must remain visible.

**Composition:** globalization and surrogate-refresh axes; compare first with all derivatives refreshed, then with a frozen work model.

## 3. Large Steps in Inverse Rendering: change coordinates, preserve the target

**Source mechanism.** Nicolet, Jacobson, and Jakob use differential coordinates u=P x with P=I+τL, then recover x by solving P x=u. Their author implementation separates this parameterization from the objective. [Paper and project](https://bnicolet.com/publications/Nicolet2021Large.html), [parameterization implementation](https://github.com/rgl-epfl/large-steps-pytorch/blob/master/largesteps/parameterize.py), [matrix construction](https://github.com/rgl-epfl/large-steps-pytorch/blob/master/largesteps/geometry.py).

For ordinary gradient descent in u, the chain rule gives the Cartesian update

\[
\Delta x=-\eta P^{-1}P^{-T}\nabla_x\Phi.
\]

This derivation distinguishes reparameterization from applying a single P⁻¹ smoothing pass.

**Proposed transfer.** On the spring graph, use a mass-whitened Laplacian or an SPD approximation of the elastic stiffness to spread local corrections. Compare diagonal mass scaling with a graph preconditioner, while evaluating the original objective and original constraints. For x,v,W, transform each derivative block consistently; edge-work variables need their own scale.

**What does not transfer.** Smooth geometric updates are useful in inverse rendering, but high-frequency vibrations may be valid mechanics. The metric must accelerate the inner solve without silently damping those modes. An invertible fixed coordinate change preserves stationary points, but a small fixed iteration budget can still change the simulated motion. Do not replace the physical mass matrix by the solver preconditioner.

**Composition:** coordinate/preconditioner axis; most informative on the chain and mass-ratio sweeps. Also applicable later to FEM element blocks: use deformation/stiffness structure, not SGD matrix reshaping.

## 4. FEPR: projection toward a fixed tentative state, with momentum structure

**Source mechanism.** FEPR projects an already computed tentative state using a quadratic metric containing

\[
\tfrac12\|x-x^*\|_M^2+\tfrac{h^2}{2}\|v-v^*\|_M^2.
\]

Its full formulation also uses momentum constraints and auxiliary scalars. The paper recommends h² for dimensional compatibility and explicitly does not guarantee improved ODE accuracy. [Dinev et al., 2018, Sections 3 and 5](https://users.cs.utah.edu/~ladislav/dinev18FEPR/dinev18FEPR.pdf). The local archived paper text was also inspected because the web PDF reader rejects its size.

**Proposed transfer.** Keep a fixed-target projection as an independent branch. Extend the current energy-only control with applicable momentum structure, then compare with the live coupling penalty

\[
\frac{1}{2\beta_v}\left\|v-\frac{x-x_n}{h}\right\|_M^2.
\]

Those two formulations have different cross derivatives and are not made equivalent by β_v=1. Evaluate the correction size, direction, phase, and actual impulse balance.

**What does not transfer.** Anchors and gravity exert external impulses and torques; global linear or angular momentum should not simply be held constant. The current energy-only implementation is **FEPR-style**, not a reproduction of the complete paper. Momentum feasibility deserves an explicit mechanism rather than blindly appending hard equalities.

**Composition:** objective/constraint branch, followed by a numerical metric sweep. Keep this separate from pure preconditioning of the original objective.

## 5. Discrete gradients: fix the work identity before tuning the optimizer

**Source mechanism.** A discrete gradient obeys

\[
\overline\nabla U(x_n,x)^T(x-x_n)=U(x)-U(x_n).
\]

Combined with midpoint kinematics and matching impulse closure, this makes the conservative kinetic and potential changes cancel. The framework and its higher-order extensions are discussed in [Gonzalez's author overview](https://web.ma.utexas.edu/users/og/numerics.html) and [Eidnes, Order theory for discrete gradient methods](https://link.springer.com/article/10.1007/s10543-022-00909-z).

**Proposed transfer.** Continue using the distance-spring factorization already derived:

\[
\overline\nabla C=\frac{x+x_n}{\|x\|+\|x_n\|},\qquad
\overline\nabla U=\frac{k(C+C_n)}{2}\overline\nabla C.
\]

This is the local anchored-edge formula; use relative endpoint positions for a moving edge. Its geometric factor makes rigid-rod chord work zero and avoids coupling stretch to rotation through an erroneous work approximation.

**What does not transfer.** Exact identities require a sufficiently accurate nonlinear solve. They do not prove the root solver converges for a large h or selects the physically useful branch. Force sampled at the midpoint configuration is generally different from a discrete gradient: rotating a spring can make its midpoint configuration appear compressed. Energy, phase, and endpoint tangency remain separate measurements.

**Composition:** mechanical closure/work-law axis. The current validated DG control provides the baseline before adding numerical accelerators.

## 6. ADMM and Projective Dynamics: explicit local variables and consensus

**Source mechanism.** Narain, Overby, and Brown split implicit elastic minimization into local constitutive problems and global compatibility using ADMM; Projective Dynamics is a special case of their construction. [Original paper, Section 3](https://www-users.cs.umn.edu/~narain/files/admm-pd.pdf).

**Proposed transfer.** Give a spring and each attached mass local copies of the work exchange, with consensus equality A z=b. For a genuine linear consensus split, use

\[
\mathcal L_\rho(z,y)=\Phi(z)+y^T(Az-b)
+\frac{\rho}{2}\|Az-b\|^2,
\qquad y^{k+1}=y^k+\rho(Az^{k+1}-b).
\]

Solve the constitutive/storage block, kinetic block, and consensus block separately. Monitor both primal consensus error and multiplier/stationarity error; a small energy defect alone is insufficient.

**What does not transfer.** Our energy and work equalities are nonlinear and the feasible set can be nonconvex. Calling an arbitrary sequence of local corrections “ADMM” does not inherit convex two-block guarantees. Increasing ρ indefinitely can recreate ill-conditioning; the actual split and subproblem solutions need to be specified.

**Composition:** variable splitting and solver schedule; first compare with the same converged joint solution, then examine finite budgets and ordering.

## 7. Filter-guided nonlinear multiblock optimization: measure feasibility and stationarity separately

**Source mechanism.** Baraldi, Leyffer, and Wild's 2024 preprint addresses nonlinear, nonconvex multiblock problems with augmented Lagrangians and a filter. It tracks constraint error η=||c|| and an optimality error ω. Its Eq. 12 accepts a candidate relative to each filter entry when

\[
\eta\leq\beta\eta^{(l)}
\quad\text{or}\quad
\omega\leq\omega^{(l)}-\gamma\eta.
\]

Sufficient decrease and restoration are additional parts of the algorithm. [Primary preprint, Sections 2–4](https://wiki.mcs.anl.gov/leyffer/images/0/00/FilterADMM.pdf).

**Proposed transfer.** Track scaled balance errors alongside the actual stationarity residual r_x (and r_v,r_W where present). A constraint-only improvement should not hide growing stationarity error. After failed local sweeps, invoke a declared restoration step or shorten the physical step; record which happened. Compare this with a single weighted least-squares merit function.

**What does not transfer.** Our mass constraints, turning points, and vanishing energy gradients can violate regularity assumptions. The cited convergence analysis is not a physical-time stability theorem. A filter may correctly converge to a frozen stationary point if the chosen objective admits it.

**Composition:** acceptance/restoration axis; compatible with joint SQP and properly defined block subproblems. This is a more directly relevant constrained-optimization transfer than copying PPO clipping literally.

## Proposed discussion agenda

1. Classify each observed failure: wrong work identity, wrong mechanical target, rank loss, or incomplete nonlinear solve.
2. Choose one fixed mathematical formulation and add a physical trust region/filter. Keep dropping curvature and dropping r_x as separate switches.
3. Compare identity, diagonal, and graph metrics at equal work. Check that converged answers agree before interpreting early-iteration improvements.
4. Consider explicit-work ADMM only once its local and joint problems have the same intended stationary conditions.
5. Leave Muon-style spectral normalization last; tiny-system SVD diagnostics should establish a useful matrix structure first.

For every proposal, retain release-from-rest, mixed rotation/stretch, mass imbalance, and chain ordering tests. Report energy, momentum/impulse closure, rod tangency, phase/amplitude, correction norm, and evaluation count together. No proposed numerical method should be credited with fixing a physical model merely because it reduces a selected residual.
