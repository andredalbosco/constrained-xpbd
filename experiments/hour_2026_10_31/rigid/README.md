# Rigid chains: the old reaction matters, and energy can move to the wrong mode

The main new finding is a failure that the single-pendulum experiments hid. The plain midpoint-potential objective is only first order on the tested double and triple rigid pendulums. In a 60-second high-mode experiment it transfers essentially all energy into other modes, while total energy remains accurate to about 7e-9. Including the old constraint reaction in the predictor restores the quarter objective and removes that transfer in this test.

The study also extends the single-rod hidden velocity mode to coupled constraints, and tests canonical-coordinate discrete gradients that keep position, tangency, energy, and, with a symmetry-aware correction, angular momentum compatible. These are small-system reference formulations and controlled experiments, not a production XPBD solver.

## The omitted-reaction experiment

Let the current physical acceleration be

\[
a_n=a_{\mathrm{free},n}+a_{\mathrm{reaction},n},
\qquad a_{\mathrm{free},n}=-M^{-1}\nabla U(x_n).
\]

The quarter objective from round two is

\[
\Phi_Q(x,v)=\frac{1}{2h^2}
\left\|x-x_n-hv_n-\frac{h^2}{4}a_n\right\|_M^2
+\frac14U(x)
+\frac{1}{2\beta}\left\|v+v_n-2\frac{x-x_n}{h}\right\|_M^2,
\]

subject to endpoint rod lengths, endpoint velocity tangency, and the initial step energy. For gravitational potential, which is linear in Cartesian positions, expanding the square gives exactly the same objective up to a constant as

\[
\boxed{\Phi_Q(x,v)=
\frac{1}{2h^2}\left\|x-x_n-hv_n-
\frac{h^2}{4}a_{\mathrm{reaction},n}\right\|_M^2
+U\!\left(\frac{x+x_n}{2}\right)
+\frac{1}{2\beta}\left\|v+v_n-2\frac{x-x_n}{h}\right\|_M^2
+\mathrm{constant}.}
\]

Thus a midpoint-potential objective whose position predictor is `x_n + h v_n` drops a specific old-reaction term. This is an algebraic identity for linear gravitational potential; it is not asserted for arbitrary nonlinear potentials. `reaction_equivalence.csv` verifies zero gradient and Hessian discrepancy in the evaluated control. Position and velocity remain independent in all these objective variants.

The measured position RMS errors over one second show the resulting order difference. References pass paired-resolution checks.

| Chain / objective | h=1/30 | h=1/60 | h=1/120 | h=1/240 |
|---|---:|---:|---:|---:|
| Double: midpoint potential, omitted old reaction | 0.029507 | 0.014206 | 0.006972 | 0.003454 |
| Double: quarter, full old reaction | 0.002138 | 0.000536 | 0.000135 | 0.0000338 |
| Triple: midpoint potential, omitted old reaction | 0.029575 | 0.013708 | 0.006611 | 0.003249 |
| Triple: quarter, full old reaction | 0.003523 | 0.000887 | 0.000224 | 0.0000564 |

The quarter variant converges at approximately second order in this grid. The plain midpoint-potential variant is approximately first order. Its KKT stationarity and energy residuals are still small, so additional iterations do not repair the mechanical target.

### Total energy hides nearly complete mode transfer

For three equal unit rods/masses, form the linearized angular mass and gravity stiffness matrices, mass-normalize their generalized eigenvectors, and initialize the highest mode at amplitude 0.02. The focused transfer test adds a lower-mode seed of 2e-7 to make sensitivity reproducible. Integrate for 60 seconds at 30 Hz. The modal-energy diagnostic is the quadratic energy of these fixed linear modes divided by the initial highest-mode quadratic energy, `E_j_linear(t)/E_high_linear(0)`. Its denominator is not instantaneous total modal energy or exact nonlinear energy. The reported conserved total energy is the separate exact nonlinear mechanical energy.

| Method | Largest fraction in other modes | Smallest retained fraction in initial mode | Max relative total-energy error |
|---|---:|---:|---:|
| Refined nonlinear reference | 2.56e-7 | 0.998489 | 8.64e-12 |
| Midpoint potential, no old reaction | 0.999995 | 1.62e-11 | 7.03e-9 |
| Quarter objective, old reaction | 2.8e-7 | 0.998489 | 3.58e-10 |
| Canonical discrete gradient | 2.3e-7 | 0.998489 | 6.33e-11 |
| Symmetric energy projection | 2.6e-7 | 0.998489 | 3.86e-11 |

The reference's small variation in the highest linear-mode energy is physical nonlinear motion. The omitted-reaction result is many orders larger. Its energy migrates into lower modes; it is not simply losing total energy.

The reaction-weight/compliance grid retains 0, 0.5, 1, 1.5, or 2 times the old reaction, with beta 0.1, 1, or 10. Omitting the reaction produces almost complete transfer for every beta. Half the reaction still produces roughly 99% transfer. Full reaction suppresses transfer to the reference scale. Overweighting the reaction also suppresses this particular high-to-low transfer, but a separate 120-run order grid shows that it loses second-order accuracy:

| Old-reaction fraction, beta=1 | 0 | 0.5 | 1 | 1.5 | 2 |
|---|---:|---:|---:|---:|---:|
| Double: observed order, h=1/120 to 1/240 | 1.013 | 1.004 | 1.994 | 1.009 | 0.998 |
| Triple: observed order, h=1/120 to 1/240 | 1.025 | 1.013 | 1.990 | 1.005 | 0.989 |

This distinguishes the derived reaction coefficient from tuning a penalty or choosing a coefficient that looks stable in one modal plot.

![Reaction, convergence order, and mode transfer](figures/reaction_and_mode_transfer.png)

![Reaction weight and beta](figures/reaction_weight_beta.png)

## Rigid-chain formulations compared

Every mass is a point at a rod endpoint; the first rod is anchored at the origin. All rods have unit length. Absolute link angles use

\[
x_i=\sum_{j\le i}\ell_j(\sin\theta_j,-\cos\theta_j),
\qquad v_i=\sum_{j\le i}\ell_j\omega_j(\cos\theta_j,\sin\theta_j).
\]

The eleven primary initial conditions include small and moderate motion in two/three links, a large base kick, high opposed angular speeds, a near-upright triple chain, a free chain without gravity, and masses (1,100), (1,0.01), or (1,20,0.05). Initial conditions are explicit in `rigid.hpp`. `double_kick` is a large initial kick, not a claim that a multi-link separatrix has been located analytically.

| Family | Position / endpoint tangency | Energy | Other property and limitation |
|---|---|---|---|
| RATTLE | Both enforced | Approximate | Time symmetric; free-chain angular momentum preserved; finite position solve can fail |
| Cartesian midpoint DG | Position enforced; tangency not enforced | Exact in converged equations | Free-chain angular momentum preserved; inherits a normal-velocity reflection mode |
| Symmetric projection of acceleration-level midpoint | Both enforced | Optional invariant row | Same multiplier before/after; explicit energy and optional angular-momentum rows |
| Ordinary postprojection | Both enforced | Invariant row | Separate control; reversibility is generally lost |
| Independent x,v quarter/midpoint/residual objectives | Both enforced | Optional hard energy | Soft midpoint kinematics; converged KKT does not establish correct dynamics |
| Canonical midpoint | Both by coordinates | Approximate | Symplectic; free-chain total angular momentum preserved |
| Canonical DG / corrected AVF | Both by coordinates | Exact in converged equations | Symmetric; requires a symmetry-aware correction to preserve free-chain angular momentum |

The symmetric projection uses fixed dimensionless Cartesian coordinates `(x/l_scale, v/v_scale)` and adds the same invariant-normal correction before and after an implicit midpoint step. The base ODE enforces acceleration-level constraints away from the manifold using the current constraint mass matrix. This follows the construction in [Hairer's symmetric projection paper](https://www.unige.ch/~hairer/preprints/symproj.pdf); the chosen metric, simultaneous energy constraint, and finite nonlinear solver are experimental choices here. [Andersen's RATTLE paper](https://doi.org/10.1016/0021-9991(83)90014-1) supplies the position/velocity constraint baseline.

The `rattle_energy` projection uses the fixed dimensionless metric. The `obj_projection_verlet_E` control instead uses the previous objective module's mass metric, with position displacement divided by h and beta=1. They are different objectives, not aliases. The `EL` suffix in objective variant names means energy plus angular momentum, and is active only without gravity; it is unrelated to the historical EL solver ordering.

### Matched trajectory measurements

Position RMS in metres, 30 Hz for ten seconds, all rows completed and references certified:

| Method | Moderate double | Moderate triple |
|---|---:|---:|
| RATTLE | 0.021158 | 0.024558 |
| Cartesian DG | 0.019212 | 0.024997 |
| Symmetric energy projection | 0.008683 | 0.014921 |
| Quarter objective | 0.020235 | 0.023049 |
| Mechanical DG residual + energy | 0.018087 | 0.023345 |
| Canonical midpoint | 0.004161 | 0.009374 |
| Canonical Gonzalez DG | 0.002462 | 0.006375 |

Over 60 seconds at 60 Hz, the moderate-double RMS values are 0.03191 for RATTLE, 0.02910 for Cartesian DG, 0.01333 for symmetric energy projection, and 0.003624 for canonical DG. The corresponding triple values are 0.02790, 0.03959, 0.02942, and 0.01369. High-speed and mass-contrast cases prevent a universal ranking: large phase errors or nonlinear failures remain visible in the CSVs. Concurrent wall time is not used as a performance ranking.

Even the successful methods have finite phase error. For the highest triple mode at 30 Hz over 60 seconds, Cartesian DG, symmetric energy projection, canonical DG, and the quarter objective lag the reference by about 2.64–2.66 radians while retaining mode amplitude. RATTLE leads by about 1.36 radians and has a roughly 1.7% total-energy excursion. Exact energy is not exact frequency.

## Canonical work that respects the rotational symmetry

In absolute angles, the mass matrix is configuration dependent:

\[
A_{ij}(\theta)=\ell_i\ell_j\!
\left(\sum_{k\ge\max(i,j)}m_k\right)\cos(\theta_i-\theta_j),
\quad p=A(\theta)\omega,
\quad H(\theta,p)=\tfrac12p^TA(\theta)^{-1}p+U(\theta).
\]

The reference solves the angle/velocity Euler–Lagrange equations. The experimental canonical methods independently solve Hamilton's equations in `(theta,p)`. With constant scale `p_scale`, define `z=(theta,p/p_scale)` and scaled Hamiltonian `H/E_scale`. The Gonzalez correction makes the gradient satisfy the exact discrete chain rule. This is the construction of [Gonzalez, Time Integration and Discrete Hamiltonian Systems](https://web.ma.utexas.edu/users/og/PUBLICATIONS/paper_DisHamSys.pdf).

For zero gravity, uniform angle rotation leaves H unchanged and total angular momentum is `sum(p)`. An arbitrary exact-energy correction can spoil this invariance. We therefore project the correction direction orthogonal to uniform angle rotation. If `d=z_next-z_old`, `P` removes the mean of d's angle components and leaves its momentum components intact, and `g` is the midpoint gradient, use

\[
\bar\nabla H=g+
\frac{H(z_{n+1})-H(z_n)-g^Td}{d^TPd}\,Pd.
\]

It satisfies the discrete chain rule and keeps the sum of angle-gradient entries zero. Coupling it to the canonical skew matrix closes both energy and free-chain angular momentum. The three-point AVF variant substitutes a three-point Gauss average of the gradient for g, then applies the same exact-work correction. Quadrature alone is not claimed exact for this non-polynomial Hamiltonian.

For the free double chain at 30 Hz over ten seconds:

| Method | Position RMS | Relative energy error | Scaled angular-momentum drift |
|---|---:|---:|---:|
| Canonical midpoint | 0.003683 | 0.00504 | 1.81e-14 |
| Ordinary Gonzalez DG | 0.010285 | 2.31e-10 | 7.18e-4 |
| Symmetry-aware Gonzalez DG | 0.012608 | 1.58e-10 | 8.93e-15 |
| Symmetry-aware corrected AVF | 0.004358 | 7.49e-11 | 8.54e-15 |
| Symmetric projection, energy only | 0.013704 | 7.71e-11 | 0.00360 |
| Symmetric projection, energy + momentum | 0.015607 | 6.76e-11 | 5.28e-12 |

Restoring momentum alone does not guarantee improved phase. The free triple chain also shows this tradeoff. The momentum suite includes h=1/120 and 1/30, retaining convergence evidence beside the invariant checks.

The ordinary Gonzalez correction depends on the chosen coordinate metric. Changing `p_scale` by factors 0.1 and 10 changes the moderate-double RMS from 0.00246 to 0.00543 and 0.00553, while all variants conserve energy. Corrected AVF is much less sensitive in that case (about 0.006019 for all three scales), because its remaining exact-work correction is smaller. This is a formulation change, not just solver preconditioning.

### A cancellation failure near coherent rotation, and a numerical repair

We also tested two/three-link free chains rotating coherently at angular speed 4, with initial angle differences of 0, 1e-8, 1e-6, 1e-4, or 1e-2, at 30 Hz and 10 Hz. Exact coherent rotation is a relative equilibrium: every angle increases linearly while the shape and canonical momenta remain fixed. It is a particularly useful case where the physical solution is known directly.

The symmetry-aware work correction divides by a projected displacement that tends to zero for this motion. The straightforward implementation subtracts two close Hamiltonian values and computes `d dot P d`; these operations can lose the very small internal-motion signal. Twelve of the forty symmetry-aware canonical DG/AVF trajectories failed their nonlinear solve. Ordinary canonical midpoint/DG remained accurate on exact coherent rotation. Symmetric energy-plus-momentum projection had three separate near-degenerate failures among its twenty trajectories.

The stable implementation retains the same mathematical method while changing three evaluations. With `p_i=A_i w_i`, it uses the identity

\[
\Delta K=\tfrac12(w_1+w_0)^T\Delta p-
\tfrac12w_1^T\Delta A\,w_0.
\]

Trigonometric difference identities evaluate `Delta A` and `Delta U` without subtracting close cosines or energies. The projected denominator is computed as `||P d||²`, and the symmetry-orthogonal gradient is dotted with `P d` to remove cancellation from the irrelevant uniform-rotation direction. No regularizing energy term or altered physical mass is added.

All forty stable follow-up trajectories complete, with maximum relative energy error 9.0e-13 and scaled angular-momentum drift 3.0e-14. On exact three-link coherent rotation, the position RMS is about 1e-10 m at 30 Hz and 3e-11 m at 10 Hz over ten seconds. The ordinary Cartesian DG control has 0.650 m RMS at 10 Hz despite preserving energy and momentum; RATTLE has 1.399 m RMS there. This is a coordinate-dependent finite-step phase difference, not a convergence defect.

The original failures remain in `equilibria_results`; the repair is in `equilibria_stable_results`, with broad free-chain follow-ups in `momentum_stable_results`. When the projected denominator is at or below 1e-22, the implementation skips the correction and falls back to the midpoint or three-point average gradient. That guarded branch no longer algebraically enforces the exact discrete chain rule. The observed accuracy is a finite-precision result, not a claim of uniform conditioning or exact energy enforcement at every degenerate state. One of twelve broader stable follow-ups still fails (`double_free_contrast`, stable canonical DG, 30 Hz, after 60 completed steps); its failed prefix is retained.

## The coupled hidden normal mode and adaptive stepping

Let the columns of `nabla C` be constraint gradients. The mass-orthogonal normal projector is

\[
P_N=M^{-1}\nabla C\left(\nabla C^TM^{-1}\nabla C\right)^{-1}\nabla C^T.
\]

Starting from exact positions with an off-tangent velocity, midpoint kinematics and endpoint position constraints force, as h tends to zero,

\[
v_{n+1}\longrightarrow v_n-2P_Nv_n.
\]

This is a reflection in the mass metric, preserving kinetic energy. Two half steps reflect twice, while one full step reflects once. Their discrepancy tends to twice the inherited normal component, rather than zero. This argument applies to multiple coupled rods.

`hidden_modes.csv` starts both from a genuine previous Cartesian-DG step and from an injected normal velocity. For the moderate double pendulum after one h=0.05 step, the predicted discrepancy limit is 0.0029699476 in the fixed scaled state norm. The measured discrepancy at h=1e-4 is 0.0029699454, with energy error 5.1e-13. At h=1e-6 it remains 0.0029698289. A separate full x,v,impulse solve verifies the reflection; it is not inferred from position-differencing alone.

All eleven ten-second Cartesian-DG step-doubling controllers hit the finite depth limit. By contrast, 30 of 33 feasible-state controllers (RATTLE, symmetric energy projection, canonical DG) finish. The common exceptional case is a mass-100 heavy tip. The controller compares one full step against two half steps in the fixed dimensionless state norm, accepts the two-half-step state, and requires discrepancy <1e-4|h|. No extrapolation is used.

The focused heavy-tip follow-up separates structural failure from accuracy-budget and nonlinear-tolerance limits. Increasing the subdivision cap from 12 to 16 lets symmetric projection and canonical DG complete one second with position RMS about 6.7e-7 and 7.3e-7 m, at 31,626 and 26,902 accepted substeps. Loosening the accuracy target to 1e-3|h| lets all three methods finish with roughly 7–9e-6 m RMS. Merely tightening Newton residuals does not fix the depth-12 limit. All controller attempts and failed prefixes are retained, including tolerance, residual, and rejection reason. State-dependent subdivision is not claimed time reversible.

![Invariant composition, hidden mode, and accuracy cost](figures/invariants_hidden_modes_and_cost.png)

## Validation, references, and scope

- The independently derived angular acceleration agrees with the Cartesian acceleration-level constraint solve to relative error 3.2e-15 on the eleven cases.
- Hamiltonian gradients and invariant gradients are checked by central differences. The worst scaled errors are 1.9e-9 and 3.3e-10.
- Eight primary integration/projection methods pass small-step geometry, tangency where intended, energy where intended, and reversibility checks. Symmetry-aware momentum methods have separate energy/momentum/reversal tests.
- Eighty-five derivative columns validate the reused independent-objective gradients, Hessians, and constraint derivatives on double and triple chains; worst scaled discrepancy is 1.55e-10.
- Stable Hamiltonian differences are checked against direct differences away from cancellation, and stable symmetry-aware methods have separate energy/momentum/reversal checks. A 176-run GCC/Clang comparison of the two one-second order suites changes no status and has maximum position-RMS difference below 4e-16 m.
- Reference trajectories use absolute-angle RK4, normally at least 3840 and 7680 substeps per simulated second. A trajectory is certified when the paired scaled-state difference is below 1e-5 and reference relative energy drift below 1e-7. Heavy-tip and near-upright triple ten-second references fail this check and are excluded from trajectory rankings. Some stress references also fail; their flags remain in the tables.
- The 60-second small/moderate references pass. Modal studies use two RK4 resolutions; the focused modal-transfer reference gap is 2.26e-8, and its energy drift is 8.64e-12. The focused one-second heavy-tip controller reference is separately refined (gap 1.62e-8).
- Relative energy errors divide by fixed initial excitation energy above the hanging minimum, floored at 1e-4; they do not divide by signed total energy. Scaled state norms use mean rod length and `sqrt(E_scale / mean mass)` as fixed position and velocity scales. Reported position RMS is in metres. Angular-momentum change is an invariant diagnostic only in gravity-free cases.
- Every failed run retains its completed-step count and prefix. Its prefix RMS is not compared as if it were a completed trajectory. Runs with zero completed steps report `nan` for unavailable position/velocity RMS, not zero error. Nonlinear methods are gated on convergence; finite output alone is not accepted.
- Wall times were collected concurrently. Residual-evaluation counts are workload diagnostics, not equal-cost comparisons between automatic differentiation, finite-difference Newton, and explicit reference methods.

The nine original trajectory suites contain 736 configurations (648 complete, 88 failed); `counts.json` gives per-suite counts and reference flags. They include repeated controls and aliases. Additional experiments comprise 120 reaction-weight/order trajectories, 160 relative-equilibrium trajectories, 40 stable near-equilibrium follow-ups, 12 stable broad-motion follow-ups, 90 modal trajectories, 19 reaction-weight/compliance transfer trajectories, 18 controller-limit probes, 20 hidden-mode probes, and derivative validation. The 176 compiler-comparison reruns are validation repeats. The files are independent additions; earlier experiment directories remain unchanged.

This study is limited to two/three-link, anchored, planar, bilateral rigid chains. It does not test collisions, friction, moving anchors, damping, or simultaneous stiff compliant and rigid elements. Generalized coordinates provide useful reference constructions here, but do not by themselves supply a local XPBD implementation for arbitrary constraint graphs.

## Exported motion comparison

[Watch the 60-second MP4](animation/modal_transfer_60s.mp4), or inspect the [representative-frame montage](animation/representative_frames.png). The four panels show the refined nonlinear reference, omitted-reaction midpoint potential, reaction-corrected quarter objective, and canonical discrete gradient. They use the saved 30 Hz positions, one fixed viewport, and normal playback speed.

The video magnifies each rod angle by **16 times**, reconstructing unit-length rods solely for display. The heading labels this transformation throughout. Mode-energy bars and total-energy residuals come from the unmodified simulation. The initial excitation is the highest linear mode with a small lower-mode seed; the video makes the erroneous transfer into the low mode visible while total energy stays nearly constant. This is an exported scientific animation and does not depend on the HTML viewer or browser playback validation.

`render_animation.py` generates the MP4 and four representative PNG frames from the CSVs using Pillow and ffmpeg. `animation/manifest.json` records input hashes, the video hash, codec/frame count/duration, and the result of decoding the entire video. The representative times are 1/30, 20, 40, and 60 seconds; the saved traces begin after their first step.

## Reproduction

Configure this directory with CMake. `COMPOSITION_ROOT`, `ROUND2_ROOT`, and `EIGEN_ROOT` are cache paths for the previous study headers and Eigen. Their standalone defaults point to the canonical constrained-xpbd repository: `experiments/composition_2026_10_31`, `experiments/round2_2026_10_31`, and `external/eigen`; override the cache paths when running from staging or another checkout. Build targets `rigid_study`, `rigid_probes`, `rigid_modes`, `rigid_modal_transfer`, and `rigid_controller_probe`; run CTest. No GPU is required.

Run the main executable as `rigid_study MODE OUTPUT_DIRECTORY`, with modes `short`, `order`, `accuracy`, `stress`, `long`, `momentum`, `metric`, `objectives`, `objectives_order`, `reaction_order`, `equilibria`, `equilibria_stable`, and `momentum_stable`. Run the other executables from this directory to write their CSVs. The `momentum_long` mode exists but is not included in the saved principal counts unless explicitly run.

`analyze.py` creates the PNG/PDF figures and `counts.json`. It needs NumPy and Matplotlib. The workstation used `/home/blacknoir/Projects/mandelbrot-diffusion/.venv/bin/python`; set `MPLCONFIGDIR` to a writable directory if needed. Main methods use Newton with finite-difference residual Jacobians, default tolerance 3e-11, fallback acceptance 2e-9, and a 32-iteration cap. Objective KKT solves use the previous module's automatic differentiation with a 48-iteration cap. The modal and control scripts preserve their explicit method parameters in the output names.

`build_viewer.py` embeds three complete experiments into `viewer.html`: moderate double/triple motion for ten seconds and the highest-mode transfer for sixty seconds. The last experiment defaults to 8x angle magnification for visibility, with an explicit display control. Its numerical data remain unchanged. The twelve displayed series pass data-shape checks and the generated JavaScript passes `node --check`. Each experiment shares a viewport across its four panels; saved-frame bounds were checked for all four angle-magnification settings. Interactive playback was not verified: the supported in-app browser reports no available browser, and the parent task's automatic approval review rejected the standalone-browser fallback. No alternative browser was launched.
