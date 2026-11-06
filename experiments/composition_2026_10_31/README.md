# Constrained XPBD: compositional experiments

31 October 2026. The runnable sources, manifests, raw results, unsuccessful runs, analysis, and animations are preserved together. Existing constrained-xpbd studies were used as controls and were not replaced.

**The strongest tested combination is midpoint kinematics + matching momentum impulses + averaged-strength discrete-gradient spring forces.** It preserves the work identities and follows the reference trajectories closely in these small systems. RATTLE is the stronger rigid-pendulum reference when endpoint velocity tangency is required. The independent-velocity objective remains a research branch: accurate energy balance alone does not make its selected motion correct.

Open [the animation viewer](viewer.html), [the numerical findings](analysis/findings.md), [the scalar objective diagnosis](analysis/objective_bias.md), or [the paper-inspired research directions](paper_ideas.md). The viewer is standalone and embeds 200 representative traces; the original CSVs retain all recorded samples.

## What was conducted

The matrix contains **20,877 configured simulation runs**: 20,804 completed the requested horizon; 73 terminated early. Completion does not imply convergence or plausible motion. The analysis flags 647 completed trajectories with less than 10% of the corresponding reference path length. These counts include repeated serial timing controls and overlapping controls between tiers; they are not 20,877 distinct mathematical formulations.

Separately, 5,000 fixed-motion work probes, 1,000 impulse-accounting probes, 480 short EL/LE diagnostic runs, and 36 dense scalar-objective searches were conducted. All 45 analytic work/impulse check groups passed. The objective searches use 65,536 points on each constant-energy ellipse followed by local scalar refinement.

| Tier | Runs | Physical horizon | Purpose |
|---|---:|---:|---|
| Core | 3,720 | 10 s | Work law × work treatment × velocity treatment × scheduling, at 60/120 Hz; rigid tangency on/off |
| Parameters | 1,100 | 3 s | All ten work laws with velocity and work penalty strengths 0.1/1/10 |
| Solver | 13,104 | 0.1 s | Curvature, retained/dropped stationarity correction, frozen/refreshed derivatives, order and 1/4/16 iterations |
| Convergence controls | 720 | 1 s | Joint solves with 64 iterations, all work laws/treatments/velocity choices; exact/dropped multiplier curvature |
| Long runs | 328 | 60 s | All modeling families, ten initial conditions; includes unsuccessful alternatives |
| Stress | 1,400 | 3 s | Timestep, stiffness up to ×100, mass ratio up to 100, near-equilibrium perturbations |
| Stress reference controls | 72 | 3 s | Joint implicit and momentum solves, including 64-iteration implicit controls |
| Recorded comparisons | 328 | 10 s | Synchronized trajectories for visual comparison |
| Serial timing | 105 | 1 s | Five repeated, single-thread runs per selected case/method |

The algorithm factorial is a short-horizon convergence experiment, not a 60-second stability claim. The complete Cartesian product of every parameter and every horizon is not claimed. All executed combinations and their exact parameters are in [manifests/](manifests/); applicability and equivalences below explain why some combinations are omitted.

## Environments

| Environment | Mechanical definition | Initial conditions |
|---|---|---|
| Scalar spring | One 1 kg mass, affine extension x−1, stiffness 100 N/m, no gravity | Release and release with initial velocity |
| Planar spring | One mass attached to a fixed origin by a radial spring | Stretch release and stretching/orbiting motion |
| Rigid pendulum | One mass, fixed pivot and length, planar gravity | Angular release and tangential kick |
| Elastic pendulum | One compliant radial spring and one mass under gravity | Mixed swing/stretch and vertical release |
| Two-spring pendulum | Fixed anchor, two serial compliant springs, two moving masses, gravity | Mixed relative velocities and large angular/stretch release |

Nominal masses and rest lengths are 1; gravity is 9.81 m/s². Stress cases rebuild the initial static extension for the changed stiffness/masses. The scalar spring remains affine even if its coordinate changes sign; the legacy EL adapter uses the older radial implementation and explicitly fails if its scalar predictor or completed sweep crosses the anchor.

## Composable choices

| Axis | Implemented choices |
|---|---|
| Spatial work estimate | Old force; endpoint force; average full force; separate averaged strength/direction; old/endpoint/averaged strength with exact discrete geometry |
| Other work estimates | Fully frozen start power; old/endpoint force dotted with average velocity over the step |
| Work variables | Prescribed work; independent transfers penalized toward the estimate; free-work/global-energy equivalence |
| Velocity | Hard endpoint reconstruction; independent velocity with soft endpoint target; independent velocity with soft midpoint target |
| Rigid geometry | Position constraint; position plus endpoint tangency |
| Local solve | Joint; overlapping edge blocks including both endpoint kinetic rows; storage/kinetic passes; symmetric passes |
| Linearization | Retain/drop primal stationarity correction; retain/drop multiplier-weighted curvature; freeze/refresh derivatives; reverse visit order |
| Mechanical controls | XPBD, raw/bounded EL and LE, momentum-closed force methods, position-projected Verlet, standard RATTLE, fixed-target FEPR-style energy projection |

Unpenalized free transfers are analytically eliminated: on the connected incidence graphs here they reduce to one global energy balance. Their equality to the explicit global branch is validated rather than counted as a different algorithm. Hard endpoint reconstruction ignores beta. Tangency is only imposed on rigid rods; a compliant spring can change length.

For the momentum-closed methods, kinetic-exact work is computed from the actual impulses and common endpoint average velocity. The separate sequential-impulse audit demonstrates exact telescoping but order-dependent allocation of cross terms. It is an accounting experiment, not an additional fully implemented sequential storage-balanced integrator.

The three discrete-geometry force-strength labels coincide for a rigid rod in the momentum solver, because the reaction strength is its own solved unknown. Frozen power and kinetic force estimates are scalar work estimates; they do not define an additional force law there. These meaningless duplications are excluded.

## Objectives and solver interpretation

The independent branch minimizes the inertial position objective plus potential energy and a velocity penalty. Gravity appears in the potential, not also in its inertial predictor. The endpoint penalty is

\[
\frac{1}{2\beta_v}\left\|v-\frac{x-x_n}{h}\right\|_M^2,
\]

and the midpoint comparator uses

\[
\frac{1}{2\beta_v}\left\|v+v_n-2\frac{x-x_n}{h}\right\|_M^2.
\]

Both have the same coefficient on v. Penalized work uses W=estimate+t with penalty gamma·||t||²/(2 E_scale). Particle kinetic and spring/gravity storage balances are retained, with identically zero prescribed rows removed. Shared vertex multipliers remain shared when adjacent edge blocks visit them.

Second-order automatic differentiation includes the work estimate's derivatives. The exact-curvature option retains multiplier-weighted constraint Hessians; the dropped option retains the objective Hessian. KKT solves use an SVD rank threshold of 1e−11, diagonal regularization 1e−10, a scaled primal increment cap, and backtracking on the active residual. Iteration budgets are explicit. A small residual certifies neither a global minimum nor a desirable time-integration branch. The [scalar objective study](analysis/objective_bias.md) independently identifies both objective bias and an example of convergence to a non-minimal stationary branch.

For the scalar release and beta=1, the moving minimum of the tested endpoint objective has initial acceleration approaching **1.5×** physical acceleration as h shrinks. The tested midpoint penalty gives **1.2×**. This is a diagnostic for this objective/constraint combination, not a criticism of independent velocity in general or of FEPR's different fixed-target objective.

In the momentum family, joint solves solve the nonlinear kinematic/momentum equations. Local schedules are minimum-norm equation-block projections, not the same objective block minimization as the independent branch. Standard local XPBD uses its conventional accumulated multipliers and approximations; curvature and retained-r_x toggles apply only to its joint control. RATTLE/SHAKE position projection and the momentum solver have no separate energy-multiplier curvature toggle. Frozen derivatives in the independent/momentum solvers span a timestep; the legacy EL frozen model spans a sweep. Interpret switches within their family.

The SHAKE-labeled control is position-projected velocity Verlet with the final tangency projection omitted. RATTLE includes that projection. Midpoint discrete-gradient rods preserve their chord-work identity but do not impose endpoint tangency. The FEPR-style control is a fixed-target energy projection, optionally with rod position/tangency; it does not reproduce the full paper's momentum constraints or auxiliary variables.

## Findings to carry into the next iteration

1. **Keep matching kinematics, impulses, and work together.** Averaged-strength discrete gradients without the appropriate mechanical closure are insufficient. Their momentum-closed midpoint composition is the useful reference.
2. **Keep RATTLE as the rigid reference.** Its small bounded energy error buys exact endpoint tangency; this differs from the energy-exact midpoint rod.
3. **Revisit the independent objective before optimizing it harder.** The scalar bias persists at its actual minimum; solver improvements cannot remove that bias. The independent variables and shared-multiplier block machinery remain useful experiments.
4. **Keep failures visible.** Raw EL can divide a finite energy residual by a nearly vanishing gradient; frozen LE can achieve zero energy error by not moving. Bounded corrections change which residuals are satisfied.
5. **Try paper-inspired numerical techniques against a fixed formulation.** [Research notes](paper_ideas.md) propose physical trust regions and feasibility/stationarity filters, graph/mass preconditioning, explicit-work splitting, and cautiously bounded spectral balancing. These new transfers from Muon, PPO/TRPO, Large Steps, and related work are not measured results yet.

## Evidence and validation

- [analysis/energy_vs_motion.png](analysis/energy_vs_motion.png): energy and trajectory errors in matched physical cohorts.
- [analysis/work_kinematics_quality.png](analysis/work_kinematics_quality.png): work/velocity compositions at a fixed four-iteration budget.
- [analysis/solver_residuals.png](analysis/solver_residuals.png): iteration/scheduling behavior.
- [analysis/rigid_energy_tangency.png](analysis/rigid_energy_tangency.png): energy versus endpoint tangency.
- [analysis/coverage.csv](analysis/coverage.csv), [analysis/failures.csv](analysis/failures.csv), and [analysis/best_combinations.csv](analysis/best_combinations.csv): machine-readable coverage, failures, and qualified candidates.
- [analysis/serial_cost.csv](analysis/serial_cost.csv): single-thread compute-only timings; these small C++ prototypes use different derivative/linear-solver implementations, so timings are not production scaling claims.
- [audits_notes.md](audits_notes.md): exact identities and the reproduced EL failure trace.

All three numerical validation executables pass. The independent derivative checks cover 720 columns with relative gradient/Hessian errors below 1.3e−10. Tests also check shared-vertex incidence signs, free/global equivalence, midpoint energy and impulse closure, RATTLE length/tangency, and second-order scalar convergence.

References use paired RK4 resolutions, with subdivisions increased for stiffness. Position/velocity comparisons require the reported refinement thresholds to pass. Some 60-second coupled trajectories are not certified; their energy/path diagnostics remain available, but they are excluded from trajectory rankings. Failed prefixes are never ranked as completed simulations. The viewer's interpolation is for playback only, and failed trajectories stop visibly at their recorded endpoints.

## Reproduce

Requires C++17, Eigen from the existing constrained-xpbd checkout, CMake, OpenMP, and Python 3. NumPy and matplotlib are required only for analysis.

```bash
python3 run_suite.py --threads 8 --output rerun_results
```

This builds, validates, regenerates the explicit manifests, and executes all tiers into a new directory. Use `--tiers core stress` for a subset and `--xpbd-source /path/to/constrained-xpbd` for another checkout. Existing nonempty output directories are rejected.

```bash
/home/blacknoir/Projects/mandelbrot-diffusion/.venv/bin/python analyze.py --results rerun_results --out rerun_analysis
python3 build_viewer.py
```

The supplied viewer generator uses the recorded `results/traces` and `results/long` folders by default; see `--help` for other paths. `study --help` describes direct manifest execution. All physical and solver parameters are in the manifests. The randomized audit seed is 42; trajectory runs are deterministic.
