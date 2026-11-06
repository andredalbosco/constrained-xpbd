# Constrained XPBD: second composition study

This study extends `composition_2026_10_31` with new objective formulations, larger coupled systems, start-only work predictions, actual local primal/dual solves, and adversarial pendulum mechanics. The prior experiment and project solver are preserved. These are research prototypes, not a merged production integrator.

Read [RESULTS.md](RESULTS.md) for the synthesis and [viewer.html](viewer.html) to compare recorded motion. Each branch has its own formulation, raw results, limitations, and reproduction instructions:

| Branch | Main question | Evidence |
|---|---|---|
| [Independent objectives](objectives/README.md) | Does the objective produce the correct acceleration, including at moving states? | AD derivatives, dense scalar searches, timestep refinement, angle/speed probes, 10/60-second trajectories |
| [Adversarial mechanics](adversarial/RESULTS.md) | What fails despite exact work/energy, and what does endpoint tangency change? | Stiff/fast/separatrix cases, nonlinear acceptance, step-doubling, reversal, symmetric and angular controls |
| [Numerical preconditioning](preconditioning/RESULTS.md) | Can the same physical equations be solved better under finite budgets? | Chains up to32 masses, modal diagnostics, equivalent-root checks, setup-inclusive isolated timings |
| [Local work](local_work/RESULTS.md) | Can useful start-only work targets and shared kinetic multipliers be iterated locally? | Prediction audits, explicit-work AL, reduced curvature, short rollouts, cross-objective combinations |

## Build and verification

Eigen and the first composition study are dependencies. Their paths can be overridden:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release \
  -DXPBD_ROOT=/home/blacknoir/Projects/constrained-xpbd \
  -DCOMPOSITION_SOURCE=/home/blacknoir/Projects/constrained-xpbd/experiments/composition_2026_10_31
cmake --build build -j 4
ctest --test-dir build --output-on-failure
```

The CMake project now defaults to Release for single-configuration generators. The four CTest suites cover mechanical/projection identities, objective derivatives, equivalent preconditioned roots, and local-work derivatives/identities. `validation_release.log` records the final optimized validation. Numerical reference checks and trajectory residuals are additional evidence, not replaced by those unit checks.

Each branch README gives the data-generation commands. Objective executables write into their current directory; run them from `objectives/`. Local-work executables accept an output directory. Adversarial executables take an output directory followed by an optional mode. The preconditioning benchmark writes into its current directory; run it from `preconditioning/`.

Python analysis needs NumPy and Matplotlib. On this workstation the existing interpreter `/home/blacknoir/Projects/mandelbrot-diffusion/.venv/bin/python` supplied these dependencies. Set `MPLCONFIGDIR` to a writable temporary folder when reproducing from a sandbox.

```bash
python3 build_viewer.py
python3 -m http.server 8769 --bind 127.0.0.1
```

The viewer embeds its trajectory data and can also be opened as a standalone HTML file. It shows recorded states, reference certification, failed prefixes, and inexact solves. The nominal objective trajectories solve at120 Hz and store frames at60 Hz; adversarial rod recordings solve/store at30 Hz. Preconditioning data contain solver/mode diagnostics rather than full animation states and are shown in that branch's figures.

## Interpretation rules

- A completed run can contain unconverged steps. Both facts are recorded separately.
- Physical ODE references are refinement-checked; uncertified long/chaotic references are excluded from trajectory rankings.
- Preconditioning references solve the same discrete DG equations at the same step size. Their errors measure incomplete-solve distortion, not temporal discretization error.
- The isolated preconditioning timing mode includes setup and factorization. Other reported timings are concurrent diagnostics and are not performance evidence.
- A true energy/work identity needs compatible momentum and kinematics. A frozen predicted-work target is an estimate for the optimized endpoint.
- Every paper transfer is identified as a specific implemented analogue. Dense SVD row balancing is not the Muon optimizer; capped Newton is not PPO; the local AL prototype is not a convex ADMM convergence result.

## Work still open

The current evidence supports further combinations of repaired objectives, endpoint tangency, mechanically closed work, numerical scaling, and careful local multiplier updates. It does not establish one best combination. Remaining gaps include larger-chain local AL with corrected objectives; sparse/block alternatives to dense SVD balancing; nonlinear branch continuation; adaptive methods that control tangency and phase together; 3D constraints; collisions/friction; damping/driving; and GPU execution under a fixed real-time budget. These are not silently counted as tested.
