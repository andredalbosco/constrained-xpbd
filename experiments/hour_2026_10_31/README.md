# Constrained XPBD: one-hour continuation, 2026-10-31

This package continues the independent-position/velocity, constrained-energy and work experiments. Numerical simulations run as native C++ executables on blacknoir. Sources, saved failures, CSV measurements, trajectories and static figures are included. The main project's solver is unchanged.

The requested research hour began at 18:04:18 UTC. The session records its scope and end in [SESSION.md](SESSION.md). Start with [RESULTS.md](RESULTS.md) for the findings and proposed combinations. [COMPOSITION_MATRIX.md](COMPOSITION_MATRIX.md) distinguishes tested combinations from proposed transfers. [INTEGRATION.md](INTEGRATION.md) maps the candidates onto the current repository interfaces and lists acceptance gates before changing a default.

| Study | Main question | Detailed report |
|---|---|---|
| Rigid double/triple pendulums | Which objective and velocity constraints preserve the intended motion? | [Rigid-chain study](rigid/README.md) |
| Local shared multipliers | Why do local work/kinetic multiplier iterations stall? | [Local solver results](local_solver/RESULTS.md) |
| Numerical coordinates | Which inexpensive conditioning and nonlinear-solve policies help? | [Preconditioning study](preconditioning/README.md) |
| Work elimination and composition | Can work multipliers be removed without discarding their target? | [Exact work reduction](work_reduction/README.md) |

Watch the [60-second pendulum comparison video](rigid/animation/modal_transfer_60s.mp4) for a direct view of the modal-transfer failure. Playback is at normal speed; rod angles are magnified 16× and labeled on screen. The video passed metadata, full-decode and representative-frame checks.

The [rigid-chain recorded-motion viewer](rigid/viewer.html) includes the modal-transfer failure. It is a standalone HTML/Canvas app, not the numerical solver. Data and JavaScript checks pass; interactive browser playback remains unverified because the supported browser connection is unavailable and automatic review rejected the isolated Chrome fallback. The [previous round's viewer](../round2_2026_10_31/viewer.html) contains the five earlier environments and adversarial single-pendulum cases.

The final [five-second larger-chain follow-up](work_reduction/longer_rollouts/RESULTS.md) extends the short-run evidence and records later solver failures and reference-certification limits.

## Build and validate

The package depends on Eigen and source headers from the saved first and second experiment rounds. The root CMake file has explicit cache paths and builds all four study directories. From the saved package:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j 4
ctest --test-dir build --output-on-failure
```

On another checkout, override `XPBD_ROOT`, `COMPOSITION_ROOT`, `ROUND2_ROOT`, and `EIGEN_ROOT`. While this package was staged in Global, validation explicitly pointed `ROUND2_ROOT` at the Global round-two directory. Saved canonical paths are under the repository's `experiments/` directory. Tests write into their build directories and preserve recorded research CSVs.

Both saved packages were rebuilt from clean directories using their canonical dependency paths. All 11 hour checks and all 4 round-two checks passed; see [canonical validation](CANONICAL_VALIDATION.md).

Python analysis uses NumPy and Matplotlib; there is no GPU requirement. The local Python environment used was `/home/blacknoir/Projects/mandelbrot-diffusion/.venv/bin/python`. Set `MPLCONFIGDIR` to a writable directory when regenerating plots. Individual study reports give exact executable modes and limitations.

## Reading the evidence

A completed finite trajectory, a converged nonlinear solve, conservation of invariants, and agreement with a physical reference are separate checks. Tables retain failed prefixes and reference-certification flags. Comparisons of physical motion exclude uncertified references; speed claims use only the explicitly isolated timing passes. Counts include controls and repeated configurations rather than distinct new methods.

No new method has been merged into the main solver. The purpose is to preserve composable evidence and identify which combinations deserve an implementation next.
