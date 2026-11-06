# CPD comparison: spring motion, spinup, and pendulum robustness

The paper-equation CPD implementation prevents the unforced angular-momentum
gain seen in our energy-only methods, but it retains nearly frozen spring
motion. Its external-work target also permits physical energy drift, and some
compliant pendulum solves fail to close their nonlinear constraints.

This is a standalone numerical experiment. It does not change the application
or earlier solvers. The primary reference is Kee et al.,
[Constrained Projective Dynamics (2021)](https://doi.org/10.1145/3450626.3459878).
We implement its coupled Eq. (21) update; this is not the authors' software.

Start with [RESULTS.md](RESULTS.md). The comparison movie is
[comparison.mp4](results/figures/comparison.mp4).
It compares full CPD, removal of only its angular-momentum row, and the archived
independent-velocity LE method. The displayed two-particle coordinate is their
relative displacement, with reduced mass one, matching the anchored LE example.
See [MODE_AUDIT.md](results/MODE_AUDIT.md) for radius and momentum measurements.
The [visual review receipt](results/figures/VISUAL_REVIEW.md) records trace and
sampled-frame inspection plus a successful full FFmpeg decode.

| Artifact | Purpose |
|---|---|
| [SOURCE_FIDELITY.md](SOURCE_FIDELITY.md) | Published equations, local choices, declared adaptations |
| [DERIVATION.md](DERIVATION.md) | Objective, constraints, approximation, coupled solve, multiplier semantics |
| [OBJECTIVE_DEGENERACY.md](OBJECTIVE_DEGENERACY.md) | Why conserved energy and momentum can leave motion selection undetermined |
| [FEASIBILITY_DIAGNOSTIC.md](FEASIBILITY_DIAGNOSTIC.md) | An elastic state with incompatible exact conservation targets |
| [GRAVITY_TARGET_DIAGNOSTIC.md](GRAVITY_TARGET_DIAGNOSTIC.md) | Initial gravity-driven energy/torque target tension |
| [VALIDATION.md](VALIDATION.md) | Independent assembly, finite differences, dense KKT checks |
| [CPP_PARITY.md](CPP_PARITY.md) | Cross-check against our unchanged earlier local C++ implementation |
| [METRICS_AUDIT.md](METRICS_AUDIT.md) | Independent checks of saved results and energy accounting |
| [results/summary.csv](results/summary.csv) | All 351 primary configurations |
| [angular_ablation/summary.csv](angular_ablation/summary.csv) | Five angular-row removal controls |
| [pd_control/summary.csv](pd_control/summary.csv) | 33 ordinary PD controls |

## Reproduce

The recorded environment is Python 3.13.7, NumPy 2.5.4, Matplotlib 3.11.2,
SciPy 1.18.1, and ffmpeg. The solver requires NumPy; SciPy is used for the standalone
feasibility diagnostic. C++ parity rebuilding requires a compiler and Eigen;
the previously measured reference CSV is included.

```bash
export OPENBLAS_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1
python3 validate.py --write-report
python3 cpp_parity.py
python3 runner.py
python3 runner.py --methods PD --policies fixed --epsilon-ratios 10000 --output pd_control
python3 runner.py --cases free_release,free_orbit,free_stiff,spring_orbit,spring_stiff --hz 60 --epsilon-ratios 100000000 --policies fixed --without-angular --output angular_ablation
python3 feasibility_diagnostic.py
python3 audit_modes.py
python3 plot_results.py
python3 compare_unforced.py
python3 metrics_audit.py --all-campaigns
```

`runner.py --resume` reuses matching saved configurations. The original campaign
reused 36 pilot configurations and computed 315 additional configurations.
The independent-velocity comparison trajectories are reused from the earlier
experiment without rerunning or modifying that method. Source/input hashes
are recorded separately from the archive's own `FILES.sha256` manifest.

Raw conservation closure, full KKT stationarity, finite completion, physical
energy conservation, and plausible motion are separate measurements throughout.
A finite step returned after an iteration cap is not called converged.
