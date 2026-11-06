# CPD energy relaxation: two block-Gauss-Seidel orders

This experiment implements the formulation agreed in the conversation and compares
**EL** (energy block then material block) with **LE** (material then energy).
The energy block solves position, rho and eta together. The material block solves
position and the unscaled material multiplier together. Both retain the full
position stationarity residual, apply the block correction immediately, and
refresh the analytical residuals and derivatives after every block.

The campaign is complete: 468 five-second configurations on the spring release,
spring orbit, gravity spring, resting pendulum, kicked pendulum, and stiff spring.
See [RESULTS.md](RESULTS.md), the [derivation](DERIVATION.md), and the
[step-level instability diagnostic](DIAGNOSTIC.md).

The results contain instability, nearly frozen converged motion, and inaccurate
orbits with very small energy error. This is a measured comparison, not a claim
that the formulation is now reliable. Full original CPD and FEPR are not
implemented by this experiment's JOINT control.

The convergence flag checks stationary equations. A separate selected-root
curvature audit distinguishes those checks from a constrained-minimum claim.

## Objective and convention

The position-only objective is

\[
\min_{x,\rho}\quad
\frac{\|x-x_n-hv_n\|_M^2}{2h^2}+U(x)+\frac\epsilon2\rho^2,
\qquad
\frac{\|x-x_n\|_M^2}{2h^2}+U(x)=(1-\rho)E_\star+\rho K_\star.
\]

For a compliant distance spring, the unscaled convention is
\(\lambda=-C(x)/\alpha\), \(r_\lambda=C(x)+\alpha\lambda\),
and the spring contribution to position stationarity is
\(-\lambda\nabla C\). It differs in scale from a timestep-scaled XPBD multiplier.
The energy multiplier eta uses the plus sign in the Lagrangian. No material
multiplier is substituted into the analytical physical-energy constraint.

Gravity is conservative and included in physical potential. Its objective-force
term remains outside the energy multiplier. A hard pendulum uses zero material
compliance and no elastic spring potential; no endpoint-tangency equality is
added. Returned velocity is always \((x-x_n)/h\).

## What is held fixed

- Cold initialization at each physical step: \(x=x_n+hv_n\), all multipliers and
  rho zero. Accumulate them through every block visit within that step.
- Energy target: the actual incoming physical energy; fixed throughout the step.
- The adapted reference energy is the potential minimum in the same gravity gauge.
- epsilon is a selected ratio times initial energy above that minimum, fixed for
  the entire run. Ratios are 1, 100, 10000, 100000000.
- Frame rates are 60 and 120 Hz, with 1, 4 and 16 sweeps per step. Each EL/LE sweep
  has one material visit and one energy visit. No early convergence exit is used.
- JOINT solves the full current approximate system once per iteration, with the
  same targets, omitted curvature and derivative refresh. XPBD omits the energy block.
- No correction cap, line search, energy rescale, hidden subdivision, or fallback
  is applied. Guard failures retain and label the offending candidate.

## Reproduce

Requires Python, NumPy, SciPy and Matplotlib. ffmpeg is needed only for the video.
The environment used here is `/tmp/xpbd-nonlinear-venv/bin/python`.

```bash
OPENBLAS_NUM_THREADS=1 /tmp/xpbd-nonlinear-venv/bin/python validate.py --write-report
OPENBLAS_NUM_THREADS=1 /tmp/xpbd-nonlinear-venv/bin/python run_experiments.py --duration 5 --output results
OPENBLAS_NUM_THREADS=1 /tmp/xpbd-nonlinear-venv/bin/python metrics_audit.py
OPENBLAS_NUM_THREADS=1 /tmp/xpbd-nonlinear-venv/bin/python diagnose.py
OPENBLAS_NUM_THREADS=1 /tmp/xpbd-nonlinear-venv/bin/python curvature_audit.py
OPENBLAS_NUM_THREADS=1 /tmp/xpbd-nonlinear-venv/bin/python summarize.py
MPLCONFIGDIR=/tmp/xpbd-block-mpl OPENBLAS_NUM_THREADS=1 /tmp/xpbd-nonlinear-venv/bin/python plot_results.py
```

The plots script also generates synchronized playback. Its selected cases and
parameters are recorded in `results/figures/selected_runs.csv`. Failed trajectories stop
at the accepted prefix and are labeled explicitly while the reference continues.

## Files

| File | Purpose |
|---|---|
| `solver.py` | Physical residuals, material and energy blocks, joint control, failure guards |
| `run_experiments.py` | Fixed-budget trajectories, diagnostics, targets, CSV and trace output |
| `reference.py` | Independent continuous DOP853 spring and angular-pendulum references |
| `validate.py`, `VALIDATION.md` | Derivative, dense block, frozen-GS, replay and edge-case checks |
| `metrics_audit.py`, `METRICS_AUDIT.md` | Independent recomputation of every summary and serialized visit |
| `DERIVATION.md` | Objective, full stationary equations, approximations and both block formulas |
| `diagnose.py`, `DIAGNOSTIC.md`, `diagnostics/failing_step.json` | Identical-incoming-state EL/LE/joint failure analysis |
| `curvature_audit.py`, `diagnostics/curvature_audit.json` | Tangent-space curvature classification of selected stationary roots |
| `results/summary.csv` | All 468 configurations; completion and convergence are separate |
| `results/trajectories/` | All physical states, continuous references, residuals and multipliers |
| `results/first_step_traces.json` | Post-block snapshots from the first step of every configuration |
| `results/reference_audit.csv`, `results/metadata.json` | Reference refinement and exact run configuration |
| `results/comparison.mp4`, `results/figures/` | Synchronized motion and scientific plots |

Root application code and earlier experiment archives are unchanged. This
directory is a standalone experiment inside constrained-xpbd.
