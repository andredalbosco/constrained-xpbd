# Nonlinear energy solves, constrained PD, and FEPR

**The explosive scalar EL correction is a linearization failure; exact nonlinear energy fixes that jump but does not fix the trajectory.** With the original endpoint velocity, exact constrained solves can stop at a turning point or freeze from rest. FEPR avoids that particular position-only restriction by projecting independent x and v, but it still has measurable phase and angular-momentum errors. The tested variants have distinct formulations and should not be merged under an “energy projection” label.

This is an isolated research implementation. The existing simulator and historical experiment code were not modified. The current comparison contains 168 five-second planar configurations, 112 scalar configurations, and nine twenty-second continuations, plus implementation validation, 30 targeted FEPR refinement runs, and the one-step full-objective oracle. Failure rows and offending raw-EL guard candidates remain in the outputs; incomplete prefixes are not successful runs. A zero-step failure has an undefined RMS rather than a misleading zero error.

## Implementations

| Family | Code | What is actually solved |
|---|---|---|
| Exact nonlinear E then material L | `scalar_exact.py`, `planar.py` | Exact scalar roots or nearest point on the true nonlinear energy contour; immediate accumulated material-multiplier update afterward. It is not full XPBD KKT. |
| Strict constrained PD-inspired | `planar.py` | Local spring projection, quadratic global objective, true nonlinear energy constraint. Hard-rod control solves rod/energy intersections. |
| Frozen-energy PD control | `planar.py` | Quadratic local/global objective and frozen-projection energy; physical residual is measured after projection refresh. |
| Inertia-only objective control | `planar.py` | Remove spring potential from the objective while keeping physical energy constrained; gravity remains in the predictor. |
| Published CPD energy-only specialization | `cpd_energy.py` | Source-matched relaxed-energy SQP/PD update with its rho variable and external-work target. The hard-rod adaptation is explicitly separate and has failed/frozen cases. |
| FEPR | `fepr.py` | Paper's fast current-iterate projection, independent positions/velocities, energy + linear/angular momentum targets, unrestricted s/t, mass/time-scaled metric. |
| FEPR hard-rod extension | `fepr.py` | RATTLE tentative state, plus simultaneous exact-length and endpoint-tangency rows in the same fast projection. These added rows are not claimed as a method tested in the FEPR paper. |
| Discrete-gradient controls | `references.py` | Compatible midpoint work/momentum with independent endpoint velocity. Spring uses Cartesian midpoint kinematics; rod uses angular midpoint kinematics. |

The primary FEPR plots use the incoming-step energy target, matching the paper. Separate `fepr` rows use a fixed run-initial target; the old-state feasibility argument does not automatically extend to that fixed target. `fepr_base` retains the base integrator without projection (exact radial implicit Euler or RATTLE). This is necessary to isolate what the FEPR projection changes.

## Observed behavior

- Exact nonlinear E→L removes the original second-step spring explosion. With one sweep it can travel to a turning point and stop; more sweeps can freeze release immediately. A small or zero energy residual is therefore insufficient.
- Strict constrained PD freezes released states with the documented nearest-old-state tie. The exact implicit objective is constant on the fixed-energy surface from rest. Removing the spring part of the objective does not repair this: it can favor the stationary point even more strongly.
- On the orbiting spring, strict constrained PD has nearly exact energy and materially wrong angular motion. This is a distribution/trajectory problem after solving energy, not evidence of an inaccurate energy solver.
- The source-matched relaxed CPD spring also nearly freezes release. Its hard-rod extension is not a reliable candidate: multiple runs fail their nonlinear cap, and the endpoint-velocity tangency defect remains.
- FEPR sustains oscillation in all standard spring/pendulum test cases. Its rod extension keeps geometry/tangency to the declared tolerance. Compared with RATTLE alone, it improves energy but can worsen phase accuracy.
- FEPR does not solve the entire reliability problem. The stiff orbiting spring can retain excellent energy while accumulating large trajectory and angular-momentum errors. Substepping is not monotonically improving over the tested finite resolutions. At 30 Hz the second stiff-spring FEPR step jumps 2.008 m with a tiny energy residual. A full fixed-center objective solve avoids that jump from some starts but finds distinct local branches, including a nearly stationary position with nonzero velocity. It is not yet a validated repair. See the targeted refinement analysis and one-step oracle.
- The discrete-gradient controls have the smallest reference errors among these conserving candidates on the regular cases, while changing the mechanical equations and using nonlinear solves. This result does not establish a scalable XPBD implementation or unconditional large-step robustness.

At 60 Hz over five seconds, the ordinary released spring's energy loss reaches essentially 100% under implicit Euler/XPBD. FEPR retains its energy and oscillation, with 0.0359 m RMS position error. The energy-exact E→L eight-sweep and strict-PD solutions stay at the initial position (about 0.307 m RMS). The discrete-gradient spring control has about 0.0119 m RMS. On the released pendulum, FEPR plus hard-rod rows has about 0.00315 m RMS; RATTLE alone has 0.000742 m, and the angular discrete-gradient control 0.000857 m. These are trajectory comparisons, not equal-cost performance benchmarks.

Twenty-second FEPR runs remain finite, but errors accumulate: released spring RMS 0.147 m, orbiting spring 0.502 m, released pendulum 0.00496 m. Corresponding discrete-gradient controls have 0.0470 m, 0.0177 m, and 0.00328 m. Long-run phase and angular momentum matter even when energy looks excellent.

## Derivation, evidence, and playback

- [Main equations and solver ordering](DERIVATION.md)
- [Scalar derivation](SCALAR_MATH.md) and [scalar results](SCALAR_RESULTS.md)
- [PD objective analysis](PD_MATH_REVIEW.md)
- [FEPR paper mapping](FEPR_DERIVATION.md) and [independent audit](fepr_independent_validation/SUMMARY.md)
- [FEPR refinement and branch failures](FEPR_REFINEMENT.md)
- [Full fixed-center FEPR one-step oracle](FEPR_FULL_ORACLE.md)
- [CPD source mapping](CPD_DERIVATION.md)
- [Planar implementation audit](planar_validation/SUMMARY.md)
- [Reference/control validation](VALIDATION_REVIEW.md)
- [Complete measured tables](results/TABLES.md), [summary CSV](results/summary.csv)
- [Motion and energy figure](results/motion_energy.png)
- [Synchronized five-second playback](results/comparison.mp4)

![Motion and energy](results/motion_energy.png)

Raw output lives in `scalar_results`, `planar_results`, `raw_recheck`, `stress_results`, `substep_results`, and `long_results`. `raw_recheck` replaces the initial raw-EL rows in the merged report: it uses the original stable free-fall gradient expression. Planar raw EL refreshes the material normal after E and matches the original refreshed-row implementation; it is distinct from the frozen-normal historical scheme outside the affine scalar case. The historical two-step scalar counterexample is checked directly by compiling the original unmodified C++ solver.

## Reproduction

Python 3.13 with the pinned `requirements.txt` was used. Plotting additionally uses the installed `ffmpeg` executable for the MP4. Build an isolated environment, then run from this directory:

```sh
python3 -m venv .venv
.venv/bin/pip install -r requirements.txt
.venv/bin/python scalar_exact.py --output scalar_results --duration 5
.venv/bin/python run_planar.py --group endpoint --hz 60,240 --duration 5
.venv/bin/python run_planar.py --group fepr --hz 60,240 --duration 5
.venv/bin/python run_planar.py --group cpd --hz 60,240 --duration 5
.venv/bin/python run_planar.py --group dg --hz 60,240 --duration 5
.venv/bin/python summarize.py
```

The archived execution used `/tmp/xpbd-nonlinear-venv/bin/python`, with `OPENBLAS_NUM_THREADS=1` and `MPLCONFIGDIR=/tmp/xpbd-nonlinear-mpl`. `summarize.py` also reads the retained targeted reruns/stress outputs; use `run_planar.py --help` to regenerate those exact case/method selections. Audits: `validation_planar.py`, `validation_fepr_independent.py`, `references.py`, `fepr.py`, and `cpd_energy.py`. The historical C++ parity checks require the sibling canonical repositories and a C++ compiler; they intentionally validate original code rather than an embedded rewritten copy.

## What to carry forward

Keep the nonlinear solve and explicit failure reporting: they identify the dangerous linear step and prevent treating an unconverged visit as an energy solution. Keep separate x/v for methods requiring both hard geometry and endpoint tangency. Treat FEPR as a verified baseline, not a universal fix or the same objective as constrained XPBD. For the next scalable local-solver step, preserve the mechanical/momentum equations alongside work balance and compare against the discrete-gradient controls; energy-only success is not enough to choose a solver.

The single-spring energy-contour code supports the component that does not enclose the anchor. It reports unsupported topology rather than fabricating a root. No claims are made for multi-pendulums, arbitrary networks, collisions, cloth, GPU speed, unit-invariant FEPR stopping, or all initial states. The new FEPR source is a paper-based implementation, not claimed as author-code parity; no public author-maintained source release was found.
