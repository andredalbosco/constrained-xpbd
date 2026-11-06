# Focused 60-second stability check

All six scenes completed 60 seconds with the recommended **LE, beta_v=1,
60 Hz, 16 sweeps, epsilon ratio 100 million** setting. There were no guard
failures, detected sustained stalls, negative kinematic factors, or approaches
to the near-zero-factor branch. Worst actual physical energy drift was
**0.01984%**, in pendulum release.

This extends one selected setting rather than rerunning the parameter grid.
No solver code, targets, guards, damping, or convergence rules were changed.
The test does not remove the previously documented turning-point velocity
changes, central-force spinup, or independent/chord velocity mismatch.

![Physical energy drift and continued accepted position travel](long_horizon/energy_and_motion.png)

## Scope and energy convention

The run contains 3,600 physical timesteps per scene. Epsilon is fixed for each
run at 100 million times that scene's initial energy above the potential floor.
The per-step target remains incoming accepted energy, fixed within its solver
visits. Endpoint velocity remains independent and is returned without an
overwrite or projection.

“Actual energy” means

\[
\frac12mv^Tv+U(x),
\]

using that returned velocity. Relative drift is divided by run-initial energy
above the common potential floor. It is not normalized by an arbitrary absolute
gravitational potential value, and it is not chord-velocity energy.

| Scene | Maximum actual energy drift | Position travel during final second | Largest accepted position increment |
|---|---:|---:|---:|
| Spring release | 0.0002617% | 1.499 m | 0.04131 m |
| Orbiting spring | 0.0000009653% | 1.270 m | 0.02134 m |
| Spring with gravity | 0.00005209% | 1.323 m | 0.03510 m |
| Pendulum release | 0.01984% | 1.201 m | 0.03329 m |
| Kicked pendulum | 0.00001254% | 1.091 m | 0.02878 m |
| Stiff orbiting spring | 0.000005707% | 3.241 m | 0.05434 m |

Spring release, orbit, gravity spring, and stiff spring end with small energy
losses. Both pendulums end with small gains. The plot shows the growth of those
drifts over the full horizon, so approximate energy conservation is not
presented as exact or indefinitely bounded conservation.

Every scene retains substantial position motion during the final second.
The runner's final-tail stall flag is false and its longest detected sustained
stall duration is zero for all six. This does not imply that short turning-point
pauses or sharp velocity restarts have disappeared.

## Velocity relation and pendulum geometry

Across all 21,600 accepted endpoint states, the minimum value of
\(1+\beta_v\eta\) is **0.4748984**. It is never negative and never below 0.05
in magnitude. The selected long runs therefore do not enter the kinematic
factor regime associated with the frozen JOINT branches in the shorter study.

The largest length errors are 0.0001063 m for pendulum release and
0.000004804 m for the kicked pendulum, on unit-length rods. No endpoint-tangency
constraint has been added: maximum radial independent velocities remain
approximately 0.03329 and 0.02487 m/s, respectively.

Independent velocity differs from the position chord. The maximum mismatch is
0.2477 m/s for release, 0.00885 for the orbit, 0.04260 for gravity spring,
0.06548 and 0.04668 for the pendulums, and 0.4346 for the stiff spring. These are
measured outputs of the chosen compliant coupling, not values projected away.
For comparison, the maximum energy drift computed from chord velocity is 2.95%
for release and 15.89% for the stiff spring, despite the much smaller actual
independent-velocity energy drift. These two kinetic-energy interpretations must
remain distinguished when assessing visible motion or downstream uses of v.

## What this check does and does not establish

The selected solver setting retains bounded positions, small physical energy
drift, and ongoing movement for the tested 60 seconds in all six scenes. This
is stronger evidence than the original five-second rollout. It is not a proof
for other initial conditions, systems with additional constraints, arbitrary
duration, or a different iteration budget.

Some fixed-budget steps do not meet the strict 1e-9 scaled convergence tolerance.
Converged fractions are 99.972% for spring release, 99.944% for gravity spring,
97.361% for pendulum release, and 96.778% for the kicked pendulum; orbit and
stiff spring converge on every step. These finite outputs are preserved and
measured rather than silently solved with extra sweeps.

Previously observed motion concerns remain. Release's largest independent
velocity increment is 0.6999 m/s in one 1/60 s step, including a sharp
turning-point restart. The unforced orbit and stiff spring still gain angular
momentum: final values are approximately 1.29085 and 3.27786 from initial 0.88.
The long check supports numerical survival and approximately conserved energy;
it does not make those existing behaviors physically faithful.

## Independent saved-data audit and reproduction

The independent metrics audit passed **694,705 checks** over all **21,606 CSV
rows**. It independently recomputed mechanical energy, chord energy, mismatch,
geometry, residuals, summary statistics, and failure bookkeeping. It also
replayed all six serialized first-step traces using independent dense block
assembly: **192 block visits**. No production rollout was rerun by the audit.

Results and provenance:

- `long_horizon/summary.csv` and `long_horizon/trajectories/`: all measured states.
- `long_horizon/metadata.json`: run settings and six scene definitions.
- `long_horizon/first_step_traces.json`: all six first-step block traces.
- `long_horizon/AUDIT.json`: independent audit results.
- `long_horizon.log`: runner output.

Reproduce the focused run in the experiment environment:

```bash
OPENBLAS_NUM_THREADS=1 python3 runner.py --duration 60 --hz 60 --sweeps 16 --beta-v 1 --epsilon-ratios 100000000 --orders LE --output long_horizon --trace-first
```

Audit without overwriting the primary campaign report:

```python
import json
from pathlib import Path
import metrics_audit

result = metrics_audit.audit_results(Path("long_horizon"))
Path("long_horizon/AUDIT.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
```
