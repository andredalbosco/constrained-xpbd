# Independent endpoint velocity: energy-constrained XPBD experiment

Independent endpoint velocity substantially improves the tested split-order
solves. The current candidate is **LE, beta_v = 1, 60 Hz, 16 sweeps, epsilon
ratio 1e8**. All six scenes keep moving through five seconds; the worst maximum
actual mechanical-energy error is **0.000706%**. Four sweeps also give useful
moving trajectories, with a worst error of 0.0924%. EL improves too: at 16
sweeps its worst error across the six scenes is 0.0793%, with ongoing motion.
Extending the selected LE setting to 60 seconds preserves motion in all six
scenes, with no guard rejection and a worst energy error of 0.01984%.

These are finite-duration measurements of approximate conservation and motion.
Sharp turning-point changes and smooth artificial circularization remain.
The joint-system control can still select nearly frozen branches. No phase,
frequency, or reference-position RMS score is used to reject a trajectory.

## Objective and velocity treatment

The experiment keeps the implicit objective and uses a positive primal velocity
compliance:

\[
\min_{x,v,\rho}\quad
\frac{\|x-x_n-hv_n\|_M^2}{2h^2}+U(x)
+\frac1{2\beta_v}\left\|v-\frac{x-x_n}{h}\right\|_M^2
+\frac\epsilon2\rho^2,
\]
subject to
\[
\frac12v^TMv+U(x)=(1-\rho)E_\star+\rho K_\star.
\]

Return the optimized velocity independently and use it in the next predictor
and physical-energy evaluation. The coupling is a primal penalty; beta_v must
be positive. A coefficient pi saved in the CSV is reconstructed as
\(M[v-(x-x_n)/h]/\beta_v\); it is not a separately iterated multiplier.
The derived-velocity control is reused from the prior archive, rather than
evaluating this primal penalty at zero compliance.

The energy target is incoming actual physical energy, held fixed throughout
each timestep solve. The adapted reference energy is minimum potential in the
same gravity gauge. Epsilon equals its selected ratio times run-initial energy
above the potential minimum and stays fixed throughout the run. Gravity is in
the potential; it is not also added to the predictor or energy target.

The material convention remains \(\lambda=-C/\alpha\), with
\(r_\lambda=C+\alpha\lambda\). Energy uses the analytical spring potential,
not the current material multiplier. Hard rods have zero material compliance,
no elastic potential, and no additional endpoint-tangency equation in this
primary comparison.

## Approximation and iteration

Keep the full position/velocity coupling Hessian and its cross blocks. Drop
material curvature and energy-multiplier curvature as in the earlier study.
Retain both stationarity residuals \(r_x,r_v\).

- Material visit: jointly solve position, velocity, and material multiplier.
- Energy visit: jointly solve position, velocity, rho, and energy multiplier.
- Immediately apply each full block, accumulate multiplier increments, and
  refresh all analytical residuals and gradients after every visit.
- EL visits energy then material; LE reverses those blocks. JOINT solves the
  full current approximate system and refreshes once per iteration.
- Cold-start each step at \(x=x_n+hv_n\), \(v=v_n\), and zero multipliers/rho.

There is no step limiter, energy rescale, line search, hidden substepping,
velocity overwrite, or fallback. A JOINT iteration does different work from an
EL/LE sweep, and these measurements are not an equal-compute performance test.

[DERIVATION.md](DERIVATION.md) gives objective, stationarity, curvature
approximation, algebraic material substitution, full system, and block formulas.
Its optional length/tangency extension is derived but was not implemented or
tested in the primary campaign.

## Measured energy and motion

The following share beta_v = 1, 60 Hz, and epsilon ratio 1e8. All entries reach
five seconds and sustain motion. Energy errors are the maximum actual change
from run-initial kinetic plus potential energy, divided by initial energy above
minimum potential, expressed in percent.

| Scene | EL, 4 sweeps | LE, 4 sweeps | EL, 16 sweeps | LE, 16 sweeps |
|---|---:|---:|---:|---:|
| Released spring | 4.28% | 0.0924% | 0.0000185% | 0.0000217% |
| Orbiting spring | 0.0000625% | 0.00000025% | 0.00000025% | 0.00000025% |
| Spring with gravity | 0.198% | 0.000231% | 0.00000306% | 0.00000301% |
| Released pendulum | 1.49% | 0.0162% | 0.0793% | 0.000705% |
| Kicked pendulum | 1.43% | 0.0160% | 0.00595% | 0.00000260% |
| Stiff orbiting spring | 1.04% | 0.00637% | 0.00000094% | 0.00000099% |

The released spring remains a moving, bounded oscillator rather than the old
large-excursion case. Both pendulums keep swinging in these split-order runs.
The spring with gravity keeps bouncing and swinging. The unforced orbit and
stiff spring remain smooth-looking in sampled paths, but radial oscillation
progressively becomes rotation; velocity disentangling has not removed that
mode transfer. The released spring also retains strong velocity changes at
turning points, rather than becoming uniformly smooth.

The full grid is six scenes, beta_v 0.1/1/10, 60/120 Hz, 1/4/16 sweeps, epsilon
ratios 1e4/1e8, and EL/LE/JOINT: **648 five-second configurations**.

| Method | Completed | Guard rejections | Nearly stationary final-second tails |
|---|---:|---:|---:|
| EL | 216/216 | 0 | 0 |
| LE | 216/216 | 0 | 0 |
| JOINT | 216/216 | 0 | 44 |

Completion means the fixed time budget was reached within the explicit guards;
it does not certify energy quality or plausible motion. The stationary-tail
flag means final-second path below 1e-4 rest lengths and is a reporting
descriptor, not a user-defined acceptance tolerance. In particular, JOINT can
retain spring energy while its position barely moves. Some joint pendulum runs
stall at an apex or restart on an unusual branch.

511 configurations include at least one step above the residual tolerance.
Convergence is reported independently from energy and motion. Lower sweep
counts can accumulate significant energy drift even without a runaway guard;
for example, the one-sweep beta-one LE released spring reaches 46% error.

## Focused 60-second stability check

The selected LE setting was extended to 60 seconds in all six scenes without
changing the objective, epsilon, or iteration budget. All six complete, none
develops a stationary tail, and the velocity-stationarity factor never becomes
negative or approaches zero (minimum about 0.4749).

| Scene | Maximum actual energy error over 60 s |
|---|---:|
| Released spring | 0.000262% |
| Orbiting spring | 0.000000965% |
| Spring with gravity | 0.0000521% |
| Released pendulum | 0.01984% |
| Kicked pendulum | 0.0000125% |
| Stiff orbiting spring | 0.00000571% |

The released pendulum accumulates more drift than in its five-second run, so
these remain approximate conservation results. The [longer-run report](LONG_HORIZON.md)
records motion, coupling, constraints, and a compact energy/path plot. Independent
reconstruction passes for all 21,606 saved states; 192 serialized first-step
visits across the six runs are additionally replayed by dense block solves.

## Why the known EL excursion becomes smaller

Free velocity changes the energy Schur denominator to

\[
h^2\left[\nabla U(x)+\frac{Mv}{h}\right]^TM^{-1}
\left[\nabla U(x)+\frac{Mv}{h}\right]
+\beta_vv^TMv+\frac{(E_\star-K_\star)^2}{\epsilon}.
\]

The additional positive \(\beta_vv^TMv\) term prevents the old combined-gradient
cancellation from eliminating this response when velocity is nonzero.

The [same-incoming-state diagnostic](DIAGNOSTIC.md) holds the state and epsilon
fixed in the old exploding spring example. At beta_v = 1, the first energy
position correction shrinks from 0.07197 m to 0.000185 m, and the following
material position correction shrinks from 27.60 m to 0.01293 m. One sweep still
has nonlinear energy error; this is evidence for a specific improved response,
not a global stability guarantee.

## Independent velocity and visible motion

Actual state energy uses \(v\). We also report the energy obtained by replacing
it with \((x-x_n)/h\), as a kinematic diagnostic. Those energies need not be
equal when endpoint velocity is independent. At beta_v = 1, 60 Hz, and four
sweeps, the largest chord-energy discrepancy is 15.9% in the stiff spring,
while its actual-energy discrepancy is 0.00637%; the largest velocity mismatch
is about 0.434 m/s. At 120 Hz and 16 sweeps the worst chord-energy discrepancy
across the six cases falls to about 4.08%, with maximum mismatch about 0.167 m/s.
These are different configurations, rather than a convergence-order claim.

At velocity stationarity,
\[
\frac{x-x_n}{h}=(1+\beta_v\eta)v.
\]
A near-zero factor can support small visible displacement despite independent
velocity. A negative factor can reverse their directions when stationarity is
satisfied. The factor is monitored alongside its residual and actual dot
product, rather than used alone as an automatic rejection rule. The selected
LE runs avoid persistent stalls; the joint control demonstrates that the
objective still admits problematic branches.

Independent rod velocities retain a small radial component because no separate
velocity-tangency constraint is imposed. The optional coupled length/tangency
formulation in the derivation remains a distinct future experiment.

## Inspect and reproduce

- [Playback: beta-one EL/LE/JOINT, 60 Hz and four sweeps](results/figures/comparison.mp4)
- [Actual energy, chord energy, mismatch, and coupling factor](results/figures/energy_coupling.png)
- [Position paths, spring lengths, and both speeds](results/figures/paths_lengths_speeds.png)
- [Recommended 16-sweep setting: paths and speeds](results/figures_16/paths_lengths_speeds.png)
- [Visual review and inspection limits](results/figures/VISUAL_REVIEW.md)
- [All 648 measured configurations](results/summary.csv)
- [All shared settings](results/figures/settings.csv)
- [60-second runs and independent audit](LONG_HORIZON.md)
- [Equation and update validation](VALIDATION.md)
- [Independent trajectory-metric audit](METRICS_AUDIT.md)

Visual inspection covers sampled frames, paths, and velocity/energy traces,
rather than continuous viewing of every animation frame. The primary campaign
does not save first-step traces; independent validation replayed 630 refreshed
block visits separately. The metric audit reconstructs all 292,248 recorded
states independently and covers accepted-state versus rejected-state semantics.
There were no rejected runtime examples in this campaign.

The experiment is standalone Python within constrained-xpbd; it is not integrated
into the native application. Existing application code and earlier experiments
are preserved. The archived prior 468-run summary is reused as a derived-velocity
baseline without rerunning it. Its target/maximum-energy conventions should be
read alongside the new runner's explicitly separate accepted/trial statistics.
Each scene has one free planar point mass and a fixed origin anchor; networks,
multiple pendulums, contacts, and arbitrary initial conditions remain outside
this experiment.

The environment used was `/tmp/xpbd-nonlinear-venv/bin/python`, with NumPy and
Matplotlib; ffmpeg is required for playback encoding. From this directory:

```bash
OPENBLAS_NUM_THREADS=1 /tmp/xpbd-nonlinear-venv/bin/python validate.py
OPENBLAS_NUM_THREADS=1 /tmp/xpbd-nonlinear-venv/bin/python runner.py
OPENBLAS_NUM_THREADS=1 /tmp/xpbd-nonlinear-venv/bin/python metrics_audit.py
OPENBLAS_NUM_THREADS=1 /tmp/xpbd-nonlinear-venv/bin/python diagnostic.py
MPLCONFIGDIR=/tmp/xpbd-independent-mpl OPENBLAS_NUM_THREADS=1 /tmp/xpbd-nonlinear-venv/bin/python plot_results.py
OPENBLAS_NUM_THREADS=1 /tmp/xpbd-nonlinear-venv/bin/python runner.py --duration 60 --hz 60 --sweeps 16 --beta-v 1 --epsilon-ratios 100000000 --orders LE --output long_horizon --trace-first
```

The canonical old solver is needed for the diagnostic and hard-coupling-limit
validation. `BASELINE_INPUTS.sha256` identifies the old solver and archived
comparison inputs. `FILES.sha256` covers the final experiment archive,
excluding itself.
