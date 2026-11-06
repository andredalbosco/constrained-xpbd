# Rejecting oversized energy updates and reusing gradients

These experiments test two proposed solver changes on the original endpoint-energy XPBD formulation: skip an energy increment whose displacement is too large, or reuse an older energy gradient. Skipping helps EL with dropped position stationarity residual on the released and radial spring. Older gradients can freeze motion or destabilize it. Neither change establishes an energy- and momentum-conserving method.

The objective, ordinary predictor, endpoint velocity reconstruction, hard initial-energy target, and curvature-dropped mass matrix are unchanged. There are no reservoirs, midpoint terms, half-force predictors, damping factors, velocity repairs, or softened energy constraints in these experiments.

## Equations and update

For one moving mass, with acceleration `a`,

\[
y=x_n+hv_n+h^2a,\qquad v(x)=(x-x_n)/h,
\]
\[
\min_x \frac{m}{2h^2}\|x-y\|^2+U(x),\qquad
H(x)=\frac m2\|v(x)\|^2+U(x)-ma^Tx-E_\star=0.
\]

The spring uses `C=|x|-L`, `U=C^2/(2 alpha)`. The pendulum uses the hard material constraint `C=0`, no elastic potential, and gravity in `a`. The fixed target `E_star` is initial mechanical energy.

Using scaled material and energy multipliers `lambda`, `eta`, define current physical residuals

\[
R=h^2r_x=m(x-y)-n\lambda-g\eta,\quad
r_C=C+\widetilde\alpha\lambda,\quad r_E=H,
\]

where `n=grad C`, `g=grad H`, and `alpha_tilde=alpha/h^2`. Each outer sweep evaluates these residuals afresh. Reusing a gradient changes only the approximate Jacobian, denoted by `g_bar` below; even `R` still uses the current analytical gradient.

After dropping curvature, eliminating positions gives

\[
\begin{bmatrix}A&B\\B&D\end{bmatrix}
\begin{bmatrix}\Delta\lambda\\\Delta\eta\end{bmatrix}
=\begin{bmatrix}b_C\\b_E\end{bmatrix},
\]
\[
A=n^Tn/m+\widetilde\alpha,\quad B=n^T\bar g/m,\quad D=\bar g^T\bar g/m,
\]
\[
d_0=-R/m\ \text{(retain }r_x\text{), or }0\ \text{(drop }r_x\text{)},
\quad b_C=-r_C-n^Td_0,\quad b_E=-H-\bar g^Td_0.
\]

Each sweep starts its increments at zero and takes one scalar Gauss-Seidel visit per row:

- **EL:** `deta=bE/D`, then `dlambda=(bC-B*deta)/A`.
- **LE:** `dlambda=bC/A`, then `deta=(bE-B*dlambda)/D`.

Positions are then updated once by `dx=d0+(n*dlambda+g_bar*deta)/m`, and both multipliers accumulate. Both scalar visits use the same linearization; the second visit uses the first multiplier increment. This is the previously tested deferred Gauss-Seidel scheme, with relinearization at the next outer sweep. It is not Jacobi and does not reevaluate gradients between the two scalar visits.

### Displacement rejection

For the energy candidate only, compute

\[
d_E=\|\bar g\Delta\eta/m\|,\qquad
d_{\rm step}=h\|v_n\|+h^2\|f_n\|/m.
\]

If `dE > gamma*dstep`, reject the energy increment: set `deta=0` and leave `eta` unchanged. The material row still runs. This tests the raw energy contribution before any cancellation by the material row. Equality is accepted; negative `gamma` disables rejection. A zero energy gradient is handled separately as a singular row, with an incompatibility counter when its right-hand side cannot be satisfied.

The spring force is evaluated from the old position; the pendulum scale uses gravity, without an estimated constraint reaction. The threshold is therefore a heuristic scale, not a universal error bound. It does not limit the objective-residual displacement or the material displacement.

### Older-gradient policies

| Policy | Gradient used in the approximate Jacobian |
|---|---|
| `fresh` | Current analytical gradient at sweep input |
| `lag_sweep` | Previous sweep input gradient, reset to predictor gradient each timestep |
| `frozen_step` | Predictor gradient for every sweep of the timestep |
| `previous_step` | Final gradient of the preceding timestep; first timestep uses fresh gradients |

All reported energies, constraint residuals, and `r_x` values remain current physical quantities.

## Experiment matrix

There are **1,350 ten-second runs**, with **225 complete or guarded-prefix traces** saved at `h=1/60`, five sweeps:

- Released spring, radial-velocity spring, tangential-velocity spring, released pendulum, moving pendulum.
- Timesteps `1/60` and `1/240`; one, five, and twenty outer sweeps.
- EL and LE; retained and dropped `r_x`.
- Fresh gradient: rejection off, or `gamma=0.25,1,4,16`.
- Each older-gradient policy: rejection off or `gamma=1`.
- Ordinary XPBD baseline for every scene/timestep/sweep count.

Spring parameters: `m=L=1`, stiffness `100`, initial position `(1.25,0)`, and velocity `(0,0)`, `(-0.7,0)`, or `(0,1)`. Pendulum parameters: `m=L=1`, gravity `9.81`, angle `0.65` radians, tangential speed zero or `1.2`. Multipliers reset each timestep; only the specified gradient history persists.

Continuous references use RK4 with 32 substeps, independently checked against 64: maximum position disagreement `1.48e-10`, maximum absolute reference energy error `2.13e-11`. Safety limits terminate divergent trajectories; they never clip accepted states. Of the runs, 932 complete without incompatible rows, 240 complete with explicitly skipped incompatible rows, and 178 hit a safety limit. Completion alone does not establish physical correctness.

## Results

At 60 Hz, five sweeps, EL, dropped `r_x`, fresh gradient, `gamma=4`:

| Case | Maximum relative energy error | Position RMS error | Path length / reference |
|---|---:|---:|---:|
| Released spring | 2.884% | 0.0563 | 0.9933 |
| Radial spring | 2.884% | 0.0732 | 0.9937 |
| Tangential spring | 0.00000171% | 1.5152 | 1.2576 |
| Released pendulum | 0.2442% | 0.3269 | 0.9354 |
| Moving pendulum | 0.1440% | 0.2542 | 0.9658 |

Energy errors are normalized by initial energy above the minimum potential: initial mechanical energy for the spring, and `E_initial+mgL` for the pendulum. Position errors use the continuous reference over all completed timesteps. Failed prefixes are not treated as successful full-duration comparisons.

**Skipping restores spring oscillation in this configuration.** The released spring makes 31 velocity reversals, matching the reference, and rejects 156 of 3,000 energy attempts. At 240 Hz its energy error is 0.1733%, but position RMS is still 0.0684. Twenty sweeps barely improve these spring results. One sweep is insufficient: at 60 Hz, `gamma=4`, the released spring reaches 224% energy error. The threshold matters: `gamma=1` gives 6.64% energy error and 0.193 position RMS at five sweeps, while `gamma=16` gives 55.5% energy error. Four is a favorable tested value, not a generally established constant.

**The pendulum favors a different threshold.** With `gamma=1`, dropped `r_x`, five sweeps, the released pendulum has 0.3128% energy error and 0.0379 position RMS; the moving pendulum has 0.2788% and 0.1015. This is better motion than `gamma=4`, despite slightly larger energy error.

**Retaining `r_x` can still explode.** Fresh EL with `gamma=4`, five sweeps, stops after 22 timesteps for the released spring and six for the radial spring. The cutoff limits the new energy contribution, while the retained base displacement can remain large because of previously accumulated multipliers. It is not a bound on the full position update.

**Old gradients do not reliably resolve the problem.** Predictor-frozen EL almost immobilizes the released spring: energy error is just 0.00249%, but path length is `2.10e-7` of the reference and position RMS is 0.307. The unbounded radial version hits a safety limit after seven steps. Adding `gamma=1` prevents that explosion but leaves very little motion. Previous-sweep and previous-timestep gradients also produce failures or frozen/chattering motion. No tested older-gradient policy provides a consistently better replacement across scenes. There are useful individual cases: predictor-frozen LE with gamma=1 and dropped r_x gives the lowest moving-pendulum position RMS at 60 Hz and five sweeps (0.03425), at the cost of 2.069% energy error.

**LE with gamma=1 still freezes the released spring.** Its mutually cancelling corrections need not be large enough to trigger rejection. A tighter gamma=0.25 changes that outcome: fresh LE oscillates with position RMS 0.08455 from rest and 0.07900 with radial velocity, at 10.38% maximum energy error. **The tangential spring remains wrong:** fresh EL at `gamma=1` or `4` rejects no energy updates in the representative case, and angular momentum error reaches about 1.533 from initial angular momentum 1.25. Conserving energy alone does not constrain the orbit. Pendulum angular momentum about the anchor is not conserved because gravity exerts torque.

## Why the cutoff catches the original overshoot

For the released scalar spring, let `s=kh^2/m` and `q` be extension. After the first EL step,

\[
q_1=q_0/(1+s),\qquad v_1=-hkq_1/m.
\]

At the second timestep's predictor,

\[
g=-skq_1,\quad H=-\tfrac32skq_1^2,\quad
\Delta q_E=-H/g=-\tfrac32q_1.
\]

The motion scale is `2s|q1|`, hence the candidate ratio is `3/(4s)`: **27 at 60 Hz, 432 at 240 Hz**. Both `gamma=1` and `4` reject it. This is the raw energy contribution; EL's subsequent material correction can partially cancel it.

Reusing an older gradient instead can select the wrong direction. After the first backward-Euler material step, the current energy gradient vanishes. Using the predictor gradient `kq0` gives an energy correction `+(s/2)q1`, toward the original position. A small energy residual can therefore coexist with frozen motion. Closing a stale linearized energy row does not guarantee that the true current energy improves.

## Reproduce and inspect

From the repository root:

```bash
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release
cmake --build build --target xpbd-gradient-policies-study xpbd-gradient-policies-tests
ctest --test-dir build --output-on-failure
./build/xpbd-gradient-policies-study results/gradient_policies
python3 experiments/gradient_policies/plot_policies.py results/gradient_policies
```

Plotting requires NumPy and Matplotlib. The run used `/home/blacknoir/Projects/mandelbrot-diffusion/.venv/bin/python` for plotting.

- [Summary CSV](../../results/gradient_policies/summary.csv): every setting, physical errors, motion metrics, update counts, rejection ratios, and failure status.
- [Time traces](../../results/gradient_policies/traces.csv): actual energy and Cartesian positions, velocities, reference trajectory, physical residuals, and cumulative counters.
- [Plot audit](../../results/gradient_policies/plots/plot_audit.json): plot list and independent endpoint-energy reconstruction audit.
- [Spring comparison](../../results/gradient_policies/plots/spring_skip_vs_old_gradient.png).
- [Fresh EL, gamma 4, dropped residual](../../results/gradient_policies/plots/fresh_EL_gamma4_drop_rx.png).
- [Same settings, retained residual](../../results/gradient_policies/plots/fresh_EL_gamma4_retain_rx.png).
- [Older-gradient comparison](../../results/gradient_policies/plots/policies_EL_gammaoff_drop_rx.png).
- [Threshold comparison](../../results/gradient_policies/plots/fresh_gamma_EL_drop_rx.png).

All 120 fresh, unbounded controls reproduce the prior deferred LE/EL study exactly for energy error, path length, completion count, and maximum physical `r_x`. Regression tests cover independent linear-system calculations, current residuals under stale Jacobians, history lifecycle, energy-only rejection, strict threshold boundaries, and the second-step overshoot. Optimized and AddressSanitizer/UndefinedBehaviorSanitizer tests pass. Every accepted bounded correction obeys its threshold; attempt counters balance, and the independent plotting audit reconstructs trace energies within `9.81e-16` scaled discrepancy.

See the [complete comparison and recommendations](COMPARISON.md) for the energy/motion tradeoffs, representative trajectory plots, and all-six-settings scorecards.
