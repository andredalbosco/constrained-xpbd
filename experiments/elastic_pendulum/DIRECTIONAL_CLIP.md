# EL with corrections opposing objective descent clipped to zero

Reject an energy correction that points uphill for the original implicit objective, then solve the material row using the accepted energy increment. This experiment tests the direction alone, without a displacement-size bound.

## Exact algorithm

Let

\[
F(x)=\frac{1}{2h^2}(x-y)^TM(x-y)+U_s(x),
\quad y=x_n+hv_n+h^2a,
\]
\[
q=\nabla F(x)=\frac{M}{h^2}(x-y)+\nabla U_s(x).
\]

The objective direction is downhill, `-M^{-1}q`. With the mass inner product, a proposed energy correction opposes that direction precisely when `q^T dx_E > 0`. For this single scalar mass, this also agrees with the ordinary Euclidean direction test.

One fresh EL sweep, with curvature dropped and `r_x` omitted, is:

\[
\Delta\eta_{\rm raw}=-H/(g^TM^{-1}g),\qquad
\Delta x_E=M^{-1}g\,\Delta\eta_{\rm raw},
\]
\[
\Delta\eta=
\begin{cases}
0,&q^T\Delta x_E>0,\\
\Delta\eta_{\rm raw},&q^T\Delta x_E\le0,
\end{cases}
\]
\[
\Delta\lambda=\frac{-r_C-B\Delta\eta}{A},\qquad
\Delta x=M^{-1}(n\Delta\lambda+g\Delta\eta).
\]

Here `g=grad H`, `n=grad C`, `A=n^T M^{-1}n+alpha/h^2`, and `B=n^T M^{-1}g`. Both multiplier visits use the same linearization. The material solve uses the accepted energy increment, and positions move once afterward. Multipliers accumulate during the timestep and reset at its beginning.

The gate uses the full **primal objective gradient**, not the multiplier-dependent stationarity residual `r_x`. The objective gradient is recomputed even though `r_x` is omitted from the increment equation. The rule checks the raw energy contribution before material response, not the combined position update. Exact zero first-order objective change is allowed. Existing zero-energy-gradient handling remains in place. There is no displacement-size limit, nonlinear objective line search, or additional curvature term in this experiment.

## Tests

The same five elastic-pendulum initial conditions, gravity, compliance, predictor, target energy, and reference trajectories are used as in [the original elastic-pendulum study](README.md). Each policy directory has 90 run records: 45 EL configurations (five cases, three timesteps, three sweep counts) plus 45 ordinary-XPBD controls. The unguarded directory repeats the same 45 baselines, so there are 135 distinct scene/method/timestep/sweep configurations across both directories.

All runs complete numerically. Independent audits verify:

- No accepted energy increment has positive first-order objective work when the gate is enabled.
- Every energy attempt is accepted, rejected by the directional rule, or classified as a singular row; displacement-size rejection is disabled.
- Unguarded controls reproduce the previous unguarded solver results exactly.
- Saved energies reconstruct from positions and velocities within `2.36e-15` scaled discrepancy.

All 18 CTest checks pass; the directional tests also pass ASan and UBSan. Regression tests check the objective gradient independently, uphill rejection, downhill and neutral acceptance, response of the material row to the accepted increment, and unchanged behavior when the gate is off. Neutral acceptance is intentional: at an objective stationary point, the first-order test has no preferred direction.

## Results

Fresh directional EL, dropped `r_x`, five sweeps, ten-second runs:

| Case | Max energy error, 60 Hz | Position RMS, 60 Hz | Max energy error, 240 Hz | Position RMS, 240 Hz |
|---|---:|---:|---:|---:|
| Vertical oscillation | 4,841.6% | 0.1600 m | 2.153% | 0.01635 m |
| Small swing | 1.927% | 0.08379 m | 0.1479% | 0.00495 m |
| Larger swing from rest | 2.572% | 0.4939 m | 0.1486% | 0.3569 m |
| Swing with initial velocity | 2.167% | 0.1266 m | 0.1911% | 0.04203 m |
| Coupled stretching and swinging | 1.597% | 0.9996 m | 0.1622% | 0.8215 m |

The useful change is **preventing vertical freezing at 240 Hz**: the mass completes 31 reversals, matching the reference, with 0.999 times its path length. At twenty sweeps, this case improves to 0.00517 m position RMS and 0.408% energy error. The 60 Hz case instead experiences a large energy spike. Its finite output and moderate ten-second position RMS must not be interpreted as stable, accurate dynamics.

The directional rule also improves the initial-velocity swing at 60 Hz compared with unguarded EL, but energy error increases. Larger swings and coupled radial/tangential motion remain inaccurate, including at 240 Hz. The multi-spring extension remains conditional on resolving these single-spring failures.

## Why a downhill direction can still fail

The test enforces only the sign of the first-order term `q^T dx_E`. It does not limit the positive quadratic term in the finite objective change. A large downhill step can cross the objective minimum and overshoot.

For the earlier released scalar spring, write `s=kh^2/m` and `q1` for extension after the first step. At the second predictor, the pathological raw energy correction is `-1.5 q1`. For `s<1`, the objective gradient is `k(1-s)q1`, so their dot product is negative: the large correction passes the requested directional test. This explains why directional clipping and a size bound address different issues.

There is a second limitation: the gate does not change the hard-energy surface or the objective's degeneracy at rest. With endpoint velocity,

\[
\nabla F=g-Mv_n/h.
\]

At rest, `grad F = g`, so the raw energy correction has first-order objective work `-H`: energy removal is allowed and energy restoration is rejected, except where the existing singular-row handling applies. The later material update still changes the state. These observations explain possible failures without adding exceptions to the rule the user requested.

## Reproduce

```bash
cmake --build build --target xpbd-elastic-pendulum-study xpbd-directional-clip-tests
ctest --test-dir build --output-on-failure
./build/xpbd-elastic-pendulum-study results/elastic_pendulum_directional directional
./build/xpbd-elastic-pendulum-study results/elastic_pendulum_unguarded unguarded
python3 experiments/elastic_pendulum/plot_directional.py results/elastic_pendulum_directional results/elastic_pendulum_unguarded --output results/elastic_pendulum_directional/plots
```

The default study mode remains the original size-rejection experiment. Shared solver options default to `objective_gate=false`, preserving all existing experiments. Set `objective_gate=true` and `gamma=-1` for this directional-only rule.

- [Directional summary](../../results/elastic_pendulum_directional/summary.csv)
- [Directional traces](../../results/elastic_pendulum_directional/traces.csv)
- [Unguarded controls](../../results/elastic_pendulum_unguarded/summary.csv)

The plotter requires NumPy and Matplotlib. The validated local interpreter is `/home/blacknoir/Projects/mandelbrot-diffusion/.venv/bin/python`.

- [Energy and motion comparison](../../results/elastic_pendulum_directional/plots/directional_energy_and_motion.png)
- [All five cases at 60 Hz](../../results/elastic_pendulum_directional/plots/directional_vs_unguarded_h60_sweeps5.png)
- [All five cases at 240 Hz](../../results/elastic_pendulum_directional/plots/directional_vs_unguarded_h240_sweeps5.png)
