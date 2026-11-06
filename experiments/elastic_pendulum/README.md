# EL on an elastic pendulum under gravity

The selected EL method works on some mild motions, but the single elastic pendulum does not yet work reliably. The usual cutoff can preserve energy while freezing a vertical oscillation. Larger coupled motions have substantial trajectory and torque-balance errors. The multi-spring extension was therefore held, following the requested condition to establish the single-spring case first.

## Method and physical setup

This study uses the unchanged `gradient_policies.hpp` update: fresh energy gradient, dropped position stationarity residual `r_x`, curvature replaced by the mass matrix, hard initial-energy target, and the ordinary gravity predictor. Each sweep solves energy first, rejects its increment if the raw energy displacement is too large, then solves the material multiplier using the accepted energy increment. The position receives their sum after the two scalar solves. There is no relinearization between those two visits.

The mass is connected to a fixed anchor by an elastic spring, and is free to stretch and swing in 2D:

\[
U_s(x)=\frac{k}{2}(\|x\|-L)^2,\quad
E(x,v)=\frac{m}{2}\|v\|^2+U_s(x)+mgy,
\]
\[
y_{\rm predictor}=x_n+hv_n+h^2a,\quad a=(0,-g),\quad v_{n+1}=(x_{n+1}-x_n)/h.
\]

Parameters are `m=1 kg`, `L=1 m`, `k=100 N/m`, `g=9.81 m/s^2`. The equilibrium radius is

\[
r_{\rm eq}=L+mg/k=1.0981\ {
m m},\qquad
V_{\min}=-mgL-\frac{(mg)^2}{2k}=-10.2911805\ {
m J}.
\]

Energy errors are normalized by `E_initial-V_min`. This includes the elastic sag contribution. The cutoff uses the net old physical force, gravity plus spring force, rather than adding their magnitudes separately.

## Tests and references

There are **270 ten-second runs**: five initial conditions, 60/120/240 Hz, one/five/twenty sweeps, and ordinary XPBD plus EL with rejection disabled or `gamma=0.25,1,4,16`. All complete without reaching numerical safety limits. Full traces cover both 60 and 240 Hz at five sweeps, with every cutoff and baseline: 60 trajectories.

Angles are measured from the downward vertical. Positive radial velocity points away from the anchor; positive tangential velocity increases the angle.

| Case | Radius | Angle | Radial speed | Tangential speed |
|---|---:|---:|---:|---:|
| Vertical oscillation | 1.2981 | 0 | 0 | 0 |
| Small swing | 1.0981 | 0.15 | 0 | 0 |
| Larger swing from rest | 1.25 | 0.65 | 0 | 0 |
| Swing with initial velocity | 1.0981 | 0.65 | 0 | 1.2 |
| Coupled stretching and swinging | 1.35 | 0.65 | -0.7 | 0.9 |

Reference dynamics include both the elastic force and gravity. RK4 with 32 substeps was checked against 64: maximum position disagreement is `1.51e-10 m`, and maximum reference energy drift is `3.66e-11 J`. The vertical reference also matches the exact solution

\[
x(t)=0,\qquad y(t)=-1.0981-0.2\cos(10t)
\]

within `1.22e-10 m`.

Independent regression tests check the full energy gradient against finite differences; net-force cutoff scaling; the scalar EL update against an independent linear-system calculation; preservation of static equilibrium; the exact backward-Euler damping formula for baseline vertical oscillations; and parity with the original zero-gravity update. Optimized and ASan/UBSan runs pass.

## Results

With the preferred **fresh EL, gamma=4, dropped r_x, five sweeps**:

| Case | Max energy error at 60 Hz | Position RMS at 60 Hz | Max energy error at 240 Hz | Position RMS at 240 Hz |
|---|---:|---:|---:|---:|
| Vertical oscillation | 0.00723% | 0.2404 m | 0.0000715% | 0.2401 m |
| Small swing | 0.000386% | 0.0869 m | 0.00000335% | 0.00413 m |
| Larger swing from rest | 0.00112% | 0.5077 m | 0.00000980% | 0.3577 m |
| Swing with initial velocity | 0.2242% | 0.2290 m | 0.00120% | 0.0439 m |
| Coupled stretching and swinging | 0.000182% | 0.9979 m | 0.000000000489% | 0.8216 m |

**The vertical case is a clear failure hidden by tiny energy error.** It travels from one turning point to the other, about `0.4 m` total, then nearly stops. The reference travels about `12.76 m` over ten seconds. At 60 Hz, 578 of 600 timesteps move less than `1e-10 m`; at 240 Hz, 2,322 of 2,400 do so. These are completed but physically incorrect runs.

**The small swing is a useful positive case at 240 Hz.** It has about 4 mm position RMS error with almost constant energy. The initial-tangential-velocity case also benefits from a tighter cutoff: with `gamma=0.25`, energy/position errors are `0.504% / 0.0563 m` at 60 Hz and `0.0645% / 0.0149 m` at 240 Hz.

**A tighter cutoff is not a general cure.** For the vertical case, `gamma=0.25` restores oscillation at 240 Hz, with 2.73% energy error and 0.0727 m position RMS, but loses almost all excitation energy at 60 Hz. It does not fix the larger or coupled swings. For coupled motion at 240 Hz, no energy increments are rejected for gamma 0.25, 1, or 4, and all three give the same incorrect trajectory.

Baseline XPBD is checked first and has the expected numerical damping. On the vertical mode its energy above equilibrium follows the exact backward-Euler factor `(1+100*h^2)^(-n)`. Completing the baseline run does not mean baseline XPBD conserves energy; it verifies that the scene and force implementation have the expected behavior.

## Angular momentum with gravity

Angular momentum about the anchor should change because gravity exerts torque. The relevant diagnostic is

\[
L_z(t)-L_z(0)-\int_0^t[-mg\,x(s)]\,ds.
\]

The study approximates the integral by the trapezoidal rule over stored timesteps and reports the same quadrature error on the continuous reference. For coupled motion at 240 Hz, EL gamma 4 has maximum torque-balance discrepancy `0.7883 kg m^2/s`, compared with reference quadrature error `4.86e-5` and baseline discrepancy `0.0366`. Thus the motion failure is not just an energy-normalization artifact. `max_angular_change` records raw change from initial angular momentum and is not a conservation error; use `max_torque_balance_error` for this comparison.

## Why gravity changes the observed behavior without fixing the underlying problem

Vertically, the physics is exactly a harmonic oscillator shifted by the spring sag. The solver nevertheless begins with an explicit gravity predictor and a separate elastic material row. That split changes the first EL correction and the numerical trajectory selected relative to the earlier zero-gravity spring experiment.

The structural objective issue remains. With `v_+=(x-x_n)/h`,

\[
F(x)=\frac12\|v_+-v_n-ha\|_M^2+U_s(x)
=E(x)-v_n^TMv_++\text{constant}.
\]

On the hard-energy surface, the objective is flat when the old velocity is zero. Otherwise it favors alignment with the old velocity. This does not uniquely select the physical reversal at a turning point. It helps explain the observed freezing and why closing the energy equation is insufficient evidence of correct dynamics; it is not a proof that every finite EL iteration must freeze.

## Reproduce

```bash
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release
cmake --build build --target xpbd-elastic-pendulum-study xpbd-elastic-pendulum-tests
ctest --test-dir build --output-on-failure
./build/xpbd-elastic-pendulum-study results/elastic_pendulum
python3 experiments/elastic_pendulum/plot_results.py results/elastic_pendulum
```

The plotting script requires NumPy and Matplotlib. The local run uses `/home/blacknoir/Projects/mandelbrot-diffusion/.venv/bin/python`.

- [Energy versus motion: failure detail](../../results/elastic_pendulum/plots/elastic_pendulum_energy_vs_motion.png)
- [All five cases at 60 Hz](../../results/elastic_pendulum/plots/elastic_pendulum_h60_sweeps5.png)
- [All five cases at 240 Hz](../../results/elastic_pendulum/plots/elastic_pendulum_h240_sweeps5.png)
- [All 270 runs](../../results/elastic_pendulum/summary.csv)
- [Energy, motion, and torque traces](../../results/elastic_pendulum/traces.csv)
- [Reference refinement checks](../../results/elastic_pendulum/reference_checks.csv)

The original solver headers remain unchanged. This is a test of the selected EL method, not a new integration rule. The multi-spring pendulum is not implemented in this study because the single-spring prerequisite remains unresolved.

The [directional-clipping follow-up](DIRECTIONAL_CLIP.md) tests EL with energy increments set to zero when they oppose descent of the full implicit objective. It has separate results and preserves the size-rejection runs above.
