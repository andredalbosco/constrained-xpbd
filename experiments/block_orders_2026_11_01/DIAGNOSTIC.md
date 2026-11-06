# Same-incoming-state diagnostic of the first large EL correction

The first large EL correction in the spring-release test occurs on physical
step 2, in the **material block after the energy block**. The energy block
first overshoots its nonlinear energy row and creates a large energy multiplier.
Refreshing the energy gradient then creates a large position-stationarity
residual. The material block responds to that residual with a 27.6 m correction.

This explanation is verified by exact remainder identities for this scalar
branch. It does not establish that every observed long-run instability has the
same mechanism, or that a different solver is a validated repair.

All numbers and per-visit states are in
[diagnostics/failing_step.json](diagnostics/failing_step.json). The diagnostic
calls the production increment functions from `solver.py`, refreshes after each
block, and uses the production finite-state/runaway guards. For every recorded
step and comparison, its output position and multipliers matched `block_step`
bit for bit. No production solver or runner was modified.

Reproduce the trace with `python3 diagnose.py` in the experiment's NumPy-capable
environment. The script resolves its output relative to its own file, so an
absolute invocation also works from another working directory. Rerunning it
reproduced the existing JSON byte for byte (SHA-256
`57c9cf363f714b2a5e62f687b2897375be1045813dfa9df7f6989c737825a497`).

## Configuration and reproducible incoming state

The scene is the one-particle spring released from rest:

\[
m=1,\quad k=100,\quad\ell=1,\quad
x_0=(1.25,0),\quad v_0=(0,0),\quad h=1/60.
\]

Gravity is zero, \(K_\star=0\), and the fixed run parameter is

\[
\epsilon=10^4\left[\frac12m\|v_0\|^2+U(x_0)-K_\star\right]=31250.
\]

Each physical step uses its incoming energy as \(E_\star\), cold-starts
\(x=x_n+hv_n\), \(\lambda=\rho=\eta=0\), and performs the specified number of
full sweeps. The first-large-visit criterion is a position increment larger than
\(0.1\ell\); it is a diagnostic threshold, not a correction cap or production
failure condition.

The first EL physical step has a no-op E visit, followed by the material visit.
It yields the exact incoming state used for **all** order comparisons below:

\[
x_n=(1.2432432432432432,0),\qquad
v_n=(-0.4054054054054079,0),
\]

\[
E_\star=3.0405405405405403.
\]

The cold predictor for step 2 is
\(x=(1.2364864864864864,0)\). Replaying EL from the initial scene therefore
finds the first large visit after only two physical steps, rather than after a
long divergence between separately evolved EL and LE trajectories.

## What the first E visit does

At this cold predictor,

\[
r_x=0,\qquad r_\rho=0,\qquad
r_\eta=-0.16207085463842308,
\]

and the analytical endpoint-energy gradient is

\[
\frac{m(x-x_n)}{h^2}+\nabla U(x)
=(-0.6756756756758371,0).
\]

The kinetic contribution and spring contribution nearly cancel. The two
denominator contributions are

\[
\frac{h^2}{m}
\left\|\frac{m(x-x_n)}{h^2}+\nabla U(x)\right\|^2
=0.00012681600519444416,
\]

\[
\frac{(E_\star-K_\star)^2}{\epsilon}
=0.00029583637691745797.
\]

Their sum is \(0.00042265238211190216\). The position-stationarity and
relaxation terms in the energy numerator are both zero at this cold start.
Consequently, the energy block makes

\[
\Delta\eta=-383.46135381655773,\quad
\Delta\rho=0.037309753344313724,\quad
\Delta x=(-0.07197097481543718,0).
\]

After the visit,

\[
x=(1.1645155116710493,0),\qquad
\frac{m\|x-x_n\|^2}{2h^2}+U(x)=12.509767972325513.
\]

The physical energy has increased despite the intended energy correction. The
relaxed energy residual is \(9.582669249385926\), not a small closure error.

## Exact nonlinear remainder and the subsequent material jump

Throughout this first E visit the particle remains on the positive radial axis.
The distance gradient is constant and the endpoint energy is quadratic, with
curvature

\[
\frac{m}{h^2}+k=3700.
\]

The E visit exactly satisfies its frozen linear position equation: its frozen
position residual is \((0,0)\). But the analytical energy gradient changes by

\[
\left(\frac{m}{h^2}+k\right)\Delta x
=(-266.2926068171176,0).
\]

Because the new multiplier is large, refreshing this gradient gives

\[
\boxed{
r_x^{\mathrm{after\ E}}
=\eta^{\mathrm{after\ E}}
\left(\frac{m}{h^2}+k\right)\Delta x
=(102112.92352143214,0).
}
\]

The identity agrees with the directly evaluated residual to
\(8.74\times10^{-11}\) in absolute norm. The relaxed energy residual has the
exact second-order remainder

\[
\boxed{
r_\eta^{\mathrm{after\ E}}
=\frac12\left(\frac{m}{h^2}+k\right)\|\Delta x\|^2
=9.582669249385942.
}
\]

The directly measured value differs by \(1.60\times10^{-14}\).

Next, the material block sees \(\lambda=0\) and
\(r_\lambda=0.16451551167104927\). Its numerator terms are

\[
-r_\lambda=-0.16451551167104927,
\qquad
\frac{h^2}{m}\nabla C^Tr_x=28.364700978175595.
\]

Its denominator is \(\alpha+h^2/m=0.010277777777777778\). Thus the huge
position-stationarity contribution dominates the material constitutive defect.
For this one-dimensional branch the same position update can be written

\[
\boxed{
\Delta x=
-\frac{r_x+r_\lambda/\alpha}{m/h^2+1/\alpha}
=(-27.6025338034052,0).
}
\]

It crosses the anchor and ends at
\(x=(-26.438018291734153,0)\), with output energy
\(1.4116086710326045\times10^6\). The distance normal flips after that crossing,
and the material residual then becomes \(52.876036583468306\).

The large post-E multiplier also makes the dropped energy-multiplier curvature
large. The material visit nevertheless uses the stipulated mass/constitutive
approximation. These observations explain this measured large push; they do not
justify silently removing \(r_x\), adding damping, or restoring curvature and
calling the resulting method unchanged.

## Changing order and budget from exactly the same incoming state

Every row below restarts from the incoming state above, uses the same fixed
\(\epsilon\), and cold-initializes the same unknowns. JOINT is the complete
approximate system solve with refreshed residuals, serving as a reference for
block partitioning. All these particular runs pass the production finite-state
guard, even when their final residual is very large. “Accepted” here means a
finite fixed-budget output, not a converged physical solution.

| Order | Sweeps | Output x-coordinate | Physical energy | Scaled residual | Largest block displacement |
|---|---:|---:|---:|---:|---:|
| EL | 1 | -26.43801829 | 1,411,608.671 | 14,116.05744 | 27.60253380 |
| EL | 4 | 0.8273815112 | 312.7836215 | 3.097432246 | 27.60253380 |
| EL | 16 | 1.2273100688 | 3.040452254 | 2.57e-9 | 27.60253380 |
| EL | 64 | 1.2273100614 | 3.040452511 | 1.11e-16 | 27.60253380 |
| LE | 1 | 1.2267226516 | 3.061431944 | 0.003372308 | 0.006391527 |
| LE | 4 | 1.2273097245 | 3.040464177 | 7.72e-5 | 0.006391527 |
| LE | 16 | 1.2273100614 | 3.040452511 | 1.15e-11 | 0.006391527 |
| LE | 64 | 1.2273100614 | 3.040452511 | 8.93e-17 | 0.006391527 |
| JOINT | 1 | 1.1613778247 | 13.36564426 | 29.04810005 | 0.07510866 |
| JOINT | 4 | 1.2243776896 | 3.157903787 | 0.002311008 | 0.07510866 |
| JOINT | 16 | 1.2273100614 | 3.040452511 | 8.98e-17 | 0.07510866 |
| JOINT | 64 | 1.2273100614 | 3.040452511 | 8.98e-17 | 0.07510866 |

The scaled residual follows `solver.scaled_residual`; it combines position
stationarity, material constitutive error, relaxation stationarity, and relaxed
energy closure with the declared scene scales. Full unscaled residual vectors
and their values before/after every visit are in the JSON trace.

For this incoming state, increased budget lets all three methods recover to the
same nonlinear residual root, despite EL's very large intermediate excursion.
That does not make the excursion harmless: a more complicated scene may cross
undefined geometry, trigger a guard, or enter a different branch before recovery.
This diagnostic does not test those hypothetical outcomes.

JOINT with one iteration also overshoots the energy row and produces
\(\eta=-376.29503178709064\). Joint coupling therefore removes the particular
split-block response but does not make the approximate nonlinear solve accurate
in one step.

The same-state comparison establishes a concrete finite-budget ordering effect.
It does not establish that LE is stable for every scene or timestep, that EL
always diverges, or that convergence to this root gives accurate long-time
motion. Those are rollout-level questions for the main experiment.
