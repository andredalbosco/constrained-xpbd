# Same incoming state: does independent velocity soften the known EL failure?

Yes, it substantially reduces the known first energy/material position
excursion. It does not make one sweep an accurate nonlinear energy solve or
establish stability for arbitrary states.

This diagnostic uses **exactly the same incoming state** for the archived
position-only solver and the three independent-velocity variants. It does not
compare states reached after their separately evolved trajectories. The only
new physical-step solves are one EL sweep and a bounded 16-sweep EL control
from that state for each formulation. A separate factor screen reads existing
saved campaign CSVs without running trajectories.

Run `python3 diagnostic.py` in the experiment's NumPy-capable environment to
regenerate `diagnostics.json`. The old solver defaults to the canonical archived
`block_orders_2026_11_01/solver.py`; `--old-source` can override that location.
The script resolves its new solver and output relative to its own file.

## Fixed incoming state and parameters

The radial spring has \(m=1\), \(k=100\), \(\ell=1\), and no gravity. The
incoming state is the one immediately preceding the old step-2 EL explosion:

\[
x_n=(1.2432432432432432,0),\qquad
v_n=(-0.4054054054054079,0),\qquad h=1/60.
\]

All variants use

\[
\epsilon=31250,\quad K_\star=0,\qquad
E_\star=\frac12m\|v_n\|^2+U(x_n)=3.0405405405405403.
\]

The fixed epsilon corresponds to ratio 10,000 times the old run-initial energy
3.125, rather than a retuned parameter for each method. Initialize
\(x=x_n+hv_n\), \(\lambda=\eta=\rho=0\), and, for independent velocity,
\(v=v_n\). The velocity penalty initially has zero residual. Physical energy
always uses the returned independent velocity in the new method and the chord
velocity in the old method.

## The energy denominator gains the predicted positive term

The cold energy numerator is approximately \(-0.1620708546\) for every variant.
Both primal-stationarity contributions and relaxation-stationarity contribution
are initially zero. The position-only analytical energy gradient nearly
cancels:

\[
\nabla U(x)+\frac{Mv_n}{h}=(-0.6756756756758371,0).
\]

For the same cold state, independent velocity adds precisely

\[
\beta_vv_n^TMv_n=0.1643535427319232\,\beta_v
\]

to the energy Schur denominator. No step-size cap or correction limiter is
introduced.

| Formulation | First E denominator | First \(\Delta\eta\) | First E position correction | Following L position correction |
|---|---:|---:|---:|---:|
| Position only | 0.0004226524 | -383.4614 | 0.0719710 m | 27.60253 m |
| Independent, \(\beta_v=0.1\) | 0.01685801 | -9.613880 | 0.00180441 m | 0.0844527 m |
| Independent, \(\beta_v=1\) | 0.1647762 | -0.9835817 | 0.000184606 m | 0.0129282 m |
| Independent, \(\beta_v=10\) | 1.643958 | -0.09858576 | 0.0000185033 m | 0.00703180 m |

The table reports increment magnitudes. Full signed increments, before/after
states, physical energy, residuals, denominator and numerator, and the separate
position responses to \(r_x\) and \(r_v\) are saved in `diagnostics.json`.
One-sweep outputs match the corresponding production `block_step` states and
multipliers bit for bit.

## Nonlinear remainders remain, mainly in velocity

For this positive-axis spring branch, the first E visit has exact remainder
identities. In the old position-only formulation,

\[
r_x^{\mathrm{after\ E}}
=\Delta\eta\left(\frac{m}{h^2}+k\right)\Delta x,
\qquad
r_\eta^{\mathrm{after\ E}}
=\frac12\left(\frac{m}{h^2}+k\right)\|\Delta x\|^2.
\]

With independent velocity they instead become

\[
\boxed{
\begin{aligned}
r_x^{\mathrm{after\ E}}&=\Delta\eta\,k\Delta x,\\
r_v^{\mathrm{after\ E}}&=\Delta\eta\,m\Delta v,\\
r_\eta^{\mathrm{after\ E}}&=\frac k2\|\Delta x\|^2+\frac m2\|\Delta v\|^2.
\end{aligned}}
\]

These identities use the zero initial multipliers/residuals and the fact that
the spring remains on its positive radial branch during the first E visit.
They are not identities for arbitrary nonlinear constraints or all later
visits. Measured residuals match them to approximately 1e-12 absolute error
for the independent variants, and energy remainder errors are below 3e-15.

For \(\beta_v=1\), the first E velocity correction is
\((-0.40982572,0)\) m/s. After refreshing, \(\|r_x\|=0.0181575\),
\(\|r_v\|=0.4030971\), and the relaxed energy residual is 0.0839803 J.
In contrast, the old method's refreshed position residual is 102112.9 and its
energy remainder is 9.58267 J. The large kinetic energy-gradient sensitivity
has moved into the independent-velocity direction, where the retained coupling
matrix and additional Schur contribution change its response.

The following material block still responds to both refreshed residuals.
For \(\beta_v=1\), its common position response includes approximately
\(-hM^{-1}r_v\), not just the spring constitutive correction. Omitting that
term would be a different method.

| Formulation | Physical energy after first E | Physical energy after one EL sweep | Physical energy after 16 EL sweeps |
|---|---:|---:|---:|
| Position only | 12.50977 | 1,411,608.671 | 3.0404523 |
| Independent, \(\beta_v=0.1\) | 3.161869 | 21.92564 | 3.0404702 |
| Independent, \(\beta_v=1\) | 3.124230 | 4.482844 | 3.0405163 |
| Independent, \(\beta_v=10\) | 3.120825 | 3.949453 | 3.0405374 |

All three independent variants satisfy the production convergence criterion at
16 sweeps for this particular incoming state. The old control narrowly misses
its 1e-9 scaled tolerance at that budget, while already recovering close energy.
Finite positive velocity compliance changes the optimization problem, so the
different converged positions and independent velocities are not expected to
equal the position-only root.

The improvement is substantial, but a single independent-velocity EL sweep
still has material energy error. A positive denominator alone does not imply
that the nonlinear update respects the energy surface.

## Existing beta-one trajectories and the kinematic factor

At velocity stationarity the defining relation is

\[
\frac{x-x_n}{h}=(1+\beta_v\eta)v.
\]

The diagnostic records this factor and its actual identity defect for the
72 existing beta-one, 60 Hz, 4/16-sweep, two-penalty, six-scene, three-order
saved runs. The identity only holds to the accuracy of \(r_v\); a factor sign
at a large-residual iterate is not enough to infer a real velocity reversal.

The useful LE beta-one, 60 Hz, 16-sweep high-penalty release trajectory retains
motion: final-second maxima are approximately 2.4925 m/s independent speed and
2.4754 m/s chord speed. Its final-second factor stays between 0.5230 and
0.9932. It does not approach the near-zero-factor freeze branch.

Other branches still do. For the beta-one JOINT release at 60 Hz, four sweeps,
and epsilon ratio 100 million, the final-second independent speed is at most
1.10e-4 m/s while chord speed is at most 7.52e-8 m/s. The factor remains between
-0.000525 and -0.000388. The visible position is nearly stationary despite
stored spring energy. JOINT with 16 sweeps also stalls, with final-second chord
speed below 7.48e-8 m/s. Kicked-pendulum JOINT at four sweeps likewise stalls
near its apex, with final-second chord speed below 1.74e-7 m/s.

Thus independent velocity makes the known cancellation-driven EL correction
much gentler and gives useful moving split-order trajectories, while leaving
possible energy-preserving freeze branches. The existing good release path
also retains sharp velocity changes at turning points; this diagnostic does
not claim that velocity disentangling removes every visible motion artifact.
