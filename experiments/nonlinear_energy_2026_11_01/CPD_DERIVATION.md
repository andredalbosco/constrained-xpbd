# What the CPD control implements

`cpd_energy.py` implements the **pinned, energy-only specialization** of Kee, Um, Jeong and Han, *Constrained Projective Dynamics* (2021), DOI [10.1145/3450626.3459878](https://doi.org/10.1145/3450626.3459878). The spring implementation matches the existing local comparison baseline; it is not a post-processing energy projection and not the full free-body energy-and-momentum method. The rigid-pendulum function is a separately labeled extension described below.

The primary equations were checked in the archived [paper text](/home/blacknoir/Projects/dynamics-research/archive/resources/papers/Constrained%20Projective%20Dynamics.md), especially Eqs. (12), (15), (19), and (21). The implementation contract is [baselines.md](/home/blacknoir/Projects/dynamics-research/docs/work_balance/baselines.md), and the executable comparison source is [baselines.cpp](/home/blacknoir/Projects/dynamics-research/experiments/work_balance/baselines.cpp). Those files remain unchanged.

## Anchored spring objective and energy target

For a planar particle attached to an origin pin, define

\[
y=x_n+h v_n+h^2 g,\qquad U(x)=\tfrac12 k(\|x\|-l)^2,
\]
\[
F(x)=\frac{m}{2h^2}\|x-y\|^2+U(x),\qquad
v(x)=\frac{x-x_n}{h}.
\]

The internal kinetic-plus-elastic energy constrained by this paper convention is

\[
H(x)=\frac{m}{2h^2}\|x-x_n\|^2+U(x),\qquad
H_*=\tfrac12 m\|v_n\|^2+U(x_n)+h m g\cdot v_n.
\]

The last term is first-order external power. Gravity potential is **not** included in `H` or its target. The objective still includes gravity through `y`. Incoming energy uses the actual incoming velocity.

The historical baseline disables both momentum constraints for pinned systems; with neither target enabled, its compatible kinetic lower bound is `K_*=0`. The specialization therefore solves

\[
\min_{x,\rho}\ F(x)+\frac\epsilon2\rho^2
\quad\text{subject to}\quad H(x)-(1-\rho)H_*=0.
\]

`rho` is an unbounded real optimization variable; it is never clamped to `[0,1]`. `epsilon=1e8` is the historical default in the current unit system, and has energy units. Finite epsilon permits nonzero relaxation even when the strict energy constraint is feasible. Setting `relax=False` removes rho and imposes strict `H=H_*`; this is a separately selectable variant.

For a central anchored spring without gravity, angular momentum about the anchor is physically conserved, but this **energy-only** specialization does not constrain it. Its orbit error must not be attributed to the full momentum-constrained free-body CPD algorithm.

## Source-matched coupled update

Let `q=(x,rho)`, `G=F+epsilon rho^2/2`, `c=H-(1-rho)H_*`, and use the plus-sign Lagrangian `G+lambda c`. Then

\[
\nabla_xF=\frac{m}{h^2}(x-y)+k(\|x\|-l)\frac{x}{\|x\|},
\]
\[
\nabla_xH=\frac{m}{h^2}(x-x_n)+k(\|x\|-l)\frac{x}{\|x\|},\quad
\nabla_\rho G=\epsilon\rho,\quad \nabla_\rho c=H_*.
\]

The anchored distance-spring PD metric is constant:

\[
A=\operatorname{diag}((m/h^2+k)I_2,\epsilon).
\]

With a row Jacobian `J`, one iteration is exactly the paper's coupled quasi-Newton SQP update:

\[
(JA^{-1}J^T)\lambda^{new}=c-JA^{-1}\nabla G,
\]
\[
q^{new}=q-A^{-1}(\nabla G+J^T\lambda^{new}).
\]

The multiplier is freshly solved, not accumulated. The entire step is accepted. There is no line search, added Schur regularization, energy projection, or velocity correction. The constant PD metric omits the nonlinear Lagrangian curvature. The initialization is `(y,0)`, matching the baseline code. The Schur solve uses symmetric row normalization and a relative rank threshold `1e-12`, matching its rank-revealing treatment.

Convergence checks both endpoint stationarity and the energy/relaxation equations using the historical scales: `h^2/(ml)` for position stationarity, `ml^2/h^2+|H_*|` for energy, and `epsilon+|H_*|` for rho stationarity. A tiny scaled energy residual alone is insufficient. `nfev` counts evaluations, while `iterations` counts coupled SQP steps.

Source mapping:

| Implementation choice | Existing source |
|---|---|
| Position predictor and initial iterate | `initialize` |
| Actual incoming-velocity energy plus first-order power | `conservation` |
| Backward kinetic energy, rho row, energy gradient | `conservationRows` |
| Constant PD metric and undamped coupled Schur step | `runProjective` |
| Backward endpoint velocity | `finalize` |

## Why the released spring nearly freezes

When `v_n=g=0`, the position objective equals `H(x)`. Strict energy makes the objective constant on the feasible energy surface, and the old state is feasible and stationary. This is a degeneracy of the objective, not a failure of the linear solver.

On the regularized branch, stationarity gives `lambda=-1` and `rho=H_*/epsilon`. Consequently

\[
H=H_*-H_*^2/\epsilon.
\]

The state is still not uniquely determined by that scalar energy surface. The local update from the old/predicted state chooses a nearly frozen branch for large epsilon. For the historical fixture `m=2,k=100,l=1,x_n=(1.3,0),v_n=0,h=.02`, the implementation returns `v_x=-3.37500194e-7`, `rho=4.5e-8`, and residual `1.18e-14`. Strict energy returns the exactly frozen state. This is an obstruction regression, not a successful physical spring release.

## Explicit hard-rod extension and its limits

The source baseline rejects zero compliance. `cpd_pendulum_step` therefore **extends** the pinned objective by removing spring storage, using `A_x=(m/h^2)I`, and adding the equality

\[
c_L(x)=\tfrac12(\|x\|^2-l^2)=0.
\]

This row is solved jointly with the energy row. Its gradient is `x`; no large finite stiffness is substituted. It keeps the same first-order-power target, relaxation variable, and backward velocity. It does not add endpoint tangency or independently free the velocity.

Exact incoming and outgoing length imply

\[
x_{n+1}\cdot v_{n+1}=\frac{\|x_{n+1}-x_n\|^2}{2h}
=\frac h2\|v_{n+1}\|^2.
\]

Thus nonzero motion necessarily has a normal endpoint-velocity component. Adding exact endpoint tangency without changing the backward velocity relation would allow only zero motion. This extension cannot match the exact position-plus-tangent-velocity state of the angular DG reference while retaining the paper's kinematics.

There is another exact degeneracy at release. With `v_n=0`, the paper target is `H_*=0`, so its energy row is `m||x-x_n||^2/(2h^2)=0`, independent of rho. The only feasible position is `x_n`; its energy gradient vanishes. At a nonvertical state the gravity tangent cannot be balanced by a finite multiplier on the length row. In floating point, a tiny displacement can supply a tiny energy gradient and a huge energy multiplier, leading to a deceptively small KKT residual. At angle `.65,h=.08`, the measured returned displacement is `1.11e-16`, energy multiplier about `4.30e14`, and residual `5.68e-17`. The row-normalized Schur rank is two, but the dimensionally scaled constraint Jacobian has rank one. This is frozen, singular behavior; a `converged` flag does not validate regularity or physical release.

`CPDResult` therefore additionally exposes the multiplier vector and constraint-Jacobian rank/smallest singular value. For this rank diagnostic, position variables are scaled by `l`, rho is dimensionless, and rows use the energy and length-squared scales above. These diagnostics do not change the update or its historical stopping rule.

For the kick fixture `theta=.65,omega=1.2,h=.08,m=l=1,g=9.81`, the extended solve converges with exact length but radial speed `0.0120` and mechanical energy change about `-0.3023 J`. The latter is consistent with the **different energy target**, not evidence that the relaxed source energy equation failed.

## Executed independent validation

`cpd_source_check.cpp` was compiled against the unchanged `baselines.cpp` and its Eigen/core headers. `cpd_source_validation.csv` records 96 cases: four scalar/planar/gravity initial states, two step sizes, epsilon `1e4` and `1e8`, strict/relaxed targets, and budgets 1, 5, and 200. Compared with the Python port:

| Check | Maximum discrepancy |
|---|---:|
| Position | 4.44e-16 m |
| Velocity | 1.11e-14 m/s |
| Rho | 1.19e-18 |
| Scaled stationarity | 2.51e-16 |
| Convergence flags / iteration counts | 0 mismatches / 0 mismatches |

`cpd_validation.json` saves the summary. `cpd_hard_rod_validation.csv` preserves all 30 extension cases, with 26 reported residual convergences and four failures. Residual convergence includes the singular frozen states discussed above and is not a mechanics pass count.

The main five-second driver results in `planar_results/cpd_summary.csv` are a separate trajectory test. At 60 Hz the spring-release path length is only `1.17e-6 m`, while its position RMS error is `0.307 m`. The gravity spring fails at 3.55 s at 60 Hz; its 240 Hz run completes with up to about 3.56% physical-energy deviation. All four hard-rod extension trajectories stop early (including a zero-step release failure at 60 Hz). Their retained prefixes must not be compared as full five-second trajectories.

These trajectory energy statistics use **mechanical energy including gravity potential**. The CPD residual uses the source's relaxed internal-energy/first-order-power target. Those are separate quantities, and should remain separate in any comparison or plot. No successful global minimum search, long-time stability, or exact hard-rod CPD algorithm is claimed.

Reproduce the native parity fixture with:

```sh
g++ -std=c++20 -O1 \
  -I/home/blacknoir/Projects/dynamics-research/experiments/work_balance \
  -I/home/blacknoir/Projects/dynamics-research/src/core \
  -I/home/blacknoir/Projects/dynamics-research/external/eigen \
  cpd_source_check.cpp \
  /home/blacknoir/Projects/dynamics-research/experiments/work_balance/baselines.cpp \
  -o /tmp/xpbd_cpd_source_check
/tmp/xpbd_cpd_source_check > cpd_source_validation.csv
/tmp/xpbd-nonlinear-venv/bin/python cpd_energy.py
```
