# Scalar release: objective bias versus solve error

The independent-velocity objective tested here can select a biased release acceleration even when its constrained minimum is accurately found. This statement concerns the scalar spring, a turning-point initial state, exact total-energy conservation, and the particular objective below. It is not a general impossibility result for independent velocities or a proof about every multidimensional objective.

## Objective actually tested

Take one affine spring with extension \(C_n=x_n-\ell>0\), mass \(m\), stiffness \(k\), no gravity, and \(v_n=0\). Define \(\omega^2=k/m\). Position and endpoint velocity are independent unknowns.

\[
\min_{x,v}\quad
F(x,v)=\frac{m}{2h^2}(x-x_n)^2
+\frac{k}{2}(x-\ell)^2
+\frac{m}{2\beta}\left[v-a\frac{x-x_n}{h}\right]^2,
\]

subject to

\[
\frac{m}{2}v^2+\frac{k}{2}(x-\ell)^2
=\frac{k}{2}C_n^2=E_n.
\]

Here \(a=1\) is the soft endpoint target. At this release state, \(a=2\) is the soft midpoint target implemented as \(v+v_n-2(x-x_n)/h\). Both penalties have unit coefficient on the unknown \(v\); beta therefore has the same normalization in the two cases. The code divides the objective by an energy scale, which does not change its minimizers.

## Eliminate the exact energy constraint locally

Let \(q=x-x_n\). On the energy-ellipse branch near the initial position,

\[
q=C_n\left[\sqrt{1-\frac{v^2}{\omega^2C_n^2}}-1\right]
=-\frac{v^2}{2\omega^2C_n}+O(v^4).
\]

Write the negative release velocity as

\[
v=-h\omega^2C_n s.
\]

The physical one-step velocity has \(s\to1\) as \(h\to0\). Exact energy conservation then gives

\[
q=-\frac12h^2\omega^2C_n s^2+O(h^4).
\]

On that same constraint, potential energy is \(E_n-mv^2/2\). Therefore

\[
F-E_n
=\frac{m}{2}\left[\left(\frac qh\right)^2-v^2\right]
+\frac{m}{2\beta}\left[v-a\frac qh\right]^2.
\]

Substitution yields

\[
F-E_n
=\frac12 m h^2\omega^4 C_n^2\,p_a(s)+O(h^4),
\qquad
p_a(s)=\frac14s^4-s^2
+\frac1\beta\left(s-\frac a2s^2\right)^2.
\]

The leading stationary equation is

\[
p_a'(s)=\frac{s}{\beta}
\left[(\beta+a^2)s^2-3as+2(1-\beta)\right]=0.
\]

## Moving branch

For the soft endpoint target, the larger positive stationary root is

\[
\boxed{s_E(\beta)=
\frac{3+\sqrt{1+8\beta^2}}{2(1+\beta)}.}
\]

For the normalized soft midpoint target,

\[
\boxed{s_M(\beta)=
\frac{3+\sqrt{1+6\beta+2\beta^2}}{\beta+4}.}
\]

These are the moving minima of the leading local expansion. In particular,

| Velocity target | beta | Predicted velocity / physical release velocity as h → 0 |
|---|---:|---:|
| Soft endpoint | 1 | 1.5 |
| Soft midpoint | 1 | 1.2 |
| Soft midpoint | 0.1 | 1.042144… |

For fixed positive beta, these moving branches retain an acceleration bias as the timestep shrinks. Soft midpoint approaches the physical value as beta tends to zero; the soft endpoint moving minimum instead approaches twice the physical release acceleration. A timestep-dependent compliance would be a separate formulation to analyze and test.

The frozen state \((x_n,0)\) also remains stationary. Along the energy manifold its second variation is

\[
\delta^2F=m\left(\frac1\beta-1\right)(\delta v)^2.
\]

Thus beta < 1 can make the frozen state a local minimum even when a lower moving minimum exists. This is a branch-selection issue in addition to the moving branch's acceleration bias.

## Dense-search evidence

`landscape.cpp` searches the entire scalar energy ellipse at 65,536 angles, then refines the best sampled interval with 100 ternary-search iterations. It compares that independently searched minimum with the exact-curvature constrained optimizer using 128 iterations. The 36 cases cover release and moving initial states, 60/120/240 Hz, both soft targets, and beta = 0.1/1/10. This is numerical evidence for the global minima in that grid, rather than a formal global proof.

For scalar release, the moving global minimum was found in the searched grid cases. Representative ratios from `results/landscape.csv` are:

| Target and beta | 120 Hz velocity / reference | 240 Hz velocity / reference | Predicted limit |
|---|---:|---:|---:|
| Endpoint, beta = 1 | 1.49396655 | 1.49848360 | 1.5 |
| Midpoint, beta = 1 | 1.19740200 | 1.19934804 | 1.2 |
| Midpoint, beta = 0.1 | 1.04117873 | 1.04190232 | 1.042144… |

The endpoint beta = 1 case at 120 Hz has normalized objective gap \(2.27\times10^{-12}\), energy residual \(9.06\times10^{-12}\), and stationarity residual \(4.26\times10^{-10}\). Its approximately 49.4% excess release velocity therefore cannot be explained by a large solve residual.

There is also a distinct optimization failure in the grid. For `scalar_mixed`, 60 Hz, soft midpoint, beta = 10, the numerical KKT point has normalized objective 1.0080857551, while the dense search finds 0.9824530531: a gap of **0.0256327020**. The KKT point has energy residual \(1.06\times10^{-12}\) and stationarity residual \(9.24\times10^{-12}\), but its endpoint velocity is −0.0617173 versus +0.4914291 at the searched minimum. Small energy and stationarity residuals do not certify global optimality.

## Consequence for the research direction

Two issues need separate experiments: whether the solver reaches the intended objective's desired branch, and whether that branch represents the desired mechanics. More iterations or exact curvature can address solve error; they do not remove the demonstrated release-acceleration bias of this objective. The evidence motivates revising the objective, its velocity/position model, or its mechanical closure while retaining independent velocity as an available variable.

Sources in this experiment suite: `independent.hpp`, `landscape.cpp`, and `results/landscape.csv`. No additional simulations were needed for this derivation.
