# What CPD fixes, and what it retains

Full CPD stops unforced self-spinup when its angular-momentum row closes.
It does not reliably restore intuitive spring oscillation or robust,
energy-conserving pendulum motion in these tests. This conclusion concerns
our source-equation implementation and declared numerical choices; it is not
a reproduction of the authors' scenes or executable.

## Experiment scope and fidelity

The primary grid contains 351 five-second configurations: 13 scenes, timesteps
1/30, 1/60, and 1/120 s, three epsilon ratios (100, 10,000, 100 million), and
three stopping protocols. Epsilon is fixed during each rollout at the ratio
times initial mechanical energy above the declared potential minimum. These
are experimental choices; no universal numeric epsilon is specified in the
inspected paper text.

The solver uses the full PD matrix, refreshed spring projections, simultaneous
momentum/energy updates, and fresh conservation multiplier values on every
global iteration. It uses backward endpoint velocity, without an independent
velocity unknown. No damping, line search, rescaling, constraint-row deletion,
or pseudoinverse fallback is applied.

| Family | Primary configurations | Full five seconds | Interpretation |
|---|---:|---:|---|
| Free finite-spring bodies with P/L/energy | 135 | 135 | Original free-body equation core |
| Central anchored finite springs and compliant pendulums | 162 | 162 | Anchored restriction: support exchanges linear momentum; angular momentum includes known external torque |
| Exact hard pendulums | 54 | 33 | Separately declared hard-distance extension; 21 singular-Schur failures |

Completion means the numerical states remained finite and passed the declared
runaway guards. Several completed rollouts contain unconverged capped steps.
Completion alone is not evidence of acceptable motion or energy conservation.
Five angular-row ablations, 33 ordinary PD controls, and twelve one-step
feasibility fixtures supplement the grid. Earlier independent-velocity LE
trajectories are reused for the movie.

The paper-style stop checks raw conservation residual below 1e-4 after at
least one update. The other protocols use fixed 16 updates or scaled full-KKT
residual below 1e-9, with cap 200. The chosen norm, minimum-one-update rule,
and cap are explicit local interpretations/settings. Full source differences,
including the physically consistent torque sign, are in
[SOURCE_FIDELITY.md](SOURCE_FIDELITY.md).

## Main matched results

The table uses 60 Hz, 16 global iterations, epsilon ratio 100 million, and five
seconds. Energy error means maximum absolute mechanical-energy change divided
by the initial energy above the declared potential floor. For unforced springs
the floor is zero. Phase/frequency accuracy is not an acceptance criterion.

| Example | Measured result | Implication for motion |
|---|---|---|
| Free spring released from 25% extension | Relative travel only 3.75e-7 m; energy error 0.000300% | Nearly frozen despite stored elastic energy |
| Free rotating spring, k=100 | Last-second radius range 1.88e-8 m around radius 1.10 m; energy error 0.000108% | Rotation continues with essentially no spring extension oscillation |
| Free rotating stiff spring, k=1000 | Last-second radius range 2.83e-8 m; angular drift below 5.5e-14; energy error 0.000264% | Prevents self-spinup while retaining nearly frozen extension |
| Free spring with opposite constant loads | Mechanical-energy error 4.2543%, despite relaxed energy-row residual near roundoff | Incoming-velocity work does not equal actual displacement work |
| Anchored compliant pendulum, release, k=10,000 | Mechanical-energy error 14.221% | Fails the requested conservation criterion |
| Anchored compliant pendulum, kick, k=10,000 | Mechanical-energy error 12.497% | Fails the requested conservation criterion |

The free two-mass springs use mass two per endpoint and reduced mass one.
Their relative coordinates match the corresponding mass-one anchored springs;
the equivalence was independently checked numerically. Opposite loads are not
uniform gravity. Pendulums in this table are finite springs under gravity,
with the anchored momentum restriction explicitly declared.

Removing only the angular row, while keeping K_star unchanged, increases the
stiff spring's angular momentum from 0.88 to 3.277856 over five seconds.
With the angular row it remains 0.88. Thus the row prevents a specific
unphysical energy redistribution. It does not repair radial motion by itself.
The controlled ablation is not labeled full original CPD.

## More iterations and stopping rules

The release freeze and fixed-radius rotation persist under the strict KKT
protocol. They are not eliminated by converging the coupled numerical system.
At 60 Hz and epsilon ratio 100 million, the paper-style stop takes one global
update in the two-particle unforced examples. Its five-second energy errors
are 1.7558% for the rotating spring and 0.33750% for the stiff version;
fixed16/strict KKT reduce those errors substantially while retaining the
unintuitive spring motion.

The three-particle free triangle is less benign numerically: at these settings,
paper/fixed16/strict-KKT mechanical-energy errors are respectively 0.084915%,
0.0078024%, and 0.000026398%. Three strict-KKT physical steps hit the cap.
These examples do not establish universal freezing of all deformable bodies.

The compliant pendulums expose a further nonlinear-solve problem. At 60 Hz,
epsilon ratio 100 million, release errors are 623.60%/14.221%/623.25% for
paper/fixed16/strict-KKT, and kick errors are 6.4262%/12.497%/137.61%.
More iteration is not a monotonic cure for these undamped solves.

The worst completed compliant setting is the kicked pendulum at 30 Hz,
paper-style stopping, epsilon ratio 100 million. At t=1/3 s its relaxed
energy residual is about 2018.686 J. Accumulated constraint closure error
accounts for 2018.762 J of the 2031.841 J mechanical-energy increase; work
approximation accounts for 13.078 J and relaxation is negligible. This is an
unconverged capped solve, not a conserved-energy result. The runner records
such finite capped states so the failure remains inspectable.

The hard-rod extension also develops severe conditioning problems. At 60 Hz
and epsilon ratio 100 million all three policies fail: release near 0.25 s
and kick near 0.50–0.52 s. These are extension results, separate from the
paper's finite-spring model.

## Why conserving these quantities is insufficient

For zero external force and an incoming rigid translation/rotation velocity
field, expansion of the actual implicit objective gives the exact identity

\[
\frac{\|x-x_n-hv_n\|_M^2}{2h^2}+U_s(x)
=\frac{\|x-x_n\|_M^2}{2h^2}+U_s(x)
-t^TP(x)-\omega L_{\mathrm{ang}}(x)+\frac12\|v_n\|_M^2.
\]

On fixed energy and momentum constraints, every term on the right is fixed.
The objective cannot rank feasible motions with those same conserved values.
At release from rest this reduces directly to objective = endpoint energy.
Adding rho leaves a scalar energy-level preference without restoring a
positional preference on that level. The metric and iteration select a branch.
This explains a structural vulnerability consistent with the observed freeze;
it is not a proof that every CPD trajectory must freeze. See
[OBJECTIVE_DEGENERACY.md](OBJECTIVE_DEGENERACY.md).

There is also an exact feasibility obstruction. A rest-length spring with pure
tangential incoming velocity has E_star=K_star, so rho loses all ability to
change the energy target. Under the angular row and backward velocity, the
minimum possible kinetic-plus-elastic energy is strictly above that target.
The strong fixture has a certified minimum gap 0.001169726 J, above the chosen
1e-4 tolerance, and fails to converge after 200 updates. This affects the
equations themselves; accurate linear solves cannot close incompatible
nonlinear constraints. The separate
[diagnostic](FEASIBILITY_DIAGNOSTIC.md) includes the proof and saved steps.

For the gravity-driven compliant pendulum initially at rest length and at
rest, the initial source energy target is zero, while the gravity torque
requires positive kinetic energy. Exact energy/angular feasibility requires
rho greater than one. At the large tested epsilon this implies an enormous
energy multiplier at stationarity, exposing a tension with the dropped
multiplier curvature. This is a feasible unbounded-rho problem, not the
infeasible unforced fixture above. The analytic
[gravity-target diagnostic](GRAVITY_TARGET_DIAGNOSTIC.md) separates that
initial target issue from the measured later capped-step failures.

## Work bookkeeping

For the prescribed constant external forces, every saved physical step obeys

\[
E_{\mathrm{mech},n+1}-E_{\mathrm{mech},n}
=r_\eta-\rho(E_\star-K_\star)
+\sum_i f_i^T\bigl[h v_{n,i}-(x_i-x_{n,i})\bigr].
\]

The three terms distinguish nonlinear energy-row defect, intentional energy
relaxation, and work approximation. The free opposite-load example's 4.2543%
error occurs with a closed relaxed row: by t=0.4 s, the work difference is
-0.0953563 J, relaxation contributes -4.42e-7 J, and accumulated constraint
defect is near roundoff. The compliant blowup above is dominated instead by
unclosed nonlinear constraints. Both mechanisms must be tracked.

## Validation and next useful ingredient

Independent dense-KKT/gradient/target validation passes 1,521 checks. Eighteen
one-step cases match our unchanged earlier C++ baseline to about 2.8e-16 in
velocity. That reference is our earlier implementation, not authors' code.
The separate saved-data audit covers all 389 primary/control runs, 128,260
trajectory rows, and 13,241 first-step global updates. Free-body CPD dense
replays agree to about 1.1e-13 in multipliers. The hard-rod extension has much
larger discrepancies in ill-conditioned systems, reaching about 0.55 in its
scaled dense-multiplier comparison; its numerical failures are not treated as
well-conditioned evidence against the original finite-spring core.
Hard-rod endpoint stationarity cannot be fully reconstructed from their CSVs'
unsigned geometry-multiplier summary; signed first-step traces are audited.

The useful ingredient to carry forward is angular-momentum control. A focused
next candidate is to add its appropriate constraint to the independent-velocity
LE formulation, then test whether radial motion survives. The present results
do not establish that combination; it has not been implemented in this study.
