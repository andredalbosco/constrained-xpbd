# Original CPD: source contract and boundaries of the comparison

The primary reference is Kee, Um, Jeong, and Han (2021),
[Constrained Projective Dynamics](https://doi.org/10.1145/3450626.3459878),
especially Eqs. (7), (11)–(21), Algorithm 1, and Sections 4.4–6. The full local
paper is
[archived here](</home/blacknoir/Projects/dynamics-research/archive/resources/papers/Constrained Projective Dynamics.md>).
This audit also inspected the existing
[baseline contract](/home/blacknoir/Projects/dynamics-research/docs/work_balance/baselines.md),
[C++ implementation](/home/blacknoir/Projects/dynamics-research/experiments/work_balance/baselines.cpp),
[options](/home/blacknoir/Projects/dynamics-research/experiments/work_balance/baselines.hpp),
and the historical
[energy-only specialization](/home/blacknoir/Projects/constrained-xpbd/experiments/nonlinear_energy_2026_11_01/CPD_DERIVATION.md).

## What must be present to call the core solver original CPD

| Component | Source-faithful contract |
|---|---|
| State unknowns | Endpoint positions and one energy-relaxation variable rho |
| Endpoint velocity | Exactly (x-x_n)/h; no independent velocity unknown |
| Objective | Original implicit PD inertia plus elastic energy, with epsilon rho²/2 |
| Elastic treatment | PD local projections refreshed per iteration; constant assembled PD matrix |
| Conservation rows | Total linear momentum, total angular momentum, relaxed internal total energy |
| Energy accounting | Kinetic plus elastic energy; gravity is external force and old-velocity power in the target |
| Global step | Coupled Eq. (21) SQP/Schur solve of all active conservation rows |
| Multiplier semantics | Fresh subproblem multiplier value, not an accumulated increment |
| Step acceptance | Full step as written; no line search, damping, rescaling, or post-projection |
| Energy slack | Unbounded optimized rho; no implicit [0,1] clamp |

In the planar version, there are two linear-momentum rows, one angular row,
and one energy row. This is the planar reduction of the paper's seven 3D rows.
Keeping zero or dependent 3D angular rows in a planar problem and blindly
inverting the Schur matrix would introduce an avoidable rank problem.

An energy-only Gauss–Newton JOINT solver is not original CPD. It differs in both
the elastic matrix and the momentum constraints. The earlier pinned energy-only
PD specialization did reproduce a restriction of Eq. (21), but its missing
angular constraint prevents using its spinup as evidence against full CPD.

## Parameter values: published versus locally selected

The paper's reported experiments use h=1/30 s and a conservation-error threshold
of 1e-4. It describes an iteration cap but does not supply one universal numeric
cap in the inspected text. The numeric epsilon value is also not supplied
there. The default ten iterations mentioned in Section 5 belongs to the PD and
implicit-midpoint comparison solvers, not a declared ten-iteration CPD cap.

The old local baseline uses epsilon=1e8, cap 200, and a dimensionless 1e-9
stationarity-plus-constraint tolerance. These are **local comparison settings**,
not established author defaults. Epsilon has energy units; a unit-dependent
value must be disclosed rather than called universally standard.

The paper stops on conservation residual, whereas the old local baseline also
checks endpoint stationarity. For a faithful study, record the chosen stopping
rule, its norm and units, the cap, and full stationarity independently. A raw
1e-4 criterion on mixed momentum/energy rows is not unit invariant. The paper
does not fully specify how its initial feasible iterate is treated by the
convergence check. The old code executes at least one update before checking;
that matters for release-from-rest fixtures whose initial conservation residual
is already zero. Preserve and declare the policy.

The paper treats its small Schur matrix as invertible. Actual degenerate
configurations require explicit rank handling. The old baseline's normalized
rank-revealing solve at relative threshold 1e-12 is a numerical robustness
addition; any incompatible residual must remain visible. Do not silently drop
momentum rows or add diagonal regularization and call that the same update.

## Gravity and the printed torque sign

Original CPD puts external force in the predictor and updates its internal
kinetic-plus-elastic energy target with h times force dot **incoming** velocity.
It does not constrain total energy including gravitational potential exactly.
Report both the source energy equation and actual mechanical energy including
gravity. Small source residual cannot establish exact mechanical-energy balance.

Equation (10) defines angular momentum using position cross momentum, but
Eq. (12) prints force cross position in the torque update. Those have opposite
signs. The existing baseline explicitly uses the physically consistent position
cross force. This is a declared correction to the printed equation, not a
verbatim reproduction of that sign.

The paper writes an inertia inverse evaluated from the incoming positions.
Its presentation does not fully specify the center-of-mass bookkeeping needed
to interpret translational and rotational kinetic energies independently.
The local baseline uses a fixed incoming COM origin per timestep and transforms
the angular target consistently. In 2D the inertia is one positive scalar;
in 3D use a pseudoinverse only with an explicit range/feasibility check.

## Attachments, exact pins, and hard rods

Section 4.4 describes finite zero-rest-length springs for attachments and
repulsion springs for collisions. Algorithm 1 updates the conservation targets
with the associated external-force accounting. It does not provide an exact
hard-pin reaction solve or an exact rod-length-plus-tangent-velocity algorithm.
The local Markdown conversion truncates part of Section 4.4's force-update
prose, so implementation details beyond the stated finite-spring model should
be checked in the original PDF or author code rather than inferred.

Exact pins generally exchange momentum with their support. Constraining the
free particles to their old total linear momentum is therefore wrong unless
the support impulse is included. The existing local baseline disables both
momentum families whenever a pin is present. That is its explicit adaptation;
it is overly broad for a single central anchor, whose angular momentum about
that anchor is still physically conserved in the absence of external torque.

For a single central anchor at the origin:

- Drop linear-momentum conservation, because the anchor exchanges linear impulse.
- Retain angular momentum about the anchor. Its support torque is zero.
- With gravity, add the known old-position gravity torque to the angular target.
- A general off-origin/multiple-support configuration needs its actual support
  torque bookkeeping; the previous simplification does not apply automatically.

Adding an exact rod-length row is a separately labeled hard-rod extension of
the CPD numerical solve. Its velocity remains a chord, so endpoint tangency is
not exact. Enforcing exact length at both endpoints, hard chord kinematics,
and exact endpoint tangency simultaneously would force zero displacement.
Changing to a finite stiff spring is another mechanical model and must be
labeled, rather than quietly substituting stiffness for a hard rod.

## Recommended original-core test scenes

The cleanest original-core comparison uses a **free two-particle spring** with
no pins. Choose masses 2 and 2, put particles at +r/2 and -r/2, and assign
velocities +u/2 and -u/2. Its reduced mass is one. Consequently its relative
coordinate has the same spring equation, internal energy, and angular momentum
as the previous mass-one anchored spring with position r and velocity u.
The full planar P/L/E rows are valid and both linear and angular momenta can
be enforced without a guessed support impulse.

Reuse the release, orbit, and stiff-orbit relative initial states. For a
gravity-relative counterpart, use known **opposite constant external loads**:
force +g on the first mass-two particle and -g on the second. The relative
acceleration is g, net external force is zero, torque is r cross g, and incoming
power is g dot u. This reproduces the old anchored relative problem while
remaining a free-body system with fully specified external forces. It is not
uniform gravity and must be labeled accordingly. Uniform gravity on the free
pair is a separate source-target sanity check; it accelerates COM and cancels
from the relative equation.

| Test family | Label and purpose |
|---|---|
| Free two-mass finite spring, P/L/E | Original CPD core; clean test of spinup, release, turns, and external-power accounting |
| Single central anchored spring, angular/energy | Physically justified anchored restriction; separate from original free-body model |
| Single central anchored spring with gravity torque target | Anchored restriction with external-power and torque approximation |
| Exact hard pendulum with angular/energy plus length | Explicit hard-rod extension; not an authored hard-rod implementation |
| Energy-only PD solve | Existing restriction/control, useful for isolating the effect of angular conservation |

## Selected experiment contract

The authorized comparison uses 13 cases: four free-pair relative fixtures
(release, orbit, stiff orbit, and opposite-load gravity counterpart), one free
triangle, four central anchored spring fixtures, two compliant pendulums with
finite stiffness 10,000, and two explicitly labeled exact hard-rod extensions.
The finite-stiffness pendulums are native elastic PD models, not exact rods.

The three selected solver-budget policies are:

| Policy | Chosen rule |
|---|---|
| Paper-style post-update | At least one global Eq. (21) update; then raw maximum conservation residual below 1e-4; cap 200 |
| Fixed budget | Exactly 16 global updates |
| Stationarity diagnostic | Full endpoint KKT and constraint residual below declared scaled 1e-9 tolerance; cap 200 |

The minimum-one-update interpretation, infinity norm, and cap 200 are explicit
local implementation choices. A separate initial-precheck release variant
exposes the zero-iteration feasible freeze; it does not silently replace the
primary policy. The strict stationarity policy is a diagnostic enhancement,
not the paper's published stopping criterion.

Epsilon values are reported as selected experiment parameters, not author
defaults. The new core solve reports eigenvalue-based rank/conditioning of its
normalized Schur matrix, using a 1e-12 relative diagnostic threshold. That rank
estimate does **not** reject a step by itself: the solver attempts a direct
scaled solve of every active row. An actual singular solve, nonfinite result,
or declared geometry/runaway guard is an explicit failure. There is no COD or
pseudoinverse least-squares fallback, conservation-row deletion, or diagonal
regularization. The old code remains a separately instrumented historical
reference rather than a claim of byte-for-byte author implementation parity.

## What the comparison can establish

Full momentum constraints should prevent unforced self-spinup **to their actual
residual accuracy**. That addresses a real missing ingredient of the earlier
energy-only experiments. It does not prove that radial motion remains natural,
that energy-only turning artifacts disappear, or that release cannot freeze.

At zero individual incoming velocities and zero external force, the original
position objective equals its constrained endpoint energy. The release
degeneracy therefore survives even when the momentum rows are restored. Large
epsilon may select a nearly frozen branch. This is a prediction to test, not a
claim that the published examples have been reproduced and found deficient.
The paper's zero-total-momenta example does not by itself specify every
particle's incoming velocity or all initialization details.

For the user's current criterion, record actual energy drift, accepted motion,
stored energy versus radial/angular kinetic energy, sustained non-equilibrium
stalls, turning-point velocity changes, and angular-momentum drift. Do not use
phase or frequency error as a rejection criterion. The paper itself disclaims
physical accuracy beyond conservation and notes limitations in periodic elastic
motion and elastic waves (Section 6).

## Read-only implementation audit

Inspection of `cpd_solver.py` matches the declared predictor, PD spring metric,
fresh Eq. (21) multiplier solve, analytical energy gradient, planar momentum
rows, central-anchor angular target, and backward velocity reconstruction.
Its normalized Schur rank estimate is diagnostic and does not reject a step
by itself. This is source-equation consistency, not author implementation
parity or an empirical validation claim.

The optional angular-row ablation holds K_star unchanged while removing only
that row, deliberately isolating the effect of enforcing angular momentum.
That ablation's retained reference is not the kinetic minimum of its reduced
constraint set. It is a controlled variant and should not be described as the
original full momentum formulation or the historical energy-only K_star=0
restriction.

The standalone `FEASIBILITY_DIAGNOSTIC.md` shows a separate exact obstruction
when an unforced spring starts at rest length with tangential velocity, making
E_star=K_star. It is not part of the primary 13-scene grid. Its proof and actual
solver results are reported separately from the main motion comparison.
