# Original CPD: source audit and independent validation plan

Primary source read: Kee et al., *Constrained Projective Dynamics* (2021),
sections 4–5, equations 10–21, and the attachment discussion in section 4.4.
The archived text is
`/home/blacknoir/Projects/dynamics-research/archive/resources/papers/Constrained Projective Dynamics.md`.
Implementation inspected read-only:
`/home/blacknoir/Projects/dynamics-research/experiments/work_balance/baselines.cpp`,
`baselines.hpp`, and `docs/work_balance/baselines.md`.

## Source-faithful core

- Unknowns are positions and one energy-relaxation scalar. Endpoint velocity is
  reconstructed from displacement; no independent velocity or velocity penalty.
- The objective is the implicit inertial term plus elastic potential and
  epsilon times relaxation squared divided by two. Gravity enters the predictor.
- Full free-body constraints are linear momentum, angular momentum, and relaxed
  internal energy. In a planar reduction there are four scalar rows; in general
  three-dimensional scenes there are seven before eliminating dependent rows.
- The position metric is the PD matrix, mass over timestep squared plus the
  assembled spring incidence stiffness. A radial spring contributes k times the
  identity, not its rank-one Gauss–Newton normal outer product.
- Equation 21 solves the conservation multiplier afresh for each quadratic
  subproblem. It is not an accumulated XPBD multiplier increment. Position and
  relaxation are then updated jointly with a full step.
- The current objective gradient is retained. The paper drops multiplier and
  nonlinear material curvature in the matrix approximation, rather than
  discarding the objective gradient.
- Paper external-energy transfer is first-order power at the incoming velocity.
  The constrained energy includes kinetic and elastic energy; gravitational
  potential belongs in the separate physical-energy diagnostic.

## Differences already present in the older local baseline

The local baseline implements the correct coupled equation-21 algebra and PD
matrix. Its following choices must be disclosed rather than attributed to the
original authors' implementation:

| Local choice | Consequence |
|---|---|
| Disable both momentum families whenever any pin exists | Pinned scenes become an energy-only adaptation. |
| Incoming center-of-mass origin | Conditions free-body angular rows; targets and inertia must use the same origin. |
| Correct torque from printed force-cross-position to position-cross-force | Matches the paper's angular-momentum definition; the printed equation 12 has the opposite sign. |
| COD pseudoinverse/least-squares Schur solve at threshold 1e-12 | Handles rank loss, but inconsistent linearized rows need not be satisfied. This is not an unconditional exact inverse. |
| Epsilon default 1e8 | Local numerical choice; no numerical epsilon default was found in the archived paper text. |
| Endpoint KKT stationarity plus scaled constraints, tolerance 1e-9 | Stronger/different stopping from the paper's conservation-error example 1e-4. |
| Reject exact zero compliance for PD/CPD | Correctly avoids claiming an unspecified hard-rod extension is original PD. |

## Anchored versus free-body physics

A single free mass attached to a fixed anchor exchanges linear momentum with the
support. Conserving that mass's linear momentum would remove its spring response.
Its angular momentum about the anchor *is* conserved for a central spring without
gravity. With gravity, torque about that anchor is position cross gravitational
force. The old blanket pin policy is therefore deliberately more restrictive than
the symmetry of this particular anchored scene.

Using the free mass's own incoming COM as angular origin makes its angular row
identically zero. An anchored angular test requires the anchor origin. A single
free particle with both linear and angular rows also has dependent rows; the
angular row is a linear combination of the two planar momentum rows.

A free two-mass spring with zero COM position/momentum is the clean full-CPD
control. Its internal forces cannot change total P or L. Planar P_x/P_y/L_z/E
rows avoid the identically zero angular rows of a two-particle 3D embedding.

The paper describes finite attachment springs and shows a pendulum, but does not
give an exact hard-distance/tangency algorithm for that scene. A stiff compliant
spring pendulum, prescribed-pin energy-only test, and exact hard rod must retain
distinct labels. Support-force target bookkeeping cannot be silently replaced
with conservation of the free subsystem's momentum.

## Independent numerical checks

1. Finite differences of elastic/objective gradients and every conservation row,
   with off-solution positions, nonzero incoming velocities, gravity, and an
   origin held fixed during differentiation.
2. Frozen spring projections give the assembled PD stiffness; refreshed
   projections give the true radial elastic gradient. Check the distinction
   against a tangential perturbation where the GN stiffness would vanish.
3. Paper Schur equation 21 versus a dense full KKT solve, including relaxation,
   objective gradient, and fresh multipliers. Independently verify both block
   residuals after the correction.
4. Translation/origin consistency, incoming-velocity targets, inertia-derived
   kinetic reference, and first-order power versus exact gravitational work.
5. Rank tests: independent planar rows, redundant single-particle momentum rows,
   and inconsistent zero energy rows. Any reduced-rank solve must identify rank
   and remaining compatibility residual. No hidden diagonal regularizer or
   unreported least-squares acceptance.
6. Release from rest: with no external force, the objective equals the constrained
   kinetic-plus-elastic energy. Full P/L constraints do not remove the stationary
   near-frozen branch. At zero momentum and kinetic reference, finite epsilon
   permits relaxation E_target/epsilon even though strict energy is feasible.
   Conservation and numerical convergence alone do not certify useful motion.

Source fidelity, algebraic correctness, and motion quality are separate outcomes.
No claim of parity with the authors' unavailable source code is made.
