# Independent validation

The audit independently assembles the dense linearized equations. It does not
reuse the production closed-form block increments as the reference.

## Coverage

- Forty-eight arbitrary spring/rod states at beta_v 0.1, 1, and 10, with nonzero
  position and velocity stationarity residuals and independently carried material,
  energy, and relaxation variables.
- Finite differences of the original implicit objective plus velocity penalty,
  physical mechanical energy, and full residual Jacobian. The latter includes the
  material/energy curvature omitted by the numerical approximation.
- Full positive coupling Hessian and its analytic inverse; material, energy, and
  joint production increments versus independently assembled dense solves.
- Frozen EL and LE iterations versus the same joint linear solution, retaining
  position and velocity stationarity and pushing both primal variables immediately.
- Explicit velocity-coupling multiplier equations, with arbitrary inconsistent
  current multipliers and the full consistency residual. Elimination reproduces
  the primal formulation, including beta_v=1.
- Six hundred thirty refreshed nonlinear visits across all six scenes, both
  orders and JOINT, all three positive beta values, and one/two/four sweeps.
  Every trace state checks refresh timing and accumulation. Returned velocity is
  compared with the final independent velocity iterate, not overwritten from
  displacement.
- Twelve beta-limit states: material position updates equal the old position-only
  equations; the energy block approaches them as beta tends to zero. The explicit
  dual formulation at beta exactly zero reproduces the hard position-only joint
  system. The primal beta-zero penalty remains undefined and is not used.
- Exact cancellation of the reconstructed position energy gradient still leaves
  the beta_v-weighted independent-velocity contribution in the energy denominator.
  Compatible and incompatible zero energy rows are handled explicitly.

## Interpretation

These checks validate the stated adapted objective, curvature approximation,
block equations, and their implementation. They do not establish stable or
physically plausible long-time motion. Frozen GS convergence is distinct from
convergence of the refreshed nonlinear iteration. Independent velocity can also
permit displacement inconsistent with physical velocity; that requires separate
trajectory diagnostics even when actual mechanical energy is accurate.

Run with a NumPy environment:

```bash
OPENBLAS_NUM_THREADS=1 python3 validate.py --write-report
```

## Numerical receipt

```json
{
  "checks": 6575,
  "details": {
    "arbitrary_states": 48,
    "beta_limit_states": 12,
    "max_frozen_GS_factor": 0.7474658547711234,
    "production_arbitrary_states": 48,
    "refreshed_replay_visits": 630,
    "returned_velocities_distinct_from_displacement": 153
  },
  "largest_relative_errors": {
    "compatible zero energy row stays at rest": 0.0,
    "coupling metric inverse": 2.8047926103327004e-14,
    "energy beta tends to zero recovers derived velocity increment": 1.4551941712004252e-09,
    "energy beta tends to zero recovers position-only block": 7.072828576234328e-10,
    "exact hard dual coupling equals position-only joint system": 2.3424017611535127e-14,
    "exact hard dual velocity increment": 1.4602703977091186e-15,
    "exact reconstructed energy-gradient cancellation": 8.881784197001252e-16,
    "frozen EL versus joint solution": 2.9969892331222632e-09,
    "frozen LE versus joint solution": 2.7407956318289737e-09,
    "frozen block rows after immediate primal push": 1.4256302878705618e-13,
    "full residual finite differences including dropped curvature": 1.0312780380044036e-10,
    "independent multiplier Schur solve": 2.961654155177619e-13,
    "independent velocity leaves positive energy Schur at cancellation": 0.0,
    "material position response equals position-only equations": 2.235668078938117e-16,
    "original implicit plus coupling objective gradient finite differences": 1.0552004898457868e-10,
    "physical energy finite differences": 1.0132964205854872e-10,
    "post-visit physical residual trace": 1.3233195198391836e-12,
    "primal versus explicit dual E block": 1.0425847478334132e-14,
    "primal versus explicit dual L block": 3.3950685574966545e-15,
    "primal versus explicit dual coupling with residual": 4.891315622907788e-15,
    "production E restores included linearized rows": 7.96650130081304e-13,
    "production E versus independently assembled dense block": 8.424755691092958e-15,
    "production J restores included linearized rows": 7.109069885348552e-14,
    "production J versus independently assembled dense block": 9.563153573581622e-16,
    "production L restores included linearized rows": 7.744915819785092e-13,
    "production L versus independently assembled dense block": 3.3546185681453584e-15,
    "production analytical residuals": 7.49962430217668e-17,
    "production compliant velocity mismatch": 0.0,
    "production full coupling Hessian": 0.0,
    "production inverse application": 7.594268044542055e-15,
    "production physical energy gradient": 0.0,
    "production physical independent-velocity energy": 4.440892098500626e-16,
    "refreshed dense replay of x v and accumulated scalars": 1.722242645760812e-14,
    "returned independent velocity matches last primal iterate": 1.722242645760812e-14
  },
  "status": "passed"
}
```
