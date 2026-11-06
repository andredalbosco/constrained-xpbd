# Independent validation

The audit uses independently assembled dense linear systems rather than a
second copy of the production closed-form increment formulas.

## Coverage

- Sixty deterministic random spring states: physical energy and potential
  derivatives, the exact residual Jacobian, elimination of the carried material
  multiplier including nonzero `r_lambda`, and frozen EL/LE versus a full joint
  solve.
- Eighty arbitrary spring/rod states, with nonzero accumulated multipliers and
  relaxation variables: production residuals, both full blocks, and the joint
  increment versus independent dense solves. Gravity contributes to the
  objective outside the energy multiplier as well as to physical energy.
- Full refreshed replay of 504 visits across all six scenes, both orders,
  epsilon 1/100/10000, and one/two/four sweeps. Every post-visit state and residual
  is compared, checking multiplier accumulation, immediate position response,
  refresh timing, block ordering, and derived velocity.
- Gravity-adapted energy floors are checked against independently derived
  potential minima. Compatible zero energy rows are flagged; incompatible rows
  fail explicitly. Equilibrium stays at rest in both orders. A zero-sweep budget
  is explicitly rejected by the public API.

## Interpretation

For a frozen spring linearization the reduced two-multiplier system has positive
diagonals and the GS factor is the squared off-diagonal divided by their product.
The audited random states all have factor below one; both orders converge to the
same dense solution. This does not imply convergence of the refreshed nonlinear
iteration or stability of the physical rollout. Finite sweep budgets can retain
nonzero residuals, and the implementation reports that separately from a finite
accepted step.

The checks validate the stated adapted equations and their implementation. They
do not establish equivalence to the original CPD authors' software, optimal
motion selection, or reliable long-time stability. The hard rod is an extension
with alpha zero, not a claim about the original CPD spring model.

Run from this directory:

```bash
OPENBLAS_NUM_THREADS=1 /tmp/xpbd-nonlinear-venv/bin/python validate.py --write-report
```

## Numerical receipt

```json
{
  "checks": 11600,
  "details": {
    "max_frozen_GS_contraction": 0.9648474997138048,
    "production_arbitrary_states": 80,
    "random_states": 60,
    "replayed_nonlinear_visits": 504
  },
  "largest_relative_errors": {
    "compatible zero energy row returns zero increment": 0.0,
    "dualization elimination including r_lambda": 2.1185341966557947e-15,
    "energy holds lambda": 0.0,
    "energy restores rho-eta equation": 1.2286172079711832e-09,
    "equilibrium stays at rest": 0.0,
    "final state versus dense refreshed block replay": 1.3680052986502741e-14,
    "frozen EL versus joint solution": 2.9727241455052434e-09,
    "frozen GS exact spectral radius": 1.1102230246251565e-16,
    "frozen LE versus joint solution": 2.946918402699208e-09,
    "frozen block equations after immediate push": 8.038875305077471e-14,
    "full residual finite differences": 5.568875674312442e-11,
    "gravity-adapted reference floor formula": 0.0,
    "independent multiplier Schur solve": 2.6391477244039755e-13,
    "material holds rho and eta": 0.0,
    "physical energy finite differences": 8.59455317284281e-11,
    "physical potential Hessian finite differences": 7.266054220523747e-11,
    "post-visit physical energy trace": 4.9966888363234e-16,
    "post-visit refreshed residual trace": 4.266955553179733e-14,
    "production E immediate included-row residual": 9.31322574625643e-10,
    "production E increment versus independent dense block": 6.62163922884977e-14,
    "production L immediate included-row residual": 7.105427357601002e-14,
    "production L increment versus independent dense block": 1.0187446402380924e-15,
    "production analytical residuals": 6.53606562910829e-17,
    "production energy gradient": 0.0,
    "production joint increment versus independent full matrix": 4.02794002382813e-16,
    "production normal": 0.0,
    "production physical energy": 2.3562917652941024e-16,
    "production potential gradient": 0.0,
    "reference floor attained at physical potential minimum": 0.0,
    "refreshed dense-block replay versus trace": 1.3680052986502741e-14,
    "velocity derived from displacement": 0.0
  },
  "status": "passed"
}
```
