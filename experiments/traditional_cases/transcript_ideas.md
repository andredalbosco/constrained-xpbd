# Ideas from the supplied 89-page transcription

Source: [Brouillon-transcription.md](/home/blacknoir/.codex/attachments/54257b07-213f-4872-b696-edd2cc0f43f9/Brouillon-transcription.md). Page numbers below refer to the transcription headings. The original handwritten PDF was not inspected in this pass. The transcription has inconsistent signs, factors and multiplier scalings; equations below are normalized derivations, not verbatim transcriptions. Personal reminders in the source are not task instructions.

## Full research agenda — pages 77–82

Page 81 explicitly lists conditioning, alpha dualization, work balance (free, fixed and moving target), the necessity of the endpoint velocity update, reservoir work and work normalization/numerical drift. The sections below distinguish those directions. The original short list focused too narrowly on immediate EL changes.

### Conditioning is distinct from changing the equations

The primal material Gauss–Newton matrix is `A_p=M/h²+J^T alpha^-1 J`. The material dual matrix is `A_d=J M^-1 J^T+alpha_tilde`. High stiffness and mass contrasts can make these difficult, but neither the condition number nor Gauss–Seidel convergence is determined by stiffness alone; constraint geometry and dependencies matter.

For a positive-diagonal dual system `S delta z=b`, diagonal equilibration uses `D=diag(S)` and solves `D^-1/2 S D^-1/2 u=D^-1/2 b`, with `delta z=D^-1/2 u`. This preserves the solution. Consistently transformed scalar GS visits also map to the same exact-arithmetic visits, so diagonal scaling alone is not a new stabilizing solver.

For a two-row material/energy block `S=[[a,b],[b,d]]`, the equilibrated off-diagonal entry is `rho=b/sqrt(a*d)`. Approaching `|rho|=1` exposes nearly dependent rows; rescaling cannot remove that dependency. A zero diagonal requires separate singular-row handling. Better blocks, relaxation, a different primal metric, or a changed constraint are separate interventions.

### Dualize the inverse compliance, not just a potential-energy expression

For positive compliance the exact identity, in the adopted scaled convention, is

```
U(x) = max_lambda [-lambda^T C(x)/h²
                  -lambda^T alpha_tilde lambda/(2h²)]
lambda = -h² alpha^-1 C(x).
```

This is the baseline XPBD dualization: stiffness appears as compliance in `C+alpha_tilde lambda=0`, rather than explicitly through `alpha^-1` in the material equation. It is different from merely rewriting the value of U in an energy constraint. When deriving work/energy gradients, the question is which additional occurrences of `alpha^-1 C` can be dualized without discarding constitutive residuals or coupling terms. Pages 60–62 develop this concern. The hard-constraint limit is a constrained model; finite-compliance multiplier-energy formulas should not be evaluated by substituting alpha=0.

### Work normalization: recoverable meaning versus a proposed reconstruction

The source explicitly says “Work normalization (numerical drift)” on page 81, but does not give a row/column normalization formula. The following is a plausible reconstruction, not a recovered original derivation.

With rows indexed by particles and columns by potential terms, work balance requires

```
W 1 = delta K
W^T 1 = -delta U.
```

Column balancing makes each potential term's total exported work equal its analytic potential decrease. Row balancing makes each particle's imported work equal its analytic kinetic increase. These are the two analytical energy choices discussed on page 77. A fractional allocation `W_ij=P_ij T_j` with `sum_i P_ij=1` distributes a chosen total work T_j once across a column. Summing actual work entries to one would erase their energy units.

Matching both prescribed marginals requires at least `sum_i delta K_i + sum_j delta U_j=0` on each connected incidence component. Sign restrictions and missing edges can add feasibility conditions. For fixed x and signed W, an alternative to division by small sums is the constrained least-squares repair

```
min_W 1/2 sum_(i,j in support) omega_ij (W_ij-W_hat_ij)²
subject to W 1=delta K, W^T 1=-delta U.
```

This is distinct from normalizing a solver matrix for conditioning. If the totals disagree, no work-only normalization can fix them while preserving both analytic marginals: positions/velocities, relaxed balances or reservoirs must participate. Work is signed, so nonnegative matrix-balancing methods are not a drop-in choice. A single moving particle and single spring without gravity has a 1-by-1 work matrix, where the row/column distinction disappears. Multiple bodies or potential terms expose the allocation question.

### Is endpoint velocity reconstruction necessary?

`v^+=(x^+-x^-)/h` is the backward-Euler kinematic relation, not a general law requiring every integrator to use it. It makes every position correction also a velocity correction of size `delta x/h`, and ties endpoint kinetic energy directly to positions. Alternatives require a stated integration/consistency equation; freeing v with no relation would underdetermine the physical update.

The notes suggest retaining x and v as independent variables and replacing the hard relation with, for example, the penalty `gamma/(2h²) ||x^+-x^- - h v^+||_M²`, or imposing separate position and tangent-velocity constraints for a rigid rod. These are model/integrator changes, rather than new GS orderings. This direction remains shelved; it has not been reactivated as an implementation task.

### Reservoirs store unsatisfied local work balances

For an actual accepted step define

```
bK_i = delta K_i - sum_j W_ij
bU_j = delta U_j + sum_i W_ij
RK_i^+ = RK_i^- - bK_i
RU_j^+ = RU_j^- - bU_j.
```

These reservoirs record signed work defects. With conservative potentials included and no unaccounted external work, summing the updates gives exact bookkeeping conservation of `sum K+sum U+sum RK+sum RU`, up to arithmetic error. This statement does not imply that physical K+U is conserved or that the trajectory is accurate.

There are two separate variants: record reservoirs only, or feed them into subsequent work targets/penalize the prospective reservoir values so that the solver attempts repayment. The latter needs a strength and an acceptance rule; ledger conservation alone does not guarantee stable repayment. This preserves the user's intended definition through work balance rather than adding a global energy reservoir to the present endpoint-energy method.

### Numerical drift: diagnose before correcting

Work quadrature error, linearization error, incomplete GS solves, geometric/tangent inconsistency and floating-point accumulation are distinct sources. Re-evaluating actual accepted-step work defects identifies what reaches the reservoir; linearized residuals alone can hide quadratic defects.

For example, for a scalar spring `U(s)=k*s²/2` with old-force work `W_old=-k*s^-*(s^+-s^-)`,

```
delta U + W_old = k/2 * (s^+-s^-)².
```

Thus a perfectly solved linearized/old-force work equation can still disagree with the analytic potential change. This is discretization error, not roundoff. Using the average spring force instead makes this particular potential-work identity exact, but does not automatically match the chosen kinetic update. Smaller timesteps reduce these local defects without selecting the correct motion branch or enforcing momentum.

These agenda items are analytical interpretations and proposed experiments. No new solver variants or comparative trajectories were run when expanding this report.

## Reconstruct positions before building lambda

Yes. Start at `y = x^- + h v^- + h² a` with cold multipliers. Write `H(x)=E(x)-E_initial`, `g=∇E`, and `alpha_tilde=alpha/h²`. One refreshed EL sweep is

```
delta nu = -H(y)/(g(y)^T M^-1 g(y))
nu^+ = nu^- + delta nu
x_E = y + M^-1 g(y) delta nu

C_E = C(x_E), J_E = J(x_E)
delta lambda = -(C_E + alpha_tilde lambda^-)
               /(J_E M^-1 J_E^T + alpha_tilde)
lambda^+ = lambda^- + delta lambda
x^+ = x_E + M^-1 J_E^T delta lambda.
```

This is the raw, single-material-row formula; multiple material rows are visited sequentially. The bounded experiment filters the whole resulting proposal afterward. Every local equation is linearized. Reconstructing positions alone, without rebuilding residuals/gradients, is the existing frozen-row EL schedule.

The refreshed comparison still drops curvature and `r_x`. At 480 Hz for 10 seconds, bounded EL's radial-spring RMS stays 11.068 cm and the released rigid pendulum stays 0.466 cm. Refreshed LE improves bounded radial-spring RMS from 0.959 to 0.500 cm but worsens the released rigid pendulum from 0.491 to 0.748 cm. Tangential spring motion remains a counterexample. See `el_le_results/summary.csv` and `el_le_plots/` for all results and raw spikes.

## Retain the stationarity residual after rebuilding — page 82

The notes explicitly retain the objective stationarity residual when deriving a material update. In our scaled convention, including the energy multiplier,

```
r_x = M(x-y) - J^T lambda - g nu
r_lambda = C + alpha_tilde lambda
M delta x - J^T delta lambda - g delta nu = -r_x.
```

If a material-only visit holds nu fixed, elimination gives

```
(J M^-1 J^T + alpha_tilde) delta lambda
    = -r_lambda + J M^-1 r_x
delta x = M^-1 (J^T delta lambda - r_x).
```

This retains stationarity while still dropping multiplier-weighted curvature. `r_x` is initially zero at the cold predictor, but generally becomes nonzero after position/gradient updates. If the energy gradient remains frozen after E and lambda is initially zero, the immediate E displacement still gives zero `r_x`; rebuilding all gradients is a separate choice that must be stated. The refreshed tests above do not include this residual-retaining variant.

## Potential-energy representations — pages 46, 49, 53–55, 62

The notes explore position, multiplier and mixed representations of elastic energy. With XPBD's scaled lambda, the corresponding quantities are

```
U_C      = C^T alpha^-1 C / 2
U_lambda = lambda^T alpha lambda / (2 h^4)
U_mix    = -lambda^T C / (2 h²)
U_dual   = -lambda^T C / h² - lambda^T alpha lambda / (2 h^4).
```

They agree when `r_lambda=C+alpha_tilde lambda=0`. In particular,

```
U_C - U_dual = r_lambda^T alpha^-1 r_lambda / 2.
```

This identity applies to soft constraints with positive compliance. Away from the constitutive manifold, replacing analytic potential by a multiplier expression changes the energy equation unless the residual terms are kept. This may be a useful solver design choice, but it is not an algebraically equivalent substitution at a general iterate. Physical-energy validation must continue to use `U_C`. This is an analytic observation; no new representation variant was run here.

## Penalize the work estimate while keeping work balances hard — pages 81–82, 85

The notes distinguish free, fixed and moving work estimates, and propose optimizing the transfers W with a penalty around an estimate. A normalized version is

```
min_(x,W) F_XPBD(x) + kappa/2 sum_(i,j) (W_ij - W_hat_ij(x))²
subject to
    K_i^+ - K_i^- = sum_j W_ij
    U_j^+ - U_j^- = -sum_i W_ij.
```

Gravity must be included as a conservative potential contribution or as an explicit external-work balance. The displayed signs follow work delivered to kinetic energy; the transcription contains a conflicting sign on page 85.

Using kinetic multipliers mu_i and potential multipliers nu_j, stationarity in W gives

```
kappa (W_ij - W_hat_ij) - mu_i + nu_j = 0
W_ij = W_hat_ij + (mu_i - nu_j)/kappa,   kappa > 0.
```

With free work (`kappa=0`), incident multipliers become equal; on a connected incidence graph this recovers the global-energy structure. Finite kappa controls deviation from the work estimate without softening the balance equations. Summing exactly satisfied balances conserves total energy. A finite number of linearized GS visits leaves balance residuals, so exact conservation is not automatic in an implementation.

If W_hat depends on x, its derivatives belong in the primal stationarity and linearization. Freezing the estimate is another explicit approximation. This formulation does not by itself eliminate frozen branches: an estimate that vanishes with displacement can still permit a frozen state. No new work-balance experiment was run in this pass.

## Other useful directions

- **Known kinetic curvature, pages 53 and 56–57.** The kinetic energy Hessian is the inexpensive constant `M/h²`. Retaining this part of multiplier-weighted curvature could affect later visits, but it vanishes from that term at the initial cold `nu=0`; it cannot prevent the first problematic linearized E solve by itself. A sufficiently large multiplier can also destroy positive definiteness.
- **Momentum rows, pages 34, 38, 42 and 48.** With old positions fixed, endpoint momentum and angular momentum are linear in new positions: `P=sum m(x^+-x^-)/h`, `L=sum (m/h) x^- cross x^+`. Their Hessians are zero. Constraining energy corrections to admissible momentum directions is relevant to the orbit failure. Anchors and external gravity require the appropriate reaction/torque targets.
- **Energy slack, pages 2, 9–11 and 85.** An auxiliary slack with a quadratic penalty gives a compliant energy row and a principled strength parameter. This relaxes conservation itself, unlike the soft work-estimate formulation. It is not a new tested method here.
- **Independent velocity, pages 81–82.** The notes consider relaxing `v=(x^+-x^-)/h`; this idea remains shelved as requested.

The immediate follow-up closest to the current EL question is the residual-retaining rebuilt visit. The work-estimate penalty is the larger structural research direction. Neither has been validated by the new refreshed-row experiment.
