# Implicit candidates with a stability requirement

This continues the spring, rigid pendulum and double spring pendulum study using the backward-Euler/XPBD objective and endpoint velocity reconstruction throughout. No midpoint kinematics, velocity rescaling or adaptive substeps are used. All multiplier solves are scalar Gauss–Seidel on linearized equations with multiplier-weighted curvature dropped. The primal residual is retained.

The best common setting in this measured sweep is **compliant local work balance + angular balance + repayment of local work reservoirs**, named `workL_s1_tau005`. It uses work-balance strength 1, work-estimate strength 1, and repayment time 0.05 seconds. It improves position RMS over ordinary XPBD in all 21 nominal case/frequency comparisons. This is an in-sweep result, not independent evidence of a universally optimal setting.

Stability is enforced on accepted states with a physical-energy ceiling. Candidates that cannot advance under that ceiling are excluded, rather than allowed to produce an unstable state. Completion is checked separately. The nominated setting completed every nominal, long-run and stiffness/perturbation run, but physical energy conservation is not yet achieved.

## Scope and measured outcome

- 798 ten-second runs: 38 configurations × 7 initial conditions × 240/480/960 Hz.
- 147 sixty-second nominal runs: 21 configurations × 7 cases at 480 Hz.
- 2,016 one-second stiffness/perturbation runs and the same 2,016 scenarios extended to ten seconds.
- Twenty-one additional sixty-second probes of individual reservoir magnitudes and cumulative bookkeeping error.

The seven initial conditions are spring release, spring orbit, rigid pendulum from rest, kicked rigid pendulum, and small/large/mixed two-spring motion. Masses, stiffness, geometry and refined RK4 references are unchanged from `traditional_cases`. There are at most two moving masses.

At 480 Hz over ten seconds:

| Case | XPBD position RMS | Nominated implicit candidate RMS | XPBD maximum physical energy error | Candidate maximum physical energy error |
|---|---:|---:|---:|---:|
| Spring release | 7.416 cm | 0.490 cm | 87.54% | 3.830% |
| Rigid pendulum | 3.765 cm | 1.556 cm | 21.99% | 9.011% |
| Kicked pendulum | 2.937 cm | 1.305 cm | 21.02% | 9.399% |
| Double spring, small swing | 1.176 cm | 0.141 cm | 41.58% | 1.184% |
| Double spring, large swing | 14.498 cm | 11.503 cm | 63.43% | 9.917% |
| Spring orbit | 5.790 cm | 0.525 cm | 47.46% | 4.562% |
| Double spring, mixed motion | 7.644 cm | 2.714 cm | 46.97% | 6.565% |

Energy is analytic physical kinetic + spring + gravitational energy, excluding reservoirs. Errors are normalized by initial energy above minimum potential. Position RMS includes both Cartesian coordinates and both masses for the double spring. The positive normalization is fixed per run.

The candidate's angular-balance defect at 480 Hz is approximately 3.9e-13 kg m²/s for the orbit, 3.2e-12 for the large double spring, 1.4e-4 for the small double spring and 2.4e-6 for the rest pendulum. The final two retain finite defects from limiting/projection and finite GS; exact conservation is not claimed.

## The implicit formulations

Use old physical state `x^-`, `v^-`, timestep h, gravity a, and ordinary predictor

```
y = x^- + h v^- + h² a
v^+ = (x^+ - x^-)/h
F(x) = (x-y)^T (M/h²) (x-y)/2 + C(x)^T alpha^-1 C(x)/2.
```

Rigid rods use a hard material equation rather than an elastic energy term. Spring energies and the objective use physical alpha. `alpha_tilde=alpha/h²` appears only after dualization. Multipliers are cold at every physical timestep and accumulate only within that step: `lambda^+=lambda^-+delta lambda`, and likewise for the other multipliers.

### Compliant global endpoint energy

The energy family minimizes

```
F(x) + gamma_E/(2 E_scale) [E(x)-E_target]²,
```

with optional angular balance. Most configurations target initial physical energy; a separate control targets the previous step's energy. This is a compliant penalty, not an exact energy constraint. After scaling the objective by h², the dual energy compliance is

```
beta_E = E_scale/(h² gamma_E)
r_nu = E(x)-E_target + beta_E nu.
```

EL and LE visit energy/material rows in the specified order. `GS4` repeats four visits to the same frozen system. `outer4` performs four relinearizations with one GS sweep each. Strengths 0.1, 1, 10 and 100, and a stiffness-gated strength `gamma_E/(1+sigma)` are tested, where `sigma=max_j h² J_j M^-1 J_j^T / alpha_j` over soft springs.

Positive compliance removes the hard energy row but does not make the nonlinear objective globally convex. Some energy-only variants become unable to advance under the energy ceiling; they are not stable candidates. Adding angular balance fixes much of the orbit distortion, but energy-only corrections still do not provide a common improvement on the coupled spring system.

### Compliant local work balance

For every supported particle/potential incidence, write the independent work as

```
W_ij = W_hat_ij(x) + t_ij
W_hat_ij(x) = f_ij^- dot (x_i-x_i^-).
```

The force estimate uses the preceding physical state. Spring forces are analytic there; gravity is constant. The rigid rod uses an estimated reaction force based on speed, configuration and gravity. That reaction is a work estimate, not an exact discrete reaction law.

Define

```
bK_i = K_i(x)-K_i^- - sum_j W_ij
bU_j = U_j(x)-U_j^- + sum_i W_ij
p = 1-exp(-h/tau), or zero without repayment.
```

The actual primal objective is

```
F(x)
 + gamma_W/(2 E_scale) sum_ij t_ij²
 + gamma_B/(2 E_scale) [sum_i (bK_i-p RK_i^-)²
                       +sum_j (bU_j-p RU_j^-)²],
```

with a hard rod equation and optional angular-balance equation. Gravity has separate conservative potential nodes. A rod has zero potential energy. Both work balances and disagreement with the force-work estimate are compliant. The actual work variables remain independent through t; changing the work estimate is not silently identified with changing physical energy.

The parameter grid varies `gamma_B=1/10/100` and `tau=0.05/0.25/1` seconds at `gamma_W=1`, with angular balance. The earlier `gamma_B=10,tau=0.25` setting is retained as `work_repay_angular`. Three further controls soften angular balance at the nominated work/repayment settings.

### Stationarity, linearization, dualization and GS

First write the primal penalty problem and its stationarity equations. Linearize each C/balance at the current iterate, drop residual-weighted curvature, and dualize the resulting quadratic penalties. In the h²-scaled system,

```
beta_B = E_scale/(h² gamma_B)
kappa_t = h² gamma_W/E_scale
r_lambda = C + alpha_tilde lambda
r_mu_i = bK_i - p RK_i^- + beta_B mu_i
r_nu_j = bU_j - p RU_j^- + beta_B nu_j.
```

The position and t stationarity residuals are retained. The remaining primal metric is diagonal: M for positions and `kappa_t I` for work deviations. Writing all row gradients as G, primal residual as `r_p`, and row residuals as r, elimination gives

```
S = G H^-1 G^T + compliance
S delta z = -r + G H^-1 r_p
delta primal = H^-1 (G^T delta z - r_p).
```

One GS sweep visits material lambda, potential nu, kinetic mu and then angular multiplier zeta. The nominated setting uses one frozen linearization and one sweep per small step. The local quadratic model has a positive primal metric and positive material/work penalty compliances, with hard angular/rod rows as specified. Finite GS is an approximate solve; no claim is made that every accepted step minimizes the full nonlinear objective.

The work estimate depends on current positions through displacement, while K/U are nonlinear in positions. All first derivatives are included using the already tested small forward-mode derivative type. The material potential in those derivatives is analytic; replacing `alpha^-1 C` with a multiplier while dropping its constitutive residual is not used.

## Angular balance and why the pendulum is special

The added angular row accounts for gravity torque about the anchor:

```
ell(x) = L(x)-L^- - h/2 [T(x)+T^-] = 0
grad ell_i = (m_i/h) tangent(x_i^-) + (h m_i/2) tangent(a).
```

It is linear in new positions under endpoint velocity reconstruction, so its curvature is exactly zero. Internal spring forces are central, and the anchor is at the origin. Global linear momentum is not constant because of gravity and anchor reactions; this study does not impose an incorrect free-system linear-momentum target.

Without the angular row, compliant work/repayment can preserve energy much better while creating a bad orbit. At 480 Hz, work repayment without angular control has 91.39 cm orbit RMS; adding angular balance with the nominated settings reduces that to 0.525 cm.

For one rigid pendulum, length and the chosen angular row provide two equations for two position coordinates. On a local branch they largely determine the next position, leaving no independent position direction for work-based energy repair. The measured result is continued physical damping stored in the reservoir: rest-pendulum energy error reaches about 44.1% at 60 seconds, despite near-exact angular balance. This conclusion is specific to this discrete angular target and endpoint kinematics; it is not an impossibility theorem for all integrators.

Angular compliance adds `beta_L zeta` to its row, with `beta_L=softness * ell_gradient^T M^-1 ell_gradient`. Softness 0.1/1/10 is tested. Softening helps the rigid pendulum but worsens orbit behavior and some coupled cases. At 480 Hz, softness 0.1 lowers rest-pendulum RMS from 1.556 to 0.808 cm, but raises orbit RMS from 0.525 to 10.889 cm. It is not a common replacement for hard angular balance.

## Work reservoirs and actual energy

After limiting and accepting a proposal, recompute the actual work balances and update

```
RK_i^+ = RK_i^- - bK_i(actual)
RU_j^+ = RU_j^- - bU_j(actual).
```

This conserves physical energy plus the sum of reservoirs as bookkeeping. It does not mean physical energy itself is conserved. Repayment targets reduce the old local reservoir through work balances; it is not a global energy-reservoir constraint.

Sixty-second probes also measure the absolute sum of individual reservoirs, since a small signed sum can hide large opposite values. The nominated setting has approximately steady local reservoir magnitude for spring release/orbit and the small/mixed double spring over the measured interval. Its small-double local absolute sum stays around 7% of initial excitation energy, while physical energy error is around 1.2%. Rigid-pendulum reservoirs continue to grow as physical motion damps; this is a remaining limitation, not hidden conservation success. Probe cumulative ledger error is below 1.4e-12 of the energy scale for all three tested settings.

## Accepted-state stability policy and screening

Every candidate first computes a material XPBD baseline. Extra motion is limited to at most 10% of `h ||v^-||_M + h² ||a^-||_M`. Rods receive a final ordinary linearized material correction. Actual displacement is checked after that correction. Angular variants may not worsen the local angular residual relative to the baseline, within a floating-point tolerance.

Accept only finite states with feasible rods and

```
E(x_candidate) <= E_initial + 1e-10 E_scale.
```

If needed, halve the extra correction toward XPBD, up to twenty trials. If no admissible correction exists, try the XPBD baseline. If that is also inadmissible, return `cannot_advance_safely`, without changing the physical state or reservoirs. No adaptive timestep subdivision is used. Blending limits the correction, so exact satisfaction of proposed multiplier equations is not implied.

For these anchored spring potentials, physical energy above the potential minimum controls lengths and velocities; rigid length is checked explicitly. The policy therefore bounds accepted physical states. It does **not** prove that every problem can advance within the solve budget. Safe completion and useful motion are separately measured.

The stress grid has 96 scenarios per selected method and duration: release/orbit/rigid/large double spring; 240/960 Hz; stiffness multipliers 1/100/10,000; light-mass ratios 1/100; original or rotated by 1e-8 radians.

- The nominated setting completed all 96 one-second and all 96 ten-second stress runs, with no committed state above initial energy. It also completed all 21 short and seven long nominal runs. It used no fallback-only steps in those runs, although it sometimes reduced the extra correction.
- Energy-only strength-10 EL completed 90/96 one-second and 86/96 ten-second stress scenarios. LE completed 81/96 and 76/96. The angular LE control also had failures. These are excluded as candidates under the stability requirement.
- Across all controls, 15/798 short runs and 2/147 long runs cannot advance safely. The one-second and ten-second stress scans record 24 and 33 such failures respectively. These failures do not commit unstable states; their accepted prefixes are not treated as successful full-run accuracy.
- All work variants selected for the stress scan completed all 96 scenarios at both durations. That does not by itself validate their trajectories: the no-angular orbit counterexample remains poor.

## Accuracy, cost and limits

The nominated setting improves short-run RMS over XPBD at all seven cases and all three frequencies. It has remaining phase error and is generally first-order in these endpoint tests. It is not as accurate as the previously investigated midpoint/discrete-gradient method at nominal parameters; that is not the candidate being pursued here.

At 60 seconds and 480 Hz, nominated RMS is 1.872 cm for spring release, 23.152 cm for the rest pendulum, 0.740 cm for the small double spring and 2.919 cm for the orbit. The mixed double spring reaches 32.58 cm, slightly worse than XPBD's 30.94 cm. Thus the all-case ten-second improvement does not extend to every long trajectory. The large-double sixty-second reference fails pointwise refinement certification; its RMS is diagnostic only. Energy and accepted-state bounds remain independently measurable.

The implementation explicitly computes a baseline plus the coupled correction. The nominated method uses 5 scalar visits per step for the spring, 6 for the rigid pendulum and 11 for the two-spring system, versus XPBD's 1/1/2. It also evaluates work derivatives and assembles the small Schur matrix. These counts are not a wall-time performance benchmark. No production-code optimization is claimed.

The most useful next question is how to preserve more physical energy while keeping the angular/work structure and bounded solve cost. The rigid-pendulum rank limitation should be addressed explicitly rather than increasing an energy multiplier and expecting another degree of freedom.

## Validation and reproduction

Unit tests check a hand-derived compliant scalar spring, angular-row finite differences, energy and work dual GS against independent positive primal Gauss–Newton systems, hard/soft angular compliance, accepted energy/cap bounds, work reservoir bookkeeping, endpoint velocity reconstruction, rod feasibility and no invalid state commit on failure. Optimized CTest and ASan/UBSan pass.

The independent plotting audit reconstructs analytic physical energy, endpoint velocities, position RMS and reservoir bookkeeping from saved states. Main traces contain every 480 Hz step; long traces contain every eighth step plus the final accepted state. All 21 ordinary-XPBD summaries reproduce the earlier baseline exactly. All short references pass 16-to-32-subdivision RK4 refinement; only the large-swing sixty-second reference is uncertified.

```bash
cmake -S experiments/implicit_candidates -B /tmp/xpbd-implicit-build -DCMAKE_BUILD_TYPE=Release
cmake --build /tmp/xpbd-implicit-build -j4
ctest --test-dir /tmp/xpbd-implicit-build --output-on-failure
/tmp/xpbd-implicit-build/xpbd-implicit-study experiments/implicit_candidates/results 10
/tmp/xpbd-implicit-build/xpbd-implicit-study experiments/implicit_candidates/results 1 --stress
/tmp/xpbd-implicit-build/xpbd-implicit-study experiments/implicit_candidates/stress_long_results 10 --stress
/tmp/xpbd-implicit-build/xpbd-implicit-study experiments/implicit_candidates/long_results 60 --long
/tmp/xpbd-implicit-build/xpbd-implicit-reservoir-probe experiments/implicit_candidates/results
/home/blacknoir/Projects/mandelbrot-diffusion/.venv/bin/python experiments/implicit_candidates/plot_results.py experiments/implicit_candidates/results --long-results experiments/implicit_candidates/long_results --stress-long experiments/implicit_candidates/stress_long_results --output experiments/implicit_candidates/plots
```

The numerical solvers are C++17 with double precision. Python is used for scientific figures and independent auditing. Production `xpbd.hpp` and all previous experiments are unchanged. The full CSVs and failure statuses are preserved; `plots/recommended_candidate.png` is the compact energy/position comparison, and `plots/all_methods.png` includes all controls.
