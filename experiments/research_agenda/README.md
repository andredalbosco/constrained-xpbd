# Small-step research agenda: spring, pendulum and double spring pendulum

This is a controlled exploration of the ideas listed on pages 77–85 of the supplied notes: conditioning, inverse-compliance dualization, local work balance, work normalization, endpoint velocity reconstruction, reservoirs and numerical drift. It extends the existing studies without modifying the production solver or earlier experiments. The midpoint discrete-gradient direction already existed in `experiments/work_balance/midpoint_gs.hpp`; here it is compared in the same 2D fixtures and reference framework as the endpoint-energy methods.

There are 630 ten-second runs: 30 configurations × 7 initial conditions × 240/480/960 Hz. A further 112 runs test 16 selected configurations for 60 seconds at 480 Hz. Only the large-swing double spring loses long-horizon reference certification. All seven short-run references pass refinement. No scene has more than two moving masses.

## What was varied

| Direction | Concrete experiment |
|---|---|
| Baselines | Ordinary XPBD, the existing bounded frozen EL and LE schedules |
| Reconstructed positions | Rebuild the row after each E/L correction, retain or drop `r_x` |
| Dualization | Residual-complete dual expression for the energy gradient; an approximation omitting the constitutive residual as a negative control |
| Work estimates | Lagged force, row matching, column matching, and an exact discrete gradient of spring potential |
| Work strength | Dimensionless penalty strengths 0.1, 1 and 10 |
| GS ordering and budget | Material → potential → kinetic, material → kinetic → potential; 1, 4 or 16 frozen GS sweeps |
| Reservoirs | Ledger only; local work-target repayment with time constants 0.25 or 1 second |
| Kinematics | Endpoint velocity reconstruction versus midpoint kinematics with discrete-gradient force |
| Midpoint solver budget | 1, 4 or 8 relinearizations, one material GS sweep per relinearization; an `r_x`-dropped control |
| Conditioning | 168 algebraic fixtures spanning mass ratios 1–10,000, stiffness factors 1–10,000, two timesteps, all seven cases |

Work and endpoint variants use the ordinary predictor, cold multipliers at every timestep, and analytic physical spring/gravity energies. Work transfer estimates use the preceding physical state, not a future converged force. Reservoirs alone persist between timesteps.

The primary cases are a radial spring release, a rigid pendulum released at 0.65 radians, and small/large swings of the two-spring pendulum. Moving holdouts add a spring orbit, a kicked pendulum and mixed double-spring velocities. Fixtures and RK4 references are copied from the traditional study in `cases.hpp`; all masses, compliances, gravity and initial states are unchanged. Spring stiffness is 100 N/m, masses are 1 kg and rest lengths are 1 m.

## Findings

At the nominal parameters, the strongest direction in this comparison is **midpoint kinematics with a discrete gradient of potential energy**. It changes the discrete equations, instead of correcting the endpoint-energy row of the original implicit objective. Every numerical solve still uses linearization, dropped force curvature and scalar Gauss–Seidel. No nonlinear energy projection is performed.

**The fixed-four-pass implementation is not recommended as a stable method.** Subsequent stiffness/perturbation stress tests found failures and energy growth. The stability-constrained candidate is the version in `stability.hpp`, which accepts only physically bounded states and retries with smaller substeps. The next section records the limits and additional cost. At the nominal short/long study parameters its checks are inactive, so its trajectories equal the four-pass results tabulated below.

## Stability requirement: fixed solve rejected, guarded variant tested

Exact implicit midpoint for a linear oscillator has amplification `(1+i*omega*h/2)/(1-i*omega*h/2)`, whose magnitude is one. Explicit midpoint instead has squared magnitude `1+(omega*h)^4/4 > 1`. Neither statement proves stability of a finite-GS solve of our nonlinear discrete-gradient equations.

The dedicated stress study has 96 scenarios for each of three methods: spring release/orbit, rest pendulum and large double spring; 240/960 Hz; stiffness multipliers 1/100/10,000; light-mass ratios 1/100; unperturbed or rotated by 1e-8 radians. Each scenario runs one simulated second. The test is an energy-boundedness screen, not a high-stiffness trajectory-accuracy benchmark.

| Method | Scenarios | Completed | >1% energy growth or failure |
|---|---:|---:|---:|
| XPBD | 96 | 96 | 0 |
| Fixed four-pass midpoint | 96 | 90 | 16 |
| Guarded four-pass midpoint | 96 | 96 | 0 |

Six fixed-solve runs hit the 100-times-excitation-energy abort threshold; the worst attempted state reached about 470 times initial excitation energy. Fixed-cost midpoint is excluded from the stable candidates, despite its excellent nominal results. Existing unguarded endpoint/work variants are also not promoted on the basis of completion alone.

The guarded variant checks, before committing any state:

```
abs(E_candidate - E_initial) / E_scale <= 1e-6
maximum local kinetic/potential work defect / E_scale <= 1e-9
primal and constitutive residuals <= 1e-9 * spatial_scale
rigid length error <= 1e-9 * rest_length
all positions/velocities finite.
```

Rejected intervals are retried as two half steps, recursively to depth 12. Accepted substeps obey the same **global** energy band around initial energy, so the bound does not accumulate once per accepted step. The primal/constitutive checks also reject an energetically admissible frozen state that does not solve the mechanical equations. A failed retry returns failure without committing an invalid state; this is an accepted-state energy bound, not a guarantee that every nonlinear problem can be advanced within the retry budget. The supported anchored spring potentials are coercive, and the rod length is bounded explicitly, so bounded physical energy controls these states. This is not a general stability theorem for arbitrary materials/contact.

All 96 guarded scenarios completed. Maximum absolute energy drift was 3.52e-8 of the initial excitation scale. Worst recursion depth was 7 and the largest average accepted-substep count was 55.57 per requested step. Thus stability can require a substantial extra solve budget. No position/velocity rescaling or nonlinear energy projection is used.

The checks accept every nominal four-pass short/long result: physical energy, work and rod defects are below their thresholds throughout. Unit tests verify unchanged states when accepted, exercise a genuinely rejected stiff step and verify no invalid state commit when the retry budget is exhausted. Stability does not imply phase accuracy; the long-horizon limitations below still apply.

Representative ten-second runs at 480 Hz:

| Case | Method | Maximum physical energy error | Cartesian position RMS |
|---|---|---:|---:|
| Spring release | XPBD | 87.54% | 7.416 cm |
| | Existing bounded EL | 2.009% | 11.068 cm |
| | Rebuilt LE, retain `r_x`, bounded | 0.1462% | 0.500 cm |
| | Work + fast reservoir repayment | 6.176% | 0.487 cm |
| | Midpoint work, 1 or 4 passes | 4.95e-10% | 0.0372 cm |
| Rigid pendulum | XPBD | 21.99% | 3.765 cm |
| | Existing bounded EL | 0.006179% | 0.466 cm |
| | Work + fast reservoir repayment | 0.7607% | 0.380 cm |
| | Midpoint work, 4 passes | 1.19e-10% | 0.00324 cm |
| Double spring, small swing | XPBD | 41.58% | 1.176 cm |
| | Existing bounded EL | 1.481% | 2.629 cm |
| | Work + fast reservoir repayment | 9.159% | 1.346 cm |
| | Midpoint work, 1 pass | 0.006080% | 0.00165 cm |
| | Midpoint work, 4 passes | 2.05e-9% | 0.00260 cm |
| Double spring, large swing | XPBD | 63.43% | 14.498 cm |
| | Existing bounded EL | 2.132% | 33.012 cm |
| | Work + fast reservoir repayment | 14.73% | 19.935 cm |
| | Midpoint work, 1 pass | 0.01193% | 0.2888 cm |
| | Midpoint work, 4 passes | 1.16e-9% | 0.2263 cm |

Energy error is divided by initial energy above the analytic minimum potential. All energy columns use physical kinetic + spring + gravitational energy, excluding reservoirs. Position RMS includes both coordinates and both masses when present. Values near roundoff are observed errors, not exact arithmetic guarantees.

### Important qualifications and negative results

- **A fixed endpoint method was not rescued by simply retaining `r_x`.** Bounded rebuilt EL remains close to the previous EL. Raw retained-residual EL aborts on spring release at all three rates, and raw EL/LE abort on some large double-spring runs. Eight of 630 ten-second configurations hit the safety bound. Accepted prefixes are not reported as successful complete trajectories in `plots/all_methods.png`.
- **Normalization alone was insufficient.** At 480 Hz on radial release, lagged work gives 70.36% maximum physical energy error and 4.697 cm RMS; column/row/secant versions are close to 70.35% and 4.696 cm. On the double spring, the choice changes behavior but gives no common winner. These are finite-sweep, capped solves, not exact solutions of every work-balance equation.
- **Strength and ordering matter but do not generalize.** Increasing work strength from 0.1 to 10 improves radial-release RMS from 7.081 to 1.605 cm, but worsens rest-pendulum RMS from 3.097 to 7.193 cm. Sixteen frozen GS sweeps can produce worse trajectories than one. Reducing a frozen linear residual does not repair the nonlinear model or select a physical branch.
- **The reservoir ledger does exactly what it should.** XPBD + ledger has bitwise-identical saved states to ordinary XPBD. Recording defects adds no physical correction. Repayment helps release/rest tests, but not the spring orbit: fast repayment yields 0.6525% energy error and 92.01 cm RMS, versus XPBD's 47.46% and 5.790 cm. Local reservoirs can have substantial opposite signs, so their absolute sum is also plotted.
- **The orbit is a useful discriminator.** Midpoint work with four passes gives 0.0139 cm RMS at 480 Hz, with angular-momentum defect about 3.3e-12 kg m²/s. Endpoint energy or local reservoir repair can maintain small energy error while seriously changing the orbit.
- **One pass need not minimize trajectory error.** On the small double spring, the one-pass and dropped-residual controls have smaller ten-second RMS than converged midpoint, despite larger work/energy defects. This is error cancellation at these resolutions, not a reason to identify lower energy error with lower phase error.
- **Midpoint retains phase error.** At 480 Hz over 60 seconds, midpoint RMS is 0.2215 cm for radial release, 0.0189 cm for the rest pendulum and 0.0157 cm for the small double spring. Large/mixed double motion can separate strongly at long times even with excellent conserved quantities. No trajectory-exact claim is made.
- **This is an anchored-system study.** Angular balance accounts for gravity torque. Global linear momentum is not constant because of anchor reactions and gravity. Internal spring impulses are equal and opposite; the method is not yet validated on a general free multi-body system, cloth, collision or friction.

## Endpoint residual and dualization experiments

With `y=x^-+h v^-+h² a`, retain the adopted scaled residuals

```
r_x = M(x-y) - J^T lambda - g nu
r_lambda = C + alpha_tilde lambda
r_nu = E - E_initial
lambda^+ = lambda^- + delta lambda
nu^+ = nu^- + delta nu.
```

At every visit, the rows are rebuilt. With multiplier-weighted Hessians dropped, a material visit holding nu fixed solves

```
(J M^-1 J^T + alpha_tilde) delta lambda
    = -r_lambda + J M^-1 r_x
delta x = M^-1 (J^T delta lambda - r_x).
```

An energy visit holding material multipliers fixed uses the analogous scalar equation with J replaced by `g^T`. The energy gradient includes analytic spring and gravitational terms. Singular energy gradients are skipped using a tolerance that includes floating-point position uncertainty. The bounded proposal uses 50% predicted deficit repair and a 10% extra-motion cap, with actual positions and multipliers blended consistently against a material baseline.

For soft constraints, the exact dual-gradient identity is

```
J^T alpha^-1 C
  = -J^T lambda/h² + J^T alpha^-1 r_lambda.
```

The residual-complete expression reproduces the analytic gradient. Conditioning fixtures test it with nonzero lambda and material residuals; maximum scaled error is 3.62e-14. Trajectory metric differences are at roundoff/accumulated rounding scale. Omitting the residual term is an explicit energy-direction approximation, not another exact dualization. That control is preserved in the all-method CSV/scorecard.

## Local work formulation

Write independent work as `W_hat(x)+t`, with one t for every supported particle/potential incidence. At each timestep solve the primal problem

```
min_(x,t) h² F_XPBD(x) + h² gamma/(2 E_scale) ||t||²
subject to
    delta K_i(x) - sum_j [W_hat_ij(x)+t_ij] - p RK_i^- = 0
    delta U_j(x) + sum_i [W_hat_ij(x)+t_ij] - p RU_j^- = 0,
```

including hard rod constraints where appropriate. Here `gamma` is the work-estimate strength and `p=1-exp(-h/tau)` for repayment, or zero without repayment. Gravity has a separate conservative potential node per particle. A rigid rod has zero potential energy and a guessed lagged reaction force in its work estimate; that estimate is not an exact discrete reaction law.

Soft spring energy is dualized with `alpha_tilde=alpha/h²`. After stationarity and linearization, drop multiplier-weighted curvature. The primal metric is diagonal: M for positions and `h² gamma/E_scale` for work deviations. Eliminate the primal increments, then solve the material/potential/kinetic multiplier rows by GS. One, four and sixteen GS sweeps here refer to the **same frozen linear system**, not extra relinearizations.

The force estimate depends on current x through displacement, and row/column/secant variants have additional x dependence. All those derivatives are included using a small forward-mode derivative type and checked independently by finite differences. This implements the implicit-work coupling explicitly rather than treating work as constant after position corrections.

Work estimates:

- Lagged: `W_hat_ij=f_ij^- dot (x_i-x_i^-)`.
- Column matched: add equal corrections to incident entries so `sum_i W_hat_ij=-delta U_j(x)`.
- Row matched: add equal corrections so `sum_j W_hat_ij=delta K_i(x)`.
- Potential exact: use the discrete spring gradient defined below, with constant gravity work. For a rod this gives zero work when both endpoint lengths are equal.

This row/column matching is a proposed interpretation of the source's “work normalization” note. It is signed additive balancing, not division by a possibly zero work sum. It differs from solver-matrix equilibration. A single spring and particle without gravity has a one-entry work matrix, so the allocation distinction disappears there.

The work proposal is capped relative to ordinary XPBD: its extra mass-weighted displacement is at most 10% of `h||v^-||_M+h²||a^-||_M`. The cap is applied to the actual proposal, with all multiplier and t changes blended. This safeguard does not impose nonlinear work balance. Actual accepted-step defects are measured afterward.

## Reservoir definition

Reservoirs store signed local work residuals, not a global endpoint-energy target:

```
bK_i = delta K_i - sum_j W_ij
bU_j = delta U_j + sum_i W_ij
RK_i^+ = RK_i^- - bK_i
RU_j^+ = RU_j^- - bU_j.
```

Consequently physical energy plus the sum of reservoirs is a bookkeeping invariant. The same definitions are used for ledger-only XPBD and feedback methods. Repayment targets `bK=p RK^-`, `bU=p RU^-`; if those equations were exactly satisfied, each reservoir would decay by `exp(-h/tau)`. Finite GS solves, nonlinear residuals and caps prevent exact repayment. The independent trace audit checks the actual augmented-energy identity; it is not used as evidence of physical energy conservation.

## Midpoint kinematics and discrete work

The integrator changes the hard kinematic relation to

```
x^+ - x^- = h (v^+ + v^-)/2
M(v^+ - v^-) = h (M a + sum_j f_bar_j).
```

For an edge `e=x_b-x_a`, radius r and spring `C=r-L`, choose

```
D_j = (e_j^+ + e_j^-)/(r_j^+ + r_j^-)
C_bar_j = (C_j^+ + C_j^-)/2
f_bar_j = -D_j C_bar_j/alpha_j.
```

These choices give `D dot delta e = delta r` and therefore

```
delta U_j = C_bar_j/alpha_j * D_j dot delta e_j
delta K = sum_i (v_i^+ + v_i^-)/2 dot M_i delta v_i.
```

Using midpoint kinematics, force work cancels the exact spring potential change; constant gravity work cancels gravitational potential change. Equal/opposite internal forces are along average edge directions, giving the matching discrete angular-momentum/torque identity. A rod instead uses the quadratic distance constraint and its exact discrete gradient, so its work vanishes when its endpoint constraint is satisfied.

The implementation solves in midpoint positions `x_m=(x^++x^-)/2`, with predictor `x_m_pred=x^-+h v^-/2+h² a/4`. Its internal material multiplier `ell` has physical-force scale `4 ell/h²`; it is one quarter of the usual h²-scaled force multiplier. Soft compliance in this internal system is `4 alpha/h²`. Every update is `ell^+=ell^-+delta ell`.

Each outer pass rebuilds material rows, retains the primal stationarity residual, drops derivatives of the discrete force directions, and performs one scalar GS sweep. The force direction is D; the left material derivative is the endpoint derivative. Their Schur product is generally nonsymmetric. No dense solve or energy root finder is used. Exact work identities explain conservation only when the relevant discrete residuals vanish; one-pass errors are measured explicitly.

Across the 21 four-pass short runs, the largest normalized physical energy error is 4.88e-11. Position errors decrease by approximately four when doubling frequency on release, rest-pendulum and orbit cases, and nearly four on the large double spring over ten seconds. This supports second-order behavior on these tests, not a general theorem for arbitrary constraints or fixed finite iterations.

## Conditioning and numerical drift

Diagonal equilibration uses `D=diag(S)` and the equivalent system `D^-1/2 S D^-1/2 u=D^-1/2 b`. Scaled/unscaled GS visits agree in reconstructed positions to 8.91e-16 relative across the nonzero-gradient fixtures. Large changes in condition number therefore do not imply different exact-arithmetic corrections. The normalized material/energy correlations and dependencies remain. Twenty-four from-rest rigid-pendulum fixtures have zero energy gradient and are flagged rather than divided by zero.

There are 456 conditioning output rows: 144 nonsingular fixtures at 1/4/16 GS sweeps plus 24 zero-gradient fixtures. These are algebraic probes, not high-stiffness trajectory benchmarks.

Numerical drift is diagnosed through physical energy, actual work defects, constitutive/primal residuals, angular balance, both particle positions, step refinement and independent references. An old-force spring work estimate has the exact defect `delta U+W_old=k(delta s)²/2` even in exact arithmetic. Conversely, conserved energy does not eliminate numerical phase error or sensitive long-time trajectory separation.

## Validation and limits

- The 63 overlapping XPBD/EL/LE short runs reproduce all four checked earlier summary metrics exactly.
- Unit tests check finite-difference work derivatives, exact row/column/discrete-potential identities, Schur reconstruction against an independent pivoted long-double oracle, diagonal-scaling invariance, nonzero-residual dualization, analytic backward-Euler and midpoint spring maps, reservoir bookkeeping/cap, and multibody energy/angular identities.
- Optimized CTest and AddressSanitizer/UndefinedBehaviorSanitizer tests pass.
- Python independently reconstructs physical energy, endpoint/midpoint kinematics, reservoir sums, angular balance and position RMS from saved states. The 112 main trace groups contain every 480 Hz step. RMS reconstruction differs by at most 8.89e-16 m; ledger-only states equal baseline states exactly. Long traces save every eighth step, while summaries use every step.
- All short references are refined from 16 to 32 RK4 subdivisions. At 60 seconds, the large double-spring references differ by about 5.2 mm and 0.032 m/s while both keep energy drift near 3e-13. Its long-run RMS is explicitly uncertified in the CSV. Physical energy/constraint diagnostics remain valid independently of that reference.
- Cost is reported as solver structure, not a wall-time benchmark: ordinary XPBD uses n scalar visits per step, midpoint with p passes uses p*n visits, and local-work schedules use `n + potential_count + particle_count` visits per frozen GS sweep plus work variables/derivatives. CSV wall times include diagnostics and trace I/O and should not be used as solver benchmarks.
- This does not test a freely optimized velocity with a soft consistency penalty, general non-diagonal preconditioners, long-horizon high-stiffness trajectories, or a general conservation method for cloth/contact. The actual kinematic alternative tested is midpoint work.

## Reproduce and inspect

```bash
cmake -S experiments/research_agenda -B /tmp/xpbd-agenda-build -DCMAKE_BUILD_TYPE=Release
cmake --build /tmp/xpbd-agenda-build -j4
ctest --test-dir /tmp/xpbd-agenda-build --output-on-failure
/tmp/xpbd-agenda-build/xpbd-agenda-study experiments/research_agenda/results 10
/tmp/xpbd-agenda-build/xpbd-agenda-conditioning experiments/research_agenda/results
/tmp/xpbd-agenda-build/xpbd-agenda-stability experiments/research_agenda/results
/tmp/xpbd-agenda-build/xpbd-agenda-study experiments/research_agenda/long_results 60 --long
/home/blacknoir/Projects/mandelbrot-diffusion/.venv/bin/python experiments/research_agenda/plot_results.py experiments/research_agenda/results --long-results experiments/research_agenda/long_results --output experiments/research_agenda/plots
/home/blacknoir/Projects/mandelbrot-diffusion/.venv/bin/python experiments/research_agenda/plot_stability.py
```

All numerical solvers are C++17 using double precision; Python with NumPy/Matplotlib is used only for independent analysis and scientific plots. `results/summary.csv` preserves all methods, failures and attempted failure energies. `plots/all_methods.png` shows every 480 Hz configuration; other figures focus on formulation, normalization, reservoirs, conditioning, convergence, moving holdouts, both masses of the double pendulum and long runs.
