# Independent CPD validation

The independent reference assembles planar spring-network PD curvature, the true
objective gradient, incoming momentum/power targets, and conservation rows. It
compares paper equation 21 with a dense KKT solve, rather than duplicating the
production Schur increment formula.

## Coverage

- Forty-eight arbitrary two/three-particle and anchored network states, nonzero
  velocities, unequal masses, gravity and nonuniform external forces, off-solution
  positions and relaxation values, and explicit angular-row ablations.
- Finite differences of the objective gradient and full P/L/energy Jacobian;
  frozen projection curvature equals the full PD matrix rather than GN normal
  outer products. Translations use a consistently translated incoming origin.
- Kinetic reference checked against an independent minimum-kinetic-energy solve
  with the enabled incoming momentum rows. The production angular-row ablation
  deliberately keeps the original kinetic reference unchanged.
- Production target values, gradients, internal versus mechanical energy, PD
  metric, global position/relaxation increment, and newly solved conservation
  multipliers compared with independent assembly and dense KKT solutions.
- Fifty-two refreshed global iterations across all thirteen scenes, including
  analytic hard-geometry gradients for two explicitly separate rod extensions.
  The trace checks newly solved multipliers, not their cumulative sum.
- Eight free-pair versus anchored comparisons: reduced mass one, equal relative
  spring state and force, and equal angular target. Relative positions, velocities,
  and energy relaxation agree to roundoff.
- One axial PD global solve agrees with analytic backward Euler. PD convergence
  uses objective stationarity, independently of inactive CPD constraints.
- Compatible dependent conservation rows and an incompatible zero energy row are
  exposed analytically. No reference step silently substitutes a pseudoinverse.
- A full-free-body release counterexample satisfies P/L/relaxed-energy feasibility
  and KKT stationarity with a nonzero relaxation value despite a feasible strict
  target. An initial conservation-only precheck can return zero iterations while
  objective stationarity remains unsatisfied.

## Numerical limitations

Random production increments and multipliers match dense KKT to approximately
5.1e-14. Finite-spring trace agreement is below 1e-10. The hard rod extensions
are substantially less well conditioned: the kick example reaches scaled Schur
condition 1.28e10 and frozen linear residual 1.67e-7. Its four-step dense replay
differs by 8.57e-7 in state and 1.14e-6 in multiplier, within the explicitly
declared 2e-6 trace comparison tolerance. The receipt exposes these values rather
than presenting that extension as a well-conditioned original-CPD test.

These are equation/implementation checks, not proof of acceptable motion. The
paper conservation stop, full KKT stationarity, physical energy under external
loads, and visual motion remain separate diagnostics. The source cross-check in
CPP_PARITY.md uses the prior unchanged local C++ baseline, not authors' code.

```bash
OPENBLAS_NUM_THREADS=1 python3 validate.py --write-report
```

## Numerical receipt

```json
{
  "checks": 1521,
  "details": {
    "explicit_rank_cases": 2,
    "free_pair_anchored_comparisons": 8,
    "full_free_release_counterexample": true,
    "global_trace_visits": 52,
    "hard_rod_extension_scenes": 2,
    "production_random_states": 48,
    "random_network_states": 48,
    "trace_error_by_scene": {
      "free_opposite_load": {
        "max_linear_residual": 3.6538132011857363e-16,
        "max_scaled_schur_condition": 5.859078875626062,
        "multiplier": 7.330361167402996e-16,
        "state": 1.8973538018496328e-19
      },
      "free_orbit": {
        "max_linear_residual": 1.4014746226210998e-15,
        "max_scaled_schur_condition": 100.45243982381561,
        "multiplier": 9.127984018123981e-15,
        "state": 2.4937018239641678e-18
      },
      "free_release": {
        "max_linear_residual": 1.429301121902406e-16,
        "max_scaled_schur_condition": 1.0000000000000009,
        "multiplier": 1.1203894404305942e-16,
        "state": 0.0
      },
      "free_stiff": {
        "max_linear_residual": 1.4269963280694033e-16,
        "max_scaled_schur_condition": 2.5894897342394785,
        "multiplier": 5.781797042956611e-16,
        "state": 1.734723475976807e-18
      },
      "free_triangle_spin": {
        "max_linear_residual": 2.753793558982719e-16,
        "max_scaled_schur_condition": 5.31155137123325,
        "multiplier": 6.120318902011822e-16,
        "state": 1.5114320952053415e-18
      },
      "pendulum_kick_compliant": {
        "max_linear_residual": 4.655955122547183e-14,
        "max_scaled_schur_condition": 56453.97245890508,
        "multiplier": 5.049843608964059e-12,
        "state": 1.0603665587198487e-14
      },
      "pendulum_kick_hard": {
        "max_linear_residual": 1.6705684346973726e-07,
        "max_scaled_schur_condition": 12796711227.655949,
        "multiplier": 1.1397763842537704e-06,
        "state": 8.561130541675175e-07
      },
      "pendulum_release_compliant": {
        "max_linear_residual": 5.187398382176023e-15,
        "max_scaled_schur_condition": 166.36590867041144,
        "multiplier": 8.680034533342176e-13,
        "state": 2.3762586247130332e-15
      },
      "pendulum_release_hard": {
        "max_linear_residual": 8.147722496687493e-10,
        "max_scaled_schur_condition": 44619858.687493876,
        "multiplier": 1.6441924084926927e-09,
        "state": 8.382607509227402e-10
      },
      "spring_gravity": {
        "max_linear_residual": 2.9230505609485893e-16,
        "max_scaled_schur_condition": 5.8590788756260626,
        "multiplier": 3.433175098891678e-16,
        "state": 6.61259742849214e-20
      },
      "spring_orbit": {
        "max_linear_residual": 6.814111089093229e-16,
        "max_scaled_schur_condition": 100.45243982381533,
        "multiplier": 7.862853733620395e-14,
        "state": 2.0186583106402866e-16
      },
      "spring_release": {
        "max_linear_residual": 1.4210854715202004e-16,
        "max_scaled_schur_condition": 1.0000000000000002,
        "multiplier": 1.1102230246251565e-16,
        "state": 4.3369444779693586e-20
      },
      "spring_stiff": {
        "max_linear_residual": 1.426996298903188e-16,
        "max_scaled_schur_condition": 2.589489734239479,
        "multiplier": 3.2373430828766563e-16,
        "state": 1.58933032037739e-18
      }
    }
  },
  "largest_relative_errors": {
    "all scene analytic gradients including exact rod extension": 3.049020276088844e-15,
    "all scene conservation finite differences including hard geometry": 1.5713131564267035e-10,
    "all scene constraints including exact rod extension": 0.0,
    "consistent origin translation preserves conservation Jacobian": 1.4540455867844627e-13,
    "consistent origin translation preserves conservation rows": 1.7837462355801318e-13,
    "consistent origin translation preserves objective gradient": 6.037503829909497e-14,
    "free pair energy relaxation versus anchored relaxation": 5.89534931288993e-19,
    "free pair relative position versus anchored position": 2.0186009369967686e-16,
    "free pair relative velocity versus anchored velocity": 6.662263269364698e-15,
    "frozen local projection gives full PD curvature": 4.533742319595099e-11,
    "full P L E relaxed release branch feasible": 1.1457501614131615e-13,
    "full P L E relaxed release branch stationary": 0.0,
    "full conservation Jacobian finite differences": 1.2813232874553246e-10,
    "global step trace NEW multiplier versus dense replay": 1.1397763842537704e-06,
    "global step trace position and relaxation versus dense replay": 8.561130541675175e-07,
    "implicit plus relaxation objective gradient finite differences": 7.801592126532932e-11,
    "inertia kinetic reference versus constrained least kinetic energy": 4.436296992744553e-16,
    "linearized conservation satisfied jointly": 4.804532694219296e-13,
    "one axial PD global solve equals analytic backward Euler": 0.0,
    "paper equation21 fresh multiplier versus dense KKT": 1.086894337244368e-14,
    "paper equation21 versus full dense KKT position and relaxation": 1.241124109476153e-14,
    "production PD metric not rank-one GN metric": 0.0,
    "production analytic P L energy gradients": 1.308461472639314e-14,
    "production equation21 fresh multipliers versus independent dense KKT": 5.02670798023094e-14,
    "production equation21 increment versus independent dense KKT": 4.8440938634002973e-14,
    "production incoming COM or anchor origin": 0.0,
    "production incoming angular and torque target": 0.0,
    "production incoming internal energy plus power target": 0.0,
    "production incoming momentum target": 0.0,
    "production inertia kinetic reference": 0.0,
    "production kinetic plus elastic constrained energy": 0.0,
    "production mechanical energy diagnostic includes force potential": 0.0,
    "production true P L energy constraint values": 2.220446049250313e-16,
    "production true elastic and inertial objective gradient": 3.1792695904862887e-15,
    "quadratic subproblem stationarity satisfied": 5.237782530015209e-12,
    "redundant momentum and energy rows can be exactly compatible": 7.850462293418876e-16,
    "refreshed global trace true constraints": 4.848042832091071e-09,
    "returned velocity reconstructed from physical displacement": 0.0,
    "zero energy gradient at unstretched incoming position": 0.0
  },
  "status": "passed"
}
```
