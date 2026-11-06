# Corrected objective crossed with start-only work

The quarter predictor / quarter potential position objective and the stiffness-aware midpoint-predicted DG work target compose usefully. In fully converged joint trajectories, the radial RMS position error falls from 0.08586 (original + frozen power), to 0.003991 (quarter + frozen power), to 0.0008308 (quarter + predicted DG). The two-mass chain changes from 0.2163 to 0.008144 to 0.001568. Each requests one second at h=1/60 with beta=gamma=1 and the same soft midpoint velocity penalty, explicit work variables, and local energy constraints.

The 24 trajectories remain in `results/objective_work_cross.csv` and `results/cross_traces/`. Six gradient/Hessian and unchanged-constraint checks pass. Local AL at rho10 and512 evaluations matches the radial/chain motion closely. The corrected chain contains4/1 mildly unconverged AL steps for frozen power/predicted DG; maximum stationarity is1.51e-7/1.26e-7 against a1e-7 threshold.

## Scalar turning-point stalls

The original joint solver with its backward-Euler initialization stalls twice on each corrected scalar trajectory, at steps18 and54. Its maximum balance defect is about0.002 and stationarity about0.01, so its scalar RMS errors near0.03 are contaminated by solver failure. The corresponding AL trajectories have RMS errors0.00337/0.00220 with only small finite residual failures, and maximum balance4.97e-7/2.19e-7.

`scalar_stalls.cpp` replays the exact four old states where the joint solves fail. It compares initializations and AL on the identical objective. Results are saved separately in `results/scalar_stalls.csv`:

| Numerical choice at the four stalled states | Result |
|---|---|
| BE initialization, budget128 | All four stall after69–112 evaluations |
| BE initialization, budget4096 | The identical stalls and iteration counts |
| Start-only linear midpoint initialization, budget128 | All four converge in1–3 evaluations |
| Local AL rho10, budget2048 | All four converge in395–453 evaluations |
| Dense feasible-circle scan followed by KKT | Agrees with the successful solution, converges in1–2 evaluations |

This scalar has no gravity, so its energy manifold is the ellipse `x=l+A cos(theta), v=-omega A sin(theta)`; actual work then follows from either balance. The independent scan checks8192 angles and refines the best interval. It supports the observed minimum but is a numerical global scan, not a certified global-optimization proof. At successful solutions the reduced Hessian is positive (0.607–0.715); at the stalled infeasible states it is near zero. Constraint Jacobians remain full row rank, with minimum singular value0.189 at the stalls and0.166 at the solutions. Thus the observed failure is a poor-initialization/KKT globalization problem near a turn, not an unavoidable singular balance constraint.

For the harmonic scalar with the midpoint-predicted DG target, the linear midpoint state satisfies both energy balances, the soft midpoint relation, and zero work-target correction. The quarter objective also has zero unconstrained gradient there. It is already a KKT solution, explaining convergence in one evaluation. The frozen-power target leaves a small work-penalty force, so the same initialization takes three evaluations.

The one-second original CSVs are intentionally preserved. A longer follow-up studies initialization and globalization rather than replacing these failures after the fact.

## Reproduce

Use the include paths from `README.md`:

```bash
g++ -std=c++17 -O2 -pthread -I"$prior_composition_dir" -I"$eigen_include_dir" cross_objective.cpp -o cross_objective
g++ -std=c++17 -O2 -pthread -I"$prior_composition_dir" -I"$eigen_include_dir" scalar_stalls.cpp -o scalar_stalls
./cross_objective results
./scalar_stalls results
python3 analyze.py
```

The exact tested absolute include paths were `/home/blacknoir/Projects/constrained-xpbd/experiments/composition_2026_10_31` and `/home/blacknoir/Projects/constrained-xpbd/external/eigen`. The source itself includes headers by name.
