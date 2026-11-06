# Solving the repaired local-work objective

This follow-up keeps the quarter predictor / quarter potential objective and the explicit local-work balances fixed while testing numerical choices. The work target is always the start-only `implicit_midpoint_dg` prediction. All work values and all vertex/storage multipliers are shared across visits.

The main result so far is that a mechanically consistent initialization removes the observed scalar turning stalls, and also makes the joint KKT solve much more reliable across the five environments. Local AL still requires substantially more work on longer chains. Data include finite but unconverged steps, rather than hiding them behind an energy plot.

## Exact formulation

With `a_n` the old physical acceleration, the dimensionless objective is

\[
F(x,v,W)=\frac{1}{2E}\|(x-x_n)/h-v_n-\tfrac14h a_n\|_M^2
+\frac{U(x)}{4E}
+\frac{1}{2\beta E}\|v+v_n-2(x-x_n)/h\|_M^2
+\frac\gamma2\sum_{ij}\left(\frac{W_{ij}-\widehat W_{ij}}E\right)^2.
\]

Every run uses beta=gamma=1. Its constraints are each vertex's kinetic balance and each spring/gravity potential's storage balance. The rod additionally enforces its fixed length and tangency. Its old acceleration includes the physical reaction; the compliant cases use total spring and gravity forces.

The old-state matrix `B_n=sum k_j grad C_j grad C_j^T` gives the predictor

\[
(M+\tfrac14h^2B_n)\delta\widehat x=hMv_n+\tfrac12h^2f_n,
\qquad\widehat v=2\delta\widehat x/h-v_n.
\]

The predicted DG work is evaluated between old and predicted positions, then frozen. For the rod the predictor instead advances the angle and maps back to the circle with tangent velocity. This is a start-only estimate, not a work identity at arbitrary solved endpoints.

The physical energy scale E is fixed at the initial total excitation above the minimum. This normalization also sets the relative physical work-penalty weight; it is not merely numerical row scaling. Changing chain size changes both the physical model here and this total scale. Consequently the N-comparison tests solver behavior on its stated models, not replication invariance. The root experiment separately investigates that issue.

## Numerical comparisons

- `kkt_BE`: exact-curvature KKT Newton initialized from the earlier backward-Euler position guess and midpoint velocity reconstruction.
- `kkt_midpoint`: identical KKT equations initialized from the start-only midpoint predictor, with W initialized to its frozen prediction.
- `kkt_midpoint_warm`: the same primal predictor while retaining the previous step's multiplier vector. E and h are unchanged, so multiplier coordinates are unchanged.
- `continuation_BE`: position objectives at blend fractions0,0.5,1, with old-force shift0.25 times the fraction and potential coefficient1-0.75 times the fraction. All stages retain the exact same balance constraints and frozen target; the final objective is the requested quarter objective. The total evaluation budget includes the intermediate stages.
- `residual_LM_BE`: KKT Newton with a damped least-squares step on the KKT residual when its line search fails. A stationary residual with nonzero constraint/stationarity defects remains a failed solve.
- `merit_BE`: a positive-definite-shifted equality-SQP direction with an L1 objective/feasibility merit line search. The shift modifies the direction only.
- Local AL: overlapping primal blocks, three passes at fixed shared duals, then one shared dual update. Initializations, warm/cold duals, rho, budgets, and block partitions are recorded separately.

These are explicit prototype globalization rules; no universal convergence theorem is claimed. All methods use the same final convergence criteria: maximum balance/position/tangency residual below 1e-8 and full stationarity norm below 1e-7. The joint systems use full-pivot LU; the earlier suite's SVD version exhibits the same scalar stalls. The harmonic scalar midpoint predictor is already a feasible zero-gradient solution for this target and objective, so it legitimately stops after one evaluation.

The original AL partition has one spring-position/work block and one vertex-velocity/work block per mass. The new `vertex_full` partition jointly updates a vertex's position, velocity, and all incident work values. `edge_full` combines those variables for both endpoints, retaining every incident work value. Overlapping values and multipliers are stored once. A block update includes its effects on neighboring balances. Counts include block visits and dual updates; larger blocks and fewer visits per pass are not treated as equal FLOPs.

For one/two-mass examples, an edge-full block can cover the entire primal system. The N8/N32 experiments are necessary to assess locality. A separate `global` primal block is explicitly a dense control. `incremental.hpp` changes only the arithmetic of the Armijo value difference, using exact squared-norm/distance difference identities. Baseline and stable-difference runs use separate output files. [TURNING_AND_LOCALITY.md](TURNING_AND_LOCALITY.md) explains the scalar trigger, the local iteration spectrum, and this arithmetic plateau.

## Evaluator and validation

`problem.hpp` analytically assembles the same objective, constraints, and full first/second derivatives. It is checked against the prior independent AD evaluator across all five environments, including rod curvature and tangency. Only the diagnostic x/v/work fields omit derivative arrays; all objective/constraint derivatives used by solvers are complete. The analytic implementation avoids measuring the previous operator-overloaded AD cost on N32. It still uses dense derivative storage and dense joint factorizations, so this is not a scalable implementation claim.

The large-chain models have N=2,8,32 planar masses; rest length 1/N per spring; stiffness 100N times the selected factor; downward gravity 9.81; and a smooth angular perturbation around statically stretched lengths. Contrasts multiply all stiffnesses by10 or alternating masses by10. All references use RK4 with 64 subdivisions, with independent64-to-128 refinement checks. These are elastic chains, not multiple rigid rods.

The short large-chain run requests 8 steps at 60 Hz. The small-environment runs request60 steps. Every trace begins at time 0 and ends at the actual completed time. Failed and finite-inexact steps remain in the CSVs and traces. Curvature diagnostics are performed after the solve and do not count against its evaluation budget; positive reduced curvature establishes a checked local minimum, not global optimality.

## Reproduce

The source includes prior experiment headers by name. Set include directories appropriate to the saved layout:

The standalone CMake file also works through the hour suite's `add_subdirectory(local_solver)`. It accepts `COMPOSITION_ROOT`, `ROUND2_ROOT`, and `EIGEN_ROOT` cache paths. Targets are prefixed `hour_local_`; CTest runs the existing analytic/AD, incremental-difference, reference, branch, and iteration-matrix checks, writing validation output under the build directory.

```bash
composition_dir=../../composition_2026_10_31
round2_dir=../../round2_2026_10_31
eigen_dir=../../../external/eigen
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release \
  -DCOMPOSITION_ROOT="$composition_dir" -DROUND2_ROOT="$round2_dir" -DEIGEN_ROOT="$eigen_dir"
cmake --build build -j2
ctest --test-dir build --output-on-failure
```

```bash
round2_local_dir=/home/blacknoir/Global/constrained-xpbd-round2-2026-10-31/local_work
composition_dir=/home/blacknoir/Projects/constrained-xpbd/experiments/composition_2026_10_31
eigen_dir=/home/blacknoir/Projects/constrained-xpbd/external/eigen
g++ -std=c++17 -O2 -pthread -I"$round2_local_dir" -I"$composition_dir" -I"$eigen_dir" study.cpp -o study
g++ -std=c++17 -O2 -pthread -I"$round2_local_dir" -I"$composition_dir" -I"$eigen_dir" branches.cpp -o branches
g++ -std=c++17 -O2 -pthread -I"$round2_local_dir" -I"$composition_dir" -I"$eigen_dir" block_bench.cpp -o block_bench
g++ -std=c++17 -O2 -pthread -I"$round2_local_dir" -I"$composition_dir" -I"$eigen_dir" long_blocks.cpp -o long_blocks
g++ -std=c++17 -O2 -pthread -I"$round2_local_dir" -I"$composition_dir" -I"$eigen_dir" spectrum.cpp -o spectrum
g++ -std=c++17 -O2 -pthread -I"$round2_local_dir" -I"$composition_dir" -I"$eigen_dir" incremental_checks.cpp -o incremental_checks
g++ -std=c++17 -O2 -pthread -I"$round2_local_dir" -I"$composition_dir" -I"$eigen_dir" rod_control.cpp -o rod_control
g++ -std=c++17 -O2 -pthread -I"$round2_local_dir" -I"$composition_dir" -I"$eigen_dir" warm_blocks.cpp -o warm_blocks
./study validate results
./study small results
./study small_stable results
./study chain_smoke results
./study chains results
./branches results
./block_bench results
./incremental_checks results
./long_blocks results
./long_blocks results stable
./spectrum results
./rod_control results
./warm_blocks results
python3 analyze.py
```

Runs use up to four workers. Wall times are diagnostic, not isolated benchmarks. The `results/*.csv` tables and `results/traces/*.csv` are the source of the generated [RESULTS.md](RESULTS.md).

The CMake build uses the locally selected Clang20 compiler; the retained research runs were compiled with `g++ -O2`. CTest emits fresh validation CSVs under `build/validation_results`, keeping research output unchanged. The direct compiler commands and CMake targets run the same sources.

Validation completed successfully: all eight CMake targets build, and all four CTests pass. Research trace integrity checks also confirm each file ends at its reported completed time, contains the correct number of frames, and contains only finite values. Finite inexact optimization steps remain labeled in their summaries.
