# Traditional XPBD: spring, pendulum, and double spring pendulum

This study returns to the ordinary inertial predictor, resets all multipliers every small timestep, and compares global endpoint-energy corrections. It does not use warm starting, a symplectic predictor, work balance, or a reservoir.

The 378 runs cover 18 methods × 7 initial conditions × 240/480/960 Hz, for 10 seconds. Each timestep performs one prescribed schedule, with no outer iteration. The energy target is the initial physical energy throughout the run. Parameters were fixed before inspecting results.

## Equations and approximations

For timestep h, old positions q and velocities v, diagonal mass matrix M, and constant gravity a:

```
y = q + h v + h² a
F(x) = (x-y)ᵀ M (x-y)/(2h²) + Σ_j C_j(x)²/(2α_j)
E(x) = (x-q)ᵀ M (x-q)/(2h²) + Σ_j C_j(x)²/(2α_j) − Σ_i m_i a·x_i
E(x) = E_initial
```

A rigid rod uses C=0 as a hard material constraint, with no elastic energy term. Every spring connects either the origin to particle 1 or particle 1 to particle 2. All spring energies use physical compliance α, not α̃. Define α̃=α/h², J=∂C/∂x, g=∇E, and scaled multipliers λ and η. Dualization of each soft spring gives the saddle objective

```
L(x,λ,η) = (x-y)ᵀM(x-y)/(2h²)
           − λᵀC(x)/h² − λᵀα̃λ/(2h²) − η(E(x)−E_initial)/h².
```

The stationarity residuals (with the x residual multiplied by h²) are

```
r_x = M(x-y) − Jᵀλ − gη
r_C = C(x) + α̃λ
r_E = E(x) − E_initial.
```

Linearize the equations, drop the multiplier-weighted constraint and energy Hessians, then drop r_x in the position update. The remaining frozen linear system is

```
[A  B] [Δλ] = −[r_C]        A = J M⁻¹ Jᵀ + α̃
[Bᵀ D] [Δη]    [r_E]        B = J M⁻¹ g, D = gᵀ M⁻¹ g
Δx = M⁻¹(JᵀΔλ + gΔη).
```

The energy gradient includes gravity and the kinetic energy reconstructed from endpoint displacement. There is no exact nonlinear energy projection. At a later nonlinear iterate r_x need not vanish; dropping it is an explicit approximation. In this experiment every schedule starts at y with λ=η=0.

## Schedules

- **XPBD:** ordinary forward per-spring Gauss–Seidel, updating positions and refreshing each spring row before its solve. This is the baseline.
- **Frozen XPBD:** one forward material GS pass on a single frozen linearization. This isolates the linearization choice from adding an energy row.
- **EL:** solve the energy row, then visit material rows using the updated energy multiplier.
- **LE:** visit material rows, then solve energy using their updated multipliers.
- **Jacobi:** solve every material and energy row from the original right-hand side, then sum position contributions.
- **Joint:** solve the full frozen block (2×2 for one mass; 3×3 for two). A normalized pivot test rejects near-rank-deficient blocks and falls back to a material pass.
- **LEL/ELE:** repeat the indicated row/block once within the same frozen system.

All energy schedules freeze gradients and right-hand sides for the entire schedule. Applying intermediate position increments with those same frozen rows is algebraically equivalent to accumulating them. These are not schedules that re-evaluate E and J after every position change.

For the double spring pendulum, L means one forward GS visit to both spring rows, not an exact material-block solve. Updating the second row can reopen the first row's linear residual. EL/LEL therefore do not satisfy both material equations exactly after a single sweep.

## Guards

- **Raw:** use the schedule directly.
- **Direction:** compare its total position correction against the frozen-XPBD correction. Reject the extra correction if ∇F·extra>0.
- **EL raw gate:** the user's original direction test, applied to the raw energy contribution before the following material solve. Reject when ∇F·Δx_E>0. This differs from the total-extra test.
- **EL raw gate + cap:** the same gate, with a norm cap on the raw energy contribution before L.
- **Repair + cap:** compare the candidate to frozen XPBD. Let H₀ be the energy residual predicted after the baseline and δH the predicted extra change. Reject unless δH reduces |H₀|; limit the predicted repair to half |H₀| and the extra position correction to 10% of the motion scale. Blend all multipliers consistently with the accepted extra correction.

The mass-weighted motion scale is `h ||v||_M + h² ||a_old||_M`, where a_old includes gravity and elastic forces (for the rigid pendulum it uses gravity). The defaults are repair fraction 0.5 and cap 0.1. These are numerical safeguards, not a derived compliant energy constraint. They bound the linearized repair and do not guarantee a bound on actual nonlinear energy error.

## Cases and references

All moving masses are 1 kg, rest lengths are 1 m, soft spring stiffness is 100 N/m, and gravity where present is 9.81 m/s². The seven cases are radial spring release, spring orbit, rigid pendulum from rest, rigid pendulum with tangential velocity, and small/large/mixed motion of the two-spring pendulum. Exact initial states are in `results/cases.csv`.

The double spring equilibrium has edge lengths 1.1962 and 1.0981 m. The small case starts from those edge lengths with angles 0.15 and −0.10 radians from downward vertical. The large case adds extensions 0.1/0.2 m with angles 0.65/−0.45 radians. The mixed case adds extensions 0.05/0.1 m and nonzero radial/tangential velocities. Both masses move in 2D; the first spring is anchored at the origin.

References use Cartesian RK4 for springs and angular RK4 for the rigid pendulum. Independent tests validate the radial baseline against its exact backward-Euler map, the pendulum's length and trajectory convergence, and a double-spring normal mode against its analytic solution. Reference substep refinement from 16 to 32 subdivisions across all cases/frequencies changes positions by at most 1.86e−10 m and velocities by 1.39e−9 m/s; maximum relative reference energy drift is 2.31e−13.

## Results

Representative 480 Hz runs, 10 seconds, one schedule per timestep:

| Case | Method | Maximum energy error | Position RMS |
|---|---|---:|---:|
| Radial spring | XPBD | 87.54% | 0.07416 m |
| | EL repair + cap | 2.009% | 0.11068 m |
| | LE repair + cap | 0.1539% | 0.009590 m |
| Rigid pendulum, rest | XPBD | 21.99% | 0.03765 m |
| | EL repair + cap | 0.006179% | 0.004660 m |
| | LE repair + cap | 0.01296% | 0.004905 m |
| Double spring, small swing | XPBD | 41.58% | 0.01176 m |
| | EL repair + cap | 1.481% | 0.02630 m |
| | LE repair + cap | 0.1542% | 0.01660 m |

Energy error is divided by initial energy above the analytic minimum potential, making gravity cases comparable without dividing by a negative or arbitrary energy zero. Position RMS is `sqrt(sum_step,sum_body ||x−xref||² / (completed_steps*body_count))`; it includes both double-pendulum masses and both coordinates.

Key observations:

1. **Raw LE and joint freeze the radial spring exactly.** Their energy error is zero, but the path length is also zero and position RMS is 0.307 m. Exact energy conservation alone does not validate the dynamics.
2. **Direction clipping alone is insufficient.** At 480 Hz the EL raw gate has an energy spike of 1.33e9 percent on radial release. A run can remain finite after a severe spike; completion is not stability. Raw Jacobi at 240 Hz is the only run to trigger the large-value abort criterion, but several other raw runs are physically unusable.
3. **The bounded correction makes LE useful for radial release.** It gives substantially better energy and motion than XPBD. EL repair preserves amplitude but has appreciable phase error: the measured period is about 0.63592 s, versus 0.62919 s for LE repair and the exact 0.62832 s.
4. **The rigid pendulum benefits.** EL repair + cap improves energy and motion by large factors, although it does not enforce angular-momentum balance. Its maximum gravity-torque balance error is 0.0662 kg m²/s versus 0.0193 for XPBD in the rest case. The maximum rod-length error is 5.1e−11 m for EL repair and 3.1e−9 m for LE repair in this rest case.
5. **Two springs expose the limits.** At 480 Hz all three double-spring cases have lower position error with ordinary XPBD than with the tested energy-corrected variants listed above. For the small case, LE repair reduces energy error from 41.58% to 0.1542%, but increases position RMS by about 41%. The maximum angular-balance defect grows from 0.00517 to 0.1608 kg m²/s.
6. **The orbit is a stronger counterexample.** LE repair reduces energy error to 0.02263%, but position RMS increases from 0.05790 to 0.91962 m. Its maximum angular-momentum defect is 0.3518 versus 0.01077 kg m²/s. Correcting total energy does not control its distribution between radial and angular motion.
7. **Smaller steps help, without eliminating the tradeoff at the tested resolutions.** In the small double case at 960 Hz, ordinary XPBD RMS is 0.00627 m and LE repair RMS is 0.00854 m. The mixed case is a limited exception at 960 Hz: LE repair slightly beats baseline position RMS (0.04493 versus 0.04742 m).
8. **Joint coupling offers little benefit here.** With the same repair + cap, its results closely match LE or LEL, while requiring a block solve. This is an accuracy observation for these tests, not a performance benchmark.

The evidence supports bounded EL for the rigid pendulum and bounded LE for the radial spring. It does not establish a general energy-corrected spring-chain method. A next experiment should address the direction of the energy correction and angular balance or local work balance, rather than judge success by total energy alone.

## Validation and artifacts

`tests.cpp` independently checks forces, finite-difference gradients, incidence Jacobians, hand-derived two-spring Schur entries, nonzero-multiplier compliant residuals, fresh/frozen baseline behavior, all frozen schedules, normalized block solves against pivoted long-double Gaussian elimination, and clipping/repair/cap bounds with consistent multiplier reconstruction. These tests pass optimized and under AddressSanitizer/UndefinedBehaviorSanitizer. The main study checks reference refinement before running each configuration. Plotting independently reconstructs energies and trajectory metrics from saved traces.

Raw output is in `results/`; `summary.csv` contains every run including failure status and attempted energy on failure. `traces.csv` stores every 480 Hz timestep for seven selected methods, with both particles, reference states, multipliers, and torque integrals. An aborted run's metrics cover its accepted prefix; never compare that prefix as a successful 10-second run.

```bash
cmake -S experiments/traditional_cases -B /tmp/xpbd-traditional-build -DCMAKE_BUILD_TYPE=Release
cmake --build /tmp/xpbd-traditional-build -j4
ctest --test-dir /tmp/xpbd-traditional-build --output-on-failure
/tmp/xpbd-traditional-build/xpbd-traditional-study experiments/traditional_cases/results 10
```

The solver and tests are C++17. Python is used only for scientific plots. Source dependencies are the existing experiment vector utilities. The core production `xpbd.hpp` and the earlier studies are unchanged.

Render the figures using any Python environment with NumPy and Matplotlib:

```bash
python3 experiments/traditional_cases/plot_results.py experiments/traditional_cases/results --output experiments/traditional_cases/plots
```

On this workstation, the tested runtime is `/home/blacknoir/Projects/mandelbrot-diffusion/.venv/bin/python`. `plots/overview.png` compares energy and position error for the three primary cases; the case figures show positions, including both masses of the double pendulum. The three frequency scorecards cover all 18 methods and all seven initial conditions. `plots/plot_audit.json` records numerical reconstruction checks, and `plots/visual_qa.json` records visual inspection. The audits cover 378 summaries and 49 trace groups, with reconstructed position RMS agreeing within 1.78e−15 m.

## Focused EL figures

`el_plots/` compares ordinary XPBD with four EL variants, one method per row: raw, the original energy-direction clip, direction clip plus a motion cap, and bounded repair plus a motion cap. Each figure shows energy error and the actual position or angle against the reference. The rigid pendulum works well with bounded EL. For radial spring release, bounded EL maintains amplitude but accumulates phase error; at 480 Hz its position RMS remains worse than ordinary XPBD. The additional tangential-velocity cases are saved alongside these primary cases.

```bash
/tmp/xpbd-traditional-build/xpbd-traditional-study experiments/traditional_cases/el_results 10 --el-focus
/home/blacknoir/Projects/mandelbrot-diffusion/.venv/bin/python experiments/traditional_cases/plot_el.py experiments/traditional_cases/el_results --output experiments/traditional_cases/el_plots
```

The focused study has 72 runs at three frequencies, including a combined-direction variant retained in the CSV. Their energy, position, angular-balance and completion metrics reproduce the corresponding original-sweep results exactly. The plot audit independently reconstructs physical energy and position RMS for all 24 saved trajectory groups. Each graph row has its own vertical range; large energy spikes use a labelled symmetric logarithmic scale.

## EL/LE with positions reconstructed between visits

`relinearized.hpp` implements the requested E → position → rebuilt material row → position sequence, and its LE counterpart. Each visit evaluates its residual and gradient at the current corrected positions and solves one linearized scalar equation. Curvature and the current `r_x` are still dropped. Reconstructing positions while retaining the old rows is already algebraically equivalent to frozen EL/LE; rebuilding the rows changes the approximation.

The 108-run comparison covers radial spring release, spring orbit, rigid pendulum release and a kicked rigid pendulum, at 240/480/960 Hz for 10 seconds. It includes XPBD, frozen/rebuilt EL and LE, both raw and with the same whole-proposal repair/cap. The bounded policy blends actual proposed positions and all multiplier increments consistently. Earlier default results are preserved; the 60 overlapping original-method summaries reproduce the earlier study exactly.

Representative bounded results at 480 Hz:

| Case | Schedule | Frozen maximum energy error | Rebuilt maximum energy error | Frozen position RMS | Rebuilt position RMS |
|---|---|---:|---:|---:|---:|
| Radial spring | EL | 2.009% | 2.009% | 11.068 cm | 11.068 cm |
| Radial spring | LE | 0.1539% | 0.1462% | 0.959 cm | 0.500 cm |
| Rigid pendulum, rest | EL | 0.006179% | 0.006179% | 0.466 cm | 0.466 cm |
| Rigid pendulum, rest | LE | 0.01296% | 0.008457% | 0.491 cm | 0.748 cm |

Rebuilding EL does not repair its raw spring spike: the problematic energy correction precedes the material visit. While the spring remains on the positive radial branch, its material constraint is affine, so rebuilding that row produces exactly the same equation. Rebuilt LE avoids the exact frozen branch on release, but its first material projection makes the energy gradient zero, requiring a singular-gradient skip. Raw rebuilt LE still has substantial phase error. None of the compared energy methods gives acceptable orbit dynamics; low energy error remains insufficient evidence.

```bash
/tmp/xpbd-traditional-build/xpbd-traditional-study experiments/traditional_cases/el_le_results 10 --el-le-focus
/home/blacknoir/Projects/mandelbrot-diffusion/.venv/bin/python experiments/traditional_cases/plot_el_le.py experiments/traditional_cases/el_le_results --output experiments/traditional_cases/el_le_plots
```

`el_le_plots/` contains nine figures, including energy and actual motion for every case and frozen/rebuilt comparisons. Its independent audit covers 36 complete saved trajectories: physical-energy reconstruction agrees within 2.58e-15 of the energy scale, and position RMS within 1.12e-15 m. The C++ tests include hand-derived refreshed pendulum sequences, affine spring equivalence, zero-gradient handling, and accepted-position/multiplier blending. Both CTests pass; the refreshed tests also pass AddressSanitizer/UndefinedBehaviorSanitizer.

Ideas extracted from the supplied handwritten-notes transcription, with their derivation and validation limits, are recorded in `transcript_ideas.md`. Retaining `r_x` after rebuilding is a distinct follow-up, not a variant tested in this comparison.
