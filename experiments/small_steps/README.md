# EL controls with small steps

Every method uses **one frozen E-then-L sweep per substep**, with 4, 8, or 16 substeps per 60 Hz frame (240, 480, or 960 Hz). Multipliers reset each substep. All comparisons use the same spring, gravity, initial conditions, physical duration, and substep size. This is the default integration schedule for this experiment.

Magnitude caps and sufficiently soft energy constraints suppress the extreme energy spikes seen with directional clipping and fixed relaxation. They do not recover accurate dynamics across all five cases. Strongly weakening energy feedback improves motion toward ordinary XPBD while allowing significant energy loss. The multi-spring extension remains conditional on solving the single-spring failures.

## Equations and controls

The original implicit objective and energy constraint are

\[
F(x)=\tfrac12(x-y)^T(M/h^2)(x-y)+U_s(x),\quad
H(x)=K((x-x_n)/h)+U_s(x)+V_g(x)-E_0,
\]
\[
y=x_n+hv_n+h^2a,\quad U_s=C^2/(2\alpha).
\]

For one spring, let `n=grad C`, `g=grad H`, `A=n^T M^-1 n+alpha/h^2`, `B=n^T M^-1 g`, and `D=g^T M^-1 g+beta/h^2`. Evaluate all rows once at the predictor. With zero initial multipliers, the position residual is exactly zero there. Curvature is dropped, as in the preceding EL experiments.

The candidate energy increment is

\[
\Delta\eta_{raw}=-(H+\tilde\beta\eta)/D,
\qquad\tilde\beta=\beta/h^2.
\]

Controls act in this order:

1. **Directional gate:** set the increment to zero when `grad F dot dx_test > 0`. Equality is allowed. Raw uses `dx_test=M^-1 g de`; effective uses `(M^-1 g-M^-1 n B/A)de`.
2. **Relaxation:** multiply the accepted increment by `omega`.
3. **Magnitude cap:** scale the increment further if `|M^-1 g de| > gamma S`, where `S=h|v_n|+h^2|a_net(x_n)|`. The net acceleration includes gravity and spring force. This scales an accepted increment; it is not the earlier oversized-update rejection.

The material visit then uses the accepted energy increment:

\[
\Delta\lambda=(-r_C-B\Delta\eta)/A,
\quad\Delta x=M^{-1}(n\Delta\lambda+g\Delta\eta).
\]

Positions move once after both multiplier visits; velocity is `(x_new-x_old)/h`. There is no reservoir, velocity projection, nonlinear line search, extra iteration, or relinearization between E and L.

For compliant energy, `beta=b*(E0-Vmin)` is a fixed physical compliance within each run and across substep-size comparisons. The dimensionless `b` controls strength: larger values mean weaker feedback. With reset eta and one sweep, compliance is a state-dependent relaxation of the hard correction. The full accumulated-multiplier term remains in the implementation and is independently tested.

## Experiment matrix

The same five elastic-pendulum cases as [the preceding study](../elastic_pendulum/README.md) are run for ten seconds: vertical oscillation, small swing, larger swing from rest, swing with tangential velocity, and coupled stretching/swinging. Mass is 1 kg, rest length 1 m, spring stiffness 100 N/m, and gravity 9.81 m/s². The energy target is the initial physical energy for the whole run, not the energy at the previous substep.

There are **465 distinct runs**: 31 configurations × five cases × three substep counts.

| Family | Parameters | Other settings |
|---|---|---|
| Ordinary XPBD | None | Material row only |
| Unguarded EL | None | No direction gate |
| Directional EL | None | Raw direction gate |
| Relaxed EL | omega = 0.1, 0.25, 0.5 | Raw gate |
| Capped EL | gamma = 0.001, 0.01, 0.1, 0.25, 0.5, 1, 2, 4 | Raw gate |
| Compliant-energy EL | b = 0.0001, 0.001, 0.01, 0.1, 1, 10, 100, 1000 | Raw gate |
| Effective-direction EL | None | Effective gate |
| Effective-direction + cap | Same eight gamma values | Effective gate |

The first sweep motivated tighter caps and weaker energy feedback, extending the initial parameter range. All tested values remain in the final matrix. Representatives use one parameter across all cases and substep counts. No configuration meets the exploratory criterion of maximum energy error below 5% everywhere. Motion-priority representatives minimize worst RMS position error normalized by the reference RMS excursion; the complete parameter scorecard exposes the energy tradeoff. Additional cap 0.1 and compliance b=10 trajectories show stronger energy feedback.

The selected trajectory pass repeats 135 of these runs, producing full substep traces at 480 Hz. It does not add independent configurations. Every repeated summary field matches the sweep exactly.

## Results

At **eight substeps per frame**, with **one EL sweep per substep**:

| Method | Vertical RMS position | Vertical max energy error | Coupled RMS position | Coupled max energy error |
|---|---:|---:|---:|---:|
| Ordinary XPBD | 0.0593 m | 87.54% | 0.1736 m | 60.90% |
| Directional EL | 0.2081 m | 1,706,019% | 0.7926 m | 3.89% |
| Directional + cap gamma=0.001 | 0.0536 m | 6.46% | 0.6020 m | 4.05% |
| Directional + cap gamma=0.1 | 0.2142 m | 5.79% | 0.7907 m | 3.89% |
| Directional + compliance b=10 | 0.1346 m | 5.92% | 0.6450 m | 3.94% |
| Directional + compliance b=100 | 0.0457 m | 11.36% | 0.5288 m | 6.71% |
| Directional + compliance b=1000 | 0.0409 m | 65.02% | 0.2434 m | 32.27% |

Energy error is `max |E-E0|/(E0-Vmin)`, where `Vmin=-10.2911805 J`. Position RMS compares with an independently refined RK4 trajectory. All 465 runs complete numerically; some contain enormous transient spikes. Completion is not a stability or accuracy criterion.

- **Fixed relaxation:** every tested relaxation factor still produces a large vertical spike somewhere in the matrix. Scaling a singularly large correction is insufficient.
- **Magnitude cap:** removes the extreme spikes in the tested configurations, but often produces substantial phase and trajectory error. Very tight caps can instead allow large energy loss.
- **Energy compliance:** sufficiently large beta suppresses spikes. Weaker feedback improves trajectories while allowing more energy error. The smallest beta values still spike.
- **Effective-direction gate:** gives exactly the same trajectories as the raw gate here, both with and without caps. At the single-spring predictor, `grad F` is parallel to `n`, and the effective slope is the raw slope multiplied by the positive factor `(alpha/h²)/A`. It adds no stabilization in this setup.
- **More small steps:** improve some configurations, but do not ensure monotone improvement. Directional-only EL still has an enormous vertical spike at 960 Hz; cap 0.001 has lower energy error but worse vertical phase error at 960 Hz than 480 Hz.

Gravity exerts torque, so angular momentum itself should change. The saved torque diagnostic is `L(t)-L(0)-integral(x cross mg)dt`, with the same trapezoidal quadrature applied to the reference. For the coupled case at 480 Hz, cap 0.001 gives a maximum torque discrepancy of 0.212 kg m²/s, compared with 0.0209 for ordinary XPBD. Closing the energy balance alone does not establish the right dynamics.

See [the derivation and independent audit](ANALYSIS.md) for why the two gates coincide and how compliance becomes adaptive relaxation in this schedule.

## Validation and cost

All 19 CTest checks pass. The new independent controls tests also pass with ASan and UBSan. Tests cover legacy parity, scalar multiplier/displacement consistency, cap behavior and order, compliant residuals with nonzero eta, singular rows, the effective-direction identity, and the exact backward-Euler phase/damping of the baseline at 4/8/16 substeps.

A zero-motion equilibrium exposes an additional distinction: unbounded EL can disturb equilibrium; a finite motion-scale cap preserves the tested equilibrium because its old net force and velocity vanish.

The sweep evaluates 2,604,000 substeps. Every run has one solver evaluation per substep. Baseline XPBD performs one scalar material solve; the EL variants perform one candidate energy solve and one material solve, except at explicitly singular energy rows. Rejected candidates still count as energy solves. At eight substeps per frame this is eight versus sixteen scalar solves, plus the direction/cap vector arithmetic. These operation counts are a cost proxy, not a wall-time benchmark.

Independent CSV audits verify counter accounting, cap bounds, nonpositive accepted gate slopes, raw/effective equivalence, and selected-run parity. Saved velocities reconstruct exactly and energies agree within 1.80e-15 scaled error. RK4 references use 16 subdivisions and are checked against 32 subdivisions at every substep.

## Reproduce

```bash
cmake --build build --target xpbd-small-steps-study xpbd-small-steps-tests
ctest --test-dir build --output-on-failure
./build/xpbd-small-steps-study results/small_steps_el/sweep
./build/xpbd-small-steps-study results/small_steps_el/selected baseline,unguarded,directional,relax_0.100000,cap_0.001000,cap_0.100000,soft_10.000000,soft_1000.000000,effective
python3 experiments/small_steps/plot_results.py results/small_steps_el/sweep results/small_steps_el/selected --output results/small_steps_el/plots
```

- [All 465 runs](../../results/small_steps_el/sweep/summary.csv)
- [Selected traces](../../results/small_steps_el/selected/traces.csv)
- [Reference refinement checks](../../results/small_steps_el/sweep/reference_checks.csv)

The plotter requires NumPy and Matplotlib; the validated local interpreter is `/home/blacknoir/Projects/mandelbrot-diffusion/.venv/bin/python`.

- [Energy and position traces for all five cases at 480 Hz](../../results/small_steps_el/plots/small_steps_energy_and_motion_480hz.png)
- [Raw, relaxed, and directional failure traces](../../results/small_steps_el/plots/difficult_variants_480hz.png)
- [All parameters at four substeps](../../results/small_steps_el/plots/all_parameters_substeps4.png)
- [All parameters at eight substeps](../../results/small_steps_el/plots/all_parameters_substeps8.png)
- [All parameters at sixteen substeps](../../results/small_steps_el/plots/all_parameters_substeps16.png)
- [Parameter tradeoffs](../../results/small_steps_el/plots/parameter_tradeoffs.png)
- [Substep convergence and scalar-solve cost](../../results/small_steps_el/plots/convergence_and_cost.png)
- [Selection rule and parameters](../../results/small_steps_el/plots/selection.json)
- [Plot audit](../../results/small_steps_el/plots/plot_audit.json)
