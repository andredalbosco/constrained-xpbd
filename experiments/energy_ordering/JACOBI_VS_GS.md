# Jacobi versus Gauss–Seidel: rerun

This rerun adds independent Jacobi corrections to the existing frozen-system
experiments. Both multiplier corrections read the same old state before their
position contributions are summed. Gauss–Seidel uses each new multiplier
correction in the following row, including when positions are applied together
at the end. These are different algorithms.

![Hard energy, released spring and pendulum](../../results/frozen_schedules_jacobi/plots/jacobi_vs_gs_hard_retain_rx_rest.png)

![Compliant energy, released spring and pendulum](../../results/frozen_schedules_jacobi/plots/jacobi_vs_gs_soft_001_retain_rx_rest.png)

## Equations tested

For one frozen, curvature-dropped linearization, elimination of position gives

\[
\begin{bmatrix}A&B\\B&D\end{bmatrix}
\begin{bmatrix}\Delta\lambda\\\Delta\eta\end{bmatrix}
=\begin{bmatrix}f_L\\f_E\end{bmatrix}.
\]

The total increments start at zero in each outer sweep. One Jacobi pass is

\[
\Delta\lambda=f_L/A,\qquad \Delta\eta=f_E/D.
\]

One material-then-energy GS pass (LE) is

\[
\Delta\lambda=f_L/A,\qquad
\Delta\eta=(f_E-B\Delta\lambda)/D.
\]

Both apply
\[
\Delta x=b+M^{-1}(J^T\Delta\lambda+g\Delta\eta),
\quad b=-M^{-1}h^2r_x
\]
when retaining the position residual; dropping it sets b=0. The right-hand
sides include b, and b is applied exactly once. The next outer sweep evaluates
the gradients and residuals again. There is no curvature term, reservoir,
velocity rescaling, damping of the Jacobi iteration, or change to the energy
target. Energy is the analytical endpoint kinetic plus potential energy;
the target remains its initial value.

## Measured behavior

The primary slice uses h=1/60 s, five outer sweeps, one frozen pass per sweep,
retained r_x, and ten seconds of requested motion. Energy error below is
max |E-E_initial|/Escale; Escale is the initial energy above the minimum
potential. Position error is RMS distance from the continuous RK4 reference.
LE means material then energy; EL means energy then material.

| Case, hard energy | Jacobi | GS LE | GS EL |
|---|---|---|---|
| Released spring | Guard after 7/600 accepted steps | Completely frozen, zero energy error | Completes with enormous energy error, 1.72e8 |
| Radially moving spring | Guard after 44/600 accepted steps | Energy error 7.29e-12, position RMS 0.313 m | Energy error 1.93e4 |
| Tangentially moving spring | Energy error 3.42e-5, position RMS 1.5075 m | Energy error 8.79e-11, position RMS 1.5075 m | Energy error 3.21e-8, position RMS 1.5075 m |
| Released pendulum | Energy error 6.87, position RMS 0.109 m | Energy error 1.00e8, position RMS 14.2 m | Energy error 0.0105, position RMS 0.0200 m |
| Moving pendulum | Energy error 25.9, position RMS 0.164 m | Energy error 1.21e4, position RMS 0.664 m | Energy error 0.00245, position RMS 0.113 m |

With energy compliance beta/Escale=0.01, all tested schedules complete without
guards. For Jacobi, released/moving pendulum energy errors are 0.00336/0.00426
and position errors are 0.0270/0.0126 m. These are better than both GS LE and EL
on this slice, but Jacobi leaves larger rod and stationarity residuals: its
maximum rod errors are 0.000415/0.000294 m, versus about 2e-9/1e-9 m for GS EL.
The compliant released spring still has incorrect motion: Jacobi travels only
0.500 m in ten seconds, versus the reference's 15.947 m, and has position RMS 0.299 m.
Compliance permits energy error; completion does not establish conservation.

Dropping r_x changes the result substantially. Hard-energy Jacobi then completes
the released spring with energy error 4084. Its released/moving pendulum energy
errors are only 0.00432/0.00199, but position errors are 0.586/0.486 m. Small
energy error alone is insufficient to judge the motion.

The first released-spring step explains the sharp contrast with LE. Initially
the energy residual is zero, so Jacobi's energy increment is zero: it takes
the ordinary XPBD spring step. At that new position the energy gradient
vanishes while energy is below the target, making the next hard linearized
energy row incompatible. LE instead accounts for the material correction
inside its energy solve, which cancels the position change at this rest state.
Neither behavior resolves the underlying direct-energy formulation problem.

## All graphs and data

Every individual method chart shows actual energy, x(t), and y(t) for all five
initial states. Gray dashed curves show the RK4 reference; light gray shows
ordinary XPBD. Runaway panels use explicitly labeled symlog axes; stopped
trajectories end at their last accepted state. Incompatible rows are skipped
and counted, not silently regularized.

| Comparison | Retain r_x | Drop r_x |
|---|---|---|
| Hard, released cases | [Graph](../../results/frozen_schedules_jacobi/plots/jacobi_vs_gs_hard_retain_rx_rest.png) | [Graph](../../results/frozen_schedules_jacobi/plots/jacobi_vs_gs_hard_drop_rx_rest.png) |
| Hard, tangentially moving cases | [Graph](../../results/frozen_schedules_jacobi/plots/jacobi_vs_gs_hard_retain_rx_moving.png) | [Graph](../../results/frozen_schedules_jacobi/plots/jacobi_vs_gs_hard_drop_rx_moving.png) |
| Compliant, released cases | [Graph](../../results/frozen_schedules_jacobi/plots/jacobi_vs_gs_soft_001_retain_rx_rest.png) | [Graph](../../results/frozen_schedules_jacobi/plots/jacobi_vs_gs_soft_001_drop_rx_rest.png) |
| Compliant, tangentially moving cases | [Graph](../../results/frozen_schedules_jacobi/plots/jacobi_vs_gs_soft_001_retain_rx_moving.png) | [Graph](../../results/frozen_schedules_jacobi/plots/jacobi_vs_gs_soft_001_drop_rx_moving.png) |

| Method | Hard, retain r_x | Hard, drop r_x | Compliant, retain r_x | Compliant, drop r_x |
|---|---|---|---|---|
| Jacobi | [Graph](../../results/frozen_schedules_jacobi/plots/jacobi_hard_retain_rx.png) | [Graph](../../results/frozen_schedules_jacobi/plots/jacobi_hard_drop_rx.png) | [Graph](../../results/frozen_schedules_jacobi/plots/jacobi_soft_001_retain_rx.png) | [Graph](../../results/frozen_schedules_jacobi/plots/jacobi_soft_001_drop_rx.png) |
| Joint | [Graph](../../results/frozen_schedules_jacobi/plots/joint_hard_retain_rx.png) | [Graph](../../results/frozen_schedules_jacobi/plots/joint_hard_drop_rx.png) | [Graph](../../results/frozen_schedules_jacobi/plots/joint_soft_001_retain_rx.png) | [Graph](../../results/frozen_schedules_jacobi/plots/joint_soft_001_drop_rx.png) |
| GS LE | [Graph](../../results/frozen_schedules_jacobi/plots/gs_le_hard_retain_rx.png) | [Graph](../../results/frozen_schedules_jacobi/plots/gs_le_hard_drop_rx.png) | [Graph](../../results/frozen_schedules_jacobi/plots/gs_le_soft_001_retain_rx.png) | [Graph](../../results/frozen_schedules_jacobi/plots/gs_le_soft_001_drop_rx.png) |
| GS EL | [Graph](../../results/frozen_schedules_jacobi/plots/gs_el_hard_retain_rx.png) | [Graph](../../results/frozen_schedules_jacobi/plots/gs_el_hard_drop_rx.png) | [Graph](../../results/frozen_schedules_jacobi/plots/gs_el_soft_001_retain_rx.png) | [Graph](../../results/frozen_schedules_jacobi/plots/gs_el_soft_001_drop_rx.png) |
| GS LEL | [Graph](../../results/frozen_schedules_jacobi/plots/gs_lel_hard_retain_rx.png) | [Graph](../../results/frozen_schedules_jacobi/plots/gs_lel_hard_drop_rx.png) | [Graph](../../results/frozen_schedules_jacobi/plots/gs_lel_soft_001_retain_rx.png) | [Graph](../../results/frozen_schedules_jacobi/plots/gs_lel_soft_001_drop_rx.png) |
| GS ELE | [Graph](../../results/frozen_schedules_jacobi/plots/gs_ele_hard_retain_rx.png) | [Graph](../../results/frozen_schedules_jacobi/plots/gs_ele_hard_drop_rx.png) | [Graph](../../results/frozen_schedules_jacobi/plots/gs_ele_soft_001_retain_rx.png) | [Graph](../../results/frozen_schedules_jacobi/plots/gs_ele_soft_001_drop_rx.png) |

[All methods, released cases](../../results/frozen_schedules_jacobi/plots/overview_hard_retain_rx_rest.png)
and [moving cases](../../results/frozen_schedules_jacobi/plots/overview_hard_retain_rx_moving.png)
include the joint solve and both three-visit GS orders.

## Coverage and reproduction

The matrix contains 1,230 runs: the previous 1,110 runs plus 120 Jacobi runs.
Five initial states, two timesteps (1/60, 1/240), and three outer sweep counts
(1, 5, 20) cover hard/compliant energy and both r_x choices. GS retains both
immediate and summed position application; Jacobi applies the sum after both
independent solves. The plotted slice contains 205 traces.

There are 1,061 ordinary completions, 116 completions with incompatible updates,
and 53 guarded stops. All stops use hard energy. Jacobi contributes 93 ordinary
completions, 13 completions with incompatible updates, and 14 guarded stops.
All preexisting CSV rows are unchanged in the rerun. The staging and repository
release executables produce byte-identical results.

All 13 CTests pass, including baseline spring, pendulum, grid, and cloth.
Independent Jacobi tests check simultaneous old increments, repeated frozen
passes, singular rows, the exact released-spring first pass, and convergence
to an independently assembled mixed-system solve. ASan and UBSan pass with
LeakSanitizer disabled for the execution environment. Plot generation checks
saved energy against positions and velocities; RK4 references are checked by
halving their substep size.

~~~sh
cmake --build build -j
ctest --test-dir build --output-on-failure
./build/xpbd-frozen-schedules-study results/frozen_schedules_jacobi
python3 experiments/energy_ordering/plot_schedules.py results/frozen_schedules_jacobi
~~~

The plotting command requires NumPy and Matplotlib.

[Summary CSV](../../results/frozen_schedules_jacobi/summary.csv) ·
[Traces CSV](../../results/frozen_schedules_jacobi/traces.csv) ·
[Plot audit](../../results/frozen_schedules_jacobi/plots/plot_audit.json)
