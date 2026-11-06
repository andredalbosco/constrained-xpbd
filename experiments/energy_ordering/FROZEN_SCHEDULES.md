# Energy and position graphs for the five frozen solve orderings

The [Jacobi versus GS rerun](JACOBI_VS_GS.md) adds independent simultaneous
multiplier corrections and regenerates the comparison. This page records the
earlier GS-only matrix.

Each chart shows actual mechanical energy and Cartesian position against time.
Gray dashed curves are the continuous RK4 reference; light gray is ordinary
XPBD; the dotted energy line is the fixed initial energy target.
The primary plots retain the position residual and enforce hard energy.
They use h=1/60 s, five outer sweeps per timestep, and one frozen GS pass per
sweep over ten seconds.

![Resting spring and pendulum](../../results/frozen_schedules/plots/overview_hard_retain_rx_rest.png)

![Moving spring and pendulum](../../results/frozen_schedules/plots/overview_hard_retain_rx_moving.png)

L means material multiplier lambda; E means energy multiplier eta.
The equations and gradients stay frozen during each pass. After completing that
pass, the next outer sweep evaluates a new linearization. Every method drops
curvature. Immediate application of position contributions does not reevaluate
gradients within the pass.

## A graph for each solve choice

Each individual chart has five rows: spring at rest, spring moving radially,
spring moving tangentially, pendulum at rest, and pendulum moving tangentially.
Its three columns are total energy, horizontal position, and vertical position.
For GS, solid curves apply the summed correction; magenta dashed curves apply
each position contribution immediately.

| Solve choice | Hard energy, retain rx | Hard energy, drop rx | Soft energy, retain rx | Soft energy, drop rx |
|---|---|---|---|---|
| Joint | [Graph](../../results/frozen_schedules/plots/joint_hard_retain_rx.png) | [Graph](../../results/frozen_schedules/plots/joint_hard_drop_rx.png) | [Graph](../../results/frozen_schedules/plots/joint_soft_001_retain_rx.png) | [Graph](../../results/frozen_schedules/plots/joint_soft_001_drop_rx.png) |
| Material → energy | [Graph](../../results/frozen_schedules/plots/gs_le_hard_retain_rx.png) | [Graph](../../results/frozen_schedules/plots/gs_le_hard_drop_rx.png) | [Graph](../../results/frozen_schedules/plots/gs_le_soft_001_retain_rx.png) | [Graph](../../results/frozen_schedules/plots/gs_le_soft_001_drop_rx.png) |
| Energy → material | [Graph](../../results/frozen_schedules/plots/gs_el_hard_retain_rx.png) | [Graph](../../results/frozen_schedules/plots/gs_el_hard_drop_rx.png) | [Graph](../../results/frozen_schedules/plots/gs_el_soft_001_retain_rx.png) | [Graph](../../results/frozen_schedules/plots/gs_el_soft_001_drop_rx.png) |
| Material → energy → material | [Graph](../../results/frozen_schedules/plots/gs_lel_hard_retain_rx.png) | [Graph](../../results/frozen_schedules/plots/gs_lel_hard_drop_rx.png) | [Graph](../../results/frozen_schedules/plots/gs_lel_soft_001_retain_rx.png) | [Graph](../../results/frozen_schedules/plots/gs_lel_soft_001_drop_rx.png) |
| Energy → material → energy | [Graph](../../results/frozen_schedules/plots/gs_ele_hard_retain_rx.png) | [Graph](../../results/frozen_schedules/plots/gs_ele_hard_drop_rx.png) | [Graph](../../results/frozen_schedules/plots/gs_ele_soft_001_retain_rx.png) | [Graph](../../results/frozen_schedules/plots/gs_ele_soft_001_drop_rx.png) |

[Baseline XPBD](../../results/frozen_schedules/plots/baseline.png) is included
as a sixth row in the overview. Soft energy uses beta/Escale=0.01, where Escale
is initial energy above the minimum potential.

The following larger charts put all six methods on separate rows for one case:

- [Spring at rest](../../results/frozen_schedules/plots/spring_rest_hard_retain_rx.png)
- [Spring with radial velocity](../../results/frozen_schedules/plots/spring_radial_hard_retain_rx.png)
- [Spring with tangential velocity](../../results/frozen_schedules/plots/spring_tangent_hard_retain_rx.png)
- [Pendulum at rest](../../results/frozen_schedules/plots/pendulum_rest_hard_retain_rx.png)
- [Pendulum with tangential velocity](../../results/frozen_schedules/plots/pendulum_tangent_hard_retain_rx.png)

## Reading the failures

Bounded panels share ranges across methods. Runaway panels explicitly use a
symmetric logarithmic scale so every accepted value remains visible.
Crosses mark guarded stops, with no continuation beyond the last accepted
timestep. Incompatible local rows are skipped and counted in chart titles;
completion alone does not imply equation closure or acceptable motion.

For the primary resting-spring case, joint, LE, and ELE freeze exactly.
LEL travels only about 1e-7 over ten seconds. EL instead develops very large
energy and position errors. In the pendulum, ordering changes energy spikes,
phase error, and length error. The plots display those outcomes rather than
ranking methods by energy error alone.

Immediate and deferred application are algebraically equivalent for frozen
gradients when all pending contributions enter subsequent right-hand sides.
Independent tests verify this on regular systems. The two full trajectories
can diverge through roundoff in unstable runs. The worst local comparison in
the wider matrix combines opposing displacements of about 4.6e17 with a
material correction of 3.2166; one summation order loses that small correction.
That case was already grossly unstable. Across compliant runs, local
application differences stay below 2.5e-16.

## Test coverage and reproduction

The C++ matrix has 1,110 runs: five initial states, two timesteps (1/60 and
1/240), three outer sweep counts (1, 5, 20), both residual choices, hard/soft
energy, and both application styles for GS. It contains 968 ordinary
completions, 103 completions with incompatible updates, and 39 guarded stops.
All guarded stops use hard energy. The plotted slice contains 185 traces.

All 13 repository CTests pass. New independent tests cover each scalar row
order, repeated visits, one-time position-residual correction, immediate versus
deferred application, singular rows, baseline parity, and convergence to an
independently assembled mixed-system solve. The new test executable also passes
ASan and UBSan, with leak detection disabled for the execution environment.
Reference trajectories are checked by halving their RK4 substep size.

Plot generation independently reconstructs physical energy from the saved
positions and velocities. The maximum scaled difference is about 1.1e-15.
The release executable and staging executable produce byte-identical CSVs.

From the repository root:

~~~sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
ctest --test-dir build --output-on-failure
./build/xpbd-frozen-schedules-study results/frozen_schedules
python3 experiments/energy_ordering/plot_schedules.py results/frozen_schedules
~~~

Plots require NumPy and Matplotlib. The plotting script defaults to the
results-directory plots subfolder.

[Raw summary](../../results/frozen_schedules/summary.csv) ·
[Position/velocity/energy traces](../../results/frozen_schedules/traces.csv) ·
[Plot audit](../../results/frozen_schedules/plots/plot_audit.json)
