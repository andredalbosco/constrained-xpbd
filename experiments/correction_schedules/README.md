# Numerical application of the energy correction

This research pass compares numerical solvers and acceptance rules for the same
endpoint-energy constrained XPBD formulation. It started at 2026-10-25 10:55:53
UTC. The user requested wrapping up at approximately 11:22 UTC, superseding the
earlier request to continue for at least an hour. No multi-spring chain was run.

The strongest simple candidate in these tests is **LE with fractional linear
repair and a displacement cap**. It improves released springs and rigid
pendulums, and performs very well on a true circular spring orbit. It is not a
general energy-and-momentum conserving solver: eccentric spring orbits and
coupled spring-pendulum motion still have large phase and torque errors.

## Fixed formulation and approximations

With predictor y=x_old+h v_old+h² a,

    F(x) = (x-y)^T M (x-y)/(2h²) + U(x)
    U(x) = C(x)^T alpha^-1 C(x)/2
    E(x) = (x-x_old)^T M (x-x_old)/(2h²) + U(x) - (M a)^T x.

The energy constraint is E(x)=E_initial. A rigid pendulum has a hard distance
row and U=0; its physical potential is gravity. Every small substep resets the
material and energy multipliers. The main study linearizes once at the
predictor, drops constraint curvature and r_x, and performs one scheduled pass.
There is no reservoir, nonlinear root solve, or separate velocity correction.

The two-row frozen system is

    [ A B ] [delta_lambda] = -[r_C]
    [ B D ] [delta_eta   ]    [r_E]
    delta_x = M^-1 (n delta_lambda + g delta_eta).

`ANALYSIS.md` gives the precise rows, schedule equations, branch degeneracy,
filter equivalences, and torque attribution. Main energy compliance is zero;
this pass varies numerical application, not physical stiffness or energy target.

## Methods actually compared

- Jacobi: compute both row increments from the original residual, then add.
- EL and LE: scalar Gauss-Seidel with the second row using the first increment.
- LEL and ELE: revisit the first row of that same frozen system.
- Alternating: EL and LE on alternating substeps.
- Joint: solve the coupled 2x2 system.
- Guarded joint: use EL when the normalized block is nearly dependent.
- Spectral: apply a capped inverse square root, or capped inverse, of the
  diagonally normalized Schur block. This is an analogy to spectral
  conditioning, not the literal Muon optimizer.

Each candidate is compared with ordinary material-only XPBD. Six acceptance
rules are tested: raw, objective-direction rejection, displacement cap,
direction plus cap, fractional linear repair, and repair plus cap.

For the default repair rule, let d0 be the material-only correction and d the
candidate correction beyond d0. With H0=r_E+g^T d0, accept only a direction
that reduces |H0| in the frozen model, and limit its predicted repair to half
|H0|. Also cap |d| at 0.1 times h|v_old|+h²|f_old|/m. Both multipliers and
their total position contribution are blended toward the material-only result.
This differs from rejecting just the raw energy-force contribution inside EL.

## Completed experiments

- **1,809 primary runs:** 67 configurations, nine cases, 240/480/960 Hz, ten
  seconds each. Two raw alternating-order runs hit the explicit safety limit.
- **270 selected reruns:** full 480 Hz position/energy/correction traces;
  90 trajectory groups, 432,000 noninitial trace steps.
- **1,728 parameter runs:** repair fraction 0.1/0.25/0.5/0.75/1 and cap
  0.001/0.01/0.1/1, including controls. All completed. Tuning these does not
  remove the eccentric-orbit or coupled-motion failure.
- **162 runs lasting 100 seconds:** six controls/candidates on the same nine
  cases and three step sizes. All completed. Low energy error persists, but
  phase error accumulates; completion is not evidence of accurate dynamics.
- **16 fine-step holdouts:** eccentric orbit and coupled motion down to
  h=1/15360. Corrected motion eventually converges roughly first order, while
  energy error falls roughly second order. Its motion-error coefficient remains
  much larger than baseline for those cases.
- **90 circular-speed/direction holdouts:** the true circle is a major success;
  radial-dominated eccentric motion is not. A separately labeled angular-neutral
  direction projection did not fix the failure and is outside the original
  multiplier-system comparison.
- **Refreshed-row exploratory runs:** re-evaluate each row after applying the
  preceding position correction, with whole-proposal or local filtering. Every
  visit still solves a linearized equation. These cost two or three row
  evaluations, versus one in the main study. See the validation note below.

The nine primary cases are a vertical gravity spring, small/large spring swings,
a spring swing with tangential velocity, coupled radial/tangential motion,
a released spring without gravity, an eccentric spring orbit, and two rigid
pendulums. Maximum scope in this pass is a single moving bob. A small chain
remains future work; no cloth or grid was added.

## Representative results

Ten seconds at 480 Hz; energy error is maximum absolute error divided by
initial energy above the minimum potential. Position is RMS distance to the
convergence-checked RK4 trajectory.

| Case | XPBD position RMS | LE repair+cap RMS | LE energy error |
|---|---:|---:|---:|
| Released spring | 0.07416 m | 0.00959 m | 0.154% |
| Vertical gravity spring | 0.05933 m | 0.00683 m | 0.146% |
| Rigid pendulum from rest | 0.03765 m | 0.00491 m | 0.0130% |
| Swing with tangential velocity | 0.07507 m | 0.00611 m | 0.0147% |
| Large spring swing from rest | 0.07897 m | 0.19531 m | 0.0470% |
| Eccentric spring orbit | 0.05790 m | 0.91962 m | 0.0226% |
| Coupled spring-pendulum motion | 0.17361 m | 0.59375 m | 0.0305% |
| True circular orbit holdout | 0.8440 m | 0.000540 m | 0.00364% |

These are not rankings by energy alone. In the eccentric orbit, angular
momentum increases from 0.88 to 1.23176 and phase leads the reference by 2.0714
radians. Energy correction transfers radial oscillation energy into rotation.
The cap is inactive there; simply tuning that cap cannot repair this behavior.

## What the analysis establishes

1. At rest, the exact endpoint energy and implicit objective differ only by a
   constant: grad E=grad F+M v_old/h. Thus the energy equality fixes the objective
   value. This is an underlying branch-selection degeneracy, not merely a bad
   Gauss-Seidel implementation. Raw LE and joint can conserve energy by freezing
   a released spring exactly.
2. Fractional repair can avoid fully taking that frozen solution. Some filtered
   schedules become algebraically identical; their overlapping plots do not
   count as independent successful methods.
3. Accurate total energy does not imply correct mode partition or angular
   momentum. Circular and eccentric orbit holdouts demonstrate the distinction.
4. Reassembling a joint position update in radial/tangential coordinates reduces
   floating-point cancellation, but does not remove physical freezing. In 648
   near-dependent/blended fixtures, worst position error against a stable
   long-double oracle fell from 1.53e-6 to 9.09e-14. Rounded large multipliers
   can still leave an absolute top-row residual despite small relative error.
5. The repair guarantee concerns the frozen linear model. Independent trace
   auditing found 80,489 accepted repairs that worsened actual energy relative
   to the material-only endpoint, despite passing the predicted bound.
6. At the reset predictor of a rigid pendulum, the implicit objective gradient
   is zero. Consequently the objective-direction gate cannot help that initial
   frozen pass.

## Validation and boundaries

`tests.cpp` checks row equations against independent values, ordering,
spectral limits, rank handling, unit/mass scaling, multiplier blending, cap
bounds, and exactly one main linearization per pass. `audit_results.py` checks
the saved trajectories and reference oscillator independently. Main tests and
the assembly probe passed optimized and Address/UndefinedBehavior Sanitizers
during the research pass. `plots/plot_audit.json` and
`plots/visual_qa.json` record data and visual checks for nine primary figures;
two further figures cover the circular-speed holdout.

The first refreshed-row dataset exposed unresolved-gradient behavior at a
material minimizer reached from rest: the exact energy gradient is zero, but
position roundoff amplified by M/h² can produce a spurious huge correction.
Keep the unprotected data distinct from any rerun with the representability
guard. Refreshed methods are exploratory and are not the basis of the main
recommendation. The guard's final status is recorded in `WRAP_UP.md`.

## Remaining work

- Decide how the energy correction should distribute work among modes while
  preserving angular momentum; an angular-neutral projection alone failed.
- Resolve turning-point branch selection without abandoning linearized,
  curvature-dropped, small-step updates.
- Finish plotting and auditing the refreshed-row follow-up, and compare its
  extra evaluation cost fairly. Its representability guard and rerun are complete.
- Then validate baseline XPBD and the selected candidates on a 2- or 3-mass
  spring pendulum. No chain result exists yet.
- A chain material-block/graph preconditioner is only a proposed next numerical
  experiment. `LARGE_STEPS.md` explains what it would and would not change.

## Reproduction

From the repository, configure the self-contained experiment. The existing root build files are unchanged. Typical runs:

```sh
cmake -S experiments/correction_schedules -B /tmp/xpbd-correction-build -DCMAKE_BUILD_TYPE=Release
cmake --build /tmp/xpbd-correction-build -j4
ctest --test-dir /tmp/xpbd-correction-build --output-on-failure
/tmp/xpbd-correction-build/xpbd-correction-schedules-study /tmp/correction-sweep
/tmp/xpbd-correction-build/xpbd-correction-schedules-study /tmp/correction-selected baseline,EL_raw,LE_raw,joint_raw,LE_repair_cap,joint_repair_cap,spectral_repair_cap,alternating_repair_cap,LE_direction_cap,EL_repair_cap
/tmp/xpbd-correction-build/xpbd-correction-followups parameters /tmp/correction-parameters.csv
/tmp/xpbd-correction-build/xpbd-correction-followups long /tmp/correction-long.csv
/tmp/xpbd-correction-build/xpbd-correction-followups refresh /tmp/correction-refresh.csv
python3 experiments/correction_schedules/plot_results.py /tmp/correction-sweep /tmp/correction-selected --output /tmp/correction-plots
```

Plotting needs NumPy, pandas, and Matplotlib. Saved data, source, analysis, and
figures are kept together so the experiment can be inspected without rerunning.

Primary inspiration sources: [Muon](https://kellerjordan.github.io/posts/muon/),
[PPO](https://arxiv.org/abs/1707.06347), and
[Large Steps in Inverse Rendering](https://rgl.epfl.ch/publications/Nicolet2021Large).
The spectral filter and fractional-repair rule are research analogies, not
implementations of those optimizers or transfers of their guarantees.
