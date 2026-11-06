# Direct energy multiplier coupled to a radial spring

This standalone experiment keeps the existing repository untouched. It tests the user's joint-coupling instability claim, including ordinary XPBD as the baseline. It finds a concrete moving-case runaway for a truncated joint update, but also distinct frozen and converged-row cases. Neither matrix conditioning nor small energy error alone establishes physical correctness.

## Reproduce

    cd /tmp/energy-coupling-spring
    c++ -std=c++17 -O2 -Wall -Wextra -pedantic coupling_study.cpp -o coupling_study
    ./coupling_study

The optional first argument is duration in seconds (default 10). Files:

- coupling_results.csv: 864 cases, 31 columns.
- coupling_traces.csv: 2,835 one-step trace rows, 20 columns.
- coupling_selected_trajectories.csv: 7,200 time-series rows, 13 columns.
- run.log: mandatory analytical check results.

Parameters: mass=1, rest length=1, stiffness=100 or 1000, h=1/240, 1/60, or 1/15, iteration limits=1, 5, 20, or 80. Initial x=(1.25,0); initial v is (0,0), (0,1), or (-0.7,1).

Both energy targets are explicitly labeled:

- initial: the fixed initial mechanical energy.
- previous_actual: mechanical energy of the previous accepted state, including its accumulated numerical error.

No gravity, damping, reservoirs, or post-step energy rescaling is used.

## Equations actually implemented

Let x be the endpoint candidate, \(v_+=(x-x_-)/h\), \(\widehat x=x_-+hv_-\), \(C=\|x\|-1\), \(J=x/\|x\|\), and \(U=kC^2/2\). The energy residual and exact gradient are

\[
H=\frac{\|x-x_-\|^2}{2h^2}+U(x)-E_{\rm target},
\qquad
g_H=\frac{x-x_-}{h^2}+kCJ.
\]

The material residual is \(r_C=C+\widetilde\alpha\lambda\), with \(\widetilde\alpha=1/(kh^2)\). The gradient always uses the actual potential derivative \(kCJ\); no off-manifold substitution involving lambda is made.

The joint curvature-dropped, position-residual-dropped update solves

\[
\begin{bmatrix}a&b\\b&d\end{bmatrix}
\begin{bmatrix}\Delta\lambda\\\Delta\eta\end{bmatrix}
=-\begin{bmatrix}r_C\\H\end{bmatrix},
\qquad
\Delta x=J\Delta\lambda+g_H\Delta\eta,
\]

where \(a=1+\widetilde\alpha\), \(b=J\cdot g_H\), and \(d=\|g_H\|^2\). Both multipliers accumulate within a timestep and reset between timesteps. All rows are re-evaluated after updates.

The sequential method applies the ordinary material projection, re-evaluates H and its true gradient, then applies \(\Delta\eta=-H/\|g_H\|^2\) and \(\Delta x=g_H\Delta\eta\).

The joint block uses the positive Schur expression

\[
s=d-\frac{b^2}{a}
=\|g_H-(J\cdot g_H)J\|^2
+\frac{\widetilde\alpha}{1+\widetilde\alpha}(J\cdot g_H)^2.
\]

No subtractive determinant is used for the solve. The raw determinant is only counted diagnostically. All 864 cases had zero nonpositive raw determinants. The matrix is positive definite for nonzero \(g_H\) and positive compliance.

Condition numbers are computed after scaling each diagonal to one; they describe relative row alignment, not the absolute size of an energy gradient. A nearly zero energy gradient can be dangerous even when this normalized condition number is modest.

The two backtracking variants use the same proposed steps. Joint backtracking checks decrease of \(r_C^2/a+H^2/d\), with scales frozen during each line search. Sequential backtracking checks decrease of the energy residual magnitude after the material update. Up to 40 halvings are allowed. Ordinary guarded methods reject an energy update when \(\|g_H\|\le10^{-12}(1+\|\Delta x\|/h^2+|kC|)\). The sequential_raw control omits this near-zero-gradient guard, retaining only an exact-zero check.

## Mandatory checks

The executable first asserts:

1. Ordinary scalar XPBD matches the analytical backward-Euler spring recurrence.
2. Its 10-second scalar trajectory remains finite, bounded, and non-increasing in energy.
3. The joint turning-point update has the predicted frozen position, lambda, and eta.
4. The sequential material update from rest produces a nearly zero energy gradient with negative energy residual.

These checks pass. The experimental sweep deliberately records failing methods rather than asserting that every case succeeds. A process exit code of zero therefore means the mandatory checks passed and the outcome tables were written, not that every experimental method worked.

## Frozen joint solution from rest

At the initial turning point \(v_-=0\), the predictor is the old position and \(H=0\). Here \(g_H=kCJ\), so the joint energy row forces \(\Delta x=0\). The compliant material row can nevertheless be satisfied by

\[
\lambda=-h^2kC,\qquad\eta=h^2.
\]

All 48 turning-point cases for each joint variant remained exactly frozen over 10 seconds. At k=100, h=1/60, normalized condition=1.3934868, energy error=0, displacement=0, speed=0, and stationarity residual=0.

This is a wrong stationary branch for a stretched spring, not a numerical blow-up. It is labeled frozen_turning_state rather than a successful physical simulation.

Sequential projection from the same state first reaches the backward-Euler energy minimum: \(g_H\) is approximately zero while \(H=-0.08445945946\). Guarded methods report singular_energy_updates. The unguarded control can take enormous steps from floating-point residue: at k=100, h=1/60, one iteration, fixed initial-energy target, the accepted maximum radius is 1.51e13 and maximum relative energy error is 4.61e28.

Thus sequential projection is not automatically safe either.

## Concrete moving-case joint runaway

Use x=(1.25,0), v=(0,1), k=1000, h=1/60, one iteration per step, and the previous_actual target. Initial energy is 31.75.

| Time | Joint E/E0 | Sequential E/E0 | Ordinary XPBD E/E0 |
|---|---:|---:|---:|
| 1 s | 10.0799 | 1.83924 | 0.0248339 |
| 2 s | 935.594 | 1.90665 | 0.0235959 |
| 5 s | 1.18047e9 | 2.15264 | 0.0205357 |
| 10 s | 1.62375e19 | 2.81043 | 0.0168833 |

The joint accepted trajectory reaches radius 8.09e8 and speed 1.94e10. This is actual accepted-state growth, not an inference from conditioning: the maximum normalized condition is only 2.65945.

Joint backtracking accepts the same steps in this case. Its combined row merit decreases even while each truncated nonlinear energy correction leaves a positive residual that is incorporated into the next target.

Changing to the fixed initial-energy target removes this explosive accumulation in this case, but one joint iteration still has maximum relative energy error 0.107893. Sequential one-iteration error is 0.765340; sequential backtracking reduces it to 0.0839105.

Five joint iterations close both rows in this moving case: maximum relative energy error is 5.15e-15 for the fixed initial target and 4.44e-14 for the previous_actual target. This does not establish full stationarity.

At the coarser h=1/15, k=1000, one iteration and fixed initial target, the joint maximum relative energy error is 3.07629; undamped sequential error is 27.4566. The evidence therefore supports a vulnerability of the truncated joint method, not a universal statement that joint is always worse than sequential.

## Row convergence is not full KKT convergence

Every accepted step records

\[
r_x=x-\widehat x-J\lambda-g_H\eta.
\]

Because the updates drop this residual, closing the material and energy rows does not force \(r_x=0\).

For the rotating k=1000, h=1/60 case with five joint iterations, maximum stationarity residual is 0.0260686 even though energy and material rows close. At k=100, h=1/60, 80 iterations, joint and sequential both close the rows, yet their maximum stationarity residuals are respectively 9.01e-4 and 1.15e-4. Their trajectories need not coincide.

Statuses distinguish completed_rows_converged, completed_unresolved, frozen_turning_state, completed_with_singular_energy_updates, and runaway_accepted. The latter means an accepted radius exceeds 10 rest lengths or speed exceeds 100. This reporting threshold does not stop a run. A separate magnitude_guard stops only above 1e100 to avoid overflowing subsequent energy calculations; two cases reached that guard. Completion is never described as full KKT convergence.

For ordinary XPBD, completed_rows_converged refers only to its material row; energy is diagnostic and is not an enforced baseline constraint.

## Validation and limits

- Optimized -O2 builds with warnings enabled and full AddressSanitizer/UndefinedBehaviorSanitizer runs pass without sanitizer findings.
- Leak detection was explicitly disabled because this execution environment uses ptrace.
- All output schemas were independently parsed and all numeric CSV fields are finite. No floating-point nonfinite failure was observed; large finite failures and the magnitude guards remain recorded.
- No solver stability claim is inferred from the positive determinant. The matrix can be well-conditioned while the nonlinear or time-discrete iteration follows a wrong branch or amplifies errors.
- The experiment is a single anchored radial spring. It does not test momentum conservation of a free system, coupled multi-spring constraints, or a solver retaining full stationarity residuals.
- These results concern direct endpoint energy constraints with backward velocity reconstruction. They do not contradict the separately tested midpoint work-balance method.
