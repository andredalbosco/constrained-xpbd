# Independent theory and numerical audit

This note concerns the original endpoint-velocity objective, a fixed hard global-energy target, fresh gradients, and one frozen linearization per small substep. It does not introduce a new work definition, energy reservoir, or a different integrator. Code inspection refers to `correction.hpp`; independent probes and audit outputs are in this directory.

## 1. Predictor system and the material baseline

For a scalar mass m, let y=x_n+h v_n+h² a, C=|y|-L, n=y/|y|, and u=∇U(y)=k C n for the elastic spring. At the predictor x=y, both multipliers reset to zero, so the primal stationarity residual R is zero. The full primal objective gradient is q=u, while the energy gradient is

    g = m v_n/h + u.

For a hard pendulum u=0, so q=0 at every predictor. This makes the common objective-direction gate completely inactive for that case. It is not legitimate to substitute a constraint force or the dual residual for q merely to make that gate nontrivial.

The frozen Schur system is

    [ A B ] [dλ] = -[C]
    [ B D ] [dη]   [H]

where A=1/m+α/h², B=n·g/m, D=g·g/m, and H=E(y,(y-x_n)/h)-E_target. Define

    dλ₀ = -C/A,
    dx₀ = n dλ₀/m,
    H₀ = H + g·dx₀ = H - BC/A,
    s = D - B²/A,
    w = (g - n B/A)/m.

The material-only baseline is x_B=y+dx₀. H₀ is the energy residual predicted at this baseline by the linearization at y; it is not its analytical energy residual. The code now records both quantities.

For a candidate multiplier pair, the accepted result is consistently blended as

    dλ_acc = dλ₀ + t(dλ_candidate-dλ₀),
    dη_acc = t dη_candidate,
    x_new = y + dx₀ + t d,
    d = [n(dλ_candidate-dλ₀)+g dη_candidate]/m.

Blending only the position or only η would define a different algorithm. The inspected code blends both multipliers correctly.

## 2. What each frozen schedule actually proposes

All entries below are total increments; a repeated row replaces its previous total increment.

| Schedule | Energy-induced difference d | Predicted fractional repair -g·d/H₀ |
|---|---|---|
| Jacobi | -(H/D) g/m | H/H₀ |
| LE | -(H₀/D) g/m | 1 |
| EL | -(H/D) w | (H/H₀)(s/D) |
| LEL | -(H₀/D) w | s/D |
| Joint | -(H₀/s) w | 1 |

For ELE, let η₁=-H/D and δλ=-Bη₁/A. Then η₂=-H₀/D+(B²/(AD))η₁ and d=(nδλ+gη₂)/m. Its final energy row is exact, so its repair ratio is also 1. Its material row need not be exact.

EL, LEL and joint end with an exact material row. Their baseline-relative differences satisfy the homogeneous material row; blending preserves that property. LE, Jacobi and ELE may leave a material residual, scaled by the blend t.

The projected direction w with denominator g·w=s is exactly Schur elimination of the coupled system. Calling that construction a new tangent-energy method would obscure an algebraic equivalence. Using denominator D instead produces a damped version of the coupled energy response.

Repeated LE/EL sweeps on a fixed positive-definite 2×2 matrix have the same asymptotic spectral radius ρ², where ρ=B/sqrt(AD). Their different finite-sweep trajectories cannot be explained by different asymptotic linear convergence rates.

## 3. Fractional repair, equivalences, and the PPO analogy

The proposed repair policy asks for

    0 ≤ -H₀ [g·(t d)] ≤ ε H₀²,   0 ≤ t ≤ 1.

For H₀ nonzero and a strictly improving direction H₀(g·d)<0, this is implemented by

    t_repair = min(1, ε |H₀| / |g·d|).

Otherwise the energy-induced difference is skipped. With ε≤1, the linear residual cannot cross zero and its magnitude is not increased. The additional motion cap limits |t d| relative to h|v_n|+h²|f_n|/m. These statements apply to the frozen linear prediction, not to the analytical energy after the displacement.

When the repair bound is active and the direction improves H₀:

* Jacobi becomes the same ε-scaled LE position direction.
* EL and LEL become the same ε-scaled joint position direction.
* LEL and joint are identical under active repair, including the multiplier pair; this is visible to numerical precision in the current experiments.
* EL/Jacobi may fail the improving-direction test when H and H₀ have opposite signs. LE, LEL and joint improve H₀ for an SPD system.

This is a useful clipping analogy to PPO: bound the effect measured in a chosen linear surrogate, preserve the original residual, and avoid overly large repairs. It is not PPO's clipped likelihood-ratio objective, and PPO does not supply an energy or trajectory guarantee here. The cap is also not an exact trust-region minimizer: it simply scales an already chosen direction.

## 4. Direction-gate equivalences

For the elastic one-constraint predictor, q=kC n. The raw energy displacement is g dη/m. Its material-compensated EL difference is w dη. Exactly,

    q·(w dη) = [(α/h²)/A] q·(g dη/m).

For α>0 the factor is positive. Therefore native pre-E and common post-candidate direction gates make the same decision for one fresh EL sweep, apart from roundoff at zero. This equivalence need not hold with multiple constraints or away from the predictor. For a hard pendulum q=0, so the full-objective direction gate is neutral for every proposal in this experiment.

At a timestep that starts at rest, the full energy and objective gradients coincide everywhere: q=g. Increasing analytical energy then necessarily goes uphill in the objective locally. This is a structural conflict, not a sign mistake in the gate.

## 5. Exact release-from-rest diagnostic

Consider a scalar harmonic spring released from q₀>0 with v₀=0, and put a_h=h²k/m. The predictor is q₀, H=0, and

    dx₀ = -a_h q₀/(1+a_h),
    H₀ = -k q₀² a_h/(1+a_h).

The first proposed position steps are exactly:

* Jacobi and EL: ordinary backward Euler, dx=dx₀.
* LE, ELE and joint: dx=0. They cancel the material displacement and freeze immediately.
* LEL: dx=-a_h² q₀/(1+a_h)², of order h⁴ rather than h².

The direction policy rejects the restorative difference in the latter cases because q=g and g·d>0, returning the material baseline.

For LE or joint without the motion cap, fraction ε selects a position between backward Euler and the frozen candidate:

    dx=(1-ε) dx₀.

Its exact first-step energy is

    E_new/E₀ = 1 - (1-ε²) a_h/(1+a_h).

Fraction ε=.5 does not exactly repair energy: the analytical energy gain from the backward-Euler baseline is quadratic in ε because that baseline is the minimum of the endpoint-energy quadratic. With the current cap=.1, the initial motion cap is tighter than ε=.5: t=min(.5,.1(1+a_h)), so the first displacement is close to .9 of the BE displacement.

This small example should accompany any aggregate result. An exact frozen energy row can select a perfectly motionless state.

## 6. Spectral filtering audit

Let R=diag(sqrt A,sqrt D), Sbar=R⁻¹ S R⁻¹=[[1,ρ],[ρ,1]], and rbar=R⁻¹(C,H). The normalizations remove the different units of the two constraint rows. Both normalized RHS components have units of length times sqrt(mass).

Using eigenvalues 1±ρ and f(σ)=min(σ^(-p),G), the proposal is

    z=-R⁻¹ f(Sbar) rbar.

The code's explicit plus/minus eigenvectors implement this correctly. p=0 gives Jacobi when G≥1; p=1 without clipping gives the joint solve. p=.5 is an inverse-square-root preconditioner. If iterated against the original residual, it changes the numerical solver, not the target equations. For 0≤p≤1 and an SPD normalized matrix, the iteration eigenvalues 1-σ f(σ) lie inside (-1,1), subject to positive finite gains. Near-zero eigenvalues still yield slow modes when the gain is capped.

A normalized spectral gain cap does not bound physical displacement. If g→0 while the normalized block tends to the identity, every p still permits a correction proportional to H/|g|. A separate zero-gradient rule or physical displacement control remains necessary.

The spectral code skips null eigenmodes rather than inventing a denominator. With an incompatible RHS this is a least-squares-style partial update; it can retain a material residual. Joint rank fallback instead returns the material baseline. Those policies are intentionally different and should not be counted as identical joint solvers at singular points.

For release from rest and small h, uncapped spectral p=.5 starts with approximately half the BE displacement. p=1 is the exact frozen joint solution and freezes. A gain cap of 4 is inactive in this first-step small-h case, since the normalized eigenvalues are near one.

The Muon connection is a preconditioning/polar-factor inspiration, not literal application of Muon to η. The polar factor of a scalar gradient is just its sign. A nontrivial orthogonalization can be attached to a matrix of normalized constraint directions; the polar factor of the SPD Schur matrix itself is the identity.

## 7. Asymptotic scale and why good energy need not mean good motion

Away from v_n=0 and gradient cancellation, H at the predictor contains an O(h) potential-work term. The material displacement cancels that term in H₀, leaving the old energy defect plus O(h²). Since g=O(1/h), repairing H₀=O(h²) produces d=O(h³) and a velocity change d/h=O(h²). This supports first-order trajectory consistency in a smooth nonsingular regime.

However, the kinetic Hessian is mI/h². The Taylor remainder between H₀ and analytical energy at the material baseline is O(h²), the same order as H₀ itself. Linear-residual contraction does not prove analytical-energy contraction. Turning points violate the uniform-gradient assumption; a raw energy displacement can become order one. The cap has a genuine role there.

For a free spring, a more specific local calculation explains the angular error. Let u=∇U(x), K=∇²U(x), S=|u|²/m, and T=vᵀKv. In the smooth small-h regime, backward Euler loses approximately .5h²(S+T) energy per step. The baseline-linear residual is

    H₀ = e_n - h²(S + .5T) + O(h³).

A fixed fractional repair 0<ε≤1 has the local energy-error recurrence

    e_(n+1) = (1-ε)e_n + h²[(ε-.5)S + .5(ε-1)T] + O(h³).

After the fast energy-error transient, the leading energy displacement is

    d ≈ h³(S+T) v / [2m|v|²],

independent of ε to this leading order. It acts like antidamping along velocity. This compensates the scalar energy loss while injecting energy into all moving modes. In a radial spring, the stiff radial oscillation is strongly damped by the underlying step, while the repair also pumps tangential motion. Angular momentum and phase can drift even with excellent global energy.

This expansion assumes smooth nonzero velocity and a nonactive cap. It is not a uniform theorem through turning points, compression singularities, or very long time intervals.

## 8. Exact angular-momentum attribution

With endpoint velocity v_new=(x_new-x_old)/h, angular momentum about the fixed anchor obeys, for zero gravity,

    L_new-L_old = (m/h) x_old × (dx₀+d).

Hence the direct energy-extra contribution is exactly

    ΔL_extra = (m/h) x_old × d.

For LE and a free central spring,

    ΔL_extra/L_old = dη [1/h² + k C/(m r_pred)].

This is a noncentral numerical impulse unless the motion is radial. Merely conserving the scalar total energy does not prohibit this torque.

With constant acceleration a and the trapezoidal external torque used in the study, the energy-extra contribution to the torque-balance defect is

    (m/h) x_old×d - .5 m h d×a.

The material-only contribution is evaluated at x_B=y+dx₀ along the same corrected trajectory. These two contributions sum to the measured defect; comparing only separate baseline and corrected trajectories would not give such an exact attribution.

At 480 Hz, the free orbit with LE_repair_cap has:

* L: .88 → 1.2317585; total error +.3517585.
* Material-only contribution along the corrected path: -.0241016.
* Energy-extra contribution: +.3758602.
* Per-step decomposition error below 3.8e-14.
* Final angle lead versus reference: 2.0714 radians; ordinary baseline lags .1020 radians.
* Last-second mean [radial kinetic, tangential kinetic, elastic] energy: [.04088,.72310,.05600], versus reference [.20582,.38291,.23127]. The same total energy is being redistributed into rotation.

For coupled motion under gravity, the corresponding final balance contributions are material +.006403 and energy extra -.245478, totaling -.239075. Evidence is in `independent_trace_audit.json` and its generating script.

## 9. High-resolution holdout: convergence versus practical quality

The independent probe adds h=1/1920,1/3840,1/7680,1/15360 s without changing the algorithm. It tests 10-second free-orbit and coupled-motion trajectories. All runs complete; the cap is inactive in these tests.

| Free orbit h | LE repair RMS position | Baseline RMS position | LE repair max angular error | LE repair max energy error |
|---|---:|---:|---:|---:|
| 1/1920 | .34072 | .01752 | .15003 | .001444% |
| 1/3840 | .18132 | .00910 | .08225 | .000362% |
| 1/7680 | .09335 | .00464 | .04304 | .0000907% |
| 1/15360 | .04733 | .00235 | .02201 | .0000227% |

The corrected trajectory approaches the reference at roughly first order, while its energy error decreases at roughly second order. It is therefore inaccurate to call this variant intrinsically inconsistent based on the 240–960 Hz plots. Its practical defect is a much larger first-order torque/phase coefficient: about twenty times the baseline position error and roughly ninety times its angular error at the finest tested h.

Coupled motion similarly converges: corrected RMS .15562,.07426,.03771,.01901 versus baseline .06391,.03427,.01775,.00904. Files: `independent_asymptotic_probe.cpp` and `independent_asymptotic.csv`.

## 10. Formulation limitation that numerical solvers cannot remove

For the endpoint relation, the implicit objective and total energy satisfy

    F(x) = E(x) - v_nᵀ M v_new + constant.

At rest, fixing E makes the objective flat. For a scalar harmonic oscillator, exact hard-energy endpoints have two roots; minimizing this objective favors the endpoint velocity aligned with the previous velocity. Near a turning point that can favor a stalled branch. A numerical filter may deliberately avoid that stationary solution and give much better trajectories, but solving the constrained objective more exactly is not itself a guarantee of better dynamics.

For a hard 2D pendulum, the two independent hard rows can already determine both position coordinates. The objective then controls dual forces rather than an additional free positional direction. Near dependence, large cancelling multipliers are possible. Spectral damping, scalar scheduling and clipping can control this numerical effect; they do not alter the underlying objective-energy geometry.

## 11. Audit and useful further holdouts

The inspected proposal equations, Schur expression, and spectral normalization are correct. Repeated-row replacement is correct. Counters named `accepted` count blend>0 even if the candidate equals the material baseline, including a singular fallback, so they should be described as acceptance decisions rather than nonzero energy corrections.

Useful tests beyond aggregate energy curves:

1. Exact release-from-rest first step and second-step gradient cancellation.
2. Baseline-linear H₀ versus analytical material-baseline energy residual.
3. Raw versus net multiplier displacement and cancellation ratio.
4. Per-step external-torque balance, plus exact material/energy attribution.
5. Radial versus tangential energy partition and unwrapped phase.
6. A true circular spring orbit, with v²=k r(r-L)/m, to separate radial-mode transfer from ordinary orbital drift.
7. Short runs starting near turning points with both signs of a small velocity.
8. Changing units or mass while preserving dimensionless dynamics, to check scaling of numerical tolerances and physical caps.
9. Repeating the same time interval at finer h before claiming inconsistency or convergence.

## 12. Circular and eccentricity holdouts

`independent_direction_holdout.csv` tests six tangential speeds at initial radius 1.1. Speed is expressed as a fraction of the true circular speed sqrt(k r(r-L)/m)=sqrt(11). The original orbit's speed .8 is only .2412 of that circular speed, so it has strong radial oscillation.

At 480 Hz:

| Speed fraction | Baseline RMS position | LE repair RMS position |
|---|---:|---:|
| .1 | .03195 | 1.24573 |
| .2412 | .05790 | .91962 |
| .5 | .22945 | .32231 |
| .75 | .51902 | .06249 |
| 1 (true circle) | .84403 | .000540 |
| 1.25 | 1.13462 | .03944 |

For the true circle, corrected max energy error is .00364% and max angular error .000144, versus baseline 26.61% and .56156. At 7680 Hz corrected position RMS is .000054. Thus a statement that the method generically fails for a 2D spring is false: it is extremely effective near the circular family and poor for radial-dominated eccentric motion.

The leading torque calculation predicts this contrast. For a true circle, S=T and v²=k r(r-L)/m; the energy repair's leading antidamping torque exactly cancels backward Euler's leading torque loss. With a radial mode present, the energy repair transfers dissipation from that mode into rotation.

The same file contains a deliberately separate diagnostic, `angular_neutral_extra`. It projects the energy direction onto displacements parallel to x_old, normalizes that projected direction to repair H₀, and uses the same fraction and cap. Therefore x_old×d=0 exactly, but this changes the original frozen multiplier system and is not one of its alternative solvers.

This diagnostic fails: for the original eccentric orbit its RMS is 1.123 at 480 Hz and .889 at 7680 Hz; radius stays near 1.10 instead of oscillating down to about .918. It keeps low energy error by suppressing radial motion. Projecting away the tangential kinetic gradient leaves a small radial gradient at a radial turning point, allowing the energy correction to cancel radial acceleration. Even the true circle becomes poor. Angular neutrality of the extra impulse alone is not a solution.

The exact nonlinear objective degeneracy, with scalar and hard-pendulum branch formulas, is detailed separately in `EXACT_DEGENERACY.md`.
