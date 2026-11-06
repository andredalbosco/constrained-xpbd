# Incremental lambda with a matching predictor

The requested expectation that vanilla XPBD quickly diverges was not confirmed.
For affine constraints, consistent primal/dual warm-starting gives exactly the
same implicit step. For nonlinear constraints, dropping the position residual
changes the method and can reduce damping; all tested vanilla cases completed.
Adding an unfiltered endpoint-energy row introduces large corrections near the
energy minimum and did fail in several cases.

## Equations and conventions

Use fixed h, scaled XPBD lambda, alpha_tilde=alpha/h², and
lambda_minus=-h² alpha^-1 C(x_minus) when a constitutive old force is available.
For a changed step size, scale a stored lambda by (h_new/h_old)². Hard-rod
initial reactions are not known here: they start at zero, then carry solved
lambda. The energy multiplier starts at zero each substep.

    y = x_minus + h v_minus + h² a_gravity
    x0 = y + M^-1 J_minus^T lambda_minus
    lambda0 = lambda_minus

Here y remains the original inertial predictor; x0 is the solver's initial
position incorporating the previous constraint force. At each linearization,

    r_C = C(x) + alpha_tilde lambda
    r_x = M(x-y) - J^T lambda - g eta
    M delta_x - J^T delta_lambda - g delta_eta = -r_x
    J delta_x + alpha_tilde delta_lambda = -r_C.

This drops curvature; dropping r_x as well is an additional approximation.
Without energy and with dropped r_x, scalar XPBD is

    delta_lambda = -(C + alpha_tilde lambda)/(J M^-1 J^T + alpha_tilde)
    delta_x = M^-1 J^T delta_lambda
    lambda += delta_lambda; x += delta_x.

Retain the TOTAL lambda in the compliant residual, not only its increment.
For constant J, M(x0-y)=J^T lambda_minus, so the warm start satisfies the top
row and produces the same affine solve as initialization at y with lambda=0.
For nonlinear J, the initial residual is

    r_x0 = (J_minus^T - J(x0)^T) lambda_minus.

Discarding it omits the force contribution from the change of constraint
direction. It is therefore not merely a different initialization of the same
exact nonlinear solve. Numerical improvements in this variant are not a proof
that it is a symplectic integrator.

## Scalar spring result

Let C(q)=q, alpha=1/k, a=h²k/m, and lambda_minus=-h²k q_minus.
The matched initialization and increment are

    q0 = (1-a)q_minus + h v_minus
    delta_lambda = -(q0 + alpha_tilde lambda_minus)/(1/m+alpha_tilde).

Substitution gives exactly

    q_plus = (q_minus+h v_minus)/(1+a)
    v_plus = (v_minus-h(k/m)q_minus)/(1+a).

The state matrix determinant is 1/(1+a); the two oscillator eigenvalues have
modulus 1/sqrt(1+a). This is the dissipative backward-Euler step, for every h.
Initializing only lambda, without the matching position shift, instead gives
q_plus=q_minus+h v_minus/(1+a), so a released spring freezes.

Taking the old-force predictor WITHOUT the incremental XPBD correction gives

    v_plus=v_minus-h(k/m)q_minus
    q_plus=q_minus+h v_plus.

This is kick-drift symplectic Euler. It has bounded oscillator motion for
0<h*sqrt(k/m)<2; it is not unconditionally stable and does not conserve the
original energy exactly. Its invariant is

    I = (m v²+k q²)/2 - (h k/2) q v.

At m=1,k=100,h=1/480 over ten seconds, consistent warm XPBD differs from cold
XPBD by at most 3.69e-14 in sqrt(delta_q²+(delta_v/10)²). Both finish at
E/E_initial=0.1245707624. Predictor-only energy stays between 0.98969072 and
1.01052632 of its initial value. A deliberately large h=0.21 (h omega=2.1)
does cause predictor-only growth, crossing the explicit 1e20 energy-ratio
limit after 35 accepted steps. This is a stability-bound diagnostic, not a
proposed small-step setting.

## Adding the energy equation

Use the same physical endpoint energy as before, including gravity:

    E(x)=(x-x_minus)^T M(x-x_minus)/(2h²)+U(x)-(M a_gravity)^T x
    H=E(x)-E_target, g=grad E(x).

For one material row, with r_x dropped, the frozen multiplier system is

    [ A B ] [delta_lambda] = -[r_C]
    [ B D ] [delta_eta   ]    [ H ]

where A=J M^-1 J^T+alpha_tilde, B=J M^-1 g, D=g^T M^-1 g.
The energy row is g^T delta_x=-H. EL, LE, and joint schedules are unchanged;
the new initialization changes r_C, H, g, and the finite-iteration path.

For a scalar spring released from rest, the old-force predictor gives

    q0=(1-a)q_minus
    H0=-(k q_minus²/2) a(1-a)
    g0=-k a q_minus.

Consequently the raw linear energy correction is

    delta_q=-H0/g0=-(1-a)q_minus/2.

It remains O(q_minus) even as h approaches zero. The warm predictor is only
O(h^4) from the endpoint-energy minimum q_BE=q_minus/(1+a); a small energy
residual therefore comes with a very small energy gradient. Dropping energy
curvature is especially poor at this point. In one scalar position degree of
freedom, LE and joint both close this same energy row and have the same total
position correction. EL's later material correction slightly changes it.

With q_minus=0.25, m=1,k=100,h=1/480, q0=0.2498914931 and the linear energy
correction is -0.1249457465. The resulting energy is about 576.75 times the
initial energy. Consistent warm-starting does not remove the exact from-rest
objective-energy degeneracy; it changes the local linearization and can trade
freezing for an enormous raw energy correction.

## Physical tests and limits

`warm_start_probe.cpp` ran 405 ten-second cases: nine existing spring/pendulum
conditions, 240/480/960 Hz, cold and warm baseline controls, retained-r_x control,
and cold/warm EL/LE/joint with raw or fractional-repair+cap acceptance. References
were checked by RK4 subdivision refinement. Twelve raw warm energy runs stopped
explicitly; all baseline variants completed. Completion is not accuracy.

At 480 Hz:

| Case | Cold baseline RMS | Warm baseline, dropped r_x RMS | Max energy error, cold -> warm |
|---|---:|---:|---:|
| Spring release | 0.074157 | 0.074157 | 87.543% -> 87.543% |
| Small spring swing | 0.006256 | 0.000703 | 17.443% -> 0.638% |
| Rigid pendulum from rest | 0.037654 | 0.001447 | 21.995% -> 0.318% |
| Spring swing with tangential velocity | 0.075066 | 0.015093 | 23.645% -> 4.484% |
| Eccentric spring orbit | 0.057904 | 0.041751 | 47.464% -> 46.534% |
| Coupled spring-pendulum motion | 0.173609 | 0.212392 | 60.902% -> 54.707% |

Retaining r_x gives results close to the original cold baseline. Raw energy
warm-starts produce the predicted release spike and other failures. Repair/cap
controls suppress these large steps, but do not generally fix orbital phase:
warm LE repair+cap eccentric-orbit RMS is 0.923417 versus 0.919622 cold.
On the small swing, warm LE repair+cap RMS is 0.001327, worse than the warm
baseline's 0.000703 despite reducing energy error to 0.0058%.

Files: `warm_start_results/scalar.csv`, `first_energy_step.csv`, `physical.csv`.
Run the `xpbd-warm-start-probe` CMake target with an output directory. Assertions
check affine equivalence, the symplectic invariant, the raw first energy step,
reference convergence, and completion of every vanilla control. This is a
single-bob study; no chain, cloth, or grid was added.
