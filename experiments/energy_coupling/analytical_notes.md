# Direct energy conservation: material/energy coupling audit

Analytical companion to the spring and pendulum experiments. This concerns the direct endpoint energy constraint, not a work estimate. The kinematics remain backward reconstruction, `v+=(x-xold)/h`. A large Newton correction, a singular linearization, nonlinear divergence, a nonunique feasible branch, and an infeasible nonlinear problem are different outcomes and should be reported separately.

## 1. Objective-only dualization and the frozen linear block

Let

    y=xold+h vold+h^2 a,
    F=||x-y||_M^2/(2h^2)+U(x),
    H=||x-xold||_M^2/(2h^2)+U(x)+Vg(x)-Eold,
    Vg(x)=-a^T Mx.

For `U=C^2/(2alpha)`, dualize only the objective U and leave analytical U(x) inside H. Using the scaled energy multiplier eta, the mixed residuals are

    rx=M(x-y)-J^T lambda-eta g,
    rlambda=C+(alpha/h^2)lambda,
    renergy=H,
    J=D C,  g=grad H.

The exact position Jacobian contains `-lambda Hess C-eta Hess H`. Dropping this curvature gives, with `at=alpha/h^2`,

    S = [[J M^-1 J^T+at, J M^-1 g],
         [g^T M^-1 J^T, g^T M^-1 g]],
    S [dlambda,deta]^T = [-rlambda,-H]^T + [J;g^T]M^-1 rx,
    dx=M^-1(J^T dlambda+g deta-rx).

Omitting rx is a second approximation. Even if the two constraint residuals become small, the actual stationary equation may remain unsatisfied.

Let `u=M^-1/2 J^T`, `z=M^-1/2 g`. Then

    det(S)=at ||z||^2+||u||^2||z||^2-(u dot z)^2.

Consequences:

- Positive compliance gives a strictly positive determinant whenever g is nonzero. Parallel material and energy gradients alone do not make the compliant 2x2 block singular.
- For a hard constraint, at=0, dependent gradients make the block singular.
- g=0 makes the energy row singular for either material type. If H also vanishes the row can be redundant; if H is nonzero, no first-order correction can satisfy it.
- Positive definiteness does not prevent bad conditioning or nonlinear failure.

For diagonal entries A,D and off-diagonal B, let `rho=B/sqrt(AD)`. The diagonally normalized matrix has condition number `(1+|rho|)/(1-|rho|)`. Scalar GS on this *same frozen system* has contraction factor rho^2. Joint solution removes that iteration error; it does not create a generic linear instability. A truncated GS solve can appear gentler solely because it has not solved a nearly singular system.

In 2D for one mass m and a radial spring, J=n is a unit vector and g=g_parallel n+g_perp. The determinant is

    ||g_perp||^2/m^2 + alpha ||g||^2/(m h^2).

Along a collinear spring state, `rho^2=kh^2/(m+kh^2)`, k=1/alpha. Large stiffness makes the normalized block ill-conditioned, but finite k and nonzero g still give an SPD block.

Root's frozen-block control was independently reviewed: its Schur complement, normalized condition, GS contraction, and residual checks are correct. With a hard row and gradient angle 1e-4, the exact linearized solution needs a displacement about 3000, while ten scalar GS sweeps move about 0.2 but leave residual about 0.3. This is underconvergence of GS, not proof that the joint solve is unstable.

## 2. A general turning-point degeneracy of the optimization itself

Direct expansion gives

    F(x)=H(x)-vold^T M(x-xold)/h+constant.

On H=0, this objective only selects displacement along the old momentum. At vold=0, F is constant on the entire energy-feasible set, regardless of U. There is no objective preference between a no-motion branch and other feasible endpoints. The physical energy multiplier nu=eta/h^2 can equal one at every regular feasible point.

This explains an actual pathology without attributing it to coupling: at rest the requested optimization is degenerate. An algorithm may choose a frozen branch, another branch, or become ill-conditioned. Energy preservation alone does not establish correct acceleration.

## 3. Compliant 2D spring released from rest

Take an anchored radial spring with no gravity, xold=(l+c0)n and vold=0. The predictor equals xold. On its radial branch, write c for the current strain and assume no crossing through the anchor. The endpoint energy residual is

    H(c)=m(c-c0)^2/(2h^2)+k(c^2-c0^2)/2.

Its two radial roots are

    c=c0,
    c=c0(m-kh^2)/(m+kh^2).

The second gives `v+=-2hk c0/(m+kh^2) n`, whose leading velocity is twice the physical initial acceleration times h. The first freezes the spring. Both have identical F because the initial velocity is zero.

At the initial state, H=0 and g=kc0 n. The compliant joint block is nonsingular for c0!=0, yet its exact first update is

    dx=0,
    dlambda=-h^2 k c0,
    deta=h^2.

It satisfies all three mixed residuals and selects the frozen branch in one update. This is not divergence or a singular frozen Schur block.

By contrast, material-first projection solves the scalar spring exactly at

    c_BE=m c0/(m+kh^2).

This point is also the minimum of the reconstructed endpoint H: g=0 and

    H(c_BE)=-k^2 h^2 c0^2/[2(m+kh^2)]<0.

The subsequent energy projection therefore has a zero denominator and a nonzero residual. A division-free implementation must report this local singularity rather than claim convergence.

For m=1, k=100, h=1/60, c0=.25, the joint block has rho^2=1/37, far from a strong-coupling limit. It nevertheless freezes. Material-first reaches c_BE=.243243243..., with H=-.084459459... and zero energy gradient. This analytically checked case refutes the blanket claim that a joint block is intrinsically unstable while scalar energy projection is intrinsically stable.

The full compliant 2D energy problem is not globally infeasible for consistent initial energy: H(xold)=-Told<=0 and the inertial energy makes H grow at infinity. At rest there is generally an energy contour of tied endpoints. Extra momentum constraints, fixed geometry, or other restrictions can change feasibility.

## 4. Hard 2D pendulum: energy-circle geometry

Let |xold|=l, U=0, mass m, and constant gravity a. The endpoint energy equality is

    ||x-xold-h^2 a||^2=h^2||vold||^2+h^4||a||^2.

It is a circle centered at e=xold+h^2 a with radius `R=h sqrt(||vold||^2+h^2||a||^2)`. The hard constraint is the circle centered at zero with radius l. With d=||e||, feasibility requires

    |l-R|<=d<=l+R.

For this consistent initial setup, the triangle inequality reduces the nontrivial restriction to `R<=l+d`. With no gravity it becomes `h||vold||<=2l`: large-step direct energy plus backward velocity can be genuinely infeasible.

At a turning point vold=0, xold is always feasible, but the predictor y=e is the energy-circle center:

    g(y)=0,
    H(y)=-m h^2||a||^2/2.

The joint material/energy block is singular immediately. Material-first projection sends y to `x=l e/||e||`; there g is parallel to the hard normal, so the joint block is again singular. Except at special tangent cases, the two nonlinear circles still intersect. The failure is of the local first-order model, not global feasibility.

Unperturbed sequential normal projections from this centerline remain on the same line and cannot reach generic off-line intersections. Thus this particular initialization can also stall or cycle under scalar nonlinear GS; its apparent behavior must not be generalized from the compliant block's SPD argument.

For initial angle theta from downward vertical and speed zero, the second branch satisfies

    tan(delta_theta/2)=-g h^2 sin(theta)/(l+g h^2 cos(theta)).

Again its leading displacement has incorrect turning-point acceleration under backward velocity reconstruction. At the upright rest state with `h^2 g=l`, e=0 and R=l: the two circles coincide, giving a continuum rather than a unique endpoint. At nonzero speed, tangency of the circles gives dependent gradients; beyond the intersection condition there is no solution.

Energy alone is also not angular momentum conservation. In the gravity-free case with initially tangent speed v and a nonzero feasible chord, backward reconstruction gives

    |Lplus|=m l v sqrt(1-(hv/(2l))^2),

less than the initial m l v. Simultaneous exact energy and anchor angular momentum with these kinematics is therefore incompatible for a nonzero step in that setup.

### Retaining rx cannot change a regular joint hard-pendulum position step

For a single hard pendulum in 2D, stack the two linearized rows as `K=[J;g^T]` and `r=[C,H]^T`. If K has full rank, they uniquely determine

    K dx=-r,  dx=-K^-1 r.

The top equation with any chosen position Hessian B is

    B dx-K^T [dlambda,deta]^T=-rx.

Retaining rx or changing B changes only the multiplier increments, because the two hard rows already fix both position coordinates. Consequently the joint position trajectory in these tests is unchanged by retaining rx when initialization and the constraint-based damping rule are held fixed. This includes large Newton position steps caused by nearly dependent rows. An acceptance rule that also uses stationarity or the objective would be a different algorithm and can change the accepted trajectory.

At rank deficiency, no choice of B or retained rx can satisfy an inconsistent linearized constraint right-hand side; g=0 with H!=0 is the simplest example. If the right-hand side is compatible, however, the top equation may select among nullspace directions. Rank deficiency alone must not be equated with an impossible linearized correction. This argument is specific to the two independent hard rows in two position dimensions: compliant material rows and higher-dimensional systems retain extra freedom, so their position steps can depend on rx and B.

## 5. Two meanings of “coupled energy multiplier and lambda”

The primary experiment dualizes objective U only. The energy row retains analytical U(x), and its actual gradient is

    g=M(x-xold)/h^2+J^T C/alpha+grad Vg.

Replacing `C/alpha` with `-lambda/h^2` away from the constitutive manifold changes that gradient:

    g_substituted=g_actual-J^T rlambda/alpha.

At small compliance even a small constitutive residual can create a large gradient error. A hard-constraint multiplier must never be substituted for grad U when U=0.

Alternatively, starting with the physical Lagrangian `F-nu H` and dualizing the *entire* coefficient `(1-nu)U` gives the material term

    -kappa C/h^2-alpha kappa^2/[2h^4(1-nu)],
    kappa=(1-nu)lambda.

Its elastic maximum representation is valid for nu<1. The pole at nu=1 is a singular parameterization boundary. The turning-point degeneracy above naturally has nu=1 with finite material lambda and kappa=0, so naive kappa/(1-nu) reconstruction is 0/0 there. This is distinct from the objective-only 2x2 coupling and must not be mixed into its interpretation.

## Reporting conclusions

Useful classifications are: regular converged state, regular but ill-conditioned linear block, local gradient singularity, nonlinear residual divergence/rejection, frozen feasible branch, and globally infeasible constraints. Compare joint and sequential methods on the same frozen system before interpreting nonlinear trajectories. Record actual H, constitutive residual, rx, step norm, multiplier size, accepted damping, gradient norms, and normalized coupling. Constraint satisfaction and small iterate changes alone do not establish accurate or stable physical motion.

## Independent implementation audit

Read-only review of the spring and pendulum study sources found the analytical material/energy derivatives and both Schur-complement formulas consistent with the equations above. Both nonlinear experiments omit rx as well as curvature; their multiplier and constraint convergence therefore need separate interpretation. The pendulum `Step.status` initializes to `ok`, and its output distinguishes constraint closure from rx closure. The spring reports `completed_rows_converged` and records its actual stationarity residual.

The spring runs both initial-energy and previous-actual-energy targets; the pendulum uses previous actual energy. A material-seeded joint pendulum variant has an extra initial material projection. Joint spring backtracking uses a two-row merit, whereas sequential spring backtracking uses energy residual after the material update. These are explicit algorithm variants, not identical frozen-system comparisons.

The independently scanned final outputs contain 864 spring rows and 320 pendulum rows. All spring numerical fields are finite. Four pendulum cases intentionally report nonfinite attempted states; their failed-trial energy fields are infinite, as are some singular condition diagnostics. Those failed trials must not be described as accepted-state blowup or finite outputs.

A useful warning from the measured spring output is the turning/sequential_raw run with k=100, h=1/240, 80 sweeps, and initial energy target. Its maximum energy residual is approximately 3.43e-12, while its maximum stationarity residual is approximately 3.53e21. It closes the two projected rows but is very far from satisfying the mixed stationary equation. Its `completed_rows_converged` label must not be shortened to optimization convergence.
