# Spring work study

Standalone derivation and experiments; the existing constrained-xpbd repository was not modified.

The useful distinction is between an exact kinetic-work identity, an exact elastic-work identity, and a numerical update satisfying both with the same work. Defining work from a completed velocity update establishes the first identity, but does not enforce the second.

## Reproduction and artifacts

Working directory: /tmp/xpbd-work-research.X6BItq

    c++ -std=c++17 -O2 -Wall -Wextra -pedantic \
      -I/home/blacknoir/Projects/constrained-xpbd/external/eigen \
      spring_study.cpp -o spring_study
    ./spring_study > spring_stdout.txt

The optional first argument changes the duration in seconds; default is 20. The program writes spring_linear.csv, spring_radial.csv, spring_degeneracy.csv, and spring_linearized.csv into its current working directory. The linearized-solver study uses at most 10 seconds. All output streams check failures. It asserts the kinetic and elastic work identities, analytical oscillator phase, predicted dissipation identities, midpoint conservation for the quadratic spring, and energy/angular conservation for the nonlinear central discrete-gradient spring. The checked build and all experiments completed successfully.

## Scalar anchored spring

Use displacement from rest q, mass m=1, stiffness k=100, initial q=0.25 and v=0. Its initial mechanical energy is 3.125; the exact solution is

\[
q(t)=0.25\cos(10t),\qquad v(t)=-2.5\sin(10t).
\]

Let \(\Delta q=q_+-q_-\), \(\Delta v=v_+-v_-\), and use arithmetic endpoint averages \(\bar q,\bar v\).

The three work definitions in the CSV are

\[
W_{\rm old}=-kq_-\Delta q,\qquad
W_{\rm impulse}=m\Delta v\,\bar v,\qquad
W_{\rm DG}=-k\bar q\,\Delta q.
\]

Their algebraic properties are

\[
\Delta T-W_{\rm impulse}=0,\qquad
\Delta U+W_{\rm DG}=0,\qquad
\Delta U+W_{\rm old}=\tfrac12k(\Delta q)^2.
\]

The last expression explains why exact old-force work balance can freeze a linear elastic strain. Changing the work label on a completed trajectory does not change that trajectory.

### Backward Euler / converged scalar XPBD

\[
q_+=q_-+hv_+,\qquad m(v_+-v_-)=-hkq_+,
\]
\[
q_+=\frac{mq_-+hmv_-}{m+kh^2},\qquad
v_+=\frac{q_+-q_-}{h}.
\]

The actual endpoint force has spatial work \(W_+=-kq_+\Delta q\), and

\[
\Delta T-W_+=-\tfrac12m(\Delta v)^2,\qquad
\Delta U+W_+=-\tfrac12k(\Delta q)^2.
\]

Consequently

\[
\Delta E=-\tfrac12m(\Delta v)^2-\tfrac12k(\Delta q)^2.
\]

Replacing the work diagnostic by \(W_{\rm impulse}\) makes the kinetic residual exactly zero but leaves \(\Delta U+W_{\rm impulse}=\Delta E\). This is an accounting identity, not an energy-conservation repair.

The exact numerical energy and phase are

\[
\frac{E_N}{E_0}=(1+\omega^2h^2)^{-N},\qquad
\theta_N=N\arctan(\omega h).
\]

### Exact elastic work with backward kinematics still dissipates

Use the discrete-gradient force but retain backward displacement:

\[
q_+=q_-+hv_+,\qquad m\Delta v=-hk\bar q.
\]

Its update is

\[
q_+=
\frac{(1-kh^2/(2m))q_-+hv_-}{1+kh^2/(2m)}.
\]

Although \(W_{\rm DG}=-\Delta U\) exactly,

\[
\Delta T-W_{\rm DG}=-\tfrac12m(\Delta v)^2,
\qquad
\Delta E=-\tfrac12m(\Delta v)^2.
\]

Thus changing the force quadrature alone leaves the backward-kinematics defect.

### Midpoint kinematics and elastic discrete gradient

\[
\Delta q=h\bar v,\qquad m\Delta v=-hk\bar q.
\]

Now

\[
W_{\rm impulse}=-hk\bar q\bar v
=-k\bar q\Delta q=W_{\rm DG},
\]

so the kinetic and elastic balances hold with the same work. The explicit scalar solve is

\[
q_+=\frac{(1-a)q_-+hv_-}{1+a},\qquad
v_+=\frac{2\Delta q}{h}-v_-,\qquad a=\frac{kh^2}{4m}.
\]

Energy is conserved in exact arithmetic. The phase is still approximate:

\[
\theta_N=2N\arctan(\omega h/2).
\]

This is second-order phase accuracy, not exact dynamics.

### An XPBD-like midpoint solve with 4 alpha / h squared scaling

The unknown midpoint position is \(y=(q_-+q_+)/2\). Set

\[
\widehat y=q_-+\tfrac h2v_-,\quad
\alpha=1/k,\quad
\widetilde\alpha=4\alpha/h^2.
\]

One scalar XPBD projection, initialized with \(\lambda=0\), is

\[
\Delta\lambda=
-\frac{y}{m^{-1}+\widetilde\alpha},\qquad
y\leftarrow y+m^{-1}\Delta\lambda.
\]

Reconstruct

\[
q_+=2y-q_-,\qquad
v_+=v_-+\frac4{mh}\Delta\lambda.
\]

This is exactly the linear implicit-midpoint update. Its mass-scaled midpoint objective is

\[
\frac12m(y-\widehat y)^2+\frac{h^2}{4}U(y).
\]

The effective half-step compliance scaling belongs to the midpoint variable y. It must accompany the corresponding position and velocity reconstruction; changing the compliance scaling alone is not this method.

The impulse-based velocity reconstruction in the experiment accumulates less cancellation error than reconstructing velocity by subtracting close endpoint positions.

## Scalar numerical evidence

All figures below cover 20 seconds, approximately 31.83 exact periods.

| Timestep | BE final relative E error | DG force + backward kinematics | Midpoint-XPBD final relative E error | Midpoint phase error, rad |
|---|---:|---:|---:|---:|
| 1/30 | approximately -1 | approximately -1 | -3.00e-15 | -1.82159 |
| 1/60 | approximately -1 | -0.999999936 | +2.66e-15 | -0.461043 |
| 1/120 | -0.999999939 | -0.999752675 | -6.00e-15 | -0.115620 |
| 1/240 | -0.999757887 | -0.984334926 | +7.55e-15 | -0.0289277 |
| 1/480 | -0.984482125 | -0.874907094 | -7.44e-15 | -0.00723333 |
| 1/960 | -0.875471455 | -0.646340744 | -1.97e-14 | -0.00180842 |

At h=1/240 the BE kinetic impulse-work residual is at most 4.17e-16, although 99.9758% of the mechanical energy is lost. Its maximum per-step elastic impulse-work residual is 0.00541594. Midpoint-XPBD reduces that elastic residual to approximately 1.87e-15. This directly demonstrates why a kinetic balance satisfied by definition is insufficient.

## Turning points and no-motion branch

If one retains \(v_+=\Delta q/h\), starts at \(v_-=0\), and imposes only exact endpoint mechanical energy, the feasible displacements are

\[
\Delta q=0
\quad\hbox{or}\quad
\Delta q=-\frac{2kq_-h^2}{m+kh^2}.
\]

The nonzero branch produces \(v_+\sim-2hkq_-/m\), twice the physical initial acceleration. At h=1/60:

| Method or branch | q next | v next | E error |
|---|---:|---:|---:|
| No-motion energy branch | 0.25 | 0 | 0 |
| Moving energy branch | 0.2364864865 | -0.8108108108 | 0 |
| Midpoint | 0.2465517241 | -0.4137931034 | approximately -4.4e-16 |
| Exact solution | 0.2465358079 | -0.4147403317 | 0 |

At this initial turning point the physical-unit backward-Euler objective

\[
F_{\rm BE}(q)=\frac{m}{2h^2}(q-q_--hv_-)^2+U(q)
\]

equals the reconstructed mechanical energy. Both energy-feasible backward branches therefore have the same objective, 3.125. A hard energy constraint does not select the physical branch.

The CSV's backward_euler_momentum_residual intentionally tests the backward-Euler momentum equation, including for midpoint and exact rows; its nonzero value for those rows is not a failed test of their own dynamics.

No division by \(\Delta q\) is used. A further checked midpoint example has \(q_-=q_+=0.25\), \(v_-=+0.2083333333\), \(v_+=-0.2083333333\): zero displacement and zero net work coexist with a nonzero momentum impulse. The polynomial discrete gradient remains finite. At the true equilibrium \(q=v=0\), midpoint-XPBD remains exactly at rest.

## Translating and rotating two-particle radial spring

For a nonlinear test use

\[
\mathbf r=\mathbf x_1-\mathbf x_2,\quad
U(\mathbf r)=\tfrac12k(\|\mathbf r\|-\ell)^2,\quad
m_1=2,\ m_2=1,\ \ell=1,\ k=100.
\]

Initial relative position is (1.25,0,0), relative velocity is (0,1,0), center of mass is (0.4,0.7,-0.1), and translation velocity is (0.3,-0.2,0.1). The internal initial energy is 3.4583333333, and the constant center-of-mass kinetic energy is 0.21. Reported relative energy errors use internal energy as denominator; absolute total-energy changes are identical.

With \(\rho_\pm=\|\mathbf r_\pm\|\), a central discrete potential gradient is

\[
\boxed{
\overline{\nabla}U
=k\left(1-\frac{2\ell}{\rho_-+\rho_+}\right)
\frac{\mathbf r_-+\mathbf r_+}{2}.
}
\]

Direct multiplication verifies

\[
\overline{\nabla}U\cdot(\mathbf r_+-\mathbf r_-)
=\frac{k}{2}(\rho_-+\rho_+-2\ell)(\rho_+-\rho_-)
=U_+-U_-.
\]

Combine it with

\[
\Delta\mathbf r=h\bar{\mathbf u},\qquad
\mu\Delta\mathbf u=-h\overline{\nabla}U,\qquad
\mu=\frac{m_1m_2}{m_1+m_2}.
\]

The opposite particle forces conserve linear momentum. Their direction is parallel to the midpoint relative position, giving zero midpoint internal torque and conserving angular momentum. The discrete-gradient identity and midpoint kinematics conserve mechanical energy. The nonlinear equations were solved by Newton iteration with the analytical Jacobian and a residual-decreasing line search.

This radial construction is the quadratic-radial-potential specialization of Gonzalez's symmetry-preserving discrete-gradient central-force scheme, equation 3.43. The source also distinguishes symmetry-preserving discrete derivatives from merely imposing a scalar energy identity. [Gonzalez, Time Integration and Discrete Hamiltonian Systems, 1996](https://web.ma.utexas.edu/users/og/PUBLICATIONS/paper_DisHamSys.pdf)

Measured over 20 seconds:

| Timestep | Ordinary midpoint max relative E error | Central DG midpoint max relative E error | Central DG max absolute angular error |
|---|---:|---:|---:|
| 1/30 | 3.47166e-3 | 7.44294e-13 | 1.92071e-14 |
| 1/60 | 8.66619e-4 | 2.77800e-12 | 2.68565e-13 |
| 1/120 | 2.16578e-4 | 4.68562e-11 | 3.89355e-13 |
| 1/240 | 5.41394e-5 | 4.51994e-12 | 8.80740e-13 |

Ordinary midpoint preserves angular momentum here but does not exactly preserve nonlinear spring energy. Backward Euler at h=1/240 loses 87.5248% of internal energy and has maximum absolute angular-momentum error 0.0716594. Linear momentum is preserved by the free center-of-mass update in all three methods.

The small non-monotone DG errors arise from the stopping tolerance and roundoff; mathematical conservation requires an exact nonlinear solve. Endpoint collapse was guarded against, not tested as a physical crossing.

## Scalar block solves with linearization and omitted curvature

The additional spring_linearized.csv experiment implements the discrete-gradient midpoint equations using scalar multiplier updates, without the Newton curvature matrix. This is a single spring block, so it tests the local solve used in a GS scheme, not convergence of multiple coupled GS blocks.

Use midpoint relative position \(y\), predictor \(\widehat y=\mathbf r_-+(h/2)\mathbf u_-\), and reconstruct \(\mathbf r_+=2y-\mathbf r_-\). Define

\[
\bar C=\tfrac12(\rho_++\rho_-)-\ell,\qquad
\bar J=\frac{\mathbf r_++\mathbf r_-}{\rho_++\rho_-},\qquad
J_{\rm end}=\frac{\mathbf r_+}{\rho_+},\qquad
\widetilde\alpha=\frac4{kh^2}.
\]

The equations and residuals are

\[
r_x=\mu(y-\widehat y)-\bar J\lambda,\qquad
r_\lambda=\bar C+\widetilde\alpha\lambda.
\]

The exact derivative of \(\bar C\) with respect to y is \(J_{\rm end}\), not \(\bar J\). Omitting the curvature term \(\lambda D_y\bar J\) gives the update

\[
\boxed{
\Delta\lambda=
\frac{-r_\lambda+J_{\rm end}\cdot r_x/\mu}
{J_{\rm end}\cdot\bar J/\mu+\widetilde\alpha},
\qquad
\Delta y=\frac{\bar J\Delta\lambda-r_x}{\mu}.
}
\]

This needs a scalar denominator and vector operations. Keeping the position residual does not require keeping curvature. A second approximation replaces \(J_{\rm end}\) by \(\bar J\), making the local denominator symmetric. Both variants have the same intended fixed point if both residuals converge.

The XPBD-style residual-dropping variant substitutes zero for \(r_x\) in these updates, while still measuring the actual \(r_x\) afterward.

Reconstruct endpoint velocity by \(\mathbf u_+=4(y-\mathbf r_-)/h-\mathbf u_-\), and define the material impulse and its work by

\[
I_\lambda=\frac4h\bar J\lambda,\qquad
W_\lambda=I_\lambda\cdot\bar{\mathbf u}.
\]

The exact accepted-state defect identities are

\[
\boxed{\mu\Delta\mathbf u-I_\lambda=\frac4h r_x,}
\]
\[
\boxed{
\Delta T-W_\lambda=\frac4h\bar{\mathbf u}\cdot r_x,\qquad
\Delta U+W_\lambda=k\,r_\lambda\bar J\cdot\Delta\mathbf r.
}
\]

Thus the material row closes elastic work, while the position/momentum row closes kinetic work. Converging only the material multiplier equation cannot establish both balances.

At h=1/60 over 10 seconds, using the same translating/rotating radial spring:

| Local approximation | Sweeps | Max relative E error | Max material residual | Max force residual \(4\|r_x\|/h^2\) |
|---|---:|---:|---:|---:|
| Retain position residual, omit curvature | 1 | 6.42737e-7 | 2.21162e-9 | 2.10817e-5 |
| Retain position residual, omit curvature | 5 | 1.43663e-13 | 1.87350e-16 | 1.07833e-12 |
| Retain position residual, omit curvature | 20 / 100 | 1.01030e-13 | 1.66533e-16 | 1.11509e-12 |
| Retain position residual, use symmetric Jbar | 1 | 1.46042e-6 | 8.32634e-9 | 2.10817e-5 |
| Retain position residual, use symmetric Jbar | 5 | 1.01030e-13 | 1.66533e-16 | 1.11509e-12 |
| Drop position residual, omit curvature | 5 / 20 / 100 | 6.54676e-7 | 2.77556e-17 | 2.10817e-5 |
| Drop position residual, use symmetric Jbar | 5 / 20 / 100 | 6.54676e-7 | 2.77556e-17 | 2.10817e-5 |

With position-residual dropping, the constitutive residual reaches roundoff but the kinetic material-work error remains approximately 3.20509e-7 per step. Further sweeps do not remove it. With the residual retained, five sweeps reach the converged discrete-gradient trajectory within approximately 9.3e-13 in final relative position.

All these single-spring variants preserve angular momentum to approximately 1.5e-13 here: every local correction stays along the midpoint relative direction, so the remaining \(r_x\) is also central. This special geometry must not be generalized to coupled constraints. Linear momentum is preserved by the free center-of-mass update.

The experiment uses no line search or Hessian in the linearized variants. It establishes convergence for this case and timestep; it does not establish convergence for arbitrary coupled systems, large timesteps, or endpoint collapse.

## Scope and limitations

- The scalar midpoint XPBD projection is an exact linear-spring equivalence. For nonlinear distance springs, simply substituting midpoint positions in ordinary XPBD is not an exact energy method.
- The nonlinear radial discrete-gradient force generally has a nonsymmetric derivative with respect to the new endpoint. It is not automatically the gradient of a scalar endpoint objective. The derived block residuals and scalar update give a concrete approximate force solve; a general variational XPBD dualization has not been completed.
- The anchored scalar spring exchanges momentum with its anchor. Momentum-conservation evidence comes from the separate free two-particle experiment.
- These experiments use no reservoirs and no post-step energy rescaling. They demonstrate compatible work and kinematics, rather than compensating an incompatible update afterward.
- Energy preservation does not imply exact phase or exact trajectories. The scalar phase errors are measured against the analytical solution; no analytical phase claim is made for the rotating nonlinear spring.
- This study does not cover hard-constraint pendulum tangent velocities, collisions, contact, arbitrary elastic potentials, or finite-sweep convergence of multiple coupled XPBD blocks.
