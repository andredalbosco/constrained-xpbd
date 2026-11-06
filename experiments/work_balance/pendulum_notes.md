# Pendulum work experiments: diagnostic references

These are reference experiments, **not the proposed final Gauss-Seidel method**. The exact Cartesian midpoint result isolates kinematic and work identities; the angular-coordinate discrete-gradient result isolates the endpoint-tangency problem. Neither establishes a new method. The existing project was not edited.

Run from this directory:

```sh
g++ -std=c++17 -O2 -Wall -Wextra -pedantic pendulum_study.cpp -o pendulum_study
./pendulum_study > pendulum_run.log
```

The standalone C++17 program uses double precision and checks result-file opening/writing, small-step length and work identities, conservative-reference energy, angular-reference tangency, projected energy loss, free midpoint angular momentum, and RK4 step-halving convergence. All checks pass. Raw results are `pendulum_summary.csv`, `pendulum_trajectories.csv`, and `pendulum_singularities.csv`.

## Setup and measurement

Mass and radius are both 1. Anchor is stationary at the origin. Simulations run for 10 s. Gravity case starts at theta=0.65 rad from downward vertical, omega=0, with g=(0,-9.81). Free-rotation case uses the same angle, omega=2, and g=0. The reported physical energy is E=|v|²/2-g·q. A fixed anchor prevents global linear-momentum conservation from being an applicable invariant. Angular momentum **about the anchor** is conserved only in the zero-gravity case; gravity exerts a torque.

Reference trajectory: angular ODE theta'=omega, omega'=-9.81 sin(theta), RK4 with step at most 1/20000 s. Halving that step changes the 10 s angle by 4.88e-15 rad and angular velocity by 1.64e-13 rad/s. Free rotation uses the exact angular trajectory. Position errors below are the maximum Euclidean error across the run against this reference, not merely endpoint error.

## 1. Exact kinetic work does not make XPBD conservative

For any velocity update, with total impulse J=v_+-v_- and average endpoint velocity vbar=(v_++v_-)/2,

    ΔT = vbar·J.

With J=hg+J_constraint,

    ΔE = vbar·J_constraint + g·(h vbar - Δq).

Ordinary hard-distance XPBD has q_predict=q_-+h v_-+h²g, projects this onto the unit circle, and reconstructs v_+=Δq/h. Its exact impulse work identity has residual below 6.1e-16 in these runs, **despite substantial mechanical damping**. Choosing this W for a ledger does not change that trajectory. Both nonzero constraint impulse work and the gravity/displacement mismatch are visible in the CSV.

For gravity, h=1/60: constraint impulse work sums to -1.67990948 J and the gravity-work mismatch sums to +0.00696797 J. Their sum is the measured -1.67294151 J energy change. For free rotation, h=1/60: constraint work alone accounts for -1.14243311 J.

## 2. Cartesian midpoint: length and energy, but endpoint tangency error

Use midpoint kinematics and a radial impulse along the discrete gradient of squared radius:

    q_+ - q_- = h (v_+ + v_-)/2,
    v_+ - v_- = h g + eta (q_+ + q_-)/2,
    |q_+|² = |q_-|² = 1.

Then radial work is zero because qbar·Δq=0, and constant-gravity work equals the potential change because Δq=h vbar. Consequently physical energy is conserved algebraically.

A closed-form solution is

    s = 2 q_- + h v_- + (h²/2)g,
    q_+ = 2(q_-·s)s/|s|² - q_-,
    v_+ = 2(q_+ - q_-)/h - v_-.

This is the exact reference solve, not a curvature-dropped iteration. It preserves the previous radius, so ordinary floating-point drift can accumulate very slightly. In all experiments the maximum length error is below 2.1e-13.

The endpoint tangency defect obeys the exact identity

    q_+·v_+ + q_-·v_- = (h/2) ΔT.

Thus two tangent endpoint velocities would require ΔT=0. Under gravity the pendulum normally exchanges kinetic and potential energy, so midpoint Cartesian kinematics, fixed length, both endpoint tangencies, and energy conservation cannot all hold simultaneously. In particular, starting from rest, fixed length plus midpoint kinematics and endpoint tangency forces no motion at all.

For free rotation, ΔT=0, and this midpoint method does preserve tangency, energy, and anchor angular momentum. Its rotation angle per step is 2 atan(h omega/2), so phase accuracy still deteriorates at large steps.

## 3. Tangent projection repairs tangency by removing energy

After the midpoint solve, project v_+ onto the tangent line:

    v_projected = v_+ - q_+(q_+·v_+)/|q_+|².

The removed kinetic energy is exactly (q_+·v_+)²/(2|q_+|²). The accumulated energy loss matches the sum of these losses within floating-point error. This velocity projection also breaks the just-completed Cartesian midpoint drift equation. It is not an exact energy-preserving fix.

## 4. Angular discrete-gradient reference: retain energy and tangency, change kinematics

Parameterize q=(sin(theta),-cos(theta)) and v=omega(cos(theta),sin(theta)). Use

    theta_+ - theta_- = h (omega_+ + omega_-)/2,
    omega_+ - omega_- = -h [V(theta_+) - V(theta_-)]/(theta_+ - theta_-),
    V(theta)=9.81(1-cos(theta)).

The potential quotient is evaluated stably as

    9.81 sin((theta_+ + theta_-)/2) sinc((theta_+ - theta_-)/2).

Eliminating omega gives one scalar equation for delta=theta_+-theta_-:

    delta - h omega_- + (h² g/2) sin(theta_-+delta/2) sinc(delta/2)=0.

We solve it with safeguarded Newton in the guaranteed bracket h omega_- ± h²g/2. At the three tested small steps the scalar equation is strictly monotone. This reference preserves circle length and endpoint tangency by parameterization and energy through its discrete chain rule. **It obeys angular midpoint kinematics, not Cartesian midpoint drift.** At h=1/60 the maximum Cartesian drift residual is 7.15820e-6 m per step, while length and tangency are at roundoff. This is exactly the kinematic restriction that the reference avoids.

## Results: gravity, 10 seconds

Numbers below are final energy change in J, maximum radial velocity in m/s, and maximum position error in m.

| h | Method | Energy change | Max radial velocity | Max position error |
|---|---|---:|---:|---:|
| 1/240 | XPBD | -0.771356 | 0.00812294 | 0.164322 |
| 1/240 | Cartesian midpoint | -2.23e-12 | 2.65023e-5 | 0.000305275 |
| 1/240 | Midpoint + tangent projection | -1.67782e-6 | 2.22e-16 | 0.000305039 |
| 1/240 | Angular discrete gradient | near roundoff | 1.11e-16 | 0.000249279 |
| 1/60 | XPBD | -1.672942 | 0.0301577 | 0.406514 |
| 1/60 | Cartesian midpoint | -2.58e-13 | 0.000424280 | 0.00488315 |
| 1/60 | Midpoint + tangent projection | -0.000107214 | 2.22e-16 | 0.00486729 |
| 1/60 | Angular discrete gradient | near roundoff | 1.11e-16 | 0.00398756 |
| 1/15 | XPBD | -1.997729 | 0.0926268 | 0.611051 |
| 1/15 | Cartesian midpoint | -4.00e-14 | 0.00685146 | 0.0776155 |
| 1/15 | Midpoint + tangent projection | -0.00668190 | 1.11e-16 | 0.0763347 |
| 1/15 | Angular discrete gradient | near roundoff | 1.11e-16 | 0.0634651 |

Angular-reference maximum energy drift is below 2e-14 J over these runs. Both conservative references show approximately second-order trajectory convergence. Conservation is not evidence of large-step accuracy: at h=0.5, Cartesian midpoint still has energy drift below 1e-14 J, but radial velocity reaches 0.733689 m/s and maximum position error reaches 1.15608 m.

## Results: free rotation, 10 seconds

At h=1/60, XPBD loses 1.14243311 J of the initial 2 J, with maximum anchor angular-momentum change 0.69044683 and maximum position error almost 2 m. Cartesian midpoint has maximum energy drift 5.74e-13 J, angular change 2.10e-13, radial velocity 1.54e-12 m/s, and position error 0.00185154 m. Angular discrete gradient is exact constant-speed rotation apart from roundoff (maximum position error below 5e-13 m).

## Singularities and large steps

The reflection formula must not silently accept the antipodal limit as a finite-force solution. If s is nonzero but q_-·s=0, the formula returns q_+=-q_-, so qbar=0. A nonzero constraint impulse then cannot equal eta qbar for any finite eta. Near this point eta becomes arbitrarily large.

Concrete probe: start at the upright point q_-=(0,1), gravity (0,-9.81), and tangential velocity (1,0). At

    h_critical = 2/sqrt(9.81) = 0.638550856814101 s,

q_-·s=0 and no finite radial multiplier solves the midpoint equations. At 0.999999 h_critical the multiplier is about -319269; at 1.000001 h_critical it is about +319282. The code rejects the singular point and records the nearby runs.

With the same upright point but zero velocity, s=0 at h_critical. The equations admit a continuum of endpoint positions rather than selecting a unique branch. The exact stationary upright solution exists, but no purely algebraic reflection division can select it at that singularity; additional continuation/branch choice is necessary. The experiment rejects this degenerate case explicitly.

The mild 0.65 rad initial state completes even the large-step sweeps, but its trajectory and endpoint tangency become poor. This is different from the deliberately constructed singular initialization above.

## What the CSV work columns mean

`sum_constraint_impulse_work` is vbar dotted with the total velocity change minus hg. For the XPBD and unprojected Cartesian-midpoint runs this is the actual discrete constraint impulse. For tangent-projected runs it also includes the projection impulse. For angular-coordinate runs it is only a diagnostic residual impulse relative to Cartesian hg; it should not be interpreted as a solved radial reaction. `max_radial_impulse_residual` is measured before optional tangent projection. The kinetic impulse identity is arithmetic, not an independent conservation guarantee.
