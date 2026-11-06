# Exact objective-energy degeneracy

This statement concerns the nonlinear formulation itself, before linearization, dropping curvature, or choosing a multiplier schedule.

Let x₋,v₋ denote the previous state, y=x₋+hv₋+h²a the predictor, and v₊=(x-x₋)/h. For constant external acceleration a,

    F(x) = m/(2h²) |x-y|² + U(x),
    E(x) = m/(2h²) |x-x₋|² + U(x) - ma·x.

Expanding the squared norms gives the exact identity

    F(x) = E(x) - m v₋·v₊ + ma·x₋ + (m/2)|v₋+ha|².

Consequently,

    ∇E(x) = ∇F(x) + m v₋/h.

From rest, v₋=0, the objective and energy differ only by a constant for every x, including with gravity. Imposing E(x)=E_target fixes the objective value. The objective cannot distinguish different feasible branches. In the primal Lagrangian F+η(E-E_target), η=-1 cancels the objective gradient everywhere. This is an exact degeneracy, not a numerical ill-conditioning claim.

For the code's scaled, negative-dual convention, the equations are

    R = m(x-y) - n λ - g η = 0,
    C + (α/h²)λ = 0,
    H = E-E_target = 0.

At rest, every feasible x has an exact KKT solution

    λ = -h² C/α,
    η = h²,

because the material force and energy multiplier together cancel the full objective gradient. A stationary spring x=x₋ is therefore a valid exact constrained optimum: its elastic and energy-multiplier forces cancel while it remains motionless. A coupled solve can find this solution immediately. The problem is not cured by solving the same equations more accurately.

## Scalar spring: the two exact endpoints

Take U(q)=kq²/2, no gravity, old velocity zero, and a_h=h²k/m. The hard energy equation is

    m(q₊-q₋)²/(2h²) + kq₊²/2 = kq₋²/2.

Its roots are

    q₊ = q₋,
    q₊ = [(1-a_h)/(1+a_h)] q₋.

The first freezes; the second has endpoint velocity

    v₊ = -[2hk/m /(1+a_h)] q₋.

The unconstrained backward-Euler endpoint q_BE=q₋/(1+a_h) lies exactly halfway between these roots and is the minimum of this endpoint-energy quadratic. Both energy-conserving roots have exactly the same implicit objective value. The nontrivial root initially gives approximately twice the physical velocity increment; exact energy alone does not identify the desired branch.

This directly supports the user's observation: the added multiplier can cancel the objective's driving force, and the remaining equations can be satisfied without producing the intended motion. It should be stated as a formulation degeneracy, not merely as a warning about solving two multipliers in the same block.

## Hard pendulum counterpart

For |x|=L, U=0, and a timestep starting at rest, the energy level is the circle

    |x-(x₋+h²a)|² = h⁴|a|².

It intersects the rod circle at x₋ and generally one other point. The second point is the reflection of x₋ about the line through the origin in direction y=x₋+h²a:

    x_other = 2 y (x₋·y)/|y|² - x₋.

For a nonvertical initial angle,

    x_other-x₋ = 2h² a_tangent + O(h⁴).

Again the stationary and moving endpoints have identical objective value; the nontrivial branch doubles the leading endpoint velocity increment. At a vertical equilibrium the two intersections coincide.

In the scaled equations, η=h² and λ=0 satisfy stationarity at every feasible endpoint from rest. For the stationary endpoint, the energy multiplier cancels gravity.

The first frozen linearization at the predictor has g=0 and generally H≠0, so its energy row is incompatible even though the nonlinear problem is feasible (x₋ itself is feasible). A skipped or incompatible frozen row must not be described as proof that the nonlinear constraints have no solution.

For a moving pendulum the energy level remains a circle centered at x₋+h²a with radius

    R² = (2h²/m)[E_target + ma·x₋ + (m/2)h²|a|²].

Its intersection with |x|=L exists exactly when the two circle radii and center distance satisfy the usual triangle inequalities. On that intersection, the objective chooses the endpoint maximizing v₋·v₊. At rest this directional preference vanishes completely.
