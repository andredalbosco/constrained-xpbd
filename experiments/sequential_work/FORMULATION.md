# Sequential local work feedback in XPBD

This audit concerns the existing frozen previous-force work estimate and signed reservoirs. Changing the GS ordering does not change that work model or turn it into exact mechanical-energy conservation.

## Residuals and scalar updates

At the start of the step freeze each spring force on endpoint a as `f_j^-`; endpoint b receives its negative. Let `F_i^-` be the sum over *all* incident springs, `dx_i=x_i-x_i^-`, and `h` the current step size. Then

    b_i = .5 m_i (|dx_i|^2/h^2 - |v_i^-|^2)
          - (F_i^- + m_i gravity) dot dx_i,
    d_j = U_j(x)-U_j(x^-) + f_j^- dot (dx_a-dx_b).

Use one accumulated material lambda per spring, one accumulated `qnu=h^2 nu` per spring, and one accumulated `qmu=h^2 mu` per particle. The rows are

    rlambda_j = C_j + (alpha_j/h^2) lambda_j,
    rnu_j     = d_j - S_j^- + (gamma/h^2) qnu_j,
    rmu_i     = b_i - R_i^- + (beta/h^2) qmu_i.

For one scalar row with balance gradient G, accumulated dual q, and compliance c, the curvature-dropped XPBD update is

    dq = -r / (G^T M^-1 G + c/h^2),
    dx = M^-1 G dq.

For the material row use G=grad C and c=alpha. For work rows the analytical gradients are

    a_i = grad b_i = m_i dx_i/h^2 - F_i^- - m_i gravity,
    z_ja = grad_a d_j = (C_j/alpha_j) grad_a C_j + f_j^-,
    z_jb = -z_ja.

For a hard material constraint alpha=0, its analytical stored energy is U=0, so `z_ja=f_j^-`. Never substitute the hard material lambda for an elastic energy gradient.

For a distance spring with inverse masses wa,wb and unit direction g, the three scalar denominators are respectively

    wa+wb+alpha/h^2,
    (wa+wb)|z_ja|^2+gamma/h^2,
    wi|a_i|^2+beta/h^2.

The requested local schedule is `lambda_j -> nu_j -> mu_a -> mu_b`. Recompute the potential balance and gradient after the accepted material update; recompute each kinetic balance and gradient immediately before its update. Positive beta,gamma make even zero-gradient work rows solvable by changing only the multiplier. With backtracking, scale and accept/revert the position and accumulated multiplier together, testing the actual nonlinear row.

## Shared particle rows and scheduling comparisons

The two endpoint updates are visits to particle-owned rows, not new per-edge kinetic constraints. A particle's qmu is shared across every incident spring and its b_i always uses the complete frozen incident-force sum F_i^-. Resetting qmu per edge or restricting F_i^- to the current edge changes the model and duplicates kinetic energy.

Local scheduling visits particle i degree(i) times per sweep. This is valid nonlinear GS rescheduling; beta must not be multiplied by degree merely to compensate. Visit isolated particles once per sweep so their external-force work row is still handled. Fixed particles can update compliant dual values without moving geometry.

`PairedLocal` and `Sequential` provide the clean scheduling comparison: both use the same endpoint visits, while only the material/potential pair changes from a 2x2 solve to two scalar solves. Comparing the old global `Paired` schedule with `Sequential` also changes kinetic-row visit frequency and when neighboring corrections are observed. None of these nonlinear schedules is order-independent at finite sweeps.

## What compliance means as a strength

The reservoir formulation can be obtained from the physical objective

    I(x) + sum_j U_j(x) + sum_i R_i^2/(2 beta) + sum_j S_j^2/(2 gamma),

subject to `b_i+R_i-R_i^-=0` and `d_j+S_j-S_j^-=0`. Eliminating R,S gives `R=beta mu`, `S=gamma nu` and the quadratic dual terms `-beta mu^2/2`, `-gamma nu^2/2`. Thus beta and gamma have units of energy; larger values mean cheaper reservoir storage and weaker positional feedback. They are genuine compliances, not dimensionless relaxation factors.

For one frozen work row, let `D=G^T M^-1 G`. The fraction of its residual correction supplied by moving positions is

    strength = D/(D+c/h^2) = h^2 D/(h^2 D+c).

The remaining fraction is supplied by changing the reservoir multiplier. A desired local strength s in (0,1) corresponds to `c=h^2 D (1-s)/s`; this is a state-dependent calibration. For a fixed parameter or UI control choose a reference state and timestep, then hold the resulting energy-valued compliance fixed. Recomputing c to hold s fixed every iteration is an adaptive-compliance algorithm and changes the objective. beta,gamma should be separately calibrated because kinetic and potential balance gradients have different scales.

As compliance grows, positional feedback weakens toward passive bookkeeping. As it approaches zero, the corresponding reservoir is increasingly penalized and the row approaches a hard work condition, with the usual compatibility and conditioning concerns. Positive compliance avoids a zero denominator but does not guarantee nonlinear stability or physical accuracy.

## Accounting and convergence boundaries

Closing accepted reservoirs algebraically,

    R_i^+ = R_i^- - b_i(x_accepted),
    S_j^+ = S_j^- - d_j(x_accepted),

conserves the sum of mechanical energy plus signed reservoirs up to arithmetic error, because the same frozen spring work appears with opposite signs in b and d, and gravity work is the negative gravitational-potential change. This accounting identity holds even with incomplete GS iteration. It does not imply preservation of mechanical energy alone, nor convergence of the feedback equations `R=beta qmu/h^2`, `S=gamma qnu/h^2`. Continue reporting feedback, material, and actual position-stationarity residuals separately.

These scalar updates retain the existing omission of the top stationarity residual in addition to dropping curvature. A converged set of scalar balance rows alone need not solve the mixed stationary equations. Local scheduling does not add missing momentum constraints: pairwise material/potential corrections have opposite endpoint impulses, but the frozen-force term need not be central at the current geometry, and individual kinetic corrections need not cancel across particles.

## Why implicit work needs more than this ordering

With current impulse work `W_ij=I_ij dot (v_i^-+v_i^+)/2`, exact particle impulse balance makes the summed work equal kinetic-energy change. With midpoint drift and an elastic discrete gradient, the potential balance factors as

    Delta C_j (Cbar_j/alpha_j + force_coefficient_j).

Those identities require a consistent current-impulse/velocity/position construction. The kinetic row can then be redundant, while the potential row becomes degenerate when Delta C_j=0; the material constitutive row must remain. The present frozen-force rows have different derivatives and no such automatic identity. Merely updating lambda,nu,mu in that order does not convert them into implicit-work rows. Doing that would require changing the unknowns/residuals and retaining the derivatives with respect to current impulse and mean velocity, with any needed impulse-balance multipliers explicitly included.

## Read-only implementation review

The staged `xpbd.hpp` implements the stated order and scalar denominators. It recomputes the potential gradient after the material correction, keeps qmu particle-owned, includes all incident frozen forces, and visits isolated particles once per sweep. Both-pinned and collapsed-edge branches still visit endpoint kinetic rows under local ordering. The original XPBD/Passive arithmetic and the original Paired update sequence are unchanged by source inspection; this statement is not a substitute for numerical parity tests.

The new scalar backtracking tests only the row currently being solved. A successful local update can therefore increase other rows or the full stationarity residual. Full rx remains omitted. The work estimate is deliberately held at the previous force to isolate scheduling: the implicit-work construction is not implemented by these nu/mu projections.
