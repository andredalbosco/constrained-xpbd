# The warm prediction is symplectic; the complete XPBD step is not

This result concerns an elastic constitutive multiplier at a fixed step size:

    lambda_minus = -h² alpha^-1 C(x_minus)
    J_minus^T lambda_minus = -h² grad U(x_minus).

The initial position is

    x0 = x_minus + h v_minus + h² a_gravity + M^-1 J_minus^T lambda_minus.

Together with v0=(x0-x_minus)/h, this is precisely kick-drift symplectic Euler:

    p0 = p_minus - h grad V(x_minus)
    x0 = x_minus + h M^-1 p0
    V(x) = U(x) - (M a_gravity)^T x, p=Mv.

The claim is about this position-momentum map, not a position initialization
in isolation. It assumes the old force is the elastic conservative force.
An arbitrary cached multiplier, an incompletely solved energy/material system,
or a hard pendulum reaction does not automatically satisfy this assumption.

## Exact proof

In canonical (x,p) coordinates the symplectic form is

    Omega = [0 I; -I 0].

The momentum kick and position drift have Jacobians

    K = [I 0; -h Hess(V) I]
    D = [I h M^-1; 0 I].

Hess(V) and M^-1 are symmetric. Direct multiplication gives
K^T Omega K=Omega and D^T Omega D=Omega, hence

    (D K)^T Omega (D K) = Omega.

Gravity adds a constant momentum kick and does not alter the Hessian. This is
an exact argument for a general smooth elastic potential, not an inference
from energy curves. Symplecticity does not mean exact conservation of the
original energy or unconditional stability.

## What the XPBD correction does

For U(q)=kq²/2 and a=h²k/m, the predictor Jacobian is

    S = [1-a, h/m; -h k, 1], det(S)=1.

After a complete scalar XPBD correction, whether cold or consistently warm,

    B = [1, h/m; -h k, 1]/(1+a), det(B)=1/(1+a).

In one degree of freedom, B^T Omega B=det(B) Omega. Therefore the completed
XPBD step contracts phase area and is not symplectic for k,h>0. Warm-starting
does not remove this contraction in the affine case. Reduced damping observed
on nonlinear constraints is not by itself evidence of symplecticity.

## Direct validation

`symplectic_probe.cpp` differentiates the map by complex steps in canonical
coordinates. Its nonlinear radial-spring double-valued map is checked against
the existing production correction routine. It also tests the previous-solution
reconstruction x0=x_prev+(prediction_new-prediction_prev) on an r_x=0 solution.

There are 216 Jacobian cases: scalar and nonlinear 2D springs, three radial
states, mass 0.3/1/4, h=1/240,1/480,1/960, with and without gravity, for the
predictor, cold XPBD, and warm XPBD. CSV rows record the input canonical state.

- Largest predictor defect max_abs(A^T Omega A-Omega): 2.78e-17.
- Largest nonlinear map parity error against the production update: 1.12e-16.
- Every tested complete XPBD map has a nonzero symplectic defect.
- At m=1,k=100,h=1/480 the scalar predictor has determinant 1; the complete
  warm XPBD map has determinant 0.9995661605206075 and defect 0.0004338394793926.
- The full 4x4 condition is checked for 2D springs; determinant one alone would
  not suffice there.

Sources and data: `symplectic_probe.cpp`, `symplectic_results/jacobians.csv`.
Build target: `xpbd-symplectic-probe`; CTest: `warm-start-symplecticity`.
The four tests in this experiment's standalone build pass. No energy constraint
or chain was added to this symplecticity test.
