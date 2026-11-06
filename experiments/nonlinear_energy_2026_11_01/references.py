"""Independent continuous references and discrete-gradient diagnostic controls.

The spring is anchored at the origin and has U=.5*k*(|x|-l)**2-m*g.x.
The pendulum angle convention is x=l*(sin(theta),-cos(theta)). A scalar g
means downward acceleration (0,-g); a two-vector specifies arbitrary gravity.

These controls change endpoint kinematics. They do NOT solve the original
endpoint-velocity EL objective. Check StepResult.converged before advancing.
No energy projection, velocity rescaling, force clamp, or step subdivision is
applied. Failed solves retain their candidate and diagnostic residual.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq, least_squares, root


@dataclass
class StepResult:
    x: np.ndarray
    v: np.ndarray
    converged: bool
    residual: float  # infinity norm of position/angle equation, divided by scale
    energy_error: float  # signed mechanical energy change, in physical units
    nfev: int
    reason: str
    method: str
    theta: Optional[float] = None
    omega: Optional[float] = None
    delta_theta: Optional[float] = None  # unwrapped increment, including rotations
    unique_root_certified: bool = False


@dataclass
class ReferenceResult:
    times: np.ndarray
    x: np.ndarray
    v: np.ndarray
    success: bool
    message: str
    nfev: int
    theta: Optional[np.ndarray] = None
    omega: Optional[np.ndarray] = None


def gravity_vector(g):
    a = np.asarray(g, dtype=float)
    if a.ndim == 0:
        a = np.array([0.0, -float(a)])
    if a.shape != (2,) or not np.all(np.isfinite(a)):
        raise ValueError("gravity must be a finite scalar or planar vector")
    return a


def _state(x, v, m, l):
    x, v = np.asarray(x, dtype=float).copy(), np.asarray(v, dtype=float).copy()
    if x.shape != (2,) or v.shape != (2,) or not np.all(np.isfinite([x, v])):
        raise ValueError("x and v must be finite planar vectors")
    if not np.isfinite(m) or m <= 0 or not np.isfinite(l) or l < 0:
        raise ValueError("require m>0 and l>=0")
    return x, v


def spring_discrete_gradient(x0, x1, k, l):
    """D.U=Delta U: D=k/2*(r0+r1-2l)*(x0+x1)/(r0+r1).

    At equal nonzero endpoints this is the usual spring gradient. When l=0
    its nonsingular expression is k*(x0+x1)/2, including the origin. When
    l>0 and both endpoints are the origin, no unique gradient is defined.
    """
    x0, x1 = np.asarray(x0, dtype=float), np.asarray(x1, dtype=float)
    if k == 0:
        return np.zeros_like(x0)
    if l == 0:
        return 0.5 * k * (x0 + x1)
    s = np.linalg.norm(x0) + np.linalg.norm(x1)
    if s == 0:
        raise ValueError("positive-rest-length spring gradient undefined at origin")
    return 0.5 * k * (s - 2 * l) * (x0 + x1) / s


def spring_energy(x, v, m, k, l, g=0.0):
    x, v = np.asarray(x), np.asarray(v)
    return float(0.5 * m * np.dot(v, v) + 0.5 * k * (np.linalg.norm(x) - l)**2
                 - m * np.dot(gravity_vector(g), x))


def spring_dg_step(x, v, h, m, k, l, g=0.0, *, rtol=1e-10, maxfev=500):
    """Solve midpoint kinematics and symmetric radial discrete-gradient force.

    x1-x0=h*(v0+v1)/2; m*(v1-v0)=-h*D_s+h*m*g.
    scipy hybr starts from the acceleration predictor. On failure a reported
    least-squares fallback solves the same residual, with no constraints or
    regularization. If it still fails, reported hybr retries use the old state
    and free-flight predictor. maxfev is per attempt, nfev counts all calls.
    A small residual certifies a root, not a unique branch.
    """
    x, v = _state(x, v, m, l)
    gv = gravity_vector(g)
    if not np.isfinite(k) or k < 0 or not np.isfinite(h):
        raise ValueError("require finite k>=0 and finite h")
    if h == 0:
        return StepResult(x, v, True, 0., 0., 0, "zero step", "identity", unique_root_certified=True)
    if rtol <= 0 or maxfev < 1:
        raise ValueError("require positive rtol and maxfev")
    r0 = np.linalg.norm(x)
    if r0 == 0 and l > 0 and k > 0:
        raise ValueError("continuous initial spring force undefined at origin")
    force0 = np.zeros(2) if k == 0 else spring_discrete_gradient(x, x, k, l)
    predictor = x + h*v + .5*h*h*(gv-force0/m)
    scale = max(l, r0, abs(h)*np.linalg.norm(v), h*h*np.linalg.norm(gv), 1e-12)
    calls = 0

    def residual(q):
        nonlocal calls
        calls += 1
        ds = np.zeros(2) if k == 0 else spring_discrete_gradient(x, q, k, l)
        return (q-x-h*v + .5*h*h*(ds/m-gv))/scale

    ans = root(residual, predictor, method="hybr", options={"xtol": 1e-11, "maxfev": maxfev})
    q = np.asarray(ans.x)
    norm = float(np.linalg.norm(residual(q), ord=np.inf))
    method, message = "hybr", str(ans.message)
    if not np.isfinite(norm) or norm > rtol:
        start = q if np.all(np.isfinite(q)) else predictor
        fallback = least_squares(residual, start, ftol=1e-13, xtol=1e-13,
                                 gtol=1e-13, max_nfev=maxfev)
        alt = float(np.linalg.norm(residual(fallback.x), ord=np.inf))
        if np.isfinite(alt) and (not np.isfinite(norm) or alt < norm):
            q, norm = fallback.x, alt
        method = "hybr+least_squares"
        message += "; " + str(fallback.message)
    if not np.isfinite(norm) or norm > rtol:
        message += f"; primary residual={norm:.6g}"
        for label, seed in [("old", x), ("free", x+h*v+.5*h*h*gv)]:
            retry = root(residual, seed, method="hybr", options={"xtol": 1e-11, "maxfev": maxfev})
            alt = float(np.linalg.norm(residual(retry.x), ord=np.inf))
            method += "+hybr_"+label
            message += f"; {label} retry residual={alt:.6g}: {retry.message}"
            if np.isfinite(alt) and (not np.isfinite(norm) or alt < norm):
                q, norm = retry.x, alt
            if np.isfinite(norm) and norm <= rtol:
                break
    v1 = 2*(q-x)/h-v
    r1 = np.linalg.norm(q)
    de = (.5*m*np.dot(v1+v, v1-v)
          + .5*k*(r1+r0-2*l)*(r1-r0)-m*np.dot(gv, q-x))
    ok = bool(np.isfinite(norm) and norm <= rtol and np.all(np.isfinite(v1)))
    # The l=0 spring and k=0 free particle have a unique linear position solve.
    return StepResult(q.copy(), v1, ok, norm, float(de), calls,
                      ("verified residual; " if ok else "residual not converged; ")+message,
                      method, unique_root_certified=bool(l == 0 or k == 0))


def _pendulum_initial(x, v, m, l, feasibility_tolerance=1e-12):
    x, v = _state(x, v, m, l)
    if l <= 0:
        raise ValueError("pendulum requires l>0")
    if abs(np.linalg.norm(x)-l) > feasibility_tolerance*l:
        raise ValueError("initial pendulum position violates rod length")
    if abs(np.dot(x, v))/l > feasibility_tolerance*max(1., np.linalg.norm(v)):
        raise ValueError("initial pendulum velocity is not endpoint tangent")
    theta = float(np.arctan2(x[0], -x[1]))
    omega = float(np.dot(np.array([np.cos(theta), np.sin(theta)]), v)/l)
    return x, v, theta, omega


def _sinc_half(delta):
    return float(np.sinc(delta/(2*np.pi)))


def pendulum_dg_step(x, v, h, m, l, g=0.0, *, rtol=1e-10, maxfev=500):
    """Angular discrete-gradient control with exact endpoint length/tangency.

    theta1-theta0=h*(omega0+omega1)/2, and
    m*l*l*(omega1-omega0)=-h*D_theta U. The angular increment is unwrapped;
    rotations are not folded back onto [-pi,pi]. This is not Cartesian
    midpoint kinematics. A bracket exists about h*omega0 with radius
    h*h*|g|/(2*l); brentq verifies one root of it. Uniqueness is certified
    only when h*h*|g|/(4*l)<1 (or gravity is zero).
    Inputs must satisfy length and tangency to relative tolerance 1e-12;
    larger defects are rejected, and output energy uses the actual input.
    """
    x, v, theta, omega = _pendulum_initial(x, v, m, l)
    gv = gravity_vector(g)
    if not np.isfinite(h) or rtol <= 0 or maxfev < 1:
        raise ValueError("require finite h and positive rtol/maxfev")
    if h == 0:
        return StepResult(x, v, True, 0., 0., 0, "zero step", "identity",
                          theta, omega, 0., True)
    inertia = m*l*l
    center = h*omega
    radius = .5*h*h*np.linalg.norm(gv)/l
    scale = max(1., abs(center), radius)
    calls = 0

    def angular_gradient(delta):
        mid = theta+.5*delta
        return -m*l*_sinc_half(delta)*np.dot(gv, [np.cos(mid), np.sin(mid)])

    def equation(delta):
        nonlocal calls
        calls += 1
        return (delta-center+.5*h*h*angular_gradient(delta)/inertia)/scale

    unique = bool(.25*h*h*np.linalg.norm(gv)/l < 1.)
    if radius == 0:
        delta, message, method = center, "zero-gravity exact rotation", "exact"
    else:
        pad = 16*np.finfo(float).eps*scale
        left, right = center-radius-pad, center+radius+pad
        try:
            delta, info = brentq(equation, left, right, xtol=2e-14,
                                rtol=4*np.finfo(float).eps, maxiter=maxfev, full_output=True)
            message, method = str(info.flag), "brentq"
        except (ValueError, RuntimeError) as exc:
            delta, message, method = center, str(exc), "brentq_failed"
    norm = float(abs(equation(delta)))
    theta1, omega1 = theta+delta, 2*delta/h-omega
    x1 = l*np.array([np.sin(theta1), -np.cos(theta1)])
    v1 = l*omega1*np.array([np.cos(theta1), np.sin(theta1)])
    # Measure the actual Cartesian input/output energy, also exposing any
    # roundoff-size initial feasibility mismatch accepted by the input check.
    de = .5*m*np.dot(v1+v, v1-v)-m*np.dot(gv, x1-x)
    ok = bool(np.isfinite(norm) and norm <= rtol and method != "brentq_failed")
    reason = ("verified residual; " if ok else "residual not converged; ")+message
    if not unique:
        reason += "; finite-step angular branch uniqueness not certified"
    return StepResult(x1, v1, ok, norm, float(de), calls, reason, method,
                      theta1, omega1, float(delta), unique)


def reference_trajectory(kind, x, v, times, m, k=0.0, l=1.0, g=0.0, *,
                         rtol=2e-12, atol=2e-14, max_step=np.inf):
    """DOP853 reference with dense output at requested nonnegative times.

    The initial condition is at t=0. On integrator failure only the valid
    requested prefix is returned, with success=False; there is no extrapolation.
    Spring ODE force is undefined at the exact origin for positive rest length.
    Pendulum reference integrates unwrapped angle/angular velocity, not a
    separately projected Cartesian ODE. Refine tolerances/max_step independently
    before treating a reference as accuracy evidence.
    """
    x, v = _state(x, v, m, l)
    gv = gravity_vector(g)
    times = np.asarray(times, dtype=float)
    if times.ndim != 1 or len(times) == 0 or not np.all(np.isfinite(times)) or np.any(times < 0) or np.any(np.diff(times) < 0):
        raise ValueError("times must be a nonempty nondecreasing nonnegative vector")
    if kind == "spring":
        if not np.isfinite(k) or k < 0:
            raise ValueError("require finite k>=0")
        initial = np.concatenate([x, v])

        def ode(t, y):
            r = np.linalg.norm(y[:2])
            if r == 0 and k > 0 and l > 0:
                raise ValueError("spring ODE reached its nondifferentiable origin")
            force = k*y[:2] if l == 0 else (k*(r-l)*y[:2]/r if k > 0 else np.zeros(2))
            return np.concatenate([y[2:], gv-force/m])

    elif kind == "pendulum":
        x, v, theta, omega = _pendulum_initial(x, v, m, l)
        initial = np.array([theta, omega])

        def ode(t, y):
            return np.array([y[1], np.dot(gv, [np.cos(y[0]), np.sin(y[0])])/l])

    else:
        raise ValueError("kind must be 'spring' or 'pendulum'")
    if times[-1] == 0:
        values = np.repeat(initial[:, None], len(times), axis=1)
        success, message, nfev = True, "zero duration", 0
    else:
        sol = solve_ivp(ode, (0., float(times[-1])), initial, method="DOP853",
                        rtol=rtol, atol=atol, dense_output=True, max_step=max_step)
        times = times[times <= sol.t[-1]]
        values = sol.sol(times)
        success, message, nfev = bool(sol.success), str(sol.message), int(sol.nfev)
    if kind == "spring":
        return ReferenceResult(times, values[:2].T, values[2:].T, success, message, nfev)
    theta, omega = values
    positions = l*np.column_stack([np.sin(theta), -np.cos(theta)])
    velocities = l*omega[:, None]*np.column_stack([np.cos(theta), np.sin(theta)])
    return ReferenceResult(times, positions, velocities, success, message, nfev, theta, omega)


def _self_test():
    rng = np.random.default_rng(731)
    max_chain = 0.
    for _ in range(256):
        a, b = rng.normal(size=(2, 2))
        k, l = 17., .8
        du = .5*k*(np.linalg.norm(b)+np.linalg.norm(a)-2*l)*(np.linalg.norm(b)-np.linalg.norm(a))
        max_chain = max(max_chain, abs(np.dot(spring_discrete_gradient(a, b, k, l), b-a)-du))
    assert max_chain < 1e-12
    assert np.array_equal(spring_discrete_gradient([0., 0.], [0., 0.], 3., 0.), [0., 0.])
    spring = (np.array([1.2, -.3]), np.array([.2, .9]))
    q, w = .65, 1.2
    rod = (np.array([np.sin(q), -np.cos(q)]), w*np.array([np.cos(q), np.sin(q)]))
    maxima = {}
    for kind, initial in [("spring", spring), ("pendulum", rod)]:
        x, v = (a.copy() for a in initial)
        worst_de, worst_res, worst_tan, worst_len, nfev = 0., 0., 0., 0., 0
        for _ in range(200):
            oldx, oldv = x.copy(), v.copy()
            ans = (spring_dg_step(x, v, .02, 1., 40., 1., 3.) if kind == "spring"
                   else pendulum_dg_step(x, v, .02, 1., 1., 9.81))
            assert ans.converged, ans.reason
            x, v = ans.x, ans.v
            worst_de, worst_res = max(worst_de, abs(ans.energy_error)), max(worst_res, ans.residual)
            nfev += ans.nfev
            if kind == "spring":
                assert np.linalg.norm(x-oldx-.01*(v+oldv)) < 1e-13
            else:
                worst_len = max(worst_len, abs(np.linalg.norm(x)-1.))
                worst_tan = max(worst_tan, abs(np.dot(x, v)))
                oldomega = oldx[0]*oldv[1]-oldx[1]*oldv[0]
                assert abs(ans.delta_theta-.01*(ans.omega+oldomega)) < 1e-13
        assert worst_de < 1e-9 and worst_len < 1e-13 and worst_tan < 1e-13
        maxima[kind] = dict(max_step_energy=worst_de, max_residual=worst_res,
                            max_length=worst_len, max_tangency=worst_tan, nfev=nfev)
        forward = (spring_dg_step(*initial, .08, 1., 40., 1., 3.) if kind == "spring"
                   else pendulum_dg_step(*initial, .08, 1., 1., 9.81))
        backward = (spring_dg_step(forward.x, forward.v, -.08, 1., 40., 1., 3.) if kind == "spring"
                    else pendulum_dg_step(forward.x, forward.v, -.08, 1., 1., 9.81))
        assert forward.converged and backward.converged
        reverse = max(np.linalg.norm(backward.x-initial[0]), np.linalg.norm(backward.v-initial[1]))
        assert reverse < 1e-10
        maxima[kind]["reversal_error"] = reverse
    try:
        pendulum_dg_step([0., -1.], [0., 1.], .1, 1., 1., 9.81)
        raise AssertionError("radial initial velocity was silently accepted")
    except ValueError:
        pass
    rotation = pendulum_dg_step([0., -1.], [10., 0.], 1., 1., 1., 0.)
    assert rotation.converged and rotation.delta_theta == 10. and rotation.omega == 10.
    failed = pendulum_dg_step(*rod, .2, 1., 1., 9.81, maxfev=1)
    assert not failed.converged and failed.residual > 1e-3
    angular = pendulum_dg_step(*rod, .08, 1., 1., 9.81)
    cartesian_midpoint_gap = np.linalg.norm(angular.x-rod[0]-.04*(angular.v+rod[1]))
    assert cartesian_midpoint_gap > 1e-5
    for kind, initial, gravity in [("spring", spring, 3.), ("pendulum", rod, 9.81)]:
        t = np.linspace(0., 2., 101)
        ref = reference_trajectory(kind, *initial, t, 1., 40., 1., gravity)
        fine = reference_trajectory(kind, *initial, t, 1., 40., 1., gravity,
                                    rtol=2e-13, atol=2e-15, max_step=.01)
        assert ref.success and fine.success
        gap = max(np.max(np.linalg.norm(ref.x-fine.x, axis=1)), np.max(np.linalg.norm(ref.v-fine.v, axis=1)))
        assert gap < 1e-8
        maxima[kind]["reference_refinement_gap"] = float(gap)
        errors = []
        for h in [.1, .05, .025]:
            xx, vv = (a.copy() for a in initial)
            for _ in range(round(2/h)):
                ans = (spring_dg_step(xx, vv, h, 1., 40., 1., gravity) if kind == "spring"
                       else pendulum_dg_step(xx, vv, h, 1., 1., gravity))
                assert ans.converged
                xx, vv = ans.x, ans.v
            errors.append(dict(h=h, position=float(np.linalg.norm(xx-fine.x[-1])),
                               velocity=float(np.linalg.norm(vv-fine.v[-1]))))
        assert errors[1]["position"]/errors[2]["position"] > 3.5
        assert errors[1]["velocity"]/errors[2]["velocity"] > 3.5
        maxima[kind]["two_second_step_refinement"] = errors
    return dict(discrete_chain_rule_max_error=max_chain, cases=maxima,
                forced_failure_residual=failed.residual,
                angular_control_cartesian_midpoint_gap=float(cartesian_midpoint_gap))


if __name__ == "__main__":
    import json
    print(json.dumps(_self_test(), indent=2))
