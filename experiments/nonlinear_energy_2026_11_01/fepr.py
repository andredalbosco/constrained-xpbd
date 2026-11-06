"""Planar FEPR Eq.8/13/14, with an explicitly labeled hard-rod extension.

See FEPR_DERIVATION.md. This is not a velocity rescale or full fixed-center SQP.
Only NumPy is needed. Run this file for derivative and numerical self-tests.
"""
from __future__ import annotations

import json
import math
import numpy as np


def _gravity(g):
    a = np.asarray(g, dtype=float)
    return np.array([0., -float(a)]) if a.ndim == 0 else a.reshape(2).copy()


def _cross(x, v):
    return float(x[0] * v[1] - x[1] * v[0])


def potential(x, m, k, l, g, rigid=False):
    x = np.asarray(x, dtype=float)
    return (0. if rigid else .5 * k * (np.linalg.norm(x) - l) ** 2) - m * _gravity(g).dot(x)


def energy(x, v, m, k, l, g, rigid=False):
    return float(.5 * m * np.dot(v, v) + potential(x, m, k, l, g, rigid))


def baseline_step(x, v, h, m, k, l, g, rigid=False):
    """Exact one-spring implicit Euler, or position/velocity RATTLE for a rod."""
    x = np.asarray(x, dtype=float).reshape(2).copy()
    v = np.asarray(v, dtype=float).reshape(2).copy()
    g = _gravity(g)
    label = "rattle" if rigid else "implicit_euler_exact_radial"
    result = dict(x=x.copy(), v=v.copy(), converged=False, residual=math.inf,
                  reason="not_started", counts={"baseline_evaluations": 1}, method=label)
    if not (h > 0 and m > 0 and l > 0 and k >= 0 and np.isfinite(np.r_[x, v, g, h, m, k, l]).all()):
        result["reason"] = "invalid_input"
        return result
    if not rigid:
        y = x + h * v + h * h * g
        radius = float(np.linalg.norm(y))
        if radius:
            direction = y / radius
        elif np.linalg.norm(x):
            direction = x / np.linalg.norm(x)
        else:
            result["reason"] = "undefined_radial_direction"
            return result
        a = m / h ** 2
        xp = (a * radius + k * l) / (a + k) * direction
        vp = (xp - x) / h
        grad = a * (xp - y) + k * (np.linalg.norm(xp) - l) * direction
        residual = float(np.linalg.norm(grad))
    else:
        if abs(np.linalg.norm(x) - l) > 1e-7 * max(1., l) or abs(x.dot(v)) > 1e-7 * max(1., l * np.linalg.norm(v)):
            result["reason"] = "incompatible_initial_rod"
            return result
        y = x + h * v + .5 * h * h * g
        a = float(x.dot(x))
        b = float(x.dot(y))
        c = float(y.dot(y) - l * l)
        discriminant = b * b - a * c
        if discriminant < 0:
            result["reason"] = "rattle_position_no_real_root"
            return result
        square = math.sqrt(discriminant)
        # Select the root continuous from alpha=0 as h -> 0. Do not clip disc.
        alpha = -c / (b + square) if b + square != 0 else (-b + square) / a
        xp = y + alpha * x
        half = (xp - x) / h
        trial = half + .5 * h * g
        vp = trial - xp * (xp.dot(trial) / xp.dot(xp))
        residual = float(max(abs(xp.dot(xp) - l * l), abs(xp.dot(vp))))
    result.update(x=xp, v=vp, converged=True, residual=residual, reason="converged")
    return result


def _constraints(q, old_x, old_v, star_x, star_v, h, m, k, l, g, rigid, target):
    x, v, s, t = q[:2], q[2:4], q[4], q[5]
    r = float(np.linalg.norm(x))
    grad = -m * g.copy()
    if not rigid and k:
        if r == 0 and l:
            raise ValueError("undefined_radial_gradient")
        if r:
            grad += k * (r - l) * x / r
    dp = m * (old_v - star_v)
    dl = m * (_cross(old_x, old_v) - _cross(star_x, star_v))
    c = np.zeros(6 if rigid else 4)
    j = np.zeros((c.size, 6))
    c[0] = energy(x, v, m, k, l, g, rigid) - target
    j[0, :2], j[0, 2:4] = grad, m * v
    c[1:3] = m * (v - star_v) - s * dp
    j[1:3, 2:4], j[1:3, 4] = m * np.eye(2), -dp
    c[3] = m * (_cross(x, v) - _cross(star_x, star_v)) - t * dl
    j[3, :2] = m * np.array([v[1], -v[0]])
    j[3, 2:4] = m * np.array([-x[1], x[0]])
    j[3, 5] = -dl
    if rigid:
        c[4], c[5] = .5 * (x.dot(x) - l * l), x.dot(v)
        j[4, :2] = x
        j[5, :2], j[5, 2:4] = v, x
    return c, j


def project_fepr(old_x, old_v, star_x, star_v, h, m, k, l, g,
                 rigid=False, *, energy_target=None, epsilon=None,
                 tolerance=1e-7, max_iterations=64, determinant_threshold=1e-12,
                 schur_regularization=1e-7, record_trace=False):
    """Direct fast FEPR projection of a supplied tentative state.

    The hard-rod rows are an extension, not part of the published FEPR method.
    No fallback state is substituted when convergence fails.
    """
    old_x = np.asarray(old_x, dtype=float).reshape(2)
    old_v = np.asarray(old_v, dtype=float).reshape(2)
    star_x = np.asarray(star_x, dtype=float).reshape(2)
    star_v = np.asarray(star_v, dtype=float).reshape(2)
    g = _gravity(g)
    if not (h > 0 and m > 0 and l > 0 and k >= 0 and tolerance > 0 and max_iterations >= 0):
        raise ValueError("invalid_projection_parameters")
    eps = 1e-3 * m * l * l if epsilon is None else float(epsilon)
    if eps <= 0:
        raise ValueError("epsilon_must_be_positive")
    target = energy(old_x, old_v, m, k, l, g, rigid) if energy_target is None else float(energy_target)
    metric = np.array([m, m, h * h * m, h * h * m, eps, eps])
    inverse = 1. / metric
    q = np.r_[star_x, star_v, 0., 0.]
    initial = q.copy()
    counts = {"iterations": 0, "constraint_evaluations": 0, "schur_solves": 0,
              "regularizations": 0, "baseline_evaluations": 0}
    reason, converged, residual = "iteration_limit", False, math.inf
    c = np.full(6 if rigid else 4, math.inf)
    trace = []
    for iteration in range(max_iterations + 1):
        try:
            c, j = _constraints(q, old_x, old_v, star_x, star_v, h, m, k, l, g, rigid, target)
            counts["constraint_evaluations"] += 1
        except (ValueError, FloatingPointError) as error:
            reason = str(error)
            break
        residual = float(np.abs(c).sum())
        if record_trace:
            trace.append({"iteration": iteration, "residual": residual, "q": q.tolist()})
        if not np.isfinite(q).all() or not np.isfinite(c).all() or not np.isfinite(j).all():
            reason = "nonfinite_projection"
            break
        if residual < tolerance:
            reason, converged = "converged", True
            break
        if iteration == max_iterations:
            break
        schur = (j * inverse) @ j.T
        sign, logdet = np.linalg.slogdet(schur)
        regularized = sign <= 0 or logdet < math.log(determinant_threshold)
        if regularized:
            schur += schur_regularization * np.eye(c.size)
            counts["regularizations"] += 1
        try:
            multipliers = np.linalg.solve(schur, c)
        except np.linalg.LinAlgError:
            reason = "schur_solve_failed"
            break
        counts["schur_solves"] += 1
        q -= inverse * (j.T @ multipliers)
        counts["iterations"] += 1
    names = ["energy", "linear_x", "linear_y", "angular"] + (["length_squared_half", "tangency"] if rigid else [])
    return {"x": q[:2].copy(), "v": q[2:4].copy(), "converged": converged,
            "residual": residual, "residuals": dict(zip(names, map(float, c))),
            "counts": counts, "reason": reason, "s": float(q[4]), "t": float(q[5]),
            "energy_target": target, "target_mode": "per_step" if energy_target is None else "fixed",
            "epsilon": eps, "tolerance": tolerance,
            "projection_objective": float(.5 * np.dot(metric * (q - initial), q - initial)),
            "star_x": star_x.copy(), "star_v": star_v.copy(),
            "method": "FEPR planar + hard-rod extension" if rigid else "FEPR planar Eq8/13/14",
            "trace": trace}


def fepr_step(x, v, h, m, k, l, g, rigid=False, *, energy_target=None,
              baseline=None, **projection_options):
    """Base step followed by paper FEPR, returning an explicit status dictionary."""
    expected = "rattle" if rigid else "implicit_euler"
    if baseline not in (None, expected):
        raise ValueError(f"supported baseline for rigid={rigid}: {expected}")
    tentative = baseline_step(x, v, h, m, k, l, g, rigid)
    if not tentative["converged"]:
        return dict(tentative, reason="baseline_" + tentative["reason"],
                    energy_target=energy(x, v, m, k, l, g, rigid) if energy_target is None else energy_target)
    result = project_fepr(x, v, tentative["x"], tentative["v"], h, m, k, l, g,
                          rigid, energy_target=energy_target, **projection_options)
    result["counts"]["baseline_evaluations"] = tentative["counts"]["baseline_evaluations"]
    result["baseline_method"] = tentative["method"]
    result["baseline_residual"] = tentative["residual"]
    return result


def self_test():
    """Independent FD derivatives, full KKT update, and physical rollout checks."""
    worst_gradient = 0.
    for rigid in (False, True):
        x, v = np.array([.3, -.96]), np.array([.7, .21])
        sx, sv = np.array([.35, -.94]), np.array([.68, .23])
        q = np.r_[sx + [.01, -.02], sv + [-.03, .04], .12, -.08]
        args = (x, v, sx, sv, .04, 1.3, 80., 1., np.array([0., -9.81]), rigid, 2.)
        c, j = _constraints(q, *args)
        fd = np.empty_like(j)
        for a in range(6):
            qp, qm = q.copy(), q.copy()
            qp[a] += 1e-6; qm[a] -= 1e-6
            fd[:, a] = (_constraints(qp, *args)[0] - _constraints(qm, *args)[0]) / 2e-6
        error = np.linalg.norm(fd - j) / (1 + np.linalg.norm(j))
        worst_gradient = max(worst_gradient, error)
        assert error < 1e-8, ("constraint_gradient", rigid, error)
    # Independently solve the full linear KKT system for one nondegenerate step.
    x, v = np.array([1.15, .15]), np.array([-.2, .7])
    base = baseline_step(x, v, .04, 1., 100., 1., [0., 0.])
    target = energy(x, v, 1., 100., 1., [0., 0.])
    q = np.r_[base["x"], base["v"], 0., 0.]
    c, j = _constraints(q, x, v, base["x"], base["v"], .04, 1., 100., 1., np.zeros(2), False, target)
    d = np.diag([1., 1., .04**2, .04**2, .001, .001])
    system = np.block([[d, j.T], [j, np.zeros((4, 4))]])
    expected = np.linalg.solve(system, np.r_[np.zeros(6), -c])[:6]
    got = project_fepr(x, v, base["x"], base["v"], .04, 1., 100., 1., [0., 0.],
                       max_iterations=1, determinant_threshold=1e-100)
    actual = np.r_[got["x"], got["v"], got["s"], got["t"]] - q
    update_error = float(np.linalg.norm(actual - expected))
    assert update_error < 1e-10, ("full_KKT_update", update_error)
    cases = [("spring_release", np.array([1.25, 0.]), np.zeros(2), False, np.zeros(2)),
             ("spring_orbit", np.array([1.1, 0.]), np.array([0., .8]), False, np.zeros(2)),
             ("rod_release", np.array([math.sin(.65), -math.cos(.65)]), np.zeros(2), True, np.array([0., -9.81])),
             ("rod_kick", np.array([math.sin(.65), -math.cos(.65)]), 1.2*np.array([math.cos(.65), math.sin(.65)]), True, np.array([0., -9.81]))]
    reports = []
    for name, x0, v0, rigid, gravity in cases:
        for h in (.01, .04):
            x, v = x0.copy(), v0.copy()
            initial = energy(x, v, 1., 100., 1., gravity, rigid)
            maximum_energy = maximum_geometry = maximum_tangency = 0.
            iterations = regularizations = 0
            for step in range(round(1/h)):
                result = fepr_step(x, v, h, 1., 100., 1., gravity, rigid,
                                   energy_target=initial, tolerance=1e-10, max_iterations=128)
                assert result["converged"], (name, h, step, result["reason"], result["residual"])
                x, v = result["x"], result["v"]
                maximum_energy = max(maximum_energy, abs(energy(x, v, 1., 100., 1., gravity, rigid)-initial))
                if rigid:
                    maximum_geometry = max(maximum_geometry, abs(np.linalg.norm(x)-1))
                    maximum_tangency = max(maximum_tangency, abs(x.dot(v)))
                iterations += result["counts"]["iterations"]
                regularizations += result["counts"]["regularizations"]
            assert maximum_energy < 1e-9
            assert np.linalg.norm(x - x0) > 1e-5, ("frozen_motion", name, h)
            assert maximum_geometry < 1e-9 and maximum_tangency < 1e-9
            reports.append(dict(case=name, h=h, maximum_energy=maximum_energy,
                                maximum_geometry=maximum_geometry, maximum_tangency=maximum_tangency,
                                iterations=iterations, regularizations=regularizations))
    # An already feasible free drift is unchanged; a too-large RATTLE step fails explicitly.
    x, v = np.array([1., 0.]), np.array([.2, .3])
    r = fepr_step(x, v, .1, 1., 0., 1., [0., 0.])
    assert r["converged"] and r["counts"]["iterations"] == 0
    assert np.linalg.norm(r["x"]-(x+.1*v)) < 1e-14
    fail = fepr_step(np.array([0., -1.]), np.array([100., 0.]), 1., 1., 100., 1., [0., -9.81], True)
    assert not fail["converged"] and fail["reason"] == "baseline_rattle_position_no_real_root"
    return dict(passed=True, worst_gradient=worst_gradient, full_KKT_update_error=update_error,
                trajectories=reports)


if __name__ == "__main__":
    print(json.dumps(self_test(), indent=2))
