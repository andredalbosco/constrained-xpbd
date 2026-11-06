"""Pinned energy-only specialization of Kee et al. CPD (2021).

Matches dynamics-research/experiments/work_balance/baselines.cpp: constant PD
metric, full undamped coupled SQP steps, a newly solved multiplier each step,
unbounded energy relaxation rho, first-order external-power energy target,
and endpoint/backward velocity. No projection or solver fallback is inserted.

cpd_pendulum_step adds an exact length equality as an EXPLICIT ADAPTATION.
Published/comparison CPD requires finite compliance and does not provide this
hard-rod method. Its backward velocity is generally not endpoint tangent.
"""
from dataclasses import dataclass, field
import numpy as np

from references import StepResult, gravity_vector


@dataclass
class CPDResult(StepResult):
    rho: float = 0.0
    energy_target: float = 0.0
    compatible_energy: float = 0.0
    energy_balance: float = 0.0
    stationarity: float = np.inf
    raw_stationarity: float = np.inf
    rho_stationarity: float = np.inf
    geometry: float = 0.0
    tangency: float = 0.0
    iterations: int = 0
    schur_rank: int = 0
    constraint_jacobian_rank: int = 0
    smallest_constraint_singular_value: float = np.nan
    multiplier: np.ndarray = field(default_factory=lambda: np.empty(0))
    residual_history: list = field(default_factory=list)


def _step(x, v, h, m, k, l, g, epsilon, max_iterations, tol, relax, rigid):
    x0, v0 = np.asarray(x, dtype=float).copy(), np.asarray(v, dtype=float).copy()
    gv = gravity_vector(g)
    if x0.shape != (2,) or v0.shape != (2,) or not np.all(np.isfinite([x0, v0])):
        raise ValueError("require finite planar x and v")
    if not all(np.isfinite(a) for a in [h, m, k, l, epsilon, tol]):
        raise ValueError("require finite parameters")
    if min(h, m, l, epsilon, tol) <= 0 or k < 0 or max_iterations < 1:
        raise ValueError("require h,m,l,epsilon,tol>0, k>=0 and max_iterations>=1")
    if not rigid and k == 0:
        raise ValueError("source spring CPD requires finite positive stiffness")
    if np.linalg.norm(x0) < 1e-14*l:
        raise ValueError("initial distance normal undefined")
    if rigid and abs(np.linalg.norm(x0)-l) > 1e-10*l:
        raise ValueError("hard-rod adaptation requires feasible initial length")
    oldU = 0. if rigid else .5*k*(np.linalg.norm(x0)-l)**2
    oldH = .5*m*np.dot(v0, v0)+oldU
    target = float(oldH+h*m*np.dot(gv, v0))
    prediction = x0+h*v0+h*h*gv
    a = m/(h*h)+(0. if rigid else k)
    inv_metric = np.array([1/a, 1/a]+([1/epsilon] if relax else []))
    q = np.concatenate([prediction, [0.]]) if relax else prediction.copy()
    nrows = 2 if rigid else 1
    multiplier = np.zeros(nrows)
    history, evaluations, rank = [], 0, 0
    reason, converged = "maximum_iterations", False
    energy_scale = m*l*l/(h*h)+abs(target)

    def evaluate(q):
        nonlocal evaluations
        evaluations += 1
        xx = q[:2]
        r = np.linalg.norm(xx)
        if not np.all(np.isfinite(q)) or (not rigid and r < 1e-14*l):
            raise FloatingPointError("nonfinite iterate or undefined distance gradient")
        elastic = 0. if rigid else .5*k*(r-l)**2
        gradU = np.zeros(2) if rigid else k*(r-l)*xx/r
        dx = xx-x0
        kinetic = .5*m*np.dot(dx, dx)/(h*h)
        rho = q[2] if relax else 0.
        grad = np.zeros(len(q))
        grad[:2] = m*(xx-prediction)/(h*h)+gradU
        if relax:
            grad[2] = epsilon*rho
        constraint = [kinetic+elastic-(1-rho)*target]
        rows = [np.concatenate([m*dx/(h*h)+gradU, [target]]) if relax
                else m*dx/(h*h)+gradU]
        if rigid:
            # This extra equality is the documented hard-rod extension.
            constraint.append(.5*(np.dot(xx, xx)-l*l))
            rows.append(np.concatenate([xx, [0.]]) if relax else xx.copy())
        return grad, np.asarray(constraint), np.asarray(rows), elastic+kinetic

    def diagnostics(q, multiplier):
        grad, c, j, hvalue = evaluate(q)
        station = grad+j.T@multiplier
        scaled_station = np.linalg.norm(station[:2], ord=np.inf)*h*h/(m*l)
        scaled_rho = abs(station[2])/(epsilon+abs(target)) if relax else 0.
        scaled_constraint = abs(c[0])/energy_scale
        if rigid:
            scaled_constraint = max(scaled_constraint, abs(c[1])/(l*l))
        return max(scaled_station, scaled_rho, scaled_constraint), scaled_station, scaled_rho, station, c, hvalue, j

    last = None
    iterations = 0
    for iteration in range(max_iterations):
        try:
            grad, c, j, _ = evaluate(q)
            inverse_gradient = inv_metric*grad
            inverse_rows = inv_metric[:, None]*j.T
            schur = j@inverse_rows
            diagonal = np.sqrt(np.maximum(np.diag(schur), 0.))
            row_scale = np.ones(nrows)
            positive = diagonal > 0
            row_scale[positive] = 1/diagonal[positive]
            scaled = row_scale[:, None]*schur*row_scale[None, :]
            rhs = row_scale*(c-j@inverse_gradient)
            sol, _, rank, _ = np.linalg.lstsq(scaled, rhs, rcond=1e-12)
            multiplier = row_scale*sol  # newly solved multiplier, not an increment
            q = q-inverse_gradient-inverse_rows@multiplier
            iterations = iteration+1
            last = diagnostics(q, multiplier)
            history.append(float(last[0]))
            if not np.isfinite(last[0]):
                reason = "nonfinite_residual"
                break
            if last[0] <= tol:
                converged, reason = True, "converged"
                break
        except (FloatingPointError, np.linalg.LinAlgError) as exc:
            reason = str(exc)
            last = None
            break
    xx, vv = q[:2].copy(), (q[:2]-x0)/h
    if last is None:
        residual = stationarity = rho_stationarity = raw_stationarity = energy_balance = energy_error = np.inf
        jacobian_rank, smallest_singular = 0, np.nan
    else:
        residual, stationarity, rho_stationarity, station, c, hvalue, j = last
        raw_stationarity = float(np.linalg.norm(station[:2], ord=np.inf))
        energy_balance = float(c[0])
        energy_error = float(hvalue-oldH-m*np.dot(gv, xx-x0))
        # Unlike the row-normalized Schur rank, this diagnostic retains the
        # vanishing energy gradient at zero-motion/zero-target degeneracy.
        column_units = np.array([l, l]+([1.] if relax else []))
        row_units = np.array([energy_scale]+([l*l] if rigid else []))
        singular = np.linalg.svd(j*column_units[None, :]/row_units[:, None], compute_uv=False)
        jacobian_rank = int(np.count_nonzero(singular > 1e-12*max(1., singular[0])))
        smallest_singular = float(singular[-1])
    return CPDResult(x=xx, v=vv, converged=converged, residual=float(residual),
                     energy_error=energy_error, nfev=evaluations, reason=reason,
                     method=("cpd_energy_hard_length_adaptation" if rigid else "cpd_pinned_energy_only")
                            +("_relaxed" if relax else "_strict"),
                     rho=float(q[2]) if relax else 0., energy_target=target,
                     energy_balance=energy_balance, stationarity=float(stationarity),
                     raw_stationarity=raw_stationarity, rho_stationarity=float(rho_stationarity),
                     geometry=float(abs(np.linalg.norm(xx)-l)) if rigid else 0.,
                     tangency=float(abs(np.dot(xx, vv))/l) if rigid else 0.,
                     iterations=iterations, schur_rank=int(rank), multiplier=multiplier,
                     constraint_jacobian_rank=jacobian_rank,
                     smallest_constraint_singular_value=smallest_singular,
                     residual_history=history)


def cpd_spring_step(x, v, h, m, k, l, g=0.0, *, epsilon=1e8,
                    max_iterations=200, tol=1e-9, relax=True):
    """Source-matched pinned finite-compliance energy-only CPD specialization."""
    return _step(x, v, h, m, k, l, g, epsilon, max_iterations, tol, relax, False)


def cpd_pendulum_step(x, v, h, m, l, g=0.0, *, epsilon=1e8,
                      max_iterations=200, tol=1e-9, relax=True):
    """EXPLICIT hard-length extension; backward velocity is not tangent.

    This is not the published CPD hard-rod algorithm (none is claimed).
    Finite compliance is not substituted for length. The row is solved
    jointly with relaxed energy, and failures/rank loss remain reported.
    """
    return _step(x, v, h, m, 0., l, g, epsilon, max_iterations, tol, relax, True)


def _self_test():
    relaxed = cpd_spring_step([1.3, 0], [0, 0], .02, 2., 100., 1.)
    strict = cpd_spring_step([1.3, 0], [0, 0], .02, 2., 100., 1., relax=False)
    assert relaxed.converged and strict.converged
    assert abs(relaxed.rho-4.5e-8) < 1e-13
    assert abs(relaxed.v[0]+3.375e-7) < 1e-10
    assert np.array_equal(strict.x, [1.3, 0.]) and np.array_equal(strict.v, [0., 0.])
    assert abs(relaxed.energy_balance) < 1e-10
    q = .65
    x = np.array([np.sin(q), -np.cos(q)])
    v = 1.2*np.array([np.cos(q), np.sin(q)])
    rod = cpd_pendulum_step(x, v, .08, 1., 1., 9.81)
    if rod.converged and np.linalg.norm(rod.x-x) > 1e-5:
        assert rod.tangency > 1e-5
    return dict(release_relaxed=dict(x=relaxed.x.tolist(), v=relaxed.v.tolist(),
                                    rho=relaxed.rho, residual=relaxed.residual,
                                    energy_error=relaxed.energy_error, iterations=relaxed.iterations),
                release_strict=dict(x=strict.x.tolist(), v=strict.v.tolist(), residual=strict.residual),
                hard_rod_adaptation=dict(converged=rod.converged, residual=rod.residual,
                                        geometry=rod.geometry, tangency=rod.tangency,
                                        energy_error=rod.energy_error))


if __name__ == "__main__":
    import json
    print(json.dumps(_self_test(), indent=2))
