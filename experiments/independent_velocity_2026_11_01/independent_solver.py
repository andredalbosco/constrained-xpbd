"""Independent-velocity energy-constrained implicit objective, refreshed block GS.

The velocity coupling is primal and beta_v must be positive. Full coupling
curvature and both primal residuals are retained. Material and energy multiplier
curvature are omitted. No velocity overwrite, energy rescale, or step damping.
"""
from __future__ import annotations

from dataclasses import dataclass, field
import math
import numpy as np


@dataclass
class Scene:
    name: str
    m: float
    k: float
    l: float
    g: np.ndarray
    x0: np.ndarray
    v0: np.ndarray
    rigid: bool = False

    @property
    def alpha(self):
        return 0. if self.rigid else 1. / self.k


def scenes():
    def spring(name, x, v, k=100., gravity=(0., 0.)):
        return Scene(name, 1., k, 1., np.array(gravity), np.array(x), np.array(v))

    def rod(name, angle, speed):
        return Scene(name, 1., 0., 1., np.array([0., -9.81]),
                     np.array([np.sin(angle), -np.cos(angle)]),
                     speed*np.array([np.cos(angle), np.sin(angle)]), True)

    return [spring('spring_release', [1.25, 0.], [0., 0.]),
            spring('spring_orbit', [1.1, 0.], [0., .8]),
            spring('spring_gravity', [.55, -1.1], [.4, .2], gravity=(0., -9.81)),
            rod('pendulum_release', .65, 0.), rod('pendulum_kick', .4, 1.2),
            spring('spring_stiff', [1.1, 0.], [0., .8], k=1000.)]


def potential(scene, x):
    elastic = 0. if scene.rigid else .5*scene.k*(np.linalg.norm(x)-scene.l)**2
    return float(elastic-scene.m*np.dot(scene.g, x))


def energy(scene, x, v):
    return .5*scene.m*float(np.dot(v, v))+potential(scene, x)


def energy_floor(scene):
    force = scene.m*float(np.linalg.norm(scene.g))
    return -force*scene.l-(0. if scene.rigid else force*force/(2*scene.k))


@dataclass
class Evaluation:
    constraint: float
    normal: np.ndarray
    potential: float
    potential_gradient: np.ndarray
    energy_gradient_x: np.ndarray
    energy_gradient_v: np.ndarray
    coupling_residual: np.ndarray
    r_x: np.ndarray
    r_v: np.ndarray
    r_lambda: float
    r_rho: float
    r_eta: float
    physical_energy: float
    target: float
    kstar: float
    alpha: float


def evaluate_state(scene, x0, v0, h, x, v, lam, rho, eta, epsilon, beta_v,
                   target=None, kstar=None):
    if h <= 0 or epsilon <= 0 or beta_v <= 0:
        raise ValueError('h, epsilon and primal beta_v must be positive')
    radius = float(np.linalg.norm(x))
    if radius <= 1e-14*scene.l:
        raise FloatingPointError('undefined_distance_normal')
    if (not np.all(np.isfinite(x)) or not np.all(np.isfinite(v))
            or not all(math.isfinite(z) for z in (lam, rho, eta))):
        raise FloatingPointError('nonfinite_iterate')
    target = energy(scene, x0, v0) if target is None else target
    kstar = energy_floor(scene) if kstar is None else kstar
    normal = x/radius
    constraint = radius-scene.l
    gradient = ((np.zeros_like(x) if scene.rigid else scene.k*constraint*normal)
                -scene.m*scene.g)
    difference = v-(x-x0)/h
    r_x = (scene.m*(x-x0-h*v0)/(h*h)-lam*normal-scene.m*scene.g
           -scene.m*difference/(beta_v*h)+eta*gradient)
    r_v = scene.m*difference/beta_v+eta*scene.m*v
    current_energy = energy(scene, x, v)
    return Evaluation(constraint, normal, potential(scene, x), gradient, gradient,
                      scene.m*v, difference, r_x, r_v,
                      constraint+scene.alpha*lam,
                      epsilon*rho+eta*(target-kstar),
                      current_energy-(1-rho)*target-rho*kstar,
                      current_energy, float(target), float(kstar), scene.alpha)


def primal_matrix(scene, h, beta_v, n=None):
    n = len(scene.x0) if n is None else n
    identity = np.eye(n)*scene.m
    return np.block([[(1+1/beta_v)*identity/(h*h), -identity/(beta_v*h)],
                     [-identity/(beta_v*h), identity/beta_v]])


def inverse_primal(scene, h, beta_v, rx, rv):
    return ((h*h*rx+h*rv)/scene.m,
            (h*rx+(1+beta_v)*rv)/scene.m)


@dataclass
class Increment:
    dx: np.ndarray
    dv: np.ndarray
    dlam: float = 0.
    drho: float = 0.
    deta: float = 0.
    denominator: float = math.nan
    compatible_rank_defect: bool = False


def material_increment(scene, h, beta_v, evaluation):
    ux, uv = inverse_primal(scene, h, beta_v, evaluation.r_x, evaluation.r_v)
    denominator = evaluation.alpha+h*h/scene.m*float(evaluation.normal @ evaluation.normal)
    dlam = (-evaluation.r_lambda+float(evaluation.normal @ ux))/denominator
    return Increment(-ux+h*h/scene.m*evaluation.normal*dlam,
                     -uv+h/scene.m*evaluation.normal*dlam,
                     dlam=dlam, denominator=denominator)


def energy_increment(scene, h, epsilon, beta_v, evaluation):
    ux, uv = inverse_primal(scene, h, beta_v, evaluation.r_x, evaluation.r_v)
    sx, sv = inverse_primal(scene, h, beta_v,
                           evaluation.energy_gradient_x, evaluation.energy_gradient_v)
    gap = evaluation.target-evaluation.kstar
    combined = h*evaluation.energy_gradient_x+evaluation.energy_gradient_v
    # Algebraically identical sum of squares avoids cancellation between
    # the position/velocity cross terms in the positive Schur quadratic.
    denominator = (float(combined @ combined)/scene.m
                   +beta_v*float(evaluation.energy_gradient_v @ evaluation.energy_gradient_v)/scene.m
                   +gap*gap/epsilon)
    numerator = (evaluation.r_eta-float(evaluation.energy_gradient_x @ ux)
                 -float(evaluation.energy_gradient_v @ uv)-gap*evaluation.r_rho/epsilon)
    if denominator == 0.:
        if numerator != 0.:
            raise FloatingPointError('incompatible_zero_energy_row')
        return Increment(-ux, -uv, drho=-evaluation.r_rho/epsilon,
                         denominator=0., compatible_rank_defect=True)
    if denominator < 0. or not math.isfinite(denominator):
        raise FloatingPointError('invalid_energy_schur_denominator')
    deta = numerator/denominator
    return Increment(-ux-sx*deta, -uv-sv*deta,
                     drho=-(evaluation.r_rho+gap*deta)/epsilon,
                     deta=deta, denominator=denominator)


def joint_increment(scene, h, epsilon, beta_v, evaluation):
    n = len(evaluation.normal)
    matrix = np.zeros((2*n+3, 2*n+3))
    matrix[:2*n, :2*n] = primal_matrix(scene, h, beta_v, n)
    matrix[:n, 2*n] = -evaluation.normal
    matrix[2*n, :n] = evaluation.normal
    matrix[2*n, 2*n] = evaluation.alpha
    gradient = np.r_[evaluation.energy_gradient_x, evaluation.energy_gradient_v]
    matrix[:2*n, 2*n+2] = matrix[2*n+2, :2*n] = gradient
    matrix[2*n+1, 2*n+1] = epsilon
    matrix[2*n+1, 2*n+2] = matrix[2*n+2, 2*n+1] = evaluation.target-evaluation.kstar
    rhs = -np.r_[evaluation.r_x, evaluation.r_v, evaluation.r_lambda,
                 evaluation.r_rho, evaluation.r_eta]
    try:
        change = np.linalg.solve(matrix, rhs)
    except np.linalg.LinAlgError as exc:
        raise FloatingPointError('singular_joint_system') from exc
    return Increment(change[:n], change[n:2*n], change[2*n],
                     change[2*n+1], change[2*n+2])


def scaled_residual(scene, h, epsilon, beta_v, evaluation):
    energy_scale = max(abs(evaluation.target-evaluation.kstar),
                       scene.k*scene.l**2 if not scene.rigid else scene.m*np.linalg.norm(scene.g)*scene.l,
                       1e-12)
    return max(float(np.linalg.norm(evaluation.r_x))*h*h/(scene.m*scene.l),
               float(np.linalg.norm(evaluation.r_v))*h/(scene.m*scene.l),
               abs(evaluation.r_lambda)/scene.l,
               abs(evaluation.r_rho)/max(epsilon, energy_scale),
               abs(evaluation.r_eta)/energy_scale)


@dataclass
class StepResult:
    x: np.ndarray
    v: np.ndarray
    lam: float
    rho: float
    eta: float
    accepted: bool
    converged: bool
    reason: str
    evaluation: Evaluation | None
    visits: int
    history: list = field(default_factory=list)
    max_block_displacement: float = 0.
    max_velocity_increment: float = 0.
    rank_defect_visits: int = 0


def block_step(scene, x0, v0, h, epsilon, beta_v=1., order='EL', sweeps=4,
               tol=1e-9, stop_on_convergence=False, trace=False):
    if order not in ('EL', 'LE', 'JOINT'):
        raise ValueError('order must be EL, LE, or JOINT')
    if h <= 0 or epsilon <= 0 or beta_v <= 0 or sweeps < 1:
        raise ValueError('h, epsilon, beta_v and sweeps must be positive')
    x0, v0 = np.asarray(x0, float), np.asarray(v0, float)
    x, v = (x0+h*v0).copy(), v0.copy()
    lam = rho = eta = 0.
    target, kstar = energy(scene, x0, v0), energy_floor(scene)
    history, visits, maxjump, maxdv, rank_visits = [], 0, 0., 0., 0
    evaluation = None
    blocks = ('J',) if order == 'JOINT' else tuple(order)
    try:
        evaluation = evaluate_state(scene, x0, v0, h, x, v, lam, rho, eta,
                                    epsilon, beta_v, target, kstar)
        for sweep in range(sweeps):
            for block in blocks:
                if block == 'L':
                    change = material_increment(scene, h, beta_v, evaluation)
                elif block == 'E':
                    change = energy_increment(scene, h, epsilon, beta_v, evaluation)
                else:
                    change = joint_increment(scene, h, epsilon, beta_v, evaluation)
                x, v = x+change.dx, v+change.dv
                lam += change.dlam
                rho += change.drho
                eta += change.deta
                visits += 1
                rank_visits += int(change.compatible_rank_defect)
                maxjump = max(maxjump, float(np.linalg.norm(change.dx)))
                maxdv = max(maxdv, float(np.linalg.norm(change.dv)))
                evaluation = evaluate_state(scene, x0, v0, h, x, v, lam, rho, eta,
                                            epsilon, beta_v, target, kstar)
                residual = scaled_residual(scene, h, epsilon, beta_v, evaluation)
                if trace:
                    history.append(dict(sweep=sweep, block=block, x=x.tolist(), v=v.tolist(),
                                        lam=lam, rho=rho, eta=eta, residual=residual,
                                        physical_energy=evaluation.physical_energy,
                                        r_lambda=evaluation.r_lambda, r_eta=evaluation.r_eta,
                                        r_rho=evaluation.r_rho, r_x=evaluation.r_x.tolist(),
                                        r_v=evaluation.r_v.tolist(),
                                        coupling_residual=evaluation.coupling_residual.tolist(),
                                        displacement=float(np.linalg.norm(change.dx)),
                                        velocity_increment=float(np.linalg.norm(change.dv))))
                if (np.linalg.norm(x) > 1e4*scene.l or np.linalg.norm(v) > 1e6*scene.l
                        or np.linalg.norm((x-x0)/h) > 1e6*scene.l):
                    raise FloatingPointError('runaway_guard')
            if stop_on_convergence and residual <= tol:
                break
        converged = scaled_residual(scene, h, epsilon, beta_v, evaluation) <= tol
        return StepResult(x, v, lam, rho, eta, True, bool(converged),
                          '' if converged else 'fixed_budget_residual', evaluation, visits,
                          history, maxjump, maxdv, rank_visits)
    except (FloatingPointError, OverflowError, ValueError) as exc:
        return StepResult(x, v, lam, rho, eta, False, False, str(exc), evaluation,
                          visits, history, maxjump, maxdv, rank_visits)
