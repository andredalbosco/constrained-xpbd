"""Refreshed overlapping block GS for the agreed position-only CPD objective.

lambda is unscaled, with C + alpha*lambda = 0. Energy is physical kinetic
plus elastic plus gravity potential. Every block retains position stationarity.
No correction cap, energy repair, line search, or fallback is applied.
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
        return 0.0 if self.rigid else 1.0 / self.k


def scenes():
    def spring(name, x, v, k=100., g=(0., 0.)):
        return Scene(name, 1., k, 1., np.array(g, float), np.array(x, float), np.array(v, float))
    def rod(name, theta, speed):
        return Scene(name, 1., 0., 1., np.array([0., -9.81]),
                     np.array([np.sin(theta), -np.cos(theta)]),
                     speed*np.array([np.cos(theta), np.sin(theta)]), True)
    return [spring('spring_release', [1.25, 0.], [0., 0.]),
            spring('spring_orbit', [1.1, 0.], [0., .8]),
            spring('spring_gravity', [.55, -1.1], [.4, .2], g=(0., -9.81)),
            rod('pendulum_release', .65, 0.), rod('pendulum_kick', .4, 1.2),
            spring('spring_stiff', [1.1, 0.], [0., .8], k=1000.)]


def potential(scene, x):
    stretch = float(np.linalg.norm(x)) - scene.l
    elastic = 0.0 if scene.rigid else .5*scene.k*stretch*stretch
    return elastic - scene.m*float(scene.g @ x)


def energy(scene, x, v):
    return .5*scene.m*float(v @ v) + potential(scene, x)


def energy_floor(scene):
    force = scene.m*float(np.linalg.norm(scene.g))
    return -force*scene.l - (0.0 if scene.rigid else force*force/(2*scene.k))


@dataclass
class Evaluation:
    constraint: float
    normal: np.ndarray
    potential: float
    potential_gradient: np.ndarray
    energy_gradient: np.ndarray
    r_x: np.ndarray
    r_lambda: float
    r_rho: float
    r_eta: float
    physical_energy: float
    target: float
    kstar: float
    alpha: float


def evaluate_state(scene, x0, v0, h, x, lam, rho, eta, epsilon, target=None, kstar=None):
    radius = float(np.linalg.norm(x))
    if radius <= 1e-14*scene.l:
        raise FloatingPointError('undefined_distance_normal')
    if not np.all(np.isfinite(x)) or not all(math.isfinite(z) for z in (lam, rho, eta)):
        raise FloatingPointError('nonfinite_iterate')
    if target is None:
        target = energy(scene, x0, v0)
    if kstar is None:
        kstar = energy_floor(scene)
    normal = x/radius
    constraint = radius-scene.l
    elastic_gradient = np.zeros_like(x) if scene.rigid else scene.k*constraint*normal
    potential_gradient = elastic_gradient-scene.m*scene.g
    dx = x-x0
    energy_gradient = scene.m*dx/(h*h)+potential_gradient
    current_potential = potential(scene, x)
    physical_energy = .5*scene.m*float(dx @ dx)/(h*h)+current_potential
    # Gravity in the objective remains outside the spring-coefficient substitution.
    r_x = scene.m*(dx-h*v0)/(h*h)-lam*normal-scene.m*scene.g+eta*energy_gradient
    return Evaluation(constraint, normal, current_potential, potential_gradient,
                      energy_gradient, r_x, constraint+scene.alpha*lam,
                      epsilon*rho+eta*(target-kstar),
                      physical_energy-(1-rho)*target-rho*kstar,
                      physical_energy, float(target), float(kstar), scene.alpha)


@dataclass
class Increment:
    dx: np.ndarray
    dlam: float = 0.0
    drho: float = 0.0
    deta: float = 0.0
    denominator: float = math.nan
    compatible_rank_defect: bool = False


def material_increment(scene, h, evaluation):
    inverse = h*h/scene.m
    denominator = evaluation.alpha+inverse*float(evaluation.normal @ evaluation.normal)
    dlam = (-evaluation.r_lambda+inverse*float(evaluation.normal @ evaluation.r_x))/denominator
    dx = -inverse*evaluation.r_x+inverse*evaluation.normal*dlam
    return Increment(dx, dlam=dlam, denominator=denominator)


def energy_increment(scene, h, epsilon, evaluation):
    inverse = h*h/scene.m
    gap = evaluation.target-evaluation.kstar
    gradient = evaluation.energy_gradient
    denominator = inverse*float(gradient @ gradient)+gap*gap/epsilon
    numerator = (evaluation.r_eta-inverse*float(gradient @ evaluation.r_x)
                 -gap*evaluation.r_rho/epsilon)
    if denominator == 0.0:
        # The linearized energy row contains no unknowns. Choose the zero
        # multiplier increment only if that row is exactly compatible.
        if numerator != 0.0:
            raise FloatingPointError('incompatible_zero_energy_row')
        return Increment(-inverse*evaluation.r_x, drho=-evaluation.r_rho/epsilon,
                         denominator=denominator, compatible_rank_defect=True)
    deta = numerator/denominator
    dx = -inverse*evaluation.r_x-inverse*gradient*deta
    drho = -(evaluation.r_rho+gap*deta)/epsilon
    return Increment(dx, drho=drho, deta=deta, denominator=denominator)


def joint_increment(scene, h, epsilon, evaluation):
    """Independent whole-system solve of the same approximate KKT matrix."""
    n = len(evaluation.normal)
    matrix = np.zeros((n+3, n+3))
    matrix[:n, :n] = (scene.m/(h*h))*np.eye(n)
    matrix[:n, n] = -evaluation.normal
    matrix[n, :n] = evaluation.normal
    matrix[n, n] = evaluation.alpha
    matrix[:n, n+2] = evaluation.energy_gradient
    matrix[n+2, :n] = evaluation.energy_gradient
    matrix[n+1, n+1] = epsilon
    matrix[n+1, n+2] = matrix[n+2, n+1] = evaluation.target-evaluation.kstar
    rhs = -np.r_[evaluation.r_x, evaluation.r_lambda, evaluation.r_rho, evaluation.r_eta]
    try:
        change = np.linalg.solve(matrix, rhs)
    except np.linalg.LinAlgError as exc:
        raise FloatingPointError('singular_joint_system') from exc
    return Increment(change[:n], change[n], change[n+1], change[n+2])


def scaled_residual(scene, h, epsilon, evaluation):
    energy_scale = max(abs(evaluation.target-evaluation.kstar),
                       scene.k*scene.l**2 if not scene.rigid else scene.m*np.linalg.norm(scene.g)*scene.l,
                       1e-12)
    return max(float(np.linalg.norm(evaluation.r_x))*h*h/(scene.m*scene.l),
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
    max_block_displacement: float = 0.0
    rank_defect_visits: int = 0


def block_step(scene, x0, v0, h, epsilon, order='EL', sweeps=4, tol=1e-9,
               stop_on_convergence=False, trace=False):
    if order not in ('EL', 'LE', 'JOINT', 'XPBD'):
        raise ValueError('order must be EL, LE, JOINT, or XPBD')
    if h <= 0 or epsilon <= 0 or sweeps < 1:
        raise ValueError('h, epsilon and sweeps must be positive')
    x0, v0 = np.asarray(x0, float), np.asarray(v0, float)
    x = (x0+h*v0).copy()
    lam = rho = eta = 0.0
    target, kstar = energy(scene, x0, v0), energy_floor(scene)
    history, visits, maxjump, rank_visits = [], 0, 0., 0
    evaluation = None
    blocks = ('J',) if order == 'JOINT' else (('L',) if order == 'XPBD' else tuple(order))
    try:
        evaluation = evaluate_state(scene, x0, v0, h, x, lam, rho, eta, epsilon, target, kstar)
        for sweep in range(sweeps):
            for block in blocks:
                if block == 'L':
                    change = material_increment(scene, h, evaluation)
                elif block == 'E':
                    change = energy_increment(scene, h, epsilon, evaluation)
                else:
                    change = joint_increment(scene, h, epsilon, evaluation)
                x = x+change.dx
                lam += change.dlam
                rho += change.drho
                eta += change.deta
                visits += 1
                rank_visits += int(change.compatible_rank_defect)
                maxjump = max(maxjump, float(np.linalg.norm(change.dx)))
                evaluation = evaluate_state(scene, x0, v0, h, x, lam, rho, eta, epsilon, target, kstar)
                residual = scaled_residual(scene, h, epsilon, evaluation)
                if trace:
                    history.append(dict(sweep=sweep, block=block, x=x.tolist(), lam=lam,
                                        rho=rho, eta=eta, residual=residual,
                                        physical_energy=evaluation.physical_energy,
                                        r_lambda=evaluation.r_lambda, r_eta=evaluation.r_eta,
                                        r_rho=evaluation.r_rho, r_x=evaluation.r_x.tolist(),
                                        displacement=float(np.linalg.norm(change.dx))))
                velocity = (x-x0)/h
                if (not np.all(np.isfinite(velocity)) or np.linalg.norm(x) > 1e4*scene.l
                        or np.linalg.norm(velocity) > 1e6*scene.l):
                    raise FloatingPointError('runaway_guard')
            if stop_on_convergence and order != 'XPBD' and residual <= tol:
                break
        if order == 'XPBD':
            converged = max(np.linalg.norm(evaluation.r_x)*h*h/(scene.m*scene.l),
                            abs(evaluation.r_lambda)/scene.l) <= tol
        else:
            converged = scaled_residual(scene, h, epsilon, evaluation) <= tol
        return StepResult(x, (x-x0)/h, lam, rho, eta, True, bool(converged),
                          '' if converged else 'fixed_budget_residual', evaluation, visits,
                          history, maxjump, rank_visits)
    except (FloatingPointError, OverflowError, ValueError) as exc:
        return StepResult(x, (x-x0)/h, lam, rho, eta, False, False, str(exc),
                          evaluation, visits, history, maxjump, rank_visits)
