"""Planar paper-equation CPD, plus explicitly labeled anchored/rod extensions.

PD local projections, constant PD metric, and fresh Eq.21 multipliers.
External work uses the paper's incoming-velocity approximation. No line search,
velocity projection, energy rescale, substep, or least-squares row repair.
"""
from __future__ import annotations

from dataclasses import dataclass, field
import numpy as np


@dataclass
class Edge:
    a: int
    b: int
    rest: float
    stiffness: float
    hard: bool = False


@dataclass
class Scene:
    name: str
    masses: np.ndarray
    x0: np.ndarray
    v0: np.ndarray
    forces: np.ndarray
    edges: list[Edge]
    family: str
    energy_floor: float = 0.
    characteristic_length: float = 1.
    reference_case: str = ''


def rotate(a):
    return np.stack((-a[..., 1], a[..., 0]), axis=-1)


def cross(a, b):
    return a[..., 0]*b[..., 1]-a[..., 1]*b[..., 0]


def scenes():
    result = []
    configs = [('spring_release', [1.25, 0.], [0., 0.], 100., [0., 0.]),
               ('spring_orbit', [1.1, 0.], [0., .8], 100., [0., 0.]),
               ('spring_gravity', [.55, -1.1], [.4, .2], 100., [0., -9.81]),
               ('spring_stiff', [1.1, 0.], [0., .8], 1000., [0., 0.])]
    for name, position, velocity, stiffness, force in configs:
        x, v, f = map(lambda z: np.array(z, float), (position, velocity, force))
        magnitude = float(np.linalg.norm(f))
        floor = -magnitude-magnitude*magnitude/(2*stiffness)
        result.append(Scene(name, np.ones(1), x[None, :], v[None, :], f[None, :],
                            [Edge(0, -1, 1., stiffness)], 'anchored', floor,
                            reference_case=name))
        # Reduced mass is one, matching the anchored unit-mass relative motion.
        free_name = {'spring_release':'free_release', 'spring_orbit':'free_orbit',
                     'spring_gravity':'free_opposite_load', 'spring_stiff':'free_stiff'}[name]
        result.append(Scene(free_name, np.full(2, 2.), np.array([x/2, -x/2]),
                            np.array([v/2, -v/2]), np.array([f, -f]),
                            [Edge(0, 1, 1., stiffness)], 'free', floor,
                            reference_case=name))
    triangle = 1.15*np.array([[-.5, -np.sqrt(3)/6], [.5, -np.sqrt(3)/6], [0., np.sqrt(3)/3]])
    result.append(Scene('free_triangle_spin', np.ones(3), triangle, .6*rotate(triangle),
                        np.zeros_like(triangle), [Edge(0,1,1.,100.),Edge(1,2,1.,100.),Edge(2,0,1.,100.)],
                        'free'))
    for name, angle, speed in [('pendulum_release', .65, 0.), ('pendulum_kick', .4, 1.2)]:
        x = np.array([[np.sin(angle), -np.cos(angle)]])
        v = speed*np.array([[np.cos(angle), np.sin(angle)]])
        for hard in (False, True):
            stiffness = 0. if hard else 10000.
            floor = -9.81-(0. if hard else 9.81**2/(2*stiffness))
            result.append(Scene(name+('_hard' if hard else '_compliant'), np.ones(1),
                                x.copy(), v.copy(), np.array([[0.,-9.81]]),
                                [Edge(0,-1,1.,stiffness,hard)],
                                'hard_rod_extension' if hard else 'anchored', floor,
                                reference_case=name))
    return result


def edge_difference(scene, x, edge):
    return x[edge.a]-(x[edge.b] if edge.b >= 0 else 0.)


@dataclass
class Elastic:
    energy: float
    gradient: np.ndarray
    projections: list


def elastic(scene, x):
    gradient = np.zeros_like(x)
    energy = 0.
    projections = []
    for edge in scene.edges:
        difference = edge_difference(scene, x, edge)
        radius = float(np.linalg.norm(difference))
        if not np.isfinite(radius) or radius < 1e-14*edge.rest:
            raise FloatingPointError('undefined_distance_normal')
        projection = edge.rest*difference/radius
        projections.append(projection)
        if not edge.hard:
            energy += .5*edge.stiffness*(radius-edge.rest)**2
            force_gradient = edge.stiffness*(difference-projection)
            gradient[edge.a] += force_gradient
            if edge.b >= 0:
                gradient[edge.b] -= force_gradient
    return Elastic(float(energy), gradient, projections)


def kinetic_energy(scene, v):
    return .5*float(np.sum(scene.masses[:,None]*v*v))


def internal_energy(scene, x, v):
    return kinetic_energy(scene,v)+elastic(scene,x).energy


def mechanical_energy(scene, x, v):
    return internal_energy(scene,x,v)-float(np.sum(scene.forces*x))


def pd_matrix(scene, h):
    matrix = np.diag(np.repeat(scene.masses/(h*h), 2))
    for edge in scene.edges:
        if edge.hard:
            continue
        for axis in (0,1):
            a = 2*edge.a+axis
            matrix[a,a] += edge.stiffness
            if edge.b >= 0:
                b = 2*edge.b+axis
                matrix[b,b] += edge.stiffness
                matrix[a,b] -= edge.stiffness
                matrix[b,a] -= edge.stiffness
    return matrix


@dataclass
class Targets:
    origin: np.ndarray
    linear: np.ndarray
    angular: float
    energy: float
    kstar: float
    inertia: float
    linear_active: bool
    angular_active: bool


def make_targets(scene, xn, vn, h, angular=True):
    linear_active = scene.family == 'free'
    origin = (np.sum(scene.masses[:,None]*xn,axis=0)/np.sum(scene.masses)
              if linear_active else np.zeros(2))
    impulse = scene.masses[:,None]*vn+h*scene.forces
    linear = np.sum(impulse,axis=0)
    # The torque order follows L=x cross mv; the printed Eq.12 reverses it.
    angular_target = float(np.sum(cross(xn-origin,impulse)))
    target_energy = internal_energy(scene,xn,vn)+h*float(np.sum(scene.forces*vn))
    inertia = float(np.sum(scene.masses*np.sum((xn-origin)**2,axis=1)))
    if inertia <= 0:
        raise FloatingPointError('undefined_incoming_inertia')
    # In the angular-row ablation Kstar is deliberately held unchanged.
    kstar = angular_target**2/(2*inertia)
    if linear_active:
        kstar += float(linear@linear)/(2*np.sum(scene.masses))
    return Targets(origin,linear,angular_target,float(target_energy),float(kstar),
                   inertia,linear_active,angular)


@dataclass
class Evaluation:
    objective_gradient: np.ndarray
    constraint_gradients: np.ndarray
    constraints: np.ndarray
    names: list[str]
    internal_energy: float
    physical_energy: float
    elastic_energy: float
    velocity: np.ndarray
    linear: np.ndarray
    angular: float
    projections: list


def evaluate(scene, xn, vn, h, x, rho, epsilon, targets):
    if not np.all(np.isfinite(x)) or not np.isfinite(rho):
        raise FloatingPointError('nonfinite_iterate')
    n = x.size
    elastic_state = elastic(scene,x)
    predictor = xn+h*vn+h*h*scene.forces/scene.masses[:,None]
    gradient = np.r_[((scene.masses[:,None]*(x-predictor)/(h*h)+elastic_state.gradient).ravel()),epsilon*rho]
    velocity = (x-xn)/h
    linear = np.sum(scene.masses[:,None]*velocity,axis=0)
    angular = float(np.sum(cross(x-targets.origin,scene.masses[:,None]*velocity)))
    internal = kinetic_energy(scene,velocity)+elastic_state.energy
    rows, values, names = [], [], []
    if targets.linear_active:
        for axis in (0,1):
            row = np.zeros(n+1)
            row[axis:n:2] = scene.masses/h
            rows.append(row); values.append(linear[axis]-targets.linear[axis]); names.append('linear_'+str(axis))
    if targets.angular_active:
        row = np.r_[(scene.masses[:,None]*rotate(xn-targets.origin)/h).ravel(),0.]
        rows.append(row); values.append(angular-targets.angular); names.append('angular')
    for index,edge in enumerate(scene.edges):
        if edge.hard:
            delta = edge_difference(scene,x,edge)
            length = float(np.linalg.norm(delta))
            row = np.zeros(n+1)
            row[2*edge.a:2*edge.a+2] = delta/length
            if edge.b >= 0:
                row[2*edge.b:2*edge.b+2] -= delta/length
            rows.append(row); values.append(length-edge.rest); names.append('geometry_'+str(index))
    row = np.r_[(scene.masses[:,None]*(x-xn)/(h*h)+elastic_state.gradient).ravel(),targets.energy-targets.kstar]
    rows.append(row); values.append(internal-(1-rho)*targets.energy-rho*targets.kstar); names.append('energy')
    return Evaluation(gradient,np.array(rows),np.array(values),names,float(internal),
                      float(internal-np.sum(scene.forces*x)),elastic_state.energy,
                      velocity,linear,angular,elastic_state.projections)


def residual_measures(scene,h,epsilon,evaluation,multipliers,targets):
    station = evaluation.objective_gradient+evaluation.constraint_gradients.T@multipliers
    length = scene.characteristic_length
    scales = []
    for name in evaluation.names:
        if name.startswith('linear'):
            scales.append(float(np.sum(scene.masses))*length/h)
        elif name == 'angular':
            scales.append(float(np.sum(scene.masses))*length*length/h)
        elif name.startswith('geometry'):
            scales.append(length)
        else:
            scales.append(max(abs(targets.energy),abs(targets.kstar),1e-12))
    position_station = np.max(np.abs(station[:-1].reshape(-1,2))*h*h/(scene.masses[:,None]*length))
    rho_station = abs(station[-1])/max(epsilon,abs(targets.energy),abs(targets.kstar),1e-12)
    constraints_scaled = float(np.max(np.abs(evaluation.constraints)/scales))
    return dict(raw_constraints=float(np.max(np.abs(evaluation.constraints))),
                scaled_constraints=constraints_scaled,
                stationarity=float(max(position_station,rho_station)),
                kkt_scaled=float(max(position_station,rho_station,constraints_scaled)),
                r_x=station[:-1].reshape(-1,2),r_rho=float(station[-1]),
                r_eta=float(evaluation.constraints[-1]))


@dataclass
class Increment:
    dx: np.ndarray
    drho: float
    multipliers: np.ndarray
    rank: int
    condition: float
    linear_residual: float


def increment(scene,h,epsilon,evaluation,inverse=None):
    n = evaluation.objective_gradient.size-1
    if inverse is None:
        inverse = np.zeros((n+1,n+1))
        inverse[:n,:n] = np.linalg.inv(pd_matrix(scene,h))
        inverse[-1,-1] = 1/epsilon
    inverse_gradient = inverse@evaluation.objective_gradient
    inverse_rows = inverse@evaluation.constraint_gradients.T
    schur = evaluation.constraint_gradients@inverse_rows
    schur = .5*(schur+schur.T)
    scales = np.sqrt(np.maximum(np.diag(schur),0.))
    scaling = np.divide(1.,scales,out=np.ones_like(scales),where=scales>0)
    normalized = scaling[:,None]*schur*scaling[None,:]
    eigenvalues = np.linalg.eigvalsh(normalized)
    cutoff = 1e-12*max(float(np.max(eigenvalues)),1.)
    rank = int(np.count_nonzero(eigenvalues>cutoff))
    condition = (float(np.max(eigenvalues)/np.min(eigenvalues))
                 if np.min(eigenvalues)>0 else float('inf'))
    rhs = evaluation.constraints-evaluation.constraint_gradients@inverse_gradient
    try:
        # A rank estimate is diagnostic. Attempt the stated direct inverse solve;
        # dependent rows are never dropped or repaired by least squares.
        multipliers = scaling*np.linalg.solve(normalized,scaling*rhs)
    except np.linalg.LinAlgError as exc:
        raise FloatingPointError('singular_conservation_schur') from exc
    change = -inverse_gradient-inverse_rows@multipliers
    frozen = np.r_[evaluation.objective_gradient+np.linalg.solve(inverse,change)
                    +evaluation.constraint_gradients.T@multipliers,
                    evaluation.constraints+evaluation.constraint_gradients@change]
    linear_residual = float(np.max(np.abs(frozen))/max(1.,np.max(np.abs(rhs)),np.max(np.abs(evaluation.objective_gradient))))
    if not np.all(np.isfinite(change)) or not np.all(np.isfinite(multipliers)):
        raise FloatingPointError('nonfinite_subproblem_solution')
    return Increment(change[:n].reshape(scene.x0.shape),float(change[-1]),multipliers,
                     rank,condition,linear_residual)


@dataclass
class StepResult:
    x: np.ndarray
    v: np.ndarray
    rho: float
    multipliers: np.ndarray
    accepted: bool
    converged: bool
    stationary: bool
    reason: str
    iterations: int
    evaluation: Evaluation|None
    measures: dict
    targets: Targets
    rank_min: int
    condition_max: float
    max_block_displacement: float
    history: list = field(default_factory=list)
    max_linear_residual: float = 0.


def step(scene,xn,vn,h,epsilon,policy='paper',iterations=16,max_iterations=200,
         tol=1e-4,angular=True,initial_precheck=False,trace=False,method='CPD'):
    if h<=0 or epsilon<=0 or iterations<1 or max_iterations<1:
        raise ValueError('h, epsilon and iteration limits must be positive')
    if policy not in ('paper','fixed','stationary') or method not in ('CPD','PD'):
        raise ValueError('invalid method or stopping policy')
    if method=='PD' and any(edge.hard for edge in scene.edges):
        raise ValueError('PD hard rod is outside the finite-spring model')
    xn,vn = np.asarray(xn,float),np.asarray(vn,float)
    x = xn+h*vn+h*h*scene.forces/scene.masses[:,None]
    rho = 0.
    targets = make_targets(scene,xn,vn,h,angular)
    evaluation = None
    history = []
    count = 0
    maxjump = condition_max = max_linear_residual = 0.
    rank_min = 1000
    multipliers = np.zeros(0)
    measures = {}
    limit = iterations if policy=='fixed' else max_iterations
    try:
        evaluation = evaluate(scene,xn,vn,h,x,rho,epsilon,targets)
        multipliers = np.zeros(len(evaluation.constraints))
        measures = residual_measures(scene,h,epsilon,evaluation,multipliers,targets)
        if method=='CPD' and initial_precheck and measures['raw_constraints']<=tol:
            return StepResult(x,(x-xn)/h,rho,multipliers,True,True,
                              measures['kkt_scaled']<=1e-9,'initial_constraint_precheck',0,
                              evaluation,measures,targets,0,0.,0.,history)
        inverse = np.zeros((x.size+1,x.size+1))
        inverse[:-1,:-1] = np.linalg.inv(pd_matrix(scene,h))
        inverse[-1,-1] = 1/epsilon
        for visit in range(limit):
            if method=='CPD':
                update = increment(scene,h,epsilon,evaluation,inverse)
                x = x+update.dx
                rho += update.drho
                # Eq.21 solves a NEW multiplier, rather than an increment.
                multipliers = update.multipliers
                rank_min = min(rank_min,update.rank)
                condition_max = max(condition_max,update.condition)
                max_linear_residual = max(max_linear_residual,update.linear_residual)
                maxjump = max(maxjump,float(np.linalg.norm(update.dx)))
            else:
                delta = -inverse[:-1,:-1]@evaluation.objective_gradient[:-1]
                x = x+delta.reshape(x.shape)
                maxjump = max(maxjump,float(np.linalg.norm(delta)))
            count = visit+1
            evaluation = evaluate(scene,xn,vn,h,x,rho,epsilon,targets)
            measures = residual_measures(scene,h,epsilon,evaluation,multipliers,targets)
            if trace:
                history.append(dict(iteration=count,x=x.tolist(),v=evaluation.velocity.tolist(),
                                    rho=rho,multipliers=multipliers.tolist(),
                                    constraint_names=evaluation.names,
                                    constraints=evaluation.constraints.tolist(),
                                    objective_gradient=evaluation.objective_gradient.tolist(),
                                    raw_constraints=measures['raw_constraints'],kkt_scaled=measures['kkt_scaled'],
                                    rank=rank_min,condition=condition_max,
                                    linear_residual=max_linear_residual,
                                    internal_energy=evaluation.internal_energy,physical_energy=evaluation.physical_energy))
            if (np.linalg.norm(x)>1e4*scene.characteristic_length
                    or np.linalg.norm(evaluation.velocity)>1e6*scene.characteristic_length):
                raise FloatingPointError('runaway_guard')
            if method=='CPD' and policy=='paper' and measures['raw_constraints']<=tol:
                break
            if method=='CPD' and policy=='stationary' and measures['kkt_scaled']<=1e-9:
                break
        if method=='PD':
            # Conservation rows are diagnostics only for ordinary PD.
            stationary = measures['stationarity']<=1e-9
            converged = stationary
        else:
            converged = (measures['raw_constraints']<=tol if policy=='paper'
                         else measures['kkt_scaled']<=1e-9)
            stationary = measures['kkt_scaled']<=1e-9
        reason = ('converged' if converged else 'fixed_budget_residual' if policy=='fixed' else 'iteration_cap')
        return StepResult(x,(x-xn)/h,rho,multipliers,True,bool(converged),bool(stationary),
                          reason,count,evaluation,measures,targets,
                          0 if rank_min==1000 else rank_min,condition_max,maxjump,history,max_linear_residual)
    except (FloatingPointError,OverflowError,ValueError,np.linalg.LinAlgError) as exc:
        return StepResult(x,(x-xn)/h,rho,multipliers,False,False,False,str(exc),count,
                          evaluation,measures,targets,0 if rank_min==1000 else rank_min,
                          condition_max,maxjump,history,max_linear_residual)
