"""Continuous DOP853 references, independent of the block-GS discretization.

Uses the same angle convention and equations as the earlier verified references
in experiments/nonlinear_energy_2026_11_01/references.py. No endpoint-energy
projection or correction is applied.
"""
import numpy as np
from scipy.integrate import solve_ivp


def reference_trajectory(scene, times, tighter=False):
    atol, rtol = ((2e-14, 3e-13) if tighter else (1e-13, 2e-12))
    if scene.rigid:
        if abs(np.linalg.norm(scene.x0)-scene.l) > 1e-12*scene.l or abs(scene.x0 @ scene.v0) > 1e-12:
            raise ValueError('continuous rod reference requires feasible initial position and velocity')
        theta = np.arctan2(scene.x0[0], -scene.x0[1])
        omega = np.array([np.cos(theta), np.sin(theta)]) @ scene.v0/scene.l
        def ode(t, state):
            angle, speed = state
            tangent = np.array([np.cos(angle), np.sin(angle)])
            return [speed, scene.g @ tangent/scene.l]
        initial = [theta, omega]
    else:
        def ode(t, state):
            position, velocity = state[:2], state[2:]
            radius = np.linalg.norm(position)
            force = scene.k*(radius-scene.l)*position/radius
            return np.r_[velocity, scene.g-force/scene.m]
        initial = np.r_[scene.x0, scene.v0]
    solution = solve_ivp(ode, [0., float(times[-1])], initial, method='DOP853',
                         t_eval=times, rtol=rtol, atol=atol)
    if not solution.success:
        raise RuntimeError(solution.message)
    if scene.rigid:
        angle, speed = solution.y
        positions = scene.l*np.column_stack([np.sin(angle), -np.cos(angle)])
        velocities = scene.l*speed[:, None]*np.column_stack([np.cos(angle), np.sin(angle)])
    else:
        positions, velocities = solution.y[:2].T, solution.y[2:].T
    return positions, velocities
