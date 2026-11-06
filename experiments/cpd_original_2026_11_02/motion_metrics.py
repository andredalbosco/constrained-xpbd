"""Array-only motion diagnostics for anchored particles or small free networks.

These functions do not integrate a state or impose a conservation target.
"""
from __future__ import annotations

import numpy as np


def body_modes(x, v, masses):
    x, v = np.asarray(x, float).reshape(-1, 2), np.asarray(v, float).reshape(-1, 2)
    masses = np.asarray(masses, float).reshape(-1)
    total_mass = masses.sum()
    center = np.sum(masses[:, None] * x, axis=0) / total_mass
    center_velocity = np.sum(masses[:, None] * v, axis=0) / total_mass
    relative_x, relative_v = x - center, v - center_velocity
    linear = np.sum(masses[:, None] * v, axis=0)
    global_angular = np.sum(masses * (x[:, 0] * v[:, 1] - x[:, 1] * v[:, 0]))
    internal_angular = np.sum(masses * (relative_x[:, 0] * relative_v[:, 1] - relative_x[:, 1] * relative_v[:, 0]))
    kinetic = .5 * np.sum(masses * np.sum(v * v, axis=1))
    if len(masses) == 1:
        displacement, velocity, reduced_mass = x[0], v[0], masses[0]
    elif len(masses) == 2:
        displacement, velocity = x[0] - x[1], v[0] - v[1]
        reduced_mass = masses[0] * masses[1] / total_mass
    else:
        # For a network, use radial/tangential internal kinetic energy about
        # its moving COM. This is not a sum of double-counted edge energies.
        radii = np.linalg.norm(relative_x, axis=1)
        radial_speeds = np.divide(np.sum(relative_x * relative_v, axis=1), radii,
                                  out=np.zeros_like(radii), where=radii > 0)
        radial_kinetic = .5 * np.sum(masses * radial_speeds**2)
        relative_kinetic = .5 * np.sum(masses[:, None] * relative_v**2)
        inertia = np.sum(masses * radii**2)
        rotational = .5 * internal_angular**2 / inertia if inertia > 0 else 0.
        return dict(center_x=float(center[0]), center_y=float(center[1]),
                    center_vx=float(center_velocity[0]), center_vy=float(center_velocity[1]),
                    linear_px=float(linear[0]), linear_py=float(linear[1]),
                    angular_global=float(global_angular), angular_internal=float(internal_angular),
                    kinetic=float(kinetic), relative_kinetic=float(relative_kinetic),
                    radial_kinetic=float(radial_kinetic),
                    tangential_kinetic=float(max(0., relative_kinetic - radial_kinetic)),
                    mode_length=float(np.sqrt(inertia / total_mass)),
                    mode_x=np.nan, mode_y=np.nan, mode_vx=np.nan, mode_vy=np.nan,
                    radial_speed=float(np.sqrt(np.sum(masses * radial_speeds**2) / total_mass)),
                    relative_speed=float(np.sqrt(2 * relative_kinetic / total_mass)),
                    rotational_min_kinetic=float(rotational),
                    shape_kinetic=float(max(0., relative_kinetic - rotational)))
    length = float(np.linalg.norm(displacement))
    radial_speed = float(displacement @ velocity) / length if length else np.nan
    radial_kinetic = .5 * reduced_mass * radial_speed**2
    relative_kinetic = .5 * reduced_mass * float(velocity @ velocity)
    return dict(center_x=float(center[0]), center_y=float(center[1]),
                center_vx=float(center_velocity[0]), center_vy=float(center_velocity[1]),
                linear_px=float(linear[0]), linear_py=float(linear[1]),
                angular_global=float(global_angular), angular_internal=float(internal_angular),
                kinetic=float(kinetic), relative_kinetic=float(relative_kinetic),
                radial_kinetic=float(radial_kinetic),
                tangential_kinetic=float(max(0., relative_kinetic - radial_kinetic)),
                mode_length=length, radial_speed=radial_speed,
                mode_x=float(displacement[0]), mode_y=float(displacement[1]),
                mode_vx=float(velocity[0]), mode_vy=float(velocity[1]),
                relative_speed=float(np.linalg.norm(velocity)),
                rotational_min_kinetic=float(max(0., relative_kinetic - radial_kinetic)),
                shape_kinetic=float(radial_kinetic))


def trajectory_motion(rows, *, h, length_scale=1.):
    accepted = [r for r in rows if r["accepted"]]
    if not accepted:
        raise ValueError("Need an initial accepted state")
    keys = [k for k in accepted[0] if k.startswith("x_") or k.startswith("y_")]
    positions = np.array([[r[k] for k in keys] for r in accepted])
    increments = np.linalg.norm(np.diff(positions, axis=0), axis=1)
    times = np.array([r["time"] for r in accepted])
    end_time = times[-1]
    tail_increments = increments[times[1:] > max(0., end_time - 1.) + 1e-10]
    threshold = 1e-4 * length_scale
    stall = longest = 0.
    for jump in increments:
        stall = stall + h if jump / h < threshold else 0.
        longest = max(longest, stall)
    result = dict(path_length=float(increments.sum()),
                tail_path_length=float(tail_increments.sum()),
                tail_stalled=int(bool(len(tail_increments)) and tail_increments.sum() < threshold),
                longest_stall_duration=longest,
                max_accepted_displacement=float(increments.max()) if len(increments) else 0.)
    if np.isfinite(accepted[0].get("mode_x", np.nan)):
        mode = np.array([[r["mode_x"], r["mode_y"]] for r in accepted])
        mode_velocity = np.array([[r["mode_vx"], r["mode_vy"]] for r in accepted])
        mode_increments = np.linalg.norm(np.diff(mode, axis=0), axis=1)
        jumps = np.linalg.norm(np.diff(mode_velocity, axis=0), axis=1)
        result.update(relative_mode_path_length=float(mode_increments.sum()),
                      relative_mode_tail_path=float(mode_increments[times[1:] > max(0., end_time - 1.) + 1e-10].sum()),
                      max_relative_velocity_jump=float(jumps.max()) if len(jumps) else 0.)
    else:
        result.update(relative_mode_path_length=np.nan, relative_mode_tail_path=np.nan,
                      max_relative_velocity_jump=np.nan)
    return result
