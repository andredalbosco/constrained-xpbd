"""Reproduce the same-state block-order diagnostic using production increments.

Run with the experiment's NumPy-capable Python environment:
    python3 diagnose.py

This writes diagnostics/failing_step.json relative to this file. The report in
DIAGNOSTIC.md explains the trace; no solver state is clamped or repaired here.
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np

import solver


ROOT = Path(__file__).resolve().parent


def native(value):
    """Convert diagnostic arrays and nonfinite sentinels to strict JSON."""
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return native(value.item())
    if isinstance(value, float) and not math.isfinite(value):
        return str(value)
    if isinstance(value, dict):
        return {key: native(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [native(item) for item in value]
    return value


def evaluation_record(scene, h, epsilon, evaluation):
    data = vars(evaluation).copy()
    inverse = h*h/scene.m
    gap = evaluation.target-evaluation.kstar
    data.update(
        energy_gradient_norm=float(np.linalg.norm(evaluation.energy_gradient)),
        r_x_norm=float(np.linalg.norm(evaluation.r_x)),
        material_denominator=evaluation.alpha+inverse*float(evaluation.normal@evaluation.normal),
        energy_denominator=inverse*float(evaluation.energy_gradient@evaluation.energy_gradient)+gap*gap/epsilon,
        energy_gradient_denominator_contribution=inverse*float(evaluation.energy_gradient@evaluation.energy_gradient),
        energy_relaxation_denominator_contribution=gap*gap/epsilon,
        energy_numerator_energy_residual=evaluation.r_eta,
        energy_numerator_stationarity=-inverse*float(evaluation.energy_gradient@evaluation.r_x),
        energy_numerator_relaxation=-gap*evaluation.r_rho/epsilon,
        material_numerator_constitutive=-evaluation.r_lambda,
        material_numerator_stationarity=inverse*float(evaluation.normal@evaluation.r_x),
        common_position_stationarity_correction=-inverse*evaluation.r_x,
        scaled_residual=solver.scaled_residual(scene, h, epsilon, evaluation),
    )
    return data


def replay_step(scene, x0, v0, h, epsilon, order, sweeps):
    """Add pre-visit diagnostics to the unchanged production block arithmetic."""
    x0, v0 = np.asarray(x0, float), np.asarray(v0, float)
    x = (x0+h*v0).copy()
    lam = rho = eta = 0.
    target = solver.energy(scene, x0, v0)
    kstar = solver.energy_floor(scene)
    records = []
    reason = ''
    accepted = True
    maximum = 0.
    blocks = ('J',) if order == 'JOINT' else tuple(order)
    evaluation = solver.evaluate_state(scene, x0, v0, h, x, lam, rho, eta, epsilon, target, kstar)
    try:
        for sweep in range(sweeps):
            for block in blocks:
                before_x = x.copy()
                before_lambda, before_rho, before_eta = lam, rho, eta
                before = evaluation
                if block == 'E':
                    change = solver.energy_increment(scene, h, epsilon, evaluation)
                elif block == 'L':
                    change = solver.material_increment(scene, h, evaluation)
                else:
                    change = solver.joint_increment(scene, h, epsilon, evaluation)
                x = x+change.dx
                lam += change.dlam
                rho += change.drho
                eta += change.deta
                maximum = max(maximum, float(np.linalg.norm(change.dx)))
                evaluation = solver.evaluate_state(scene, x0, v0, h, x, lam, rho, eta, epsilon, target, kstar)
                records.append(dict(
                    sweep=sweep,
                    block=block,
                    before=dict(x=before_x, lam=before_lambda, rho=before_rho, eta=before_eta,
                                evaluation=evaluation_record(scene, h, epsilon, before)),
                    increment=vars(change).copy(),
                    after=dict(x=x.copy(), lam=lam, rho=rho, eta=eta,
                               evaluation=evaluation_record(scene, h, epsilon, evaluation)),
                ))
                velocity = (x-x0)/h
                if (not np.all(np.isfinite(velocity)) or np.linalg.norm(x) > 1e4*scene.l
                        or np.linalg.norm(velocity) > 1e6*scene.l):
                    raise FloatingPointError('runaway_guard')
    except (FloatingPointError, OverflowError, ValueError) as exc:
        accepted = False
        reason = str(exc)
    result = dict(
        order=order, sweeps=sweeps, incoming_x=x0, incoming_v=v0,
        energy_target=target, energy_floor=kstar,
        output_x=x, output_v=(x-x0)/h, lam=lam, rho=rho, eta=eta,
        accepted=accepted, reason=reason, visits=len(records), max_block_displacement=maximum,
        final_evaluation=evaluation_record(scene, h, epsilon, evaluation), trace=records,
    )
    public = solver.block_step(scene, x0, v0, h, epsilon, order, sweeps, trace=True)
    assert public.accepted == accepted
    assert len(public.history) == len(records)
    assert np.array_equal(public.x, x)
    assert public.lam == lam and public.rho == rho and public.eta == eta
    return result


def generate():
    scene = solver.scenes()[0]
    assert scene.name == 'spring_release'
    h = 1/60
    epsilon = 1e4*(solver.energy(scene, scene.x0, scene.v0)-solver.energy_floor(scene))
    x, v = scene.x0.copy(), scene.v0.copy()
    rollout = []
    first_event = None
    for step in range(120):
        result = replay_step(scene, x, v, h, epsilon, 'EL', 1)
        result['physical_step'] = step+1
        result['physical_interval'] = [step*h, (step+1)*h]
        rollout.append(result)
        for visit_index, visit in enumerate(result['trace']):
            if np.linalg.norm(visit['increment']['dx']) > .1*scene.l:
                first_event = dict(
                    physical_step=step+1, visit_index=visit_index, block=visit['block'],
                    incoming_x=x.copy(), incoming_v=v.copy(), before=visit['before'],
                    increment=visit['increment'], after=visit['after'],
                )
                break
        if first_event is not None or not result['accepted']:
            break
        x, v = result['output_x'], result['output_v']
    assert first_event is not None, 'No large visit found within the declared replay interval'
    comparisons = [
        replay_step(scene, x, v, h, epsilon, order, sweeps)
        for order in ('EL', 'LE', 'JOINT') for sweeps in (1, 4, 16, 64)
    ]
    energy_visit = comparisons[0]['trace'][0]
    before = energy_visit['before']['evaluation']
    change = energy_visit['increment']
    after = energy_visit['after']['evaluation']
    frozen_position_residual = before['r_x']+scene.m*change['dx']/(h*h)+before['energy_gradient']*change['deta']
    energy_gradient_change = (scene.m/(h*h)+scene.k)*change['dx']
    refreshed_remainder = energy_visit['after']['eta']*energy_gradient_change
    identity = dict(
        frozen_position_residual=frozen_position_residual,
        actual_refreshed_position_residual=after['r_x'],
        analytical_gradient_change=energy_gradient_change,
        measured_gradient_change=after['energy_gradient']-before['energy_gradient'],
        eta_times_energy_gradient_change=refreshed_remainder,
        equality_error=float(np.linalg.norm(after['r_x']-refreshed_remainder)),
    )
    energy_remainder = .5*(scene.m/h**2+scene.k)*float(change['dx']@change['dx'])
    identity.update(
        exact_energy_second_order_remainder=energy_remainder,
        actual_after_energy_residual=after['r_eta'],
        energy_remainder_identity_error=abs(after['r_eta']-energy_remainder),
    )
    return dict(
        scene=vars(scene), h=h, epsilon=epsilon, epsilon_ratio=1e4,
        targets='Incoming accepted energy per physical step, fixed during the step',
        refresh_policy='After every block',
        initialization='x=x_n+h*v_n; lambda=rho=eta=0 every physical step',
        large_visit_threshold=.1*scene.l,
        solver_api_parity='Bitwise output state/multiplier equality for every recorded replay',
        first_large_visit=first_event, rollout_to_first_large_visit=rollout,
        comparisons_from_identical_incoming_state=comparisons,
        first_energy_visit_remainder_identity=identity,
    )


def main():
    data = generate()
    path = ROOT/'diagnostics'/'failing_step.json'
    path.parent.mkdir(exist_ok=True)
    path.write_text(json.dumps(native(data), indent=2, allow_nan=False)+'\n')
    event = data['first_large_visit']
    print(f'Wrote {path}: first large visit is step {event["physical_step"]}, block {event["block"]}.')


if __name__ == '__main__':
    main()
