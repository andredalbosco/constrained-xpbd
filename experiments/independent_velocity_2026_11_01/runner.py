"""Measure independent-velocity EL/LE/JOINT without reference phase scoring.

The solver must be validated separately before executing this campaign.
All kinetic-energy diagnostics use the returned independent velocity unless
explicitly labeled chord/reconstructed. Failed trial states remain in CSV.
"""
from __future__ import annotations

import argparse
import csv
import json
import platform
from pathlib import Path
import sys
import time

import numpy as np

import independent_solver as solver


ROOT = Path(__file__).resolve().parent
BASELINE = Path("/home/blacknoir/Projects/constrained-xpbd/experiments/block_orders_2026_11_01/results")


def write_csv(path, rows):
    if not rows:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def load_csv(path):
    with path.open(newline="") as f:
        return list(csv.DictReader(f))


def row_at(scene, step, h, x, v, chord, energy0, energy_scale,
           epsilon, beta_v, out=None, dx=0., dv=0., chord_dv=0.):
    radius = float(np.linalg.norm(x))
    radial_speed = float(x @ v) / radius if radius else np.nan
    kinetic = .5 * scene.m * float(v @ v)
    mechanical = solver.energy(scene, x, v)
    chord_kinetic = .5 * scene.m * float(chord @ chord)
    chord_energy = chord_kinetic + solver.potential(scene, x)
    evaluation = None if out is None else out.evaluation
    residual = 0. if evaluation is None else solver.scaled_residual(scene, h, epsilon, beta_v, evaluation)
    mismatch = float(np.linalg.norm(v - chord))
    angular = scene.m * float(x[0] * v[1] - x[1] * v[0])
    return dict(step=step, time=step * h, x=float(x[0]), y=float(x[1]),
                vx=float(v[0]), vy=float(v[1]), chord_vx=float(chord[0]), chord_vy=float(chord[1]),
                energy=mechanical, rel_energy=(mechanical - energy0) / energy_scale,
                potential=solver.potential(scene, x), kinetic=kinetic,
                chord_energy=chord_energy, chord_rel_energy=(chord_energy - energy0) / energy_scale,
                chord_kinetic=chord_kinetic, velocity_chord_mismatch=mismatch,
                displacement_velocity_mismatch=h * mismatch,
                pi_x=scene.m * float(v[0] - chord[0]) / beta_v,
                pi_y=scene.m * float(v[1] - chord[1]) / beta_v,
                velocity_chord_dot=float(v @ chord),
                speed=float(np.linalg.norm(v)), chord_speed=float(np.linalg.norm(chord)),
                radius=radius, radial_kinetic=.5 * scene.m * radial_speed**2,
                tangential_kinetic=max(0., kinetic - .5 * scene.m * radial_speed**2),
                angular_momentum=angular,
                geometry=abs(radius - scene.l) if scene.rigid else 0.,
                tangency=abs(float(x @ v)) / scene.l if scene.rigid else 0.,
                chord_tangency=abs(float(x @ chord)) / scene.l if scene.rigid else 0.,
                displacement=float(dx), velocity_jump=float(dv), chord_velocity_jump=float(chord_dv),
                **{"lambda": 0. if out is None else out.lam},
                rho=0. if out is None else out.rho, eta=0. if out is None else out.eta,
                residual=residual, accepted=1 if out is None else int(out.accepted),
                r_x_norm=0. if evaluation is None else float(np.linalg.norm(evaluation.r_x)),
                r_v_norm=0. if evaluation is None else float(np.linalg.norm(evaluation.r_v)),
                r_lambda=0. if evaluation is None else evaluation.r_lambda,
                r_rho=0. if evaluation is None else evaluation.r_rho,
                r_eta=0. if evaluation is None else evaluation.r_eta,
                converged=1 if out is None else int(out.converged))


def finite_max(rows, field, default=0.):
    values = np.asarray([r[field] for r in rows], float)
    return float(np.max(values)) if len(values) else default


def run_one(scene, order, hz, sweeps, beta_v, ratio, duration, trace_first=False):
    h = 1. / hz
    x, v = scene.x0.copy(), scene.v0.copy()
    chord = v.copy()
    energy0 = solver.energy(scene, x, v)
    floor = solver.energy_floor(scene)
    scale = max(energy0 - floor, 1e-12)
    epsilon = ratio * scale
    run_id = f"{scene.name}_{order}_beta{beta_v:g}_hz{hz}_s{sweeps}_eps{ratio:g}"
    rows = [row_at(scene, 0, h, x, v, chord, energy0, scale, epsilon, beta_v)]
    total = round(duration * hz)
    accepted_steps = converged = visits = rank_defects = 0
    block_jump = velocity_block_jump = 0.
    reason = ""
    first_trace = []
    started = time.perf_counter()
    for step in range(1, total + 1):
        old_x, old_v, old_chord = x.copy(), v.copy(), chord.copy()
        with np.errstate(over="raise", invalid="raise", divide="raise"):
            out = solver.block_step(scene, x, v, h, epsilon, beta_v=beta_v,
                                    order=order, sweeps=sweeps, trace=(step == 1 and trace_first))
        visits += out.visits
        rank_defects += getattr(out, "rank_defect_visits", 0)
        block_jump = max(block_jump, out.max_block_displacement)
        velocity_block_jump = max(velocity_block_jump, getattr(out, "max_velocity_increment", 0.))
        if step == 1 and trace_first:
            first_trace = out.history
        chord = (out.x - old_x) / h
        row = row_at(scene, step, h, out.x, out.v, chord, energy0, scale, epsilon, beta_v, out,
                     np.linalg.norm(out.x - old_x), np.linalg.norm(out.v - old_v),
                     np.linalg.norm(chord - old_chord))
        rows.append(row)
        if not out.accepted:
            reason = out.reason
            break
        x, v = out.x, out.v
        accepted_steps = step
        converged += int(out.converged)
    accepted = [r for r in rows if r["accepted"]]
    tail = [r for r in accepted[1:] if r["time"] > max(0., accepted_steps * h - 1.) + 1e-10]
    initial_angular = rows[0]["angular_momentum"]
    actual_errors = [abs(r["rel_energy"]) for r in accepted]
    chord_errors = [abs(r["chord_rel_energy"]) for r in accepted]
    path = sum(r["displacement"] for r in accepted)
    velocity_travel = h * sum(r["speed"] for r in accepted[1:])
    stall_threshold = 1e-4 * scene.l
    current_stall = longest_stall = 0.
    for row in accepted[1:]:
        current_stall = current_stall + h if row["chord_speed"] < stall_threshold else 0.
        longest_stall = max(longest_stall, current_stall)
    last = accepted[-1]
    summary = dict(run_id=run_id, case=scene.name, order=order, beta_v=beta_v,
                   hz=hz, sweeps=sweeps, epsilon_ratio=ratio, epsilon=epsilon, duration=duration,
                   completed=int(accepted_steps == total and not reason), accepted_time=accepted_steps * h,
                   reason=reason, max_rel_energy=max(abs(r["rel_energy"]) for r in rows),
                   max_accepted_rel_energy=max(actual_errors), final_rel_energy=last["rel_energy"],
                   max_accepted_chord_rel_energy=max(chord_errors), final_chord_rel_energy=last["chord_rel_energy"],
                   failed_trial_rel_energy=rows[-1]["rel_energy"] if not rows[-1]["accepted"] else np.nan,
                   path_length=path, velocity_travel=velocity_travel,
                   position_velocity_travel_ratio=path / max(velocity_travel, 1e-30),
                   tail_path_length=sum(r["displacement"] for r in tail),
                   tail_velocity_travel=h * sum(r["speed"] for r in tail),
                   tail_stalled=int(bool(tail) and sum(r["displacement"] for r in tail) < 1e-4 * scene.l),
                   longest_stall_duration=longest_stall,
                   max_accepted_displacement=finite_max(accepted, "displacement"),
                   max_accepted_velocity_jump=finite_max(accepted, "velocity_jump"),
                   max_accepted_chord_velocity_jump=finite_max(accepted, "chord_velocity_jump"),
                   max_velocity_chord_mismatch=finite_max(accepted, "velocity_chord_mismatch"),
                   rms_velocity_chord_mismatch=float(np.sqrt(np.mean([r["velocity_chord_mismatch"]**2 for r in accepted]))),
                   max_displacement_velocity_mismatch=finite_max(accepted, "displacement_velocity_mismatch"),
                   reversed_velocity_fraction=float(np.mean([r["velocity_chord_dot"] < -1e-10 for r in accepted[1:]])) if accepted_steps else np.nan,
                   max_geometry_error=finite_max(accepted, "geometry"),
                   max_tangency_error=finite_max(accepted, "tangency"),
                   max_chord_tangency_error=finite_max(accepted, "chord_tangency"),
                   max_residual=finite_max(accepted, "residual"),
                   converged_fraction=converged / max(accepted_steps, 1),
                   max_block_displacement=block_jump, max_block_velocity_increment=velocity_block_jump,
                   max_abs_rho=max(abs(r["rho"]) for r in accepted),
                   max_abs_eta=max(abs(r["eta"]) for r in accepted),
                   initial_angular_momentum=initial_angular, final_angular_momentum=last["angular_momentum"],
                   max_angular_error=max(abs(r["angular_momentum"] - initial_angular) for r in accepted),
                   angular_conservation_expected=int(not scene.rigid and np.linalg.norm(scene.g) == 0.),
                   min_radius=min(r["radius"] for r in accepted), max_radius=finite_max(accepted, "radius"),
                   final_radial_kinetic=last["radial_kinetic"], final_tangential_kinetic=last["tangential_kinetic"],
                   seconds=time.perf_counter() - started, visits=visits, rank_defect_visits=rank_defects,
                   energy_scale=scale, initial_energy=energy0, energy_floor=floor,
                   target_policy="incoming physical energy using independent velocity")
    return summary, rows, first_trace


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--duration", type=float, default=5.)
    parser.add_argument("--hz", default="60,120")
    parser.add_argument("--sweeps", default="1,4,16")
    parser.add_argument("--beta-v", default="0.1,1,10")
    parser.add_argument("--epsilon-ratios", default="10000,100000000")
    parser.add_argument("--orders", default="EL,LE,JOINT")
    parser.add_argument("--cases", default="")
    parser.add_argument("--output", type=Path, default=ROOT / "results")
    parser.add_argument("--trace-first", action="store_true")
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    scenes = [s for s in solver.scenes() if not args.cases or s.name in args.cases.split(",")]
    summaries, traces = [], {}
    for scene in scenes:
        for hz in map(int, args.hz.split(",")):
            for sweeps in map(int, args.sweeps.split(",")):
                for order in args.orders.split(","):
                    for beta_v in map(float, args.beta_v.split(",")):
                        if beta_v <= 0:
                            raise ValueError("Primal independent-velocity experiment requires beta_v > 0")
                        for ratio in map(float, args.epsilon_ratios.split(",")):
                            summary, rows, trace = run_one(scene, order, hz, sweeps, beta_v,
                                                           ratio, args.duration, args.trace_first)
                            summaries.append(summary)
                            write_csv(args.output / "trajectories" / (summary["run_id"] + ".csv"), rows)
                            if trace:
                                traces[summary["run_id"]] = trace
                            print(f"{summary['run_id']} complete={summary['completed']} "
                                  f"time={summary['accepted_time']:.3f} "
                                  f"E={summary['max_accepted_rel_energy']:.4g} "
                                  f"chordE={summary['max_accepted_chord_rel_energy']:.4g} "
                                  f"mismatch={summary['max_velocity_chord_mismatch']:.4g} "
                                  f"tail={summary['tail_path_length']:.4g} {summary['reason']}", flush=True)
                            if len(summaries) % 12 == 0:
                                write_csv(args.output / "summary.csv", summaries)
    write_csv(args.output / "summary.csv", summaries)
    (args.output / "first_step_traces.json").write_text(json.dumps(traces, indent=2))
    metadata = dict(command=sys.argv, python=platform.python_version(), numpy=np.__version__,
                    actual_energy=".5*m*dot(returned_v,returned_v) + U(returned_x)",
                    chord_energy=".5*m*norm((x_new-x_old)/h)^2 + U(x_new)",
                    mismatch="returned_v - (x_new-x_old)/h; never overwritten after a step",
                    target="incoming actual physical energy fixed during each step",
                    refresh="after every full block", stationarity="retained in every block",
                    epsilon="epsilon_ratio times initial energy above potential floor; fixed per run",
                    failure_policy="stop after solver rejection; offending candidate retained with accepted=0",
                    baseline_source=str(BASELINE), baseline_reused_without_solver_reruns=True,
                    primary_tangency="diagnostic only; no extra spring tangency constraint",
                    stall_threshold="chord speed below 1e-4*l m/s; tail is final accepted second",
                    scenes=[dict(name=s.name, m=s.m, k=s.k, l=s.l, g=s.g.tolist(),
                                 x0=s.x0.tolist(), v0=s.v0.tolist(), rigid=s.rigid) for s in scenes])
    (args.output / "metadata.json").write_text(json.dumps(metadata, indent=2))
    (args.output / "derived_velocity_baseline_summary.csv").write_bytes((BASELINE / "summary.csv").read_bytes())
    print(f"Wrote {len(summaries)} independent-velocity runs to {args.output}", flush=True)


if __name__ == "__main__":
    main()
