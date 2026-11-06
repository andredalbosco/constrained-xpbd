"""Measure paper-equation CPD, with anchored/rod extensions kept explicit.

Mechanical energy, the paper's internal-energy budget, external-power error,
and energy relaxation are recorded separately. No phase/reference score.
Execute only after independent validation of cpd_solver.py has passed.
"""
from __future__ import annotations

import argparse
import csv
from dataclasses import asdict
import hashlib
import json
from pathlib import Path
import platform
import sys
import time

import numpy as np

import cpd_solver as solver
from motion_metrics import body_modes, trajectory_motion


ROOT = Path(__file__).resolve().parent
LEGACY = Path("/home/blacknoir/Projects/constrained-xpbd/experiments/independent_velocity_2026_11_01/results")


def write_csv(path, rows):
    if not rows:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def energies(scene, x, v):
    # Scalar energies remain meaningful even at a rejected coincident spring,
    # where the solver's spring gradient is undefined.
    elastic = sum(.5 * edge.stiffness * (np.linalg.norm(solver.edge_difference(scene, x, edge)) - edge.rest)**2
                  for edge in scene.edges if not edge.hard)
    kinetic = .5 * np.sum(scene.masses[:, None] * v * v)
    internal = float(kinetic + elastic)
    return float(kinetic), float(elastic), internal, internal - float(np.sum(scene.forces * x))


def maximum(rows, field, absolute=False):
    values = np.asarray([r[field] for r in rows], float)
    values = np.abs(values) if absolute else values
    return float(np.nanmax(values)) if len(values) and not np.isnan(values).all() else np.nan


def run_name(scene,hz,ratio,policy,iterations,max_iterations,method,angular,precheck):
    suffix=f"{policy}{iterations if policy == 'fixed' else max_iterations}"
    return f"{scene.name}_{method}_{suffix}_hz{hz}_eps{ratio:g}" + ("_noL" if not angular else "") + ("_precheck" if precheck else "")


def trajectory_row(scene, step, h, x, v, energy0, internal0, scale, cumulative,
                   out=None, step_paper_work=0., step_exact_work=0.):
    kin, elastic, internal, mechanical = energies(scene, x, v)
    modes = body_modes(x, v, scene.masses)
    row = dict(step=step, time=step*h)
    for i in range(len(x)):
        row.update({f"x_{i}": float(x[i, 0]), f"y_{i}": float(x[i, 1]),
                    f"vx_{i}": float(v[i, 0]), f"vy_{i}": float(v[i, 1])})
    names = [] if out is None or out.evaluation is None else out.evaluation.names
    multipliers = {} if out is None else dict(zip(names, out.multipliers))
    target = internal0 if out is None else out.targets.energy
    kstar = 0. if out is None else out.targets.kstar
    rho = 0. if out is None else out.rho
    relaxed_target = (1-rho)*target + rho*kstar
    energy_residual = 0. if out is None else internal - relaxed_target
    paper_budget = internal0 + cumulative["paper_work"]
    relaxed_budget = paper_budget + cumulative["relaxation"]
    balance_prediction = (cumulative["paper_work"] - cumulative["exact_work"]
                          + cumulative["relaxation"] + cumulative["constraint_defect"])
    edge_lengths = [float(np.linalg.norm(solver.edge_difference(scene, x, edge))) for edge in scene.edges]
    hard_geometry = max([abs(length-edge.rest) for length, edge in zip(edge_lengths, scene.edges) if edge.hard] or [0.])
    tangent = 0.
    for edge in scene.edges:
        if edge.hard:
            dx = solver.edge_difference(scene, x, edge)
            dv = v[edge.a]-(v[edge.b] if edge.b >= 0 else 0.)
            tangent = max(tangent, abs(float(dx @ dv))/edge.rest)
    measures = {} if out is None else out.measures
    row.update(kinetic=kin, elastic=elastic, internal_energy=internal,
               energy=mechanical, rel_energy=(mechanical-energy0)/scale,
               algorithm_target=target, kstar=kstar, relaxed_target=relaxed_target,
               algorithm_target_defect=internal-target, relaxed_target_defect=energy_residual,
               paper_work=step_paper_work, exact_external_work=step_exact_work,
               cumulative_paper_work=cumulative["paper_work"], cumulative_exact_work=cumulative["exact_work"],
               cumulative_relaxation=cumulative["relaxation"],
               cumulative_constraint_defect=cumulative["constraint_defect"],
               algorithm_budget=paper_budget, relaxed_algorithm_budget=relaxed_budget,
               algorithm_budget_defect=internal-paper_budget, relaxed_budget_defect=internal-relaxed_budget,
               work_approximation_defect=cumulative["paper_work"]-cumulative["exact_work"],
               mechanical_budget_identity_error=mechanical-energy0-balance_prediction,
               rho=rho, eta=float(multipliers.get("energy", 0.)),
               angular_multiplier=float(multipliers.get("angular", 0.)),
               linear_multiplier_x=float(multipliers.get("linear_0", 0.)),
               linear_multiplier_y=float(multipliers.get("linear_1", 0.)),
               geometry_multiplier_max=max([abs(float(z)) for name,z in multipliers.items() if name.startswith("geometry")] or [0.]),
               raw_constraints=measures.get("raw_constraints", 0.),
               scaled_constraints=measures.get("scaled_constraints", 0.),
               stationarity=measures.get("stationarity", 0.), kkt_scaled=measures.get("kkt_scaled", 0.),
               linear_target_error=0. if out is None or out.evaluation is None else float(np.linalg.norm(out.evaluation.linear-out.targets.linear)),
               angular_target_error=0. if out is None or out.evaluation is None else abs(out.evaluation.angular-out.targets.angular),
               geometry=hard_geometry, tangency=tangent,
               min_edge_length=min(edge_lengths), max_edge_length=max(edge_lengths),
               nominal_length_error=max(abs(length-edge.rest) for length,edge in zip(edge_lengths,scene.edges)),
               accepted=1 if out is None else int(out.accepted),
               converged=1 if out is None else int(out.converged),
               stationary=1 if out is None else int(out.stationary),
               iterations=0 if out is None else out.iterations,
               schur_rank=0 if out is None else out.rank_min,
               conservation_rows=len(names), condition=0. if out is None else out.condition_max,
               linear_solve_residual=0. if out is None else getattr(out,"max_linear_residual",np.nan),
               max_block_displacement=0. if out is None else out.max_block_displacement)
    row.update(modes)
    return row


def run_one(scene, hz, ratio, policy, iterations, duration, max_iterations=200,
            angular=True, method="CPD", tol=1e-4, initial_precheck=False, trace_first=True):
    h = 1./hz
    x, v = scene.x0.copy(), scene.v0.copy()
    _, _, internal0, energy0 = energies(scene, x, v)
    scale = max(energy0-scene.energy_floor, 1e-12)
    epsilon = ratio*scale
    run_id = run_name(scene,hz,ratio,policy,iterations,max_iterations,method,angular,initial_precheck)
    cumulative = dict(paper_work=0., exact_work=0., relaxation=0., constraint_defect=0.)
    rows = [trajectory_row(scene, 0, h, x, v, energy0, internal0, scale, cumulative)]
    total = round(duration*hz)
    accepted_steps = count_converged = count_stationary = visits = rank_defects = caps = 0
    rank_min, condition, block_jump = 1000, 0., 0.
    reason = ""
    first_trace = []
    started = time.perf_counter()
    for index in range(1, total+1):
        old_x, old_v = x.copy(), v.copy()
        with np.errstate(over="raise", invalid="raise", divide="raise"):
            out = solver.step(scene, x, v, h, epsilon, policy=policy, iterations=iterations,
                              max_iterations=max_iterations, tol=tol, angular=angular,
                              initial_precheck=initial_precheck, trace=(index == 1 and trace_first), method=method)
        visits += out.iterations
        rank_min = min(rank_min, out.rank_min)
        condition = max(condition, out.condition_max)
        block_jump = max(block_jump, out.max_block_displacement)
        caps += int(out.reason == "iteration_cap")
        if method == "CPD" and out.evaluation is not None:
            rank_defects += int(out.rank_min < len(out.evaluation.names))
        if index == 1:
            first_trace = out.history
        paper_work = h*float(np.sum(scene.forces*old_v))
        exact_work = float(np.sum(scene.forces*(out.x-old_x)))
        _, _, current_internal, _ = energies(scene, out.x, out.v)
        relaxation = -out.rho*(out.targets.energy-out.targets.kstar)
        defect = current_internal-((1-out.rho)*out.targets.energy+out.rho*out.targets.kstar)
        trial_cumulative = {"paper_work": cumulative["paper_work"]+paper_work,
                            "exact_work": cumulative["exact_work"]+exact_work,
                            "relaxation": cumulative["relaxation"]+relaxation,
                            "constraint_defect": cumulative["constraint_defect"]+defect}
        row = trajectory_row(scene,index,h,out.x,out.v,energy0,internal0,scale,trial_cumulative,
                             out,paper_work,exact_work)
        rows.append(row)
        if not out.accepted:
            reason = out.reason
            break
        cumulative = trial_cumulative
        x,v = out.x,out.v
        accepted_steps = index
        count_converged += int(out.converged)
        count_stationary += int(out.stationary)
    accepted = [r for r in rows if r["accepted"]]
    motion = trajectory_motion(rows,h=h,length_scale=scene.characteristic_length)
    initial_angular = accepted[0]["angular_global"]
    final = accepted[-1]
    fields = dict(run_id=run_id,case=scene.name,family=scene.family,reference_case=scene.reference_case,
                  method=method,policy=policy,iterations_budget=iterations,max_iterations=max_iterations,
                  raw_constraint_tolerance=tol,
                  hz=hz,epsilon_ratio=ratio,epsilon=epsilon,angular_active=int(angular and method == "CPD"),
                  linear_active=int(scene.family == "free" and method == "CPD"),initial_precheck=int(initial_precheck),
                  duration=duration,completed=int(accepted_steps == total and not reason),
                  accepted_time=accepted_steps*h,reason=reason,
                  max_rel_energy=maximum(rows,"rel_energy",True),
                  max_accepted_rel_energy=maximum(accepted,"rel_energy",True),
                  final_rel_energy=final["rel_energy"],
                  max_algorithm_target_defect=maximum(accepted,"algorithm_target_defect",True)/scale,
                  max_relaxed_target_defect=maximum(accepted,"relaxed_target_defect",True)/scale,
                  max_algorithm_budget_defect=maximum(accepted,"algorithm_budget_defect",True)/scale,
                  max_relaxed_budget_defect=maximum(accepted,"relaxed_budget_defect",True)/scale,
                  max_work_approximation_defect=maximum(accepted,"work_approximation_defect",True)/scale,
                  max_mechanical_budget_identity_error=maximum(accepted,"mechanical_budget_identity_error",True)/scale,
                  max_abs_rho=maximum(accepted,"rho",True),min_rho=min(r["rho"] for r in accepted),max_rho=maximum(accepted,"rho"),
                  max_abs_eta=maximum(accepted,"eta",True),
                  max_raw_constraints=maximum(accepted,"raw_constraints"),
                  max_scaled_constraints=maximum(accepted,"scaled_constraints"),
                  max_stationarity=maximum(accepted,"stationarity"),max_kkt_scaled=maximum(accepted,"kkt_scaled"),
                  converged_fraction=count_converged/max(accepted_steps,1),
                  stationary_fraction=count_stationary/max(accepted_steps,1),
                  iteration_cap_steps=caps,solver_iterations=visits,
                  mean_iterations=visits/max(len(rows)-1,1),min_schur_rank=rank_min,
                  rank_defect_steps=rank_defects,max_condition=condition,
                  max_linear_solve_residual=maximum(accepted,"linear_solve_residual"),
                  max_block_displacement=block_jump,
                  max_linear_target_error=maximum(accepted,"linear_target_error"),
                  max_angular_target_error=maximum(accepted,"angular_target_error"),
                  max_angular_drift=max(abs(r["angular_global"]-initial_angular) for r in accepted),
                  initial_angular=initial_angular,final_angular=final["angular_global"],
                  final_radial_kinetic=final["radial_kinetic"],final_tangential_kinetic=final["tangential_kinetic"],
                  final_rotational_kinetic=final["rotational_min_kinetic"],final_shape_kinetic=final["shape_kinetic"],
                  max_geometry_error=maximum(accepted,"geometry"),max_tangency_error=maximum(accepted,"tangency"),
                  max_nominal_length_error=maximum(accepted,"nominal_length_error"),
                  energy_scale=scale,initial_energy=energy0,initial_internal_energy=internal0,
                  energy_floor=scene.energy_floor,seconds=time.perf_counter()-started)
    fields.update(motion)
    return fields,rows,first_trace


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--duration",type=float,default=5.)
    parser.add_argument("--hz",default="30,60,120")
    parser.add_argument("--epsilon-ratios",default="100,10000,100000000")
    parser.add_argument("--policies",default="paper,fixed,stationary")
    parser.add_argument("--iterations",type=int,default=16)
    parser.add_argument("--max-iterations",type=int,default=200)
    parser.add_argument("--tol",type=float,default=1e-4)
    parser.add_argument("--methods",default="CPD")
    parser.add_argument("--cases",default="")
    parser.add_argument("--without-angular",action="store_true")
    parser.add_argument("--initial-precheck",action="store_true")
    parser.add_argument("--no-trace-first",action="store_true")
    parser.add_argument("--resume",action="store_true")
    parser.add_argument("--output",type=Path,default=ROOT/"results")
    args = parser.parse_args()
    args.output.mkdir(parents=True,exist_ok=True)
    scene_list = [s for s in solver.scenes() if not args.cases or s.name in args.cases.split(",")]
    summaries,traces,skips = [],{},[]
    solver_hash=hashlib.sha256(Path(solver.__file__).read_bytes()).hexdigest()
    if args.resume and (args.output/"summary.csv").exists():
        with (args.output/"summary.csv").open(newline="") as f:
            summaries=list(csv.DictReader(f))
        old_metadata=json.loads((args.output/"metadata.json").read_text())
        if old_metadata.get("solver_sha256") != solver_hash:
            raise ValueError("Solver source changed; do not mix checkpoint trajectories")
        if any(float(r["duration"]) != args.duration or float(r["raw_constraint_tolerance"]) != args.tol for r in summaries):
            raise ValueError("Checkpoint duration/tolerance mismatch")
        traces=json.loads((args.output/"first_step_traces.json").read_text())
    done={r["run_id"] for r in summaries}
    for scene in scene_list:
        for method in args.methods.split(","):
            if method == "PD" and any(edge.hard for edge in scene.edges):
                skips.append(dict(case=scene.name,method=method,reason="finite-spring PD has no exact hard-rod row"))
                continue
            for hz in map(int,args.hz.split(",")):
                for ratio in map(float,args.epsilon_ratios.split(",")):
                    for policy in args.policies.split(","):
                        run_id=run_name(scene,hz,ratio,policy,args.iterations,args.max_iterations,method,
                                        not args.without_angular,args.initial_precheck)
                        if run_id in done:
                            continue
                        summary,rows,trace = run_one(scene,hz,ratio,policy,args.iterations,args.duration,
                                                     args.max_iterations,not args.without_angular,
                                                     method,args.tol,args.initial_precheck,not args.no_trace_first)
                        summaries.append(summary)
                        done.add(run_id)
                        write_csv(args.output/"trajectories"/(summary["run_id"]+".csv"),rows)
                        traces[summary["run_id"]] = trace
                        print(f"{summary['run_id']} full={summary['completed']} t={summary['accepted_time']:.3f} "
                              f"actualE={summary['max_accepted_rel_energy']:.4g} "
                              f"target={summary['max_relaxed_target_defect']:.4g} "
                              f"path={summary['path_length']:.4g} tail={summary['tail_path_length']:.4g} "
                              f"stationary={summary['stationary_fraction']:.3g} {summary['reason']}",flush=True)
                        if len(summaries)%9 == 0:
                            write_csv(args.output/"summary.csv",summaries)
    write_csv(args.output/"summary.csv",summaries)
    (args.output/"first_step_traces.json").write_text(json.dumps(traces,indent=2))
    (args.output/"metadata.json").write_text(json.dumps(dict(
        command=sys.argv,python=platform.python_version(),numpy=np.__version__,solver_sha256=solver_hash,
        actual_energy="K + elastic - sum external_constant_force dot x",
        algorithm_energy="H = K + elastic",target="Hstar = incoming H + h sum f dot incoming v",
        relaxation="H = (1-rho) Hstar + rho Kstar; rho unbounded",
        epsilon="epsilon_ratio times initial mechanical energy above potential floor, fixed per run; local experimental choices",
        velocity="derived from endpoint difference; never independent",
        stopping="paper: raw allconstraint max<=1e-4 after at least one global step; fixed:16; stationary:fullscaled KKT<=1e-9, max200",
        anchored="P row disabled, L about fixed origin retained; explicit anchored adaptation",
        hard_rod="explicit global hard-distance row extension, not paper finite-spring model",
        pd_control="unconstrained finite-spring PD; conservation rows/targets are diagnostic only; no Schur rank defects are counted",
        angular_ablation="row disabled while Kstar held unchanged" if args.without_angular else "none",
        baseline_source=str(LEGACY),baseline_reused_without_rerun=True,
        scenes=[asdict(s) for s in solver.scenes() if s.name in {r["case"] for r in summaries}],skipped=skips),indent=2,
        default=lambda z:z.tolist() if isinstance(z,np.ndarray) else str(z)))
    print(f"Wrote {len(summaries)} measured runs to {args.output}",flush=True)


if __name__ == "__main__":
    main()
