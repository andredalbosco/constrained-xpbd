"""Run EL/LE block orders with explicit failure prefixes and motion references."""
from __future__ import annotations
import argparse
import csv
import json
from pathlib import Path
import platform
import sys
import time
import numpy as np
import scipy
from solver import scenes, energy, energy_floor, block_step, scaled_residual
from reference import reference_trajectory


def write_csv(path, rows):
    if not rows:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('w', newline='') as output:
        writer = csv.DictWriter(output, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def trajectory_row(scene, step, h, x, v, refx, refv, initial_energy, energy_scale,
                   out=None, epsilon=1.):
    physical = energy(scene, x, v)
    geometry = abs(np.linalg.norm(x)-scene.l) if scene.rigid else 0.
    tangency = abs(float(x @ v))/scene.l if scene.rigid else 0.
    residual = 0. if out is None or out.evaluation is None else scaled_residual(scene,h,epsilon,out.evaluation)
    return dict(step=step,time=step*h,x=x[0],y=x[1],vx=v[0],vy=v[1],
                energy=physical,rel_energy=(physical-initial_energy)/energy_scale,
                ref_x=refx[0],ref_y=refx[1],ref_vx=refv[0],ref_vy=refv[1],
                **{'lambda':0. if out is None else out.lam},
                rho=0. if out is None else out.rho,eta=0. if out is None else out.eta,
                residual=residual,geometry=geometry,tangency=tangency,
                accepted=1 if out is None else int(out.accepted),
                converged=1 if out is None else int(out.converged))


def run_one(scene, order, hz, sweeps, ratio, duration, reference, target_policy='incoming'):
    x,v=scene.x0.copy(),scene.v0.copy()
    h=1./hz
    initial_energy=energy(scene,x,v)
    energy_scale=max(initial_energy-energy_floor(scene),1e-12)
    # This dimensional epsilon stays fixed for the entire run.
    epsilon=(ratio if order!='XPBD' else 1.)*energy_scale
    run_id=f'{scene.name}_{order}_hz{hz}_s{sweeps}_eps{ratio:g}'
    positions,velocities=reference
    rows=[trajectory_row(scene,0,h,x,v,positions[0],velocities[0],initial_energy,energy_scale)]
    maxe=geometry=tangency=residual=blockjump=path=sumsq=angular=0.
    converged=visits=accepted_steps=rank_visits=0
    maxrho=maxeta=0.
    initial_angular=scene.m*(x[0]*v[1]-x[1]*v[0])
    reason=''
    first_trace=[]
    start=time.perf_counter()
    total=round(duration*hz)
    for step in range(1,total+1):
        old=x.copy()
        with np.errstate(over='raise',invalid='raise',divide='raise'):
            out=block_step(scene,x,v,h,epsilon,order,sweeps,trace=(step==1))
        visits+=out.visits
        rank_visits+=out.rank_defect_visits
        blockjump=max(blockjump,out.max_block_displacement)
        if step==1:
            first_trace=out.history
        row=trajectory_row(scene,step,h,out.x,out.v,positions[step],velocities[step],
                           initial_energy,energy_scale,out,epsilon)
        rows.append(row)
        maxe=max(maxe,abs(row['rel_energy']))
        geometry=max(geometry,row['geometry'])
        tangency=max(tangency,row['tangency'])
        residual=max(residual,row['residual'])
        maxrho=max(maxrho,abs(out.rho))
        maxeta=max(maxeta,abs(out.eta))
        if not out.accepted:
            reason=out.reason
            break
        x,v=out.x,out.v
        accepted_steps=step
        converged+=int(out.converged)
        path+=float(np.linalg.norm(x-old))
        sumsq+=float(np.sum((x-positions[step])**2))
        if not scene.rigid and np.linalg.norm(scene.g)==0:
            current_angular=scene.m*(x[0]*v[1]-x[1]*v[0])
            angular=max(angular,abs(current_angular-initial_angular))
    summary=dict(run_id=run_id,case=scene.name,order=order,hz=hz,sweeps=sweeps,
                 epsilon_ratio=ratio,epsilon=epsilon,duration=duration,
                 completed=int(accepted_steps==total and not reason),accepted_time=accepted_steps*h,
                 reason=reason,max_rel_energy=maxe,final_rel_energy=rows[-1]['rel_energy'],
                 rms_position=np.sqrt(sumsq/accepted_steps) if accepted_steps else np.nan,
                 path_length=path,max_geometry_error=geometry,max_tangency_error=tangency,
                 max_residual=residual,converged_fraction=converged/max(accepted_steps,1),
                 max_block_displacement=blockjump,max_abs_rho=maxrho,max_abs_eta=maxeta,
                 max_angular_error=angular,seconds=time.perf_counter()-start,
                 visits=visits,rank_defect_visits=rank_visits,
                 energy_scale=energy_scale,initial_energy=initial_energy,
                 energy_floor=energy_floor(scene),target_policy=target_policy)
    return summary,rows,first_trace


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--duration',type=float,default=5.)
    parser.add_argument('--hz',default='60,120')
    parser.add_argument('--sweeps',default='1,4,16')
    parser.add_argument('--epsilon-ratios',default='1,100,10000,100000000')
    parser.add_argument('--orders',default='EL,LE,JOINT,XPBD')
    parser.add_argument('--cases',default=','.join(scene.name for scene in scenes()))
    parser.add_argument('--output',default='results')
    args=parser.parse_args()
    destination=Path(args.output)
    destination.mkdir(parents=True,exist_ok=True)
    summaries=[]
    traces={}
    references={}
    reference_audit=[]
    scene_list=[scene for scene in scenes() if scene.name in args.cases.split(',')]
    for scene in scene_list:
        for hz in map(int,args.hz.split(',')):
            times=np.arange(round(args.duration*hz)+1)/hz
            reference=reference_trajectory(scene,times)
            tighter=reference_trajectory(scene,times,tighter=True)
            delta=float(np.max(np.linalg.norm(reference[0]-tighter[0],axis=1)))
            energy0=energy(scene,scene.x0,scene.v0)
            energy_scale=max(energy0-energy_floor(scene),1e-12)
            drift=max(abs(energy(scene,xx,vv)-energy0)/energy_scale
                      for xx,vv in zip(*reference))
            if delta>1e-8 or drift>1e-8:
                raise RuntimeError(f'reference audit failed: {scene.name}, {delta}, {drift}')
            reference_audit.append(dict(case=scene.name,hz=hz,max_position_refinement=delta,
                                        max_relative_energy=drift))
            for sweeps in map(int,args.sweeps.split(',')):
                for order in args.orders.split(','):
                    ratios=[0.] if order=='XPBD' else map(float,args.epsilon_ratios.split(','))
                    for ratio in ratios:
                        summary,rows,trace=run_one(scene,order,hz,sweeps,ratio,args.duration,reference)
                        summaries.append(summary)
                        traces[summary['run_id']]=trace
                        write_csv(destination/'trajectories'/f"{summary['run_id']}.csv",rows)
                        print(f"{summary['run_id']} complete={summary['completed']} "
                              f"time={summary['accepted_time']:.3f} E={summary['max_rel_energy']:.4g} "
                              f"rms={summary['rms_position']:.4g} conv={summary['converged_fraction']:.3f} "
                              f"{summary['reason']}",flush=True)
                        if len(summaries)%10==0:
                            write_csv(destination/'summary.csv',summaries)
    write_csv(destination/'summary.csv',summaries)
    write_csv(destination/'reference_audit.csv',reference_audit)
    (destination/'first_step_traces.json').write_text(json.dumps(traces,indent=2))
    metadata=dict(command=sys.argv,python=platform.python_version(),numpy=np.__version__,
                  scipy=scipy.__version__,target='incoming total physical energy, fixed during each step',
                  refresh='after every full block',stationarity='retained in every block',
                  epsilon='epsilon_ratio times initial energy above potential floor; fixed per run',
                  failure_policy='stop after guard/rank failure; offending candidate retained with accepted=0',
                  scenes=[dict(name=s.name,m=s.m,k=s.k,l=s.l,g=s.g.tolist(),
                               x0=s.x0.tolist(),v0=s.v0.tolist(),rigid=s.rigid) for s in scene_list])
    (destination/'metadata.json').write_text(json.dumps(metadata,indent=2))
    print(f'Wrote {len(summaries)} runs to {destination}',flush=True)


if __name__=='__main__':
    main()
