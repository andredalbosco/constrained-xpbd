"""Recompute saved rollout metrics, replay traces, and compare old references.

Run after run_experiments.py has finished writing this experiment's results.
No result files are modified. numpy and scipy are required.
"""
import csv
import importlib.util
import json
from pathlib import Path
import sys
import numpy as np

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))
import solver
import reference
from validate import independent

summaries = list(csv.DictReader((ROOT/'results/summary.csv').open()))
traces = json.loads((ROOT/'results/first_step_traces.json').read_text())
scene_map = {scene.name: scene for scene in solver.scenes()}
maxima = {}
checks = 0
failed = completed = finite_nonconverged = 0
trajectory_rows = trace_visits = 0

def close(name, a, b, tolerance=2e-9):
    global checks
    checks += 1
    if np.isnan(a) and np.isnan(b):
        return
    error=abs(a-b)/max(1.0, abs(b))
    maxima[name]=max(maxima.get(name,0), error)
    assert np.isfinite(error) and error<=tolerance, (name,a,b,error)

for summary in summaries:
    run=summary['run_id']
    scene=scene_map[summary['case']]
    hz,sweeps=int(summary['hz']),int(summary['sweeps'])
    epsilon=float(summary['epsilon'])
    h=1/hz
    rows=list(csv.DictReader((ROOT/'results/trajectories'/f'{run}.csv').open()))
    accepted=[row for row in rows[1:] if int(row['accepted'])]
    trajectory_rows+=len(rows)
    rejected=[row for row in rows[1:] if not int(row['accepted'])]
    assert int(rows[0]['accepted']) == 1
    assert len(rejected)<=1 and (not rejected or rejected[0] is rows[-1])
    assert len(rows)==len(accepted)+len(rejected)+1
    initial=float(summary['initial_energy'])
    floor=float(summary['energy_floor'])
    scale=float(summary['energy_scale'])
    close('initial physical energy',initial,solver.energy(scene,scene.x0,scene.v0))
    close('normalization above potential floor',scale,max(initial-floor,1e-12))
    close('fixed dimensional epsilon',epsilon,(1 if summary['order']=='XPBD' else float(summary['epsilon_ratio']))*scale)
    close('accepted prefix time',float(summary['accepted_time']),len(accepted)*h)
    expected_completed=int(len(accepted)==round(float(summary['duration'])*hz) and not rejected)
    assert int(summary['completed']) == expected_completed
    completed+=expected_completed
    failed+=len(rejected)
    if expected_completed and float(summary['converged_fraction'])<1:
        finite_nonconverged+=1
    close('accepted prefix convergence fraction',float(summary['converged_fraction']),
          sum(int(row['converged']) for row in accepted)/max(len(accepted),1))
    squared=0.0
    path=0.0
    old=scene.x0.copy()
    initial_angular=scene.m*(scene.x0[0]*scene.v0[1]-scene.x0[1]*scene.v0[0])
    angular=0.0
    for row in accepted:
        x=np.array([float(row['x']),float(row['y'])])
        v=np.array([float(row['vx']),float(row['vy'])])
        refx=np.array([float(row['ref_x']),float(row['ref_y'])])
        squared+=float(np.sum((x-refx)**2))
        path+=float(np.linalg.norm(x-old))
        close('returned displacement velocity',float(np.linalg.norm(v-(x-old)/h)),0.0,2e-7)
        old=x
        if not scene.rigid and np.linalg.norm(scene.g)==0:
            angular=max(angular,abs(scene.m*(x[0]*v[1]-x[1]*v[0])-initial_angular))
    close('accepted prefix RMS',float(summary['rms_position']),np.sqrt(squared/len(accepted)) if accepted else np.nan)
    close('accepted prefix path length',float(summary['path_length']),path)
    close('accepted prefix angular error',float(summary['max_angular_error']),angular)
    keys={'max_rel_energy':('rel_energy',True),
          'max_geometry_error':('geometry',False),
          'max_tangency_error':('tangency',False),
          'max_residual':('residual',False),
          'max_abs_rho':('rho',True),
          'max_abs_eta':('eta',True)}
    for out_key,(in_key,absolute) in keys.items():
        values=[float(row[in_key]) for row in rows]
        expected=max(abs(value) for value in values) if absolute else max(values)
        close(f'all-row {out_key}',float(summary[out_key]),expected)
    close('last candidate final energy',float(summary['final_rel_energy']),float(rows[-1]['rel_energy']))
    for row in rows:
        x=np.array([float(row['x']),float(row['y'])])
        v=np.array([float(row['vx']),float(row['vy'])])
        physical=solver.energy(scene,x,v)
        close('CSV physical energy',float(row['energy']),physical)
        close('CSV initial-relative energy',float(row['rel_energy']),(physical-initial)/scale)
        close('CSV geometry',float(row['geometry']),abs(np.linalg.norm(x)-scene.l) if scene.rigid else 0.0)
        close('CSV tangency',float(row['tangency']),abs(float(x@v))/scene.l if scene.rigid else 0.0)
    if expected_completed:
        expected_visits=round(float(summary['duration'])*hz)*sweeps*(2 if summary['order'] in ('EL','LE') else 1)
        assert int(summary['visits'])==expected_visits
    if rejected:
        assert summary['reason']
    else:
        assert not summary['reason']
    # Independently replay the stored first step, using dense block systems and
    # fresh analytical derivatives. This additionally audits serialization.
    point=np.r_[scene.x0+h*scene.v0,0.0,0.0,0.0]
    blocks={'EL':'EL','LE':'LE','JOINT':'J','XPBD':'L'}[summary['order']]
    target=solver.energy(scene,scene.x0,scene.v0)
    args=(scene.x0,scene.v0,h,scene.m,scene.k,scene.l,scene.g,epsilon,target,floor)
    for visit,snapshot in enumerate(traces[run]):
        block=blocks[visit%len(blocks)]
        assert snapshot['block']==block and snapshot['sweep']==visit//len(blocks)
        evaluation=independent(point[:2],point[2],point[3],point[4],*args,rigid=scene.rigid)
        ids=[0,1,2] if block=='L' else ([0,1,3,4] if block=='E' else [0,1,2,3,4])
        delta=np.linalg.solve(evaluation['matrix'][np.ix_(ids,ids)],-evaluation['residual'][ids])
        point[ids]+=delta
        actual=np.r_[snapshot['x'],snapshot['lam'],snapshot['rho'],snapshot['eta']]
        for a,b in zip(actual,point):
            close('serialized first-step dense replay',float(a),float(b),2e-6)
        updated=independent(point[:2],point[2],point[3],point[4],*args,rigid=scene.rigid)
        for a,b in zip(np.r_[snapshot['r_x'],snapshot['r_lambda'],snapshot['r_rho'],snapshot['r_eta']],updated['residual']):
            close('serialized trace physical residuals',float(a),float(b),2e-5)
        trace_visits+=1

source_candidates = (
    ROOT.parent/'nonlinear_energy_2026_11_01'/'references.py',
    Path('/home/blacknoir/Projects/constrained-xpbd/experiments/nonlinear_energy_2026_11_01/references.py'),
    ROOT.parent/'constrained-xpbd-nonlinear-2026-11-01'/'references.py',
)
old_source = next((candidate for candidate in source_candidates if candidate.is_file()), None)
if old_source is None:
    raise FileNotFoundError('Earlier verified references.py unavailable; tried: '
                            + ', '.join(str(candidate) for candidate in source_candidates))
spec = importlib.util.spec_from_file_location('earlier_verified_reference', old_source)
old_reference = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = old_reference
spec.loader.exec_module(old_reference)
reference_comparisons = []
for scene in solver.scenes():
    for hz in (60, 120):
        times = np.arange(int(.5*hz)+1)/hz
        current = reference.reference_trajectory(scene, times)
        previous = old_reference.reference_trajectory(
            'pendulum' if scene.rigid else 'spring', scene.x0, scene.v0,
            times, scene.m, scene.k, scene.l, scene.g)
        assert previous.success
        position = float(np.max(np.linalg.norm(current[0]-previous.x, axis=1)))
        velocity = float(np.max(np.linalg.norm(current[1]-previous.v, axis=1)))
        close('prior verified reference position', position, 0.0, 1e-10)
        close('prior verified reference velocity', velocity, 0.0, 1e-9)
        reference_comparisons.append(dict(scene=scene.name, hz=hz,
                                         max_position_difference=position,
                                         max_velocity_difference=velocity))

result=dict(checks=checks,runs=len(summaries),trajectory_rows=trajectory_rows,
            completed=completed,failed=failed,completed_with_nonconverged_steps=finite_nonconverged,
            serialized_trace_visits=trace_visits,maxima=maxima,
            prior_reference_source=str(old_source),reference_comparisons=reference_comparisons)
print(json.dumps(result,indent=2,sort_keys=True))
