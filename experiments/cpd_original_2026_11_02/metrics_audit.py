"""Reconstruct saved CPD budgets and diagnostics without running a solver."""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import csv
import json
from pathlib import Path

import numpy as np

from validate import Case, preparation, evaluate

ROOT = Path(__file__).resolve().parent


class Audit:
    def __init__(self):
        self.count = 0
        self.maxima = {}

    def close(self, name, actual, expected, tolerance=3e-7):
        self.count += np.size(expected)
        a, b = np.asarray(actual, float), np.asarray(expected, float)
        with np.errstate(invalid='ignore',divide='ignore'):
            error = np.where((np.isnan(a)&np.isnan(b)) | (a==b), 0., abs(a-b)/np.maximum(1., abs(b)))
        value = float(np.max(error)) if error.size else 0.
        self.maxima[name] = max(self.maxima.get(name, 0.), value)
        assert np.isfinite(value) and value <= tolerance, (name, value, tolerance)

    def true(self, name, condition):
        self.count += 1
        assert condition, name


def csv_rows(path):
    with path.open(newline='') as stream:
        return list(csv.DictReader(stream))


def reference_case(scene, xn, vn, angular=True):
    n = len(scene['masses'])
    anchor = any(edge['b'] < 0 for edge in scene['edges'])
    x = np.vstack([xn, np.zeros(2)]) if anchor else xn.copy()
    v = np.vstack([vn, np.zeros(2)]) if anchor else vn.copy()
    m = np.r_[scene['masses'], 1.] if anchor else np.array(scene['masses'])
    free = np.r_[np.ones(n, bool), False] if anchor else np.ones(n, bool)
    forces = np.vstack([scene['forces'], np.zeros(2)]) if anchor else np.array(scene['forces'])
    edges = [(edge['a'], edge['b'] if edge['b'] >= 0 else n,
              0. if edge['hard'] else edge['stiffness'], edge['rest']) for edge in scene['edges']]
    return Case(x, v, m, free, edges, np.zeros(2), scene['family']=='free', angular,
                None if scene['family']=='free' else np.zeros(2), forces)


def trace_reference(scene, case, prepared, point):
    data = evaluate(case, prepared, point)
    # Angular-row ablations retain the original momentum kinetic reference.
    rows, values = [], []
    x = point[:-1].reshape(-1, 2)
    for edge in scene['edges']:
        if edge['hard']:
            delta = x[edge['a']]-(x[edge['b']] if edge['b'] >= 0 else 0.)
            radius = np.linalg.norm(delta)
            row = np.zeros(len(point))
            row[2*edge['a']:2*edge['a']+2] = delta/radius
            if edge['b'] >= 0:
                row[2*edge['b']:2*edge['b']+2] -= delta/radius
            rows.append(row); values.append(radius-edge['rest'])
    if rows:
        data['jacobian'] = np.vstack([data['jacobian'][:-1], rows, data['jacobian'][-1:]])
        data['constraints'] = np.r_[data['constraints'][:-1], values, data['constraints'][-1]]
    return data


def audit_campaign(source, audit):
    metadata = json.loads((source/'metadata.json').read_text())
    scenes = {scene['name']: scene for scene in metadata['scenes']}
    summaries = csv_rows(source/'summary.csv')
    traces = json.loads((source/'first_step_traces.json').read_text())
    counts = Counter(); families = defaultdict(Counter)
    total_rows = trace_visits = hard_unsigned = 0
    trace_diagnostics = defaultdict(lambda: dict(visits=0, max_dense_state_error=0.,
        max_dense_multiplier_error=0., max_constraint_linear_residual=0.,
        max_stationarity_linear_residual=0., ill_conditioned_visits=0, dense_singular_visits=0))
    worst = None
    for summary in summaries:
        scene = scenes[summary['case']]
        mass, forces = np.array(scene['masses']), np.array(scene['forces'])
        n, h, epsilon = len(mass), 1/int(summary['hz']), float(summary['epsilon'])
        scale = float(summary['energy_scale'])
        data = csv_rows(source/'trajectories'/(summary['run_id']+'.csv'))
        c = {key: np.array([float(row[key]) for row in data]) for key in data[0]}
        x = np.stack([np.column_stack([c[f'x_{i}'], c[f'y_{i}']]) for i in range(n)], axis=1)
        v = np.stack([np.column_stack([c[f'vx_{i}'], c[f'vy_{i}']]) for i in range(n)], axis=1)
        accepted = c['accepted'].astype(bool)
        rejected = np.flatnonzero(~accepted)
        audit.true('single retained rejected last candidate', len(rejected)<=1 and
                   (not len(rejected) or rejected[0]==len(data)-1))
        naccepted, attempted = int(accepted.sum())-1, len(data)-1
        total_rows += len(data)
        counts['runs'] += 1; counts['completed'] += int(summary['completed'])
        counts['failed'] += bool(len(rejected))
        families[scene['family']]['runs'] += 1
        families[scene['family']]['completed'] += int(summary['completed'])
        families[scene['family']]['failed'] += bool(len(rejected))
        kinetic = .5*np.sum(mass[None, :, None]*v*v, axis=(1,2))
        elastic = np.zeros(len(data)); gu = np.zeros_like(x)
        hard_geometry = np.zeros(len(data)); tangent = np.zeros(len(data))
        edge_lengths = []
        hard_gradient = np.zeros_like(x)
        for edge in scene['edges']:
            delta = x[:, edge['a']]-(x[:, edge['b']] if edge['b']>=0 else 0.)
            radius = np.linalg.norm(delta, axis=1); edge_lengths.append(radius)
            normal = delta/radius[:, None]
            if not edge['hard']:
                elastic += .5*edge['stiffness']*(radius-edge['rest'])**2
                g = edge['stiffness']*(radius-edge['rest'])[:, None]*normal
                gu[:, edge['a']] += g
                if edge['b']>=0: gu[:, edge['b']] -= g
            else:
                hard_geometry = np.maximum(hard_geometry, abs(radius-edge['rest']))
                dv = v[:,edge['a']]-(v[:,edge['b']] if edge['b']>=0 else 0.)
                tangent = np.maximum(tangent, abs(np.sum(delta*dv,axis=1))/edge['rest'])
                hard_gradient[:,edge['a']] += normal
                if edge['b']>=0: hard_gradient[:,edge['b']] -= normal
        internal = kinetic+elastic
        physical = internal-np.sum(forces[None,:,:]*x,axis=(1,2))
        actual_rel = (physical-physical[0])/scale
        paper_work = np.r_[0., h*np.sum(forces[None,:,:]*v[:-1],axis=(1,2))]
        exact_work = np.r_[0., np.sum(forces[None,:,:]*np.diff(x,axis=0),axis=(1,2))]
        origin = (np.sum(mass[None,:,None]*x[:-1],axis=1)/mass.sum()
                  if scene['family']=='free' else np.zeros((attempted,2)))
        old_relative = x[:-1]-origin[:,None,:]
        impulse = mass[None,:,None]*v[:-1]+h*forces[None,:,:]
        pstar = impulse.sum(axis=1)
        lstar = np.sum(old_relative[:,:,0]*impulse[:,:,1]-old_relative[:,:,1]*impulse[:,:,0],axis=1)
        inertia = np.sum(mass[None,:]*np.sum(old_relative**2,axis=2),axis=1)
        kstar = lstar*lstar/(2*inertia)
        if scene['family']=='free': kstar += np.sum(pstar*pstar,axis=1)/(2*mass.sum())
        target = internal[:-1]+paper_work[1:]
        rho = c['rho'][1:]
        relaxed = (1-rho)*target+rho*kstar
        defect = internal[1:]-relaxed
        relaxation = -rho*(target-kstar)
        cum_paper, cum_exact = np.cumsum(paper_work), np.cumsum(exact_work)
        cum_relax = np.r_[0., np.cumsum(relaxation)]
        cum_defect = np.r_[0., np.cumsum(defect)]
        budget = internal[0]+cum_paper
        relaxed_budget = budget+cum_relax
        computed = dict(kinetic=kinetic,elastic=elastic,internal_energy=internal,energy=physical,
            rel_energy=actual_rel,algorithm_target=np.r_[internal[0],target],kstar=np.r_[0.,kstar],
            relaxed_target=np.r_[internal[0],relaxed],algorithm_target_defect=np.r_[0.,internal[1:]-target],
            relaxed_target_defect=np.r_[0.,defect],paper_work=paper_work,exact_external_work=exact_work,
            cumulative_paper_work=cum_paper,cumulative_exact_work=cum_exact,cumulative_relaxation=cum_relax,
            cumulative_constraint_defect=cum_defect,algorithm_budget=budget,
            relaxed_algorithm_budget=relaxed_budget,algorithm_budget_defect=internal-budget,
            relaxed_budget_defect=internal-relaxed_budget,work_approximation_defect=cum_paper-cum_exact,
            mechanical_budget_identity_error=physical-physical[0]-(cum_paper-cum_exact+cum_relax+cum_defect),
            geometry=hard_geometry,tangency=tangent,min_edge_length=np.min(edge_lengths,axis=0),
            max_edge_length=np.max(edge_lengths,axis=0),
            nominal_length_error=np.max([abs(z-edge['rest']) for z,edge in zip(edge_lengths,scene['edges'])],axis=0))
        for field, values in computed.items(): audit.close('trajectory '+field,c[field],values)
        center=np.sum(mass[None,:,None]*x,axis=1)/mass.sum()
        center_v=np.sum(mass[None,:,None]*v,axis=1)/mass.sum()
        relative_x,relative_v=x-center[:,None,:],v-center_v[:,None,:]
        angular_global=np.sum(mass[None,:]*(x[:,:,0]*v[:,:,1]-x[:,:,1]*v[:,:,0]),axis=1)
        angular_internal=np.sum(mass[None,:]*(relative_x[:,:,0]*relative_v[:,:,1]-relative_x[:,:,1]*relative_v[:,:,0]),axis=1)
        if n<=2:
            mode_x=x[:,0] if n==1 else x[:,0]-x[:,1]
            mode_v=v[:,0] if n==1 else v[:,0]-v[:,1]
            reduced=mass[0] if n==1 else mass[0]*mass[1]/mass.sum()
            mode_length=np.linalg.norm(mode_x,axis=1)
            radial_speed=np.sum(mode_x*mode_v,axis=1)/mode_length
            radial_k=.5*reduced*radial_speed**2
            relative_k=.5*reduced*np.sum(mode_v*mode_v,axis=1)
            tangent_k=np.maximum(0.,relative_k-radial_k)
            rotational=tangent_k;shape=radial_k
        else:
            radii=np.linalg.norm(relative_x,axis=2)
            rs=np.divide(np.sum(relative_x*relative_v,axis=2),radii,out=np.zeros_like(radii),where=radii>0)
            radial_k=.5*np.sum(mass[None,:]*rs*rs,axis=1)
            relative_k=.5*np.sum(mass[None,:,None]*relative_v*relative_v,axis=(1,2))
            tangent_k=np.maximum(0.,relative_k-radial_k)
            inertia_internal=np.sum(mass[None,:]*radii*radii,axis=1)
            rotational=.5*angular_internal**2/inertia_internal
            shape=np.maximum(0.,relative_k-rotational)
            mode_length=np.sqrt(inertia_internal/mass.sum())
        modes=dict(angular_global=angular_global,angular_internal=angular_internal,
                   center_x=center[:,0],center_y=center[:,1],center_vx=center_v[:,0],center_vy=center_v[:,1],
                   linear_px=np.sum(mass[None,:]*v[:,:,0],axis=1),linear_py=np.sum(mass[None,:]*v[:,:,1],axis=1),
                   relative_kinetic=relative_k,radial_kinetic=radial_k,tangential_kinetic=tangent_k,
                   rotational_min_kinetic=rotational,shape_kinetic=shape,mode_length=mode_length)
        for field,values in modes.items():audit.close('trajectory body mode '+field,c[field],values)
        audit.close('derived backward velocity',v[1:],np.diff(x,axis=0)/h)
        audit.close('initial mechanical energy',float(summary['initial_energy']),physical[0])
        audit.close('initial internal energy',float(summary['initial_internal_energy']),internal[0])
        audit.close('initial mechanical energy scale',scale,max(physical[0]-scene['energy_floor'],1e-12))
        audit.close('fixed dimensional epsilon',epsilon,float(summary['epsilon_ratio'])*scale)

        # Fresh endpoint multiplier stationarity, with the hard multiplier sign
        # limitation kept explicit instead of inferring a signed multiplier.
        row_constraints = []
        row_scales = []
        if scene['family']=='free':
            pactual = np.sum(mass[None,:,None]*v[1:],axis=1)
            row_constraints.extend([(pactual-pstar)[:,0],(pactual-pstar)[:,1]])
            row_scales += [np.full(attempted,mass.sum()*scene['characteristic_length']/h)]*2
        lactual = np.sum((x[1:,:,0]-origin[:,None,0])*mass[None,:]*v[1:,:,1]
                        -(x[1:,:,1]-origin[:,None,1])*mass[None,:]*v[1:,:,0],axis=1)
        angular_on = '_noL' not in summary['run_id']
        if angular_on:
            row_constraints.append(lactual-lstar)
            row_scales.append(np.full(attempted,mass.sum()*scene['characteristic_length']**2/h))
        for edge,radius in zip(scene['edges'],edge_lengths):
            if edge['hard']:
                row_constraints.append(radius[1:]-edge['rest'])
                row_scales.append(np.full(attempted,scene['characteristic_length']))
        row_constraints.append(defect)
        row_scales.append(np.maximum(np.maximum(abs(target),abs(kstar)),1e-12))
        raw_constraints=np.max(abs(np.array(row_constraints)),axis=0)
        scaled_constraints=np.max(abs(np.array(row_constraints))/np.array(row_scales),axis=0)
        predictor=x[:-1]+h*v[:-1]+h*h*forces[None,:,:]/mass[None,:,None]
        rx=mass[None,:,None]*(x[1:]-predictor)/h**2+gu[1:]
        gx=mass[None,:,None]*(x[1:]-x[:-1])/h**2+gu[1:]
        rx += c['eta'][1:,None,None]*gx
        rx += c['angular_multiplier'][1:,None,None]*mass[None,:,None]/h*np.stack([-old_relative[:,:,1],old_relative[:,:,0]],axis=2)
        rx += mass[None,:,None]/h*np.stack([c['linear_multiplier_x'][1:],c['linear_multiplier_y'][1:]],axis=1)[:,None,:]
        rrho=epsilon*rho+c['eta'][1:]*(target-kstar)
        rho_station=abs(rrho)/np.maximum(np.maximum(epsilon,abs(target)),np.maximum(abs(kstar),1e-12))
        position_station=np.max(abs(rx)*h*h/(mass[None,:,None]*scene['characteristic_length']),axis=(1,2))
        stationarity=np.maximum(position_station,rho_station)
        if any(edge['hard'] for edge in scene['edges']):
            hard_unsigned += attempted
            possible=[]
            for sign in (-1,1):
                trial=rx+sign*c['geometry_multiplier_max'][1:,None,None]*hard_gradient[1:]
                possible.append(np.maximum(rho_station,np.max(abs(trial)*h*h/(mass[None,:,None]*scene['characteristic_length']),axis=(1,2))))
            distances=np.array([abs(z-c['stationarity'][1:]) for z in possible])
            stationarity=np.where(distances[0]<=distances[1],possible[0],possible[1])
            audit.close('unsigned hard-multiplier endpoint stationarity consistency only',c['stationarity'][1:],stationarity)
        else:
            audit.close('independent fresh endpoint multiplier stationarity',c['stationarity'][1:],stationarity)
        audit.close('independent endpoint raw constraints',c['raw_constraints'][1:],raw_constraints)
        audit.close('independent endpoint scaled constraints',c['scaled_constraints'][1:],scaled_constraints)
        audit.close('independent endpoint KKT scale',c['kkt_scaled'][1:],np.maximum(stationarity,scaled_constraints))
        audit.close('incoming-origin angular target error',c['angular_target_error'][1:],abs(lactual-lstar))
        pactual=np.sum(mass[None,:,None]*v[1:],axis=1)
        audit.close('incoming linear target error',c['linear_target_error'][1:],np.linalg.norm(pactual-pstar,axis=1))

        last=np.flatnonzero(accepted)[-1]
        jump=np.linalg.norm(np.diff(x[accepted].reshape(accepted.sum(),-1),axis=0),axis=1)
        tail=c['time'][accepted][1:]>max(0.,c['time'][last]-1.)+1e-10
        stall=longest=0.
        for z in jump:
            stall=stall+h if z/h<1e-4*scene['characteristic_length'] else 0.
            longest=max(longest,stall)
        expected=dict(accepted_time=naccepted*h,
            completed=int(naccepted==round(float(summary['duration'])/h) and not len(rejected)),
            max_rel_energy=np.max(abs(actual_rel)),max_accepted_rel_energy=np.max(abs(actual_rel[accepted])),
            final_rel_energy=actual_rel[last],path_length=jump.sum(),tail_path_length=jump[tail].sum(),
            tail_stalled=int(bool(len(jump[tail])) and jump[tail].sum()<1e-4*scene['characteristic_length']),
            longest_stall_duration=longest,max_accepted_displacement=jump.max() if len(jump) else 0.,
            converged_fraction=c['converged'][1:][accepted[1:]].sum()/max(naccepted,1),
            stationary_fraction=c['stationary'][1:][accepted[1:]].sum()/max(naccepted,1),
            solver_iterations=c['iterations'][1:].sum(),mean_iterations=c['iterations'][1:].sum()/max(attempted,1),
            min_schur_rank=np.min(c['schur_rank'][1:]),max_condition=np.max(c['condition'][1:]),
            rank_defect_steps=int(np.count_nonzero(c['schur_rank'][1:]<c['conservation_rows'][1:])) if summary['method']=='CPD' else 0,
            max_abs_rho=abs(c['rho'][accepted]).max(),min_rho=c['rho'][accepted].min(),max_rho=c['rho'][accepted].max(),
            max_abs_eta=abs(c['eta'][accepted]).max(),
            max_geometry_error=hard_geometry[accepted].max(),max_tangency_error=tangent[accepted].max(),
            initial_angular=angular_global[0],final_angular=angular_global[last],
            max_angular_drift=np.max(abs(angular_global[accepted]-angular_global[0])),
            final_radial_kinetic=radial_k[last],final_tangential_kinetic=tangent_k[last],
            final_rotational_kinetic=rotational[last],final_shape_kinetic=shape[last])
        for output,field in [('max_algorithm_target_defect','algorithm_target_defect'),
            ('max_relaxed_target_defect','relaxed_target_defect'),('max_algorithm_budget_defect','algorithm_budget_defect'),
            ('max_relaxed_budget_defect','relaxed_budget_defect'),('max_work_approximation_defect','work_approximation_defect'),
            ('max_mechanical_budget_identity_error','mechanical_budget_identity_error')]:
            expected[output]=np.max(abs(computed[field][accepted]))/scale
        for output,field in [('max_raw_constraints','raw_constraints'),('max_scaled_constraints','scaled_constraints'),
                             ('max_stationarity','stationarity'),('max_kkt_scaled','kkt_scaled'),
                             ('max_linear_target_error','linear_target_error'),('max_angular_target_error','angular_target_error')]:
            expected[output]=c[field][accepted].max()
        for field,value in expected.items(): audit.close('summary '+field,float(summary[field]),value)
        audit.true('failure reason matches rejected candidate',bool(summary['reason'])==bool(len(rejected)))

        # Per-visit dense KKT verification uses the matching previous saved state,
        # preventing ill-conditioned nonlinear trajectory drift from obscuring
        # whether each fresh multiplier solves its actual quadratic subproblem.
        point=np.r_[(np.array(scene['x0'])+h*np.array(scene['v0'])+h*h*forces/mass[:,None]).ravel(),0.]
        case=reference_case(scene,np.array(scene['x0']),np.array(scene['v0']),angular=angular_on)
        kcase=reference_case(scene,np.array(scene['x0']),np.array(scene['v0']),angular=True)
        prepared=preparation(kcase,h,epsilon)
        for snapshot in traces.get(summary['run_id'],[]):
            reference=trace_reference(scene,case,prepared,point)
            a,j=prepared['metric'],reference['jacobian']
            newpoint=np.r_[np.array(snapshot['x']).ravel(),snapshot['rho']]
            delta=newpoint-point
            mu=np.array(snapshot['multipliers'])
            scales=np.sqrt(np.maximum(np.diag(j@np.linalg.solve(a,j.T)),0.))
            scaling=np.divide(1.,scales,out=np.ones_like(scales),where=scales>0)
            s=scaling[:,None]*(j@np.linalg.solve(a,j.T))*scaling[None,:]
            condition=np.linalg.cond(s)
            group=trace_diagnostics[scene['family']+'/'+summary['method']]
            group['visits']+=1
            group['ill_conditioned_visits']+=int(condition>1e10)
            norm=max(1.,np.max(abs(reference['gradient'])),np.max(abs(reference['constraints'])))
            cg=float(np.max(abs(reference['constraints']+j@delta))/norm)
            sg=float(np.max(abs(reference['gradient']+a@delta+j.T@mu))/norm)
            group['max_constraint_linear_residual']=max(group['max_constraint_linear_residual'],cg)
            group['max_stationarity_linear_residual']=max(group['max_stationarity_linear_residual'],sg)
            audit.close('serialized fresh objective gradient',snapshot['objective_gradient'],trace_reference(scene,case,prepared,newpoint)['gradient'])
            audit.close('serialized refreshed true conservation constraints',snapshot['constraints'],trace_reference(scene,case,prepared,newpoint)['constraints'])
            try:
                if summary['method']=='PD':
                    dense=np.r_[-np.linalg.solve(a,reference['gradient']),np.zeros(len(j))]
                else:
                    dense=np.linalg.solve(np.block([[a,j.T],[j,np.zeros((len(j),len(j)))]]),
                                          -np.r_[reference['gradient'],reference['constraints']])
                de=float(np.linalg.norm(delta-dense[:len(point)])/max(1.,np.linalg.norm(dense[:len(point)])))
                me=float(np.linalg.norm(mu-dense[len(point):])/max(1.,np.linalg.norm(dense[len(point):])))
                group['max_dense_state_error']=max(group['max_dense_state_error'],de)
                group['max_dense_multiplier_error']=max(group['max_dense_multiplier_error'],me)
                if condition<1e10:
                    audit.close('conditioned serialized dense global correction',delta,dense[:len(point)],3e-5)
                    audit.close('conditioned serialized NEW multiplier',mu,dense[len(point):],3e-5)
            except np.linalg.LinAlgError:
                group['dense_singular_visits']+=1
            point=newpoint; trace_visits+=1

        if '_compliant' in summary['case'] and (worst is None or float(summary['max_accepted_rel_energy'])>worst['max_accepted_rel_energy']):
            imax=np.flatnonzero(accepted)[np.argmax(abs(actual_rel[accepted]))]
            worst=dict(run_id=summary['run_id'],completed=int(summary['completed']),
                time=float(c['time'][imax]),max_accepted_rel_energy=float(abs(actual_rel[imax])),
                mechanical_energy_change=float(physical[imax]-physical[0]),energy_scale=scale,
                work_approximation=float((cum_paper-cum_exact)[imax]),
                relaxation=float(cum_relax[imax]),constraint_defect=float(cum_defect[imax]),
                identity_error=float(computed['mechanical_budget_identity_error'][imax]))
    return dict(source=str(source.resolve()),counts={k:int(v) for k,v in counts.items()},
        families={f:{k:int(v) for k,v in cc.items()} for f,cc in families.items()},trajectory_rows=total_rows,
        trace_visits=trace_visits,trace_diagnostics=dict(trace_diagnostics),
        hard_unsigned_endpoint_rows=hard_unsigned,worst_accepted_compliant_pendulum=worst)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--results',type=Path,default=ROOT/'results')
    parser.add_argument('--all-campaigns',action='store_true')
    args=parser.parse_args()
    audit=Audit()
    sources=[args.results]
    if args.all_campaigns:
        sources.extend(path for path in (ROOT/'pd_control',ROOT/'angular_ablation') if (path/'metadata.json').exists())
    campaigns=[audit_campaign(source,audit) for source in sources]
    receipt=dict(status='passed',checks=int(audit.count),campaigns=campaigns,
                 largest_scaled_errors=audit.maxima,numerical_solver_executed=False)
    (ROOT/'METRICS_AUDIT.json').write_text(json.dumps(receipt,indent=2,sort_keys=True))
    (ROOT/'METRICS_AUDIT.md').write_text('''# Independent saved-data audit

Mechanical energy, paper internal energy, incoming-force work targets, kinetic
references, relaxed energy rows, and their cumulative budget decomposition are
reconstructed from saved state arrays without invoking a production solver.
The physical-energy change is start-point work minus exact external work, plus
energy relaxation and accumulated constraint defects. These terms are kept
separate; a small relaxed-row residual is not automatically physical-energy
conservation.

Counts are split by family. Full-free-body finite-spring tests are distinct from
anchored adaptations and exact hard-rod extensions. Finite completion is distinct
from paper conservation convergence and full KKT stationarity. Failed candidates
remain in all-trial energy maxima; accepted-prefix maxima and final accepted
energies exclude them. Mean iterations use attempted steps, including rejection.

Every finite-spring endpoint stationarity norm is reconstructed using the signed
fresh conservation multipliers. Hard-rod endpoint CSVs retain only the magnitude
of their geometry multiplier: both possible signs are checked for consistency,
without claiming a uniquely reconstructed signed stationarity vector. Serialized
first-step traces retain full signed multipliers and are checked directly.

Dense first-step checks use each actual previous saved state. Ill-conditioned
Schur systems above condition 1e10 retain their direct dense discrepancies as
diagnostics rather than treating rounded dense solutions as an exact oracle.
No pseudoinverse, row repair, damping, or numerical rollout is used by this audit.

```bash
OPENBLAS_NUM_THREADS=1 python3 metrics_audit.py --all-campaigns
```

The JSON receipt includes per-family linear residuals and dense comparisons,
all primary/control counts, and the worst accepted compliant-pendulum budget.
Internal block maxima, rank-detection thresholds, and CPU timings cannot be
fully reconstructed from endpoint CSVs and remain instrumented diagnostics.

```json
'''+json.dumps(receipt,indent=2,sort_keys=True)+'\n```\n')
    print(json.dumps(receipt,indent=2,sort_keys=True))


if __name__=='__main__': main()
