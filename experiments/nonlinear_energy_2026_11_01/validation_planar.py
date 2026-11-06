#!/usr/bin/env python3
"""Independent numerical/formulation checks for planar.py; no source edits."""
import csv
import math
from pathlib import Path
import numpy as np
import planar as p
import scalar_exact as s

OUT = Path(__file__).parent / 'planar_validation'
OUT.mkdir(exist_ok=True)
rng = np.random.default_rng(20261101)
rows = []
metrics = {}


def record(name, error, tolerance, details=''):
    error, tolerance = float(error), float(tolerance)
    rows.append(dict(check=name, error=error, tolerance=tolerance,
                     passed=bool(error <= tolerance), details=details))
    metrics[name] = max(metrics.get(name, 0.0), error)


def randomized_surface_checks():
    for case in range(64):
        angle = rng.uniform(-math.pi, math.pi)
        old = rng.uniform(.7, 1.5) * np.array([math.cos(angle), math.sin(angle)])
        v = rng.normal(size=2) * rng.uniform(0, 2)
        h = 1 / float(rng.choice([30, 60, 120, 480]))
        m, k, l = rng.uniform(.5, 2), rng.uniform(10, 300), rng.uniform(.6, 1.2)
        g = np.array([0., -9.81]) if case % 2 else np.zeros(2)
        target = p.energy(old, v, m, k, l, g)
        scale = max(1., abs(target), m * np.dot(v, v), k * (np.linalg.norm(old) - l)**2)
        curve = p.EnergyCurve(old, h, m, k, l, g, target)
        for phi in np.linspace(0, 2 * math.pi, 33):
            x = curve.point(float(phi))
            error = abs(p.endpoint_energy(x, old, h, m, k, l, g) - target) / scale
            record('radial_contour_true_energy', error, 2e-10, str(case))
        if case < 32:
            queries = [old, old + h * v + h*h*g,
                       curve.rmin * curve.axis,
                       old + h * rng.normal(size=2),
                       old + .3 * rng.normal(size=2)]
            for query in queries:
                x32, _ = curve.closest(query, old, 32)
                x128, _ = curve.closest(query, old, 128)
                f32 = m * float(np.sum((x32 - query)**2))
                f128 = m * float(np.sum((x128 - query)**2))
                record('grid32_vs128_objective', max(0, f32-f128), 2e-10*max(1.,f128), str(case))
                record('closest_true_energy', abs(p.endpoint_energy(x32,old,h,m,k,l,g)-target)/scale, 2e-10, str(case))


def rod_checks():
    for case in range(64):
        theta, omega = rng.uniform(-2, 2), (0. if case % 4 == 0 else rng.uniform(-3, 3))
        l, m = rng.uniform(.6, 1.6), rng.uniform(.5, 2)
        old = l * np.array([math.sin(theta), -math.cos(theta)])
        v = l * omega * np.array([math.cos(theta), math.sin(theta)])
        h = 1 / float(rng.choice([30, 60, 120, 480]))
        g = np.array([0., -9.81])
        target = p.energy(old,v,m,0,l,g)
        x, _ = p.rod_joint(old,v,h,m,l,g,target)
        record('rod_joint_length', abs(np.linalg.norm(x)-l), 2e-11, str(case))
        record('rod_joint_energy', abs(p.endpoint_energy(x,old,h,m,0,l,g)-target)/max(1.,abs(target)), 2e-10, str(case))
        if omega == 0:
            record('rod_rest_tie', np.linalg.norm(x-old), 1e-13, str(case))
        else:
            # Reflect the chosen circle intersection across the energy-circle
            # axis to obtain the only other feasible endpoint independently.
            c = old + h*h*g
            axis = c / np.linalg.norm(c)
            other = 2*np.dot(x,axis)*axis-x
            predictor = old+h*v+h*h*g
            delta = np.sum((x-predictor)**2)-np.sum((other-predictor)**2)
            record('rod_joint_objective', max(0,delta), 1e-11, str(case))


def pd_checks():
    for case in range(24):
        old=np.array([1.15, .1]) + .05*rng.normal(size=2)
        v=rng.normal(size=2)
        h=1 / float(rng.choice([60,120,480]))
        m,k,l=1.,100.,1.
        g=np.array([0.,-9.81]) if case%2 else np.zeros(2)
        target=p.energy(old,v,m,k,l,g)
        y=old+h*v+h*h*g
        projection=l*y/np.linalg.norm(y)
        A=m/h**2+k
        center=(m/h**2*y+k*projection)/A
        pd=p.step(old,v,h,m,k,l,g,target,'pd',1)
        # For a single radial spring, the implicit minimizer is analytic.
        analytic=((m/h**2*np.linalg.norm(y)+k*l)/A)*y/np.linalg.norm(y)
        record('pd_unconstrained_implicit_minimum',np.linalg.norm(pd.x-analytic),1e-12,str(case))
        exact=p.step(old,v,h,m,k,l,g,target,'pd_true',1)
        if not exact.converged:
            record('pd_true_feasible',1,0,exact.reason)
            continue
        curve=p.EnergyCurve(old,h,m,k,l,g,target)
        pts=np.array([curve.point(float(t)) for t in np.linspace(0,2*math.pi,1024,endpoint=False)])
        # Independently evaluate the frozen PD objective and compare all dense
        # feasible samples. The global solve uses the physical energy contour.
        objective=lambda x: .5*m*np.sum((x-y)**2)/h**2+.5*k*np.sum((x-projection)**2)
        dense=min(objective(x) for x in pts)
        record('pd_true_global_objective',max(0,objective(exact.x)-dense),2e-9,str(case))
        record('pd_true_feasible',abs(p.endpoint_energy(exact.x,old,h,m,k,l,g)-target),2e-9,str(case))
        surrogate=p.step(old,v,h,m,k,l,g,target,'pd_surrogate',1)
        if surrogate.converged:
            e_surrogate=.5*m*np.sum(((surrogate.x-old)/h)**2)+.5*k*np.sum((surrogate.x-projection)**2)-m*np.dot(g,surrogate.x)
            record('pd_surrogate_frozen_energy',abs(e_surrogate-target),2e-9,str(case))
            true_e=p.endpoint_energy(surrogate.x,old,h,m,k,l,g)
            record('pd_surrogate_overestimates_true',max(0,true_e-e_surrogate),2e-10,str(case))
            c=(m/h**2*old+m*g+k*projection)/A
            radius=np.linalg.norm(surrogate.x-c)
            pts=c+radius*np.stack([np.cos(np.linspace(0,2*math.pi,1024)),np.sin(np.linspace(0,2*math.pi,1024))],axis=1)
            dense=min(objective(x) for x in pts)
            record('pd_surrogate_global_objective',max(0,objective(surrogate.x)-dense),2e-9,str(case))


def scalar_equivalence():
    for hz in [120,480]:
        h=1/hz
        for initial_v in [0.,.8]:
            for target_mode in ['initial','previous_step']:
                for sweeps in [1,4,16]:
                    q,v=.25,initial_v
                    x=np.array([1.25,0.]); vv=np.array([initial_v,0.])
                    initial_energy=s.energy(q,v,1,100)
                    for step in range(2):
                        target=initial_energy if target_mode=='initial' else s.energy(q,v,1,100)
                        oldq,oldv=q,v
                        predictor=q+h*v
                        q=predictor
                        lam=0.
                        for _ in range(sweeps):
                            roots=s.energy_roots(oldq,h,target,1,100)
                            q,_,_=s.nearest_root(roots,q,predictor)
                            q,lam=s.material_visit(q,lam,h,1,100)
                        v=(q-oldq)/h
                        out=p.step(x,vv,h,1,100,1,np.zeros(2),target,'nonlinear_EL',sweeps)
                        label=f'hz={hz},v0={initial_v},target={target_mode},sweeps={sweeps},step={step+1}'
                        record('scalar_nonlinear_position',np.linalg.norm(out.x-np.array([1+q,0.])),2e-6,label)
                        record('scalar_nonlinear_velocity',np.linalg.norm(out.v-np.array([v,0.])),1e-3,label)
                        x,vv=out.x,out.v
    original=list(csv.DictReader((Path(__file__).parent/'scalar_results/historical_check.csv').open()))
    for hz in [60,120,240,480]:
        x,v=np.array([1.25,0.]),np.zeros(2)
        for step in [1,2]:
            out=p.step(x,v,1/hz,1,100,1,np.zeros(2),3.125,'raw_EL',1)
            row=next(r for r in original if int(r['hz'])==hz and int(r['step'])==step)
            record('raw_planar_historical_x',abs(out.x[0]-float(row['original_x'])),1e-7,f'{hz}/{step}')
            record('raw_planar_historical_v',abs(out.v[0]-float(row['original_v'])),1e-4,f'{hz}/{step}')
            x,v=out.x,out.v


def rod_historical_control():
    source=OUT/'rod_historical_check.csv'
    if not source.exists():
        raise RuntimeError('Compile and run rod_historical_check.cpp before validation')
    for row in csv.DictReader(source.open()):
        hz,step=int(row['hz']),int(row['step'])
        theta,omega,m,l=[float(row[key]) for key in ['theta','omega','mass','length']]
        if step==1:
            x=l*np.array([math.sin(theta),-math.cos(theta)])
            v=l*omega*np.array([math.cos(theta),math.sin(theta)])
            target=p.energy(x,v,m,0,l,np.array([0.,-9.81]))
            if omega==0:
                predictor=x+np.array([0.,-9.81])/hz**2
                record('rod_predictor_zero_energy_gradient',np.linalg.norm(p.gradient(predictor,x,1/hz,m,0,l,np.array([0.,-9.81]))),1e-14,str(row))
        out=p.step(x,v,1/hz,m,0,l,np.array([0.,-9.81]),target,'raw_EL',1,rigid=True)
        record('raw_rod_historical_x',np.linalg.norm(out.x-np.array([float(row['x']),float(row['y'])])),1e-9,str(row))
        record('raw_rod_historical_v',np.linalg.norm(out.v-np.array([float(row['vx']),float(row['vy'])])),1e-6,str(row))
        x,v=out.x,out.v


if __name__=='__main__':
    randomized_surface_checks()
    rod_checks()
    pd_checks()
    scalar_equivalence()
    rod_historical_control()
    with (OUT/'checks.csv').open('w',newline='') as f:
        writer=csv.DictWriter(f,['check','error','tolerance','passed','details'])
        writer.writeheader();writer.writerows(rows)
    failures=[r for r in rows if not r['passed']]
    print('checks',len(rows),'failures',len(failures))
    for key,value in metrics.items():print(key,value)
    for failure in failures[:20]:print('FAIL',failure)
    (OUT/'SUMMARY.md').write_text('# Independent planar validation\n\n'+f'{len(rows)} checks; {len(failures)} failures.\n\n'+
      'The nonlinear E visit is a numerical global closest-point search on the exact physical energy contour, using the single-particle mass metric. It is not a root along a frozen energy-gradient line.\n\n'+
      '| Check | Maximum measured error |\n|---|---:|\n'+''.join(f'| {key} | {value:.8g} |\n' for key,value in metrics.items())+
      '\nAll randomized cases are within the explicitly supported contour component that does not enclose the origin. Grid32 versus grid128 compares attained objective values, not a general proof of global minimization. The scalar comparison permits representation and angular search error and uses only the first two timesteps.\n\n'+
      'The planar raw EL implementation refreshes the material row after the energy displacement. Its rod control therefore compiles the original traditional_cases/relinearized.hpp implementation, not the frozen traditional_cases/solver.hpp schedule. These differ for a changing rod normal. For the first two positive-branch scalar spring steps, the material row is affine and the frozen and refreshed forms agree.\n')
    raise SystemExit(bool(failures))
