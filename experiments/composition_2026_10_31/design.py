#!/usr/bin/env python3
"""Explicit compositional design, with aliases and applicability recorded separately."""
import argparse, csv, itertools, json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
WORK = ['start','end','force_average','separate_average','dg_start','dg_end','dg_average','frozen_power','kinetic_start','kinetic_end']
KIN = ['hard_endpoint','soft_endpoint','soft_midpoint']
SCHEDULE = ['joint','constraint_block','separate','symmetric']
MAIN = ['scalar_release','radial_orbit','rod_release','elastic_mixed','chain_mixed']
ALL = ['scalar_release','scalar_mixed','radial_orbit','radial_release','rod_release','rod_kick','elastic_mixed','elastic_release','chain_mixed','chain_large']
FIELDS = ['id','tier','case','family','work','treatment','kinematics','schedule','curvature','beta','gamma','iterations','tangency','keep_rx','relinearize','reverse','hz','duration','stiffness_factor','mass_ratio','perturbation','trace']
DEFAULT = dict(family='independent',work='dg_average',treatment='penalized',kinematics='soft_endpoint',schedule='joint',curvature='drop',beta=1,gamma=1,iterations=4,tangency=0,keep_rx=1,relinearize=1,reverse=0,hz=120,duration=10,stiffness_factor=1,mass_ratio=1,perturbation=1,trace=0)

def independent_core():
    for law,treatment,kin,schedule in itertools.product(WORK,['prescribed','penalized'],KIN,SCHEDULE):
        yield dict(work=law,treatment=treatment,kinematics=kin,schedule=schedule)
    for kin in KIN:
        yield dict(work='none',treatment='global',kinematics=kin)

def controls():
    for it in [1,4,16]:
        yield dict(family='xpbd',work='none',treatment='none',kinematics='hard_endpoint',schedule='constraint_block',iterations=it)
    yield dict(family='xpbd',work='none',treatment='none',kinematics='hard_endpoint',schedule='joint',iterations=16,curvature='exact')
    for law,schedule,it in itertools.product(WORK[:7],SCHEDULE,[4,16]):
        yield dict(family='force_momentum',work=law,treatment='impulse_closed',kinematics='hard_midpoint',schedule=schedule,iterations=it)
    for kin in KIN:
        yield dict(family='independent',work='none',treatment='global',kinematics=kin,iterations=32,curvature='exact')
    for beta in [.1,1,10]:
        yield dict(family='fepr_projection',work='none',treatment='global',kinematics='fixed_target',beta=beta,iterations=16)
    for family in ['shake','rattle']:
        yield dict(family=family,work='none',treatment='none',kinematics='verlet',iterations=16,tangency=int(family=='rattle'))
    for order,cap,rx,fresh in itertools.product(['EL','LE'],['raw','bounded'],[0,1],[0,1]):
        yield dict(family='EL',work='none',treatment=cap,kinematics='hard_endpoint',schedule=order,iterations=1,keep_rx=rx,relinearize=fresh)

def selected():
    """Every modeling family represented in long and stress probes, not only winners."""
    yield dict(family='xpbd',work='none',treatment='none',kinematics='hard_endpoint',schedule='constraint_block',iterations=4)
    for law in ['start','end','force_average','separate_average','dg_start','dg_end','dg_average']:
        yield dict(family='force_momentum',work=law,treatment='impulse_closed',kinematics='hard_midpoint',iterations=16)
    for kin in KIN:
        yield dict(work='none',treatment='global',kinematics=kin,iterations=16)
    for law,treat,kin in itertools.product(['start','dg_average','kinetic_end'],['prescribed','penalized'],['soft_endpoint','soft_midpoint']):
        yield dict(work=law,treatment=treat,kinematics=kin,iterations=16)
    for schedule in SCHEDULE[1:]:
        yield dict(work='dg_average',treatment='penalized',kinematics='soft_midpoint',schedule=schedule,iterations=16)
    yield dict(family='fepr_projection',work='none',treatment='global',kinematics='fixed_target',iterations=16)
    for family in ['shake','rattle']:
        yield dict(family=family,work='none',treatment='none',kinematics='verlet',iterations=16,tangency=int(family=='rattle'))
    for order in ['EL','LE']:
        yield dict(family='EL',work='none',treatment='bounded',kinematics='hard_endpoint',schedule=order,iterations=1,keep_rx=1)

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--tier',choices=['smoke','core','parameters','solver','long','stress','traces','convergence','stress_reference','serial_cost'],required=True);ap.add_argument('--out',type=Path);args=ap.parse_args();tier=args.tier
    rows=[];seen=set()
    def add(case,cfg,**extra):
        row={**DEFAULT,**cfg,**extra,'tier':tier,'case':case}
        if row['family'] in ['shake','rattle'] and not case.startswith('rod'): return
        if not case.startswith('rod'):row['tangency']=0
        # Rigid reaction has its own strength unknown. Discrete geometry strength variants coincide.
        if row['family']=='force_momentum' and case.startswith('rod') and row['work'] in ['dg_start','dg_end']:return
        if row['kinematics']=='hard_endpoint':row['beta']=1
        if row['treatment']!='penalized':row['gamma']=1
        key=tuple((k,str(row[k])) for k in FIELDS if k not in ['id','trace'])
        if key in seen and tier!='serial_cost':return
        seen.add(key);row['id']=f'{tier}_{len(rows):06d}';rows.append(row)
    if tier=='smoke':
        for case in MAIN:
            for cfg in itertools.chain(selected(),controls()):add(case,cfg,duration=.1,trace=1)
    elif tier=='core':
        for case,hz,cfg in itertools.product(MAIN,[60,120],list(independent_core())+list(controls())):
            add(case,cfg,hz=hz)
            if case.startswith('rod') and cfg.get('family','independent') in ['independent','fepr_projection']:add(case,{**cfg,'tangency':1},hz=hz)
    elif tier=='parameters':
        for case,law,kin,beta,gamma in itertools.product(MAIN,WORK,['soft_endpoint','soft_midpoint'],[.1,1,10],[.1,1,10]):
            add(case,dict(work=law,kinematics=kin,beta=beta,gamma=gamma),duration=3)
        for case,law,kin,beta in itertools.product(MAIN,WORK,['soft_endpoint','soft_midpoint'],[.1,10]):
            add(case,dict(work=law,treatment='prescribed',kinematics=kin,beta=beta),duration=3)
    elif tier=='solver':
        # Full algorithm factorial on the same local-work objectives, including both moving endpoints.
        for case,law,treat,kin,schedule,it,curvature,rx,fresh,reverse in itertools.product(MAIN,['start','dg_average','kinetic_end'],['prescribed','penalized'],['soft_endpoint','soft_midpoint'],SCHEDULE,[1,4,16],['drop','exact'],[0,1],[0,1],[0,1]):
            add(case,dict(work=law,treatment=treat,kinematics=kin,schedule=schedule,iterations=it,curvature=curvature,keep_rx=rx,relinearize=fresh,reverse=reverse),duration=.1,tangency=int(case.startswith('rod')))
        for case,law,schedule,it,fresh,reverse in itertools.product(MAIN,WORK[:7],SCHEDULE,[1,4,16],[0,1],[0,1]):
            add(case,dict(family='force_momentum',work=law,treatment='impulse_closed',kinematics='hard_midpoint',schedule=schedule,iterations=it,relinearize=fresh,reverse=reverse),duration=.1)
    elif tier=='long':
        for case,cfg in itertools.product(ALL,list(selected())):
            add(case,cfg,duration=60,hz=120,trace=1)
            if case.startswith('rod') and cfg.get('family','independent') in ['independent','fepr_projection']:add(case,{**cfg,'tangency':1},duration=60,hz=120,trace=1)
    elif tier=='stress':
        for case,cfg in itertools.product(MAIN,list(selected())):
            for hz,kfactor,mratio,perturb in [(60,1,1,1),(240,1,1,1),(60,10,1,1),(60,100,1,1),(60,1,10,1),(60,1,100,1),(60,10,10,1),(60,100,100,1),(60,1,1,.001),(60,100,10,.001)]:
                if case.startswith('rod') and kfactor!=1:continue
                add(case,cfg,duration=3,hz=hz,stiffness_factor=kfactor,mass_ratio=mratio,perturbation=perturb)
                if case.startswith('rod') and cfg.get('family','independent') in ['independent','fepr_projection']:add(case,{**cfg,'tangency':1},duration=3,hz=hz,stiffness_factor=kfactor,mass_ratio=mratio,perturbation=perturb)
    elif tier=='traces':
        for case,cfg in itertools.product(ALL,list(selected())):
            add(case,cfg,duration=10,hz=120,trace=1)
            if case.startswith('rod') and cfg.get('family','independent') in ['independent','fepr_projection']:add(case,{**cfg,'tangency':1},duration=10,hz=120,trace=1)
    elif tier=='convergence':
        for case,law,treat,kin,curv in itertools.product(MAIN,WORK,['prescribed','penalized'],KIN,['drop','exact']):
            for tangent in ([0,1] if case.startswith('rod') else [0]):
                add(case,dict(work=law,treatment=treat,kinematics=kin,curvature=curv,tangency=tangent,iterations=64),duration=1)
    elif tier=='stress_reference':
        for case in MAIN:
            for kf,mr in [(1,1),(100,1),(1,100),(100,100)]:
                if case.startswith('rod') and kf!=1:continue
                for family,it in [('xpbd',4),('xpbd',64),('force_momentum',4),('force_momentum',16)]:
                    add(case,dict(family=family,work='dg_average' if family=='force_momentum' else 'none',treatment='impulse_closed' if family=='force_momentum' else 'none',kinematics='hard_midpoint' if family=='force_momentum' else 'hard_endpoint',curvature='exact',iterations=it),stiffness_factor=kf,mass_ratio=mr,duration=3,hz=60)
    elif tier=='serial_cost':
        for case in MAIN:
            for family in ['xpbd','force_momentum','independent','fepr_projection','rattle']:
                if family=='rattle' and not case.startswith('rod'):continue
                cfg=next(c for c in selected() if c.get('family','independent')==family and (family!='force_momentum' or c['work']=='dg_average') and (family!='independent' or c['treatment']=='global' and c['kinematics']=='soft_endpoint'))
                for rep in range(5):add(case,cfg,duration=1,trace=0)
    target=args.out or ROOT/'manifests'/f'{tier}.csv';target.parent.mkdir(parents=True,exist_ok=True)
    with target.open('w') as f:w=csv.DictWriter(f,fieldnames=FIELDS);w.writeheader();w.writerows(rows)
    print(json.dumps({'tier':tier,'runs':len(rows),'manifest':str(target),'simulated_steps_requested':sum(int(r['duration']*r['hz']) for r in rows)}))
if __name__=='__main__':main()
