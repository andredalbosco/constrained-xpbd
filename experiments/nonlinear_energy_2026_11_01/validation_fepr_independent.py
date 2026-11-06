#!/usr/bin/env python3
"""Independent source/equation checks of FEPR fast projection and rod extension."""
import csv
import math
from pathlib import Path
import numpy as np
import fepr as f

OUT=Path(__file__).parent/'fepr_independent_validation'
OUT.mkdir(exist_ok=True)
rng=np.random.default_rng(801314)
checks=[]


def check(name,error,tolerance,detail=''):
    checks.append(dict(check=name,error=float(error),tolerance=float(tolerance),
                       passed=bool(error<=tolerance),details=detail))


def constraints_independent(q,old,ov,star,sv,m,k,l,g,target,rigid):
    x,v,s,t=q[:2],q[2:4],q[4],q[5]
    cross=lambda a,b:a[0]*b[1]-a[1]*b[0]
    u=-m*float(g@x)+(0 if rigid else .5*k*(np.linalg.norm(x)-l)**2)
    c=[.5*m*(v@v)+u-target,
       m*(v[0]-sv[0])-s*m*(ov[0]-sv[0]),
       m*(v[1]-sv[1])-s*m*(ov[1]-sv[1]),
       m*(cross(x,v)-cross(star,sv))-t*m*(cross(old,ov)-cross(star,sv))]
    if rigid:c.extend([.5*(x@x-l*l),x@v])
    return np.asarray(c)


def derivative_and_rows():
    for rigid in [False,True]:
        for case in range(24):
            theta=rng.uniform(-2,2)
            l,m,k=1.1,1.3,80.
            old=l*np.array([math.sin(theta),-math.cos(theta)])
            ov=rng.uniform(-2,2)*np.array([math.cos(theta),math.sin(theta)])
            star,sv=old+.1*rng.normal(size=2),ov+.2*rng.normal(size=2)
            q=np.r_[star+.02*rng.normal(size=2),sv+.1*rng.normal(size=2),rng.uniform(-2,2),rng.uniform(-2,2)]
            g=np.array([0.,-9.81])
            target=f.energy(old,ov,m,k,l,g,rigid)
            args=(old,ov,star,sv,.04,m,k,l,g,rigid,target)
            c,j=f._constraints(q,*args)
            direct=constraints_independent(q,old,ov,star,sv,m,k,l,g,target,rigid)
            check('independent_constraint_rows',np.max(np.abs(c-direct)),2e-13)
            fd=np.empty_like(j)
            for coordinate in range(6):
                dq=np.zeros(6);dq[coordinate]=2e-6*max(1.,abs(q[coordinate]))
                fd[:,coordinate]=(constraints_independent(q+dq,old,ov,star,sv,m,k,l,g,target,rigid)-constraints_independent(q-dq,old,ov,star,sv,m,k,l,g,target,rigid))/(2*dq[coordinate])
            check('finite_difference_jacobian',np.linalg.norm(j-fd)/(1+np.linalg.norm(j)),2e-8,str(rigid))
            old_feasible=np.r_[old,ov,1.,1.]
            check('old_state_per_step_feasibility',np.linalg.norm(f._constraints(old_feasible,*args)[0]),3e-14,str(rigid))
            changed=list(args);changed[-1]=target+.3
            fixed=f._constraints(old_feasible,*changed)[0]
            check('fixed_target_old_state_not_automatically_feasible',abs(fixed[0]+.3)+np.linalg.norm(fixed[1:]),3e-14,str(rigid))


def fast_update_vs_kkt():
    for rigid in [False,True]:
        for case in range(12):
            theta=rng.uniform(-1,1)
            old=np.array([math.sin(theta),-math.cos(theta)])
            ov=.7*np.array([math.cos(theta),math.sin(theta)])
            star,sv=old+.04*rng.normal(size=2),ov+.15*rng.normal(size=2)
            h,m,k,l,eps=.07,1.3,80.,1.,.003
            g=np.array([0.,-9.81])
            target=f.energy(old,ov,m,k,l,g,rigid)
            metric=np.diag([m,m,m*h*h,m*h*h,eps,eps])
            args=(old,ov,star,sv,h,m,k,l,g,rigid,target)
            initial=np.r_[star,sv,0.,0.]
            q=initial.copy()
            for iteration in range(2):
                c,j=f._constraints(q,*args)
                system=np.block([[metric,j.T],[j,np.zeros((len(c),len(c)))]])
                delta=np.linalg.solve(system,np.r_[np.zeros(6),-c])[:6]
                q+=delta
                actual=f.project_fepr(old,ov,star,sv,h,m,k,l,g,rigid,energy_target=target,epsilon=eps,tolerance=1e-14,max_iterations=iteration+1,determinant_threshold=1e-250,schur_regularization=0,record_trace=True)
                got=np.r_[actual['x'],actual['v'],actual['s'],actual['t']]
                check('fast_update_full_increment_KKT',np.linalg.norm(got-q)/(1+np.linalg.norm(q)),2e-8,f'rigid={rigid},iteration={iteration+1}')
                check('iteration_count_no_hidden_extra_solve',abs(actual['counts']['iterations']-(iteration+1)),0)
                check('default_metric_override_epsilon',abs(actual['epsilon']-eps),0)
            check('finite_iteration_limit_explicit',0 if actual['converged'] or actual['reason']=='iteration_limit' else 1,0)
    # Demonstrate second-iteration FEPR is not fixed-center SQP.
    old,ov=np.array([1.2,.2]),np.array([-.3,.8])
    base=f.baseline_step(old,ov,.1,1.,100.,1.,np.zeros(2))
    sx,sv=base['x'],base['v']
    target=f.energy(old,ov,1.,100.,1.,np.zeros(2))
    d=np.diag([1.,1.,.01,.01,.001,.001])
    start=np.r_[sx,sv,0.,0.]
    first=f.project_fepr(old,ov,sx,sv,.1,1.,100.,1.,np.zeros(2),energy_target=target,max_iterations=1,determinant_threshold=1e-250,schur_regularization=0)
    q1=np.r_[first['x'],first['v'],first['s'],first['t']]
    c,j=f._constraints(q1,old,ov,sx,sv,.1,1.,100.,1.,np.zeros(2),False,target)
    system=np.block([[d,j.T],[j,np.zeros((4,4))]])
    fast_delta=np.linalg.solve(system,np.r_[np.zeros(6),-c])[:6]
    fixed_delta=np.linalg.solve(system,np.r_[-d@(q1-start),-c])[:6]
    difference=np.linalg.norm(fast_delta-fixed_delta)
    check('fast_distinct_from_fixed_center_SQP',0 if difference>1e-6 else 1,0,f'actual update difference={difference:.9g}')


def status_and_no_clipping():
    old,ov=np.array([1.25,0.]),np.zeros(2)
    result=f.fepr_step(old,ov,.04,1.,100.,1.,np.zeros(2),tolerance=1e-11,max_iterations=128)
    check('spring_release_converges',0 if result['converged'] else 1,0,result['reason'])
    check('auxiliary_s_not_clipped',0 if result['s']<0 else 1,0,f"s={result['s']:.9g}")
    check('anchored_linear_momentum_not_forced_old',0 if np.linalg.norm(result['v']-ov)>1e-3 else 1,0)
    check('velocity_independent_after_FEPR',0 if np.linalg.norm(result['v']-(result['x']-old)/.04)>1e-4 else 1,0)
    check('paper_raw_L1_stopping',max(0,result['residual']-1e-11),0)
    # Impossible target and zero gradient: regularize but do not invent motion.
    singular=f.project_fepr(np.array([1.,0.]),np.zeros(2),np.array([1.,0.]),np.zeros(2),.01,1.,100.,1.,np.zeros(2),energy_target=1.,max_iterations=3,record_trace=True)
    check('singular_explicit_failure',0 if not singular['converged'] and singular['reason']=='iteration_limit' else 1,0)
    check('singular_regularization_count',abs(singular['counts']['regularizations']-3),0)
    check('singular_return_not_hidden_fallback',np.linalg.norm(singular['x']-np.array([1.,0.]))+np.linalg.norm(singular['v']),0)
    # Per-step target is feasible, but a zero iteration cap must not claim success.
    base=f.baseline_step(old,ov,.04,1.,100.,1.,np.zeros(2))
    capped=f.project_fepr(old,ov,base['x'],base['v'],.04,1.,100.,1.,np.zeros(2),max_iterations=0)
    check('zero_iteration_cap_explicit',0 if not capped['converged'] and capped['reason']=='iteration_limit' and capped['counts']['iterations']==0 else 1,0)
    # A free drift already satisfies every row and must be left untouched.
    ov=np.array([.3,.4])
    drift=f.fepr_step(old,ov,.04,1.,0.,1.,np.zeros(2))
    check('feasible_free_drift_unchanged',np.linalg.norm(drift['x']-(old+.04*ov))+np.linalg.norm(drift['v']-ov),2e-13)
    check('feasible_free_drift_zero_projection',drift['counts']['iterations'],0)


def rattle_checks():
    for case in range(24):
        theta=rng.uniform(-1.5,1.5)
        h=1 / float(rng.choice([60,120,480]))
        m,l=1.3,1.2
        x=l*np.array([math.sin(theta),-math.cos(theta)])
        v=l*rng.uniform(-2,2)*np.array([math.cos(theta),math.sin(theta)])
        g=np.array([0.,-9.81])
        base=f.baseline_step(x,v,h,m,0,l,g,True)
        check('rattle_base_length',abs(np.linalg.norm(base['x'])-l),2e-13)
        check('rattle_base_tangency',abs(base['x']@base['v']),2e-13)
        # Reverse h by reversing velocity: time-reversibility of the base step.
        reverse=f.baseline_step(base['x'],-base['v'],h,m,0,l,g,True)
        check('rattle_base_reversibility',np.linalg.norm(reverse['x']-x)+np.linalg.norm(reverse['v']+v),2e-10)
        out=f.fepr_step(x,v,h,m,0,l,g,True,tolerance=1e-10,max_iterations=128)
        check('hard_rod_extension_reports_identity',0 if 'extension' in out['method'] else 1,0)
        if out['converged']:
            check('hard_rod_extension_length',abs(np.linalg.norm(out['x'])-l),2e-10)
            check('hard_rod_extension_tangency',abs(out['x']@out['v']),2e-10)
            check('hard_rod_extension_energy',abs(f.energy(out['x'],out['v'],m,0,l,g,True)-f.energy(x,v,m,0,l,g,True)),2e-10)
        else:
            check('hard_rod_failure_explicit',0 if out['reason']!='converged' else 1,0,out['reason'])


if __name__=='__main__':
    derivative_and_rows();fast_update_vs_kkt();status_and_no_clipping();rattle_checks()
    with (OUT/'checks.csv').open('w',newline='') as out:
        writer=csv.DictWriter(out,['check','error','tolerance','passed','details'])
        writer.writeheader();writer.writerows(checks)
    failed=[x for x in checks if not x['passed']]
    worst={}
    for row in checks:worst[row['check']]=max(worst.get(row['check'],0),row['error'])
    print('checks',len(checks),'failures',len(failed))
    for name,error in worst.items():print(name,error)
    for row in failed:print('FAIL',row)
    (OUT/'SUMMARY.md').write_text('# Independent FEPR audit\n\n'+f'{len(checks)} checks; {len(failed)} failures.\n\n'+
      'Read the primary-paper text through Equations 1, 8, 12, 13, and 14. The implementation matches the paper fast projection in planar form: four energy/momentum rows, independent x and v, unrestricted auxiliary scalars, D=diag(mI,h²mI,epsilon,epsilon), and a current-iterate distance center. It does not implement fixed-center full SQP.\n\n'+
      '| Check | Maximum measured error |\n|---|---:|\n'+''.join(f'| {name} | {error:.9g} |\n' for name,error in worst.items())+
      '\nThe additional rod length and tangency rows, together with the RATTLE base step, are an explicitly named extension. Passing this audit validates the row formulas and numerical update, not a claim that the paper evaluated hard rods. Finite-difference tests include nonzero s/t, and the KKT comparison checks two sequential updates so fixed-center SQP cannot pass accidentally. Raw dimensional residual and determinant thresholds match the declared policy; they are not unit-invariant.\n')
    raise SystemExit(bool(failed))
