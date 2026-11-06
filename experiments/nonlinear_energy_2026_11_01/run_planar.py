"""Reproducible spring/pendulum comparison; failed prefixes are not successes."""
import argparse,csv,json,time
from pathlib import Path
import numpy as np
from planar import energy,step
from references import reference_trajectory,spring_dg_step,pendulum_dg_step

ROOT=Path(__file__).resolve().parent

def cases():
    def rod(name,theta,omega):
        return dict(name=name,rigid=True,m=1.,k=0.,l=1.,g=np.array([0.,-9.81]),x=np.array([np.sin(theta),-np.cos(theta)]),v=omega*np.array([np.cos(theta),np.sin(theta)]))
    base=dict(rigid=False,m=1.,k=100.,l=1.,g=np.zeros(2))
    return [dict(base,name='spring_release',x=np.array([1.25,0.]),v=np.zeros(2)),
            dict(base,name='spring_orbit',x=np.array([1.1,0.]),v=np.array([0.,.8])),
            dict(base,name='spring_gravity',g=np.array([0.,-9.81]),x=np.array([.55,-1.1]),v=np.array([.4,.2])),
            rod('pendulum_release',.65,0),rod('pendulum_kick',.4,1.2),rod('pendulum_large',2.4,.5),
            dict(base,name='spring_stiff',k=1000.,x=np.array([1.1,0.]),v=np.array([0.,.8]))]


def serialize_state(case,method,sweeps,hz,j,x,v,E,E0,scale,refx,refv,post_E=0):
    return dict(case=case,method=method,sweeps=sweeps,hz=hz,step=j,time=j/hz,x=x[0],y=x[1],vx=v[0],vy=v[1],energy=E,rel_energy=(E-E0)/scale,
                ref_x=refx[0],ref_y=refx[1],ref_vx=refv[0],ref_vy=refv[1],post_E=post_E)


def run(c,method,sweeps,hz,duration,ref):
    x,v=c['x'].copy(),c['v'].copy();m,k,l,g=c['m'],c['k'],c['l'],c['g'];h=1/hz
    E0=energy(x,v,m,k,l,g);emin=-m*np.linalg.norm(g)*l-(.5*(m*np.linalg.norm(g))**2/k if k else 0)
    scale=max(E0-emin,1e-8);N=round(duration*hz)
    rows=[serialize_state(c['name'],method,sweeps,hz,0,x,v,E0,E0,scale,ref.x[0],ref.v[0])]
    reason='';maxe=maxrod=maxtan=maxang=maxpostE=maxveldef=0.;sumsq=path=0.;counts=ties=0
    start=time.perf_counter();j=0;lastres=0
    for it in range(1,N+1):
        previous=x.copy()
        if method=='dg':
            if c['rigid']:out=pendulum_dg_step(x,v,h,m,l,g)
            else:out=spring_dg_step(x,v,h,m,k,l,g)
            ok=out.converged;reason=out.reason if not ok else '';xx,vv=out.x,out.v
            nfev=out.nfev;pe=0.;nt=0;lastres=out.residual
        elif method.startswith('fepr_sub'):
            from fepr import fepr_step
            subs=int(method.removeprefix('fepr_sub'));xx,vv=x.copy(),v.copy();ok=True;nfev=0;reason=''
            for sub in range(subs):
                out=fepr_step(xx,vv,h/subs,m,k,l,g,rigid=c['rigid'])
                nfev+=out.get('counts',{}).get('constraint_evaluations',0)+out.get('counts',{}).get('baseline_evaluations',0)
                if not out['converged']:ok=False;reason=out['reason'];break
                xx,vv=out['x'],out['v']
            pe=0.;nt=0;lastres=out['residual']
        elif method=='cpd':
            from cpd_energy import cpd_spring_step,cpd_pendulum_step
            out=cpd_pendulum_step(x,v,h,m,l,g) if c['rigid'] else cpd_spring_step(x,v,h,m,k,l,g)
            ok=out.converged;reason=out.reason if not ok else '';xx,vv=out.x,out.v
            nfev=out.nfev;pe=0.;nt=0;lastres=out.residual
        elif method in ('fepr','fepr_step_target','fepr_base'):
            from fepr import fepr_step,baseline_step
            kwargs=dict(energy_target=E0 if method=='fepr' else None)
            out=baseline_step(x,v,h,m,k,l,g,rigid=c['rigid']) if method=='fepr_base' else fepr_step(x,v,h,m,k,l,g,rigid=c['rigid'],**kwargs)
            ok=out['converged'];reason=out.get('reason','') if not ok else '';xx,vv=np.array(out['x']),np.array(out['v'])
            nfev=out.get('counts',{}).get('constraint_evaluations',0)+out.get('counts',{}).get('baseline_evaluations',0);pe=0.;nt=0;lastres=out['residual']
        else:
            meth='pd_joint' if c['rigid'] and method=='pd_true' else method
            out=step(x,v,h,m,k,l,g,E0,meth,sweeps,c['rigid'])
            ok=out.converged;reason=out.reason if not ok else '';xx,vv=out.x,out.v
            nfev=out.nfev;pe=out.post_E;nt=out.ties;lastres=out.residual
        counts+=nfev;ties+=nt
        if not ok:break
        e=energy(xx,vv,m,k,l,g);rel=abs(e-E0)/scale
        if not np.all(np.isfinite(xx)) or not np.all(np.isfinite(vv)) or not np.isfinite(e):reason='nonfinite';break
        # Include the offending candidate in diagnostics/CSV before stopping.
        x,v=xx,vv;j=it
        rows.append(serialize_state(c['name'],method,sweeps,hz,j,x,v,e,E0,scale,ref.x[j],ref.v[j],pe))
        maxe=max(maxe,rel);sumsq+=float(np.sum((x-ref.x[j])**2));path+=float(np.linalg.norm(x-previous))
        if c['rigid']:
            maxrod=max(maxrod,abs(np.linalg.norm(x)-l));maxtan=max(maxtan,abs(np.dot(x,v))/l)
        if not c['rigid'] and np.linalg.norm(g)==0:
            a=m*(x[0]*v[1]-x[1]*v[0]);a0=m*(c['x'][0]*c['v'][1]-c['x'][1]*c['v'][0]);maxang=max(maxang,abs(a-a0))
        if np.isfinite(pe):maxpostE=max(maxpostE,abs(pe)/scale)
        maxveldef=max(maxveldef,float(np.linalg.norm(v-(x-previous)/h)))
        if rel>100 or np.linalg.norm(x)>1e4 or np.linalg.norm(v)>1e6:reason='runaway_guard';break
    summary=dict(case=c['name'],method=method,sweeps=sweeps,hz=hz,duration=duration,completed=int(j==N and not reason),accepted_time=j*h,reason=reason,
       max_rel_energy=maxe,final_rel_energy=rows[-1]['rel_energy'],rms_position=(sumsq/j)**.5 if j else float('nan'),path_length=path,max_rod_error=maxrod,max_tangency_error=maxtan,
       max_angular_error=maxang,max_post_E_rel=maxpostE,max_endpoint_velocity_defect=maxveldef,nfev=counts,ties=ties,last_residual=lastres,seconds=time.perf_counter()-start)
    return summary,rows


def write_csv(path,rows):
    if not rows:return
    with path.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--group',default='endpoint',choices=['endpoint','fepr','dg','cpd','all']);ap.add_argument('--hz',default='60,240');ap.add_argument('--duration',type=float,default=5);ap.add_argument('--cases',default='spring_release,spring_orbit,spring_gravity,pendulum_release,pendulum_kick');ap.add_argument('--output',default='planar_results');ap.add_argument('--methods',default='');args=ap.parse_args()
    out=ROOT/args.output;out.mkdir(exist_ok=True);summaries=[];states=[]
    methods=[]
    if args.group in ('endpoint','all'):methods += [('xpbd',1),('raw_EL',1),('nonlinear_EL',1),('nonlinear_EL',8),('pd_true',1),('pd_true',8),('pd_surrogate',1),('pd_surrogate',8),('pd_inertia',1)]
    if args.group in ('cpd','all'):methods += [('cpd',1)]
    if args.group in ('dg','all'):methods += [('dg',1)]
    if args.group in ('fepr','all'):methods += [('fepr',1),('fepr_step_target',1),('fepr_base',1)]
    if args.methods: methods=[(name,int(sweeps)) for name,sweeps in (item.split(':') for item in args.methods.split(','))]
    for c in cases():
        if c['name'] not in args.cases.split(','):continue
        for hz in map(int,args.hz.split(',')):
            times=np.arange(round(args.duration*hz)+1)/hz
            ref=reference_trajectory('pendulum' if c['rigid'] else 'spring',c['x'],c['v'],times,c['m'],c['k'],c['l'],c['g'])
            if not ref.success:raise RuntimeError(ref.message)
            for method,sweeps in methods:
                if c['rigid'] and method in ('pd_surrogate','pd_inertia'):continue
                s,r=run(c,method,sweeps,hz,args.duration,ref);summaries.append(s);states.extend(r)
                print(c['name'],hz,method,sweeps,'done',s['completed'],'E',f"{s['max_rel_energy']:.3g}",'rms',f"{s['rms_position']:.3g}",s['reason'],flush=True)
                write_csv(out/(args.group+'_summary.csv'),summaries)
    write_csv(out/(args.group+'_states.csv'),states)

if __name__=='__main__':main()
