"""One-step diagnostic: full fixed-center FEPR objective versus fast projection.
This is an oracle experiment, not a change to the published fast implementation.
"""
import csv,json
from pathlib import Path
import numpy as np
from scipy.optimize import minimize
import fepr
root=Path(__file__).resolve().parent
rows=[r for r in csv.DictReader(open(root/'stress_results/fepr_states.csv')) if r['case']=='spring_stiff' and r['method']=='fepr_step_target' and r['hz']=='30']
oldrow=rows[1]
x=np.array([float(oldrow['x']),float(oldrow['y'])]);v=np.array([float(oldrow['vx']),float(oldrow['vy'])]);h=1/30;m=1.;k=1000.;l=1.;g=np.zeros(2)
base=fepr.baseline_step(x,v,h,m,k,l,g)
fast=fepr.fepr_step(x,v,h,m,k,l,g)
D=np.array([m,m,h*h*m,h*h*m,.001,.001]);sqrtD=np.sqrt(D)
star=np.r_[base['x'],base['v'],0.,0.];target=fepr.energy(x,v,m,k,l,g)
args=(x,v,base['x'],base['v'],h,m,k,l,g,False,target)
# Row scaling is a nonsingular change of the constraints, not their objective.
scales=np.array([max(1,abs(target)),1.,1.,1.])
def q_from(z):return star+z/sqrtD

def con(z):return fepr._constraints(q_from(z),*args)[0]/scales

def jac(z):return fepr._constraints(q_from(z),*args)[1]/scales[:,None]/sqrtD

seeds={'tentative':star,'old_feasible':np.r_[x,v,1.,1.],'fast':np.r_[fast['x'],fast['v'],fast['s'],fast['t']]}
results=[]
for name,q0 in seeds.items():
 z0=(q0-star)*sqrtD
 opt=minimize(lambda z:.5*np.dot(z,z),z0,jac=lambda z:z,method='SLSQP',constraints=[{'type':'eq','fun':con,'jac':jac}],options={'ftol':1e-13,'maxiter':1000})
 q=q_from(opt.x);c,J=fepr._constraints(q,*args);Jw=J/sqrtD
 mult=np.linalg.lstsq(Jw.T,-opt.x,rcond=1e-12)[0]
 stationarity=np.linalg.norm(opt.x+Jw.T@mult)
 results.append(dict(seed=name,success=bool(opt.success),message=opt.message,iterations=int(opt.nit),objective=float(opt.fun),constraint_l1=float(np.abs(c).sum()),whitened_stationarity=float(stationarity),jump=float(np.linalg.norm(q[:2]-x)),q=q.tolist()))
qfast=seeds['fast'];zfast=(qfast-star)*sqrtD
fast_report=dict(objective=float(.5*np.dot(zfast,zfast)),constraint_l1=fast['residual'],jump=float(np.linalg.norm(fast['x']-x)),q=qfast.tolist())
report=dict(case='spring_stiff second step at30Hz',old_x=x.tolist(),old_v=v.tolist(),fast=fast_report,full_objective=results)
(root/'fepr_full_oracle.json').write_text(json.dumps(report,indent=2))
lines=['# One-step full FEPR objective diagnostic','','This solves the fixed-center Eq1 objective on the stiff spring second step at30Hz. It does not replace the source-faithful fast FEPR algorithm. The state is read from the saved per-step-target trajectory. Variables are whitened by the exact D metric, equality rows are nonsingularly scaled, and SLSQP uses analytic Jacobians. Three starts are tested; this is local minimum evidence, not a proof of the global minimum.','', '| Solver/seed | Success | Objective | Raw constraint L1 | Whitened stationarity | Position jump (m) |','|---|---|---:|---:|---:|---:|',f"| Fast FEPR | yes | {fast_report['objective']:.9g} | {fast_report['constraint_l1']:.4g} | not enforced | {fast_report['jump']:.9g} |"]
for r in results:lines.append(f"| Full / {r['seed']} | {r['success']} | {r['objective']:.9g} | {r['constraint_l1']:.4g} | {r['whitened_stationarity']:.4g} | {r['jump']:.9g} |")
lines += ['','This diagnostic distinguishes reaching the energy/momentum manifold from minimizing distance to the tentative state. No full-objective long trajectory or scalable local implementation is claimed. The observed local solution can motivate a safeguarded/fixed-center fallback, but this single-step result is not a validated general repair.']
(root/'FEPR_FULL_ORACLE.md').write_text('\n'.join(lines))
print(json.dumps(report,indent=2))
