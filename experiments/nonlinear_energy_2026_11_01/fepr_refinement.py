"""Reproduce the 30 scalar FEPR diagnostic runs; see FEPR_REFINEMENT.md."""
import sys,csv,math,json
import numpy as np
from pathlib import Path
root=Path(__file__).resolve().parent
sys.path.insert(0,str(root))
from fepr import fepr_step,energy
out=root/'results'/'fepr_refinement.csv';out.parent.mkdir(exist_ok=True)
fields=['family','target_mode','hz','h','tolerance','epsilon','completed_time','failed','reason','rms_position','rms_velocity','final_phase_error','max_phase_error','final_energy_error','max_energy_error','iterations','regularizations','max_residual','max_abs_s','max_position_correction','max_velocity_correction','first_position_error','first_velocity_error','first_q','first_v','mean_phase_rate','crossings','small_extension_visits']
records=[]
def run(hz,tol,target_mode,eps=None,family='paper_default'):
 h=1/hz;x=np.array([1.25,0.]);v=np.zeros(2);E=3.125;qs=[.25];vs=[0.];ts=[0.];px=pv=maxe=0.;it=reg=0;maxres=maxs=maxdx=maxdv=0;first=None;fail=False;reason='';maxstep=int(round(5*hz));small=0
 for step in range(maxstep):
  r=fepr_step(x,v,h,1.,100.,1.,[0.,0.],energy_target=E if target_mode=='fixed' else None,tolerance=tol,max_iterations=128,epsilon=eps)
  if not r['converged']:fail=True;reason=r['reason'];break
  x,v=r['x'],r['v'];t=(step+1)*h;q=x[0]-1
  refq=.25*math.cos(10*t);refv=-2.5*math.sin(10*t)
  px+=(q-refq)**2;pv+=(v[0]-refv)**2
  maxe=max(maxe,abs(energy(x,v,1.,100.,1.,[0.,0.])-E));it+=r['counts']['iterations'];reg+=r['counts']['regularizations'];maxres=max(maxres,r['residual']);maxs=max(maxs,abs(r['s']));maxdx=max(maxdx,float(np.linalg.norm(x-r['star_x'])));maxdv=max(maxdv,float(np.linalg.norm(v-r['star_v'])))
  if first is None:first=(abs(q-refq),abs(v[0]-refv),q,v[0])
  qs.append(q);vs.append(v[0]);ts.append(t);small+=abs(q)<.01
 phase=np.unwrap(np.arctan2(-np.array(vs)/10,np.array(qs)));perr=phase-10*np.array(ts);n=len(ts)-1
 row=dict(family=family,target_mode=target_mode,hz=hz,h=h,tolerance=tol,epsilon=.001 if eps is None else eps,completed_time=ts[-1],failed=fail,reason=reason,rms_position=math.sqrt(px/n) if n else math.nan,rms_velocity=math.sqrt(pv/n) if n else math.nan,final_phase_error=float(perr[-1]),max_phase_error=float(np.max(abs(perr))),final_energy_error=energy(x,v,1.,100.,1.,[0.,0.])-E,max_energy_error=maxe,iterations=it,regularizations=reg,max_residual=maxres,max_abs_s=maxs,max_position_correction=maxdx,max_velocity_correction=maxdv,mean_phase_rate=float(phase[-1]/ts[-1]) if n else math.nan,crossings=int(np.sum(np.array(qs[1:])*np.array(qs[:-1])<0)),small_extension_visits=small)
 row.update(zip(['first_position_error','first_velocity_error','first_q','first_v'],first or [math.nan]*4));records.append(row)
for target in ['fixed','per_step']:
 for tol in [1e-7,1e-10]:
  for hz in [30,60,120,240,480]:run(hz,tol,target)
# Existing-parameter ablations diagnose the auxiliary penalty; they are not silently adopted fixes.
for hz in [30,60,120,240,480]:
 for eps in [1e-9,.001*(1/hz)**2]:run(hz,1e-10,'fixed',eps,'epsilon_control')
with out.open('w') as f:
 w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(records)
for r in records:print(r['family'],r['target_mode'],r['hz'],r['tolerance'],r['epsilon'],'RMS %.6g phase %.6g E %.3g firstx %.3g firstv %.3g maxs %.3g'%(r['rms_position'],r['final_phase_error'],r['max_energy_error'],r['first_position_error'],r['first_velocity_error'],r['max_abs_s']),r['reason'])
