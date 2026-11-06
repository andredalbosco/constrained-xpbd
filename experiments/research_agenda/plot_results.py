#!/usr/bin/env python3
"""Independent physical-energy/trajectory audit and scientific comparison plots."""
import argparse, csv, json, math, os
from array import array
from collections import defaultdict
from pathlib import Path
os.environ.setdefault('MPLCONFIGDIR','/tmp/xpbd-agenda-matplotlib')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.colors import LogNorm
import numpy as np

NAMES={'spring_release':'Spring release','pendulum_rest':'Rigid pendulum','double_small':'Double spring: small swing','double_large':'Double spring: large swing','spring_orbit':'Spring with tangential velocity','pendulum_kick':'Pendulum with tangential velocity','double_mixed':'Double spring: mixed motion'}
LABEL={'xpbd':'XPBD','xpbd_ledger':'XPBD + ledger only','EL_bounded':'EL + repair/cap','LE_bounded':'LE + repair/cap','EL_rx_bounded':'EL, retain r_x','LE_rx_bounded':'LE, retain r_x','work_lagged':'Local work, lagged force','work_column':'Column-matched work','work_row':'Row-matched work','work_secant':'Potential-exact work','work_lagged_GS16':'Local work, 16 GS sweeps','reservoir_fast':'Work + repayment, τ=0.25 s','reservoir_slow':'Work + repayment, τ=1 s','midpoint_1':'Midpoint work, 1 pass','midpoint_4':'Midpoint work, 4 passes','midpoint_drop_rx':'Midpoint work, drop r_x'}
COLOR={'xpbd':'#8b969e','xpbd_ledger':'#8b969e','EL_bounded':'#d07b20','LE_bounded':'#2585af','EL_rx_bounded':'#bc6821','LE_rx_bounded':'#2585af','work_lagged':'#9566af','work_column':'#3475ba','work_row':'#20a0b1','work_secant':'#d49231','work_lagged_GS16':'#b54a59','reservoir_fast':'#cd4e66','reservoir_slow':'#aa7a55','midpoint_1':'#71ad70','midpoint_4':'#176d44','midpoint_drop_rx':'#74a32c'}
PRIMARY=['spring_release','pendulum_rest','double_large']
MAIN=['xpbd','EL_bounded','LE_rx_bounded','work_lagged','reservoir_fast','midpoint_1','midpoint_4']
def read(path):
 with path.open() as f:return list(csv.DictReader(f))
def summaries(path):
 return {(r['case'],r['method'],int(r['hz'])):{k:v if k in ('case','method','status') else float(v) for k,v in r.items()} for r in read(path/'summary.csv')}
def traces(path):
 groups=defaultdict(lambda:array('d'))
 with (path/'traces.csv').open() as f:
  r=csv.reader(f);header=next(r);cols={k:i for i,k in enumerate(header[3:])};width=len(cols)
  for row in r:groups[row[0],row[1],int(row[2])].extend(map(float,row[3:]))
 return {k:np.asarray(v).reshape(-1,width) for k,v in groups.items()},cols
def audit(groups,col,summary,models,full):
 out={'groups':len(groups),'rows':0,'max_energy_reconstruction_error':0.,'max_rms_reconstruction_error':0.,'max_cumulative_ledger_error':0.,'max_kinematic_error':0.}
 for key,a in groups.items():
  s=summary[key];m=models[key[0]];n=int(m['count']);h=1/key[2];at=lambda name:a[:,col[name]]
  assert len(a)>1 and at('step')[0]==0
  assert at('step')[-1]==s['completed'],key
  assert np.all(np.diff(at('step'))>0)
  E=np.zeros(len(a));squared=np.zeros(len(a));kin=0.
  for i in range(1,n+1):
   x,y,vx,vy=[at(k+str(i)) for k in ['x','y','vx','vy']]
   E+=m['m'+str(i)]*(.5*(vx*vx+vy*vy)-m['gx']*x-m['gy']*y)
   ex=x-(at('x'+str(i-1)) if i>1 else 0);ey=y-(at('y'+str(i-1)) if i>1 else 0)
   alpha=m['alpha'+str(i)]
   if alpha>0:E+=(np.hypot(ex,ey)-m['l'+str(i)])**2/(2*alpha)
   squared+=(x-at('ref_x'+str(i)))**2+(y-at('ref_y'+str(i)))**2
   if full:
    for position,velocity in [(x,vx),(y,vy)]:
     reconstructed=2*np.diff(position)/h-velocity[:-1] if key[1].startswith('midpoint') else np.diff(position)/h
     kin=max(kin,float(np.max(np.abs(reconstructed-velocity[1:]))))
  err=float(np.max(np.abs(E-at('energy'))))/s['energy_scale'];assert err<1e-8,(key,err)
  out['max_energy_reconstruction_error']=max(out['max_energy_reconstruction_error'],err)
  if full:
   assert np.all(np.diff(at('step'))==1)
   rms=math.sqrt(float(np.mean(squared[1:]/n)));delta=abs(rms-s['rms_position']);assert delta<1e-10,(key,delta)
   emax=float(np.max(np.abs(E[1:]-s['initial_energy'])))/s['energy_scale'];assert math.isclose(emax,s['max_energy_rel'],rel_tol=1e-9,abs_tol=1e-10),(key,emax)
   out['max_rms_reconstruction_error']=max(out['max_rms_reconstruction_error'],delta)
  if key[1].startswith(('work_','reservoir_')) or key[1]=='xpbd_ledger':
   R=sum(at(name) for name in ['RK1','RK2','RU1','RU2','RU3','RU4']);assert np.max(np.abs(R-at('reservoir')))<1e-8*s['energy_scale']
   ledger=float(np.max(np.abs(E+R-s['initial_energy'])))/s['energy_scale'];assert ledger<1e-8,(key,ledger)
   out['max_cumulative_ledger_error']=max(out['max_cumulative_ledger_error'],ledger)
  out['max_kinematic_error']=max(out['max_kinematic_error'],kin);out['rows']+=len(a)
 return out
def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('results',type=Path);p.add_argument('--long-results',type=Path);p.add_argument('--output',type=Path,required=True);args=p.parse_args();args.output.mkdir(parents=True,exist_ok=True)
 models={r['case']:{k:float(v) for k,v in r.items() if k!='case'} for r in read(args.results/'cases.csv')};summary=summaries(args.results);groups,col=traces(args.results)
 info={'main':audit(groups,col,summary,models,True),'runs':len(summary),'aborted':[{k:r[k] for k in ['case','method','hz','completed','steps','status','attempted_energy']} for r in summary.values() if r['status']!='completed']}
 print(json.dumps(info,indent=2),flush=True)
 plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.right':False,'axes.spines.top':False})
 def style(method):return {'color':COLOR[method],'ls':'--' if method in ('midpoint_1','midpoint_drop_rx','EL_rx_bounded','LE_rx_bounded') else '-','lw':1.7 if method=='midpoint_4' else 1.1}
 def data(case,method,metric,source=groups,index=col):
  a=source[case,method,480];s=summary[case,method,480];at=lambda k:a[:,index[k]];n=int(s['count'])
  if metric=='energy':v=100*np.abs(at('energy')-s['initial_energy'])/s['energy_scale']
  elif metric=='signed':v=100*(at('energy')-s['initial_energy'])/s['energy_scale']
  elif metric=='error':v=np.sqrt(sum((at('x'+str(i))-at('ref_x'+str(i)))**2+(at('y'+str(i))-at('ref_y'+str(i)))**2 for i in range(1,n+1))/n)
  elif metric=='reservoir':v=100*sum(np.abs(at(name)) for name in ['RK1','RK2','RU1','RU2','RU3','RU4'])/s['energy_scale']
  elif case.startswith('pendulum'):v=np.arctan2(at('x1'),-at('y1'))
  else:v=at('x1')-(models[case]['l1'] if case=='spring_release' else 0)
  return at('time'),v
 def graph(ax,case,methods,metric,source=groups,index=col):
  for m in methods:
   t,v=data(case,m,metric,source,index);ax.plot(t,v,**style(m))
  ax.grid(alpha=.2);ax.set_xlabel('Time (s)')
  if metric=='energy':ax.set_yscale('symlog',linthresh=1e-8);ax.set_ylabel('Physical energy error (%) · symlog');ax.set_ylim(bottom=0)
  elif metric=='signed':ax.set_ylabel('Physical energy change (%)')
  elif metric=='error':ax.set_yscale('symlog',linthresh=1e-6);ax.set_ylim(bottom=0);ax.set_ylabel('Position error vs reference (m)')
  elif metric=='reservoir':ax.set_ylabel('Sum of |local reservoirs| / E scale (%)')
  else:
   a=source[case,'xpbd',480];at=lambda k:a[:,index[k]]
   v=np.arctan2(at('ref_x1'),-at('ref_y1')) if case.startswith('pendulum') else at('ref_x1')-(models[case]['l1'] if case=='spring_release' else 0)
   ax.plot(at('time'),v,':',color='#111',lw=1.2)
   ax.set_ylabel('Angle from vertical (rad)' if case.startswith('pendulum') else 'x − rest length (m)' if case=='spring_release' else 'First mass x (m)')
 def decorate(fig,title,methods,foot,reference=False):
  fig.subplots_adjust(top=.84,bottom=.09,left=.07,right=.98,hspace=.43,wspace=.28);fig.suptitle(title,x=.07,y=.98,ha='left',fontsize=19,fontweight='bold')
  handles=[Line2D([],[],label=LABEL[m],**style(m)) for m in methods]
  if reference:handles+=[Line2D([],[],color='#111',ls=':',label='Reference motion')]
  fig.legend(handles=handles,loc='upper left',bbox_to_anchor=(.065,.945),ncol=4,frameon=False,fontsize=9)
  fig.text(.07,.015,foot,fontsize=9,color='#333')
 def panel(name,title,cases,methods,metrics,foot):
  fig,axs=plt.subplots(len(cases),len(metrics),figsize=(17,12),squeeze=False);decorate(fig,title,methods,foot,'motion' in metrics)
  for i,c in enumerate(cases):
   for j,metric in enumerate(metrics):graph(axs[i,j],c,methods,metric);axs[i,j].set_title(NAMES[c],loc='left',fontweight='bold')
  fig.savefig(args.output/name,dpi=150,facecolor='white');plt.close(fig)
 panel('overview.png','Changing the solve versus changing the integration rule',PRIMARY,MAIN,['energy','motion','error'],'480 Hz · 10 seconds · all plots use analytic physical energy; reservoirs are excluded from energy curves.\nEvery solver uses scalar Gauss–Seidel and linearized rows. Midpoint changes the kinematic relation and force discretization.')
 original_label=LABEL['midpoint_4'];LABEL['midpoint_4']='Guarded midpoint work, 4 passes'
 panel('stable_candidates.png','Candidates after the stability screen: XPBD and guarded midpoint work',PRIMARY,['xpbd','midpoint_4'],['energy','motion','error'],'480 Hz · 10 seconds · physical energy and actual positions. Black dotted curves are the refined reference.\nThe guard is inactive on these nominal trajectories; stiffness stress tests require rejection and smaller substeps. Fixed-cost unguarded midpoint is excluded.')
 LABEL['midpoint_4']=original_label
 work=['xpbd','work_lagged','work_column','work_row','work_secant','work_lagged_GS16']
 panel('work_normalization.png','Local work: lagged, row-matched, column-matched and potential-exact',PRIMARY,work,['energy','error','reservoir'],'All work proposals use the same 10% extra-motion cap. Strength = 1; the GS16 variant changes the frozen solve budget.\nMatching work sums or reducing a linear residual does not guarantee accurate motion or physical energy conservation.')
 reserves=['xpbd_ledger','work_lagged','reservoir_fast','reservoir_slow']
 panel('reservoirs.png','Recording work defects versus attempting repayment',PRIMARY,reserves,['signed','error','reservoir'],'Ledger only leaves XPBD positions unchanged. Repayment acts on local kinetic/potential work targets.\nPhysical + reservoir energy passes the independent ledger audit; that identity alone does not validate the physical trajectory.')
 endpoint=['xpbd','EL_bounded','LE_bounded','EL_rx_bounded','LE_rx_bounded']
 panel('retained_residual.png','Endpoint energy: retain r_x after reconstructing positions',PRIMARY,endpoint,['energy','motion','error'],'All corrected endpoint methods use the existing 50% repair / 10% motion bounds. Rebuilt rows are linearized.\nSolid EL/LE are the original frozen schedules; dashed variants rebuild every visited row and retain r_x.')
 panel('holdouts.png','Moving initial conditions expose energy-only corrections',['spring_orbit','pendulum_kick','double_mixed'],['xpbd','EL_bounded','LE_rx_bounded','reservoir_fast','midpoint_4'],['energy','motion','error'],'480 Hz · 10 seconds · validated refined references · physical energy and actual motion shown together.')
 # Both masses, rather than judging the double system from the first mass alone.
 fig,axs=plt.subplots(2,2,figsize=(16,10));methods=['xpbd','EL_bounded','reservoir_fast','midpoint_4'];decorate(fig,'Double spring pendulum: positions of both moving masses',methods,'480 Hz · 10 seconds · large-swing case. Solid curves are mass 1; dashed curves are mass 2.\nBlack dotted curves are the refined physical reference.',True)
 for axis,name in enumerate(['x','y']):
  for method in methods:
   a=groups['double_large',method,480]
   for i in [1,2]:axs[0,axis].plot(a[:,col['time']],a[:,col[name+str(i)]],color=COLOR[method],ls='-' if i==1 else '--',lw=1)
  a=groups['double_large','xpbd',480]
  for i in [1,2]:axs[0,axis].plot(a[:,col['time']],a[:,col['ref_'+name+str(i)]],':',color='#111',lw=1)
  axs[0,axis].set_xlabel('Time (s)');axs[0,axis].set_ylabel(name+' position (m)');axs[0,axis].grid(alpha=.2)
 for j,c in enumerate(['double_small','double_large']):
  for method in methods:
   a=groups[c,method,480]
   for i in [1,2]:axs[1,j].plot(a[:,col['x'+str(i)]],a[:,col['y'+str(i)]],color=COLOR[method],ls='-' if i==1 else '--',lw=1)
  a=groups[c,'xpbd',480]
  for i in [1,2]:axs[1,j].plot(a[:,col['ref_x'+str(i)]],a[:,col['ref_y'+str(i)]],':',color='#111',lw=1)
  axs[1,j].set_title(NAMES[c]);axs[1,j].set_xlabel('x (m)');axs[1,j].set_ylabel('y (m)');axs[1,j].set_aspect('equal',adjustable='datalim');axs[1,j].grid(alpha=.2)
 fig.savefig(args.output/'double_positions.png',dpi=150,facecolor='white');plt.close(fig)
 methods=['xpbd','EL_bounded','LE_rx_bounded','midpoint_1','midpoint_4','midpoint_drop_rx'];fig,axs=plt.subplots(3,3,figsize=(17,12));decorate(fig,'Small-step refinement: accuracy, energy and angular balance',methods,'10-second runs at 240 / 480 / 960 Hz. Angular balance includes the external gravity-torque integral.\nEnergy at roundoff is not an exact arithmetic guarantee; midpoint one-pass defects remain visible on the coupled spring system.')
 for i,c in enumerate(PRIMARY):
  for j,(metric,label,factor) in enumerate([('max_energy_rel','Maximum physical energy error (%)',100),('rms_position','Position RMS (m)',1),('max_angular_balance','Angular balance defect (kg m²/s)',1)]):
   for method in methods:
    rr=[summary[c,method,hz] for hz in [240,480,960]];axs[i,j].plot([240,480,960],[max(1e-15,factor*r[metric]) for r in rr],marker='o',**style(method))
   axs[i,j].set_xscale('log',base=2);axs[i,j].set_yscale('log');axs[i,j].set_xticks([240,480,960],['240','480','960']);axs[i,j].set_xlabel('Small steps per second');axs[i,j].set_ylabel(label);axs[i,j].set_title(NAMES[c],loc='left',fontweight='bold');axs[i,j].grid(alpha=.2)
 fig.savefig(args.output/'small_step_convergence.png',dpi=150,facecolor='white');plt.close(fig)
 cond=read(args.results/'conditioning.csv');valid=[r for r in cond if r['status']=='valid'];fig,axs=plt.subplots(2,2,figsize=(14,9));fig.subplots_adjust(top=.88,bottom=.14,hspace=.4,wspace=.24);fig.suptitle('Conditioning: scaling changes units, not the exact GS visits',fontsize=18,fontweight='bold')
 for ratio in [1,10,100,10000]:
  rr=[r for r in valid if r['case']=='double_mixed' and r['hz']=='240' and float(r['mass_ratio'])==ratio and r['gs_sweeps']=='16'];rr.sort(key=lambda r:float(r['stiffness_factor']))
  for ax,key in zip(axs.flat,['raw_condition','scaled_condition','scaled_residual','scaling_position_difference']):
   ax.loglog([float(r['stiffness_factor']) for r in rr],[max(1e-18,float(r[key])) for r in rr],marker='o',label=f'Mass ratio {ratio:g}');ax.set_xlabel('Spring stiffness multiplier');ax.grid(alpha=.2)
 for ax,label in zip(axs.flat,['Unscaled condition number (unit dependent)','Diagonal-scaled condition number','Scaled residual after 16 GS sweeps','Relative position difference: scaled vs unscaled']):ax.set_ylabel(label)
 axs[0,0].legend(frameon=False);fig.text(.08,.035,'Double spring, mixed state · 240 Hz · zero-gradient energy rows are flagged separately, without division by zero.\nThe residual-complete dual-gradient identity is also checked at nonzero material residuals across all conditioning fixtures.',fontsize=9)
 fig.savefig(args.output/'conditioning.png',dpi=150,facecolor='white');plt.close(fig)
 # Include all schedules and negative controls, including aborted runs.
 methods=list(dict.fromkeys(k[1] for k in summary));cases=list(NAMES);fig,axs=plt.subplots(1,2,figsize=(19,14));fig.subplots_adjust(top=.90,bottom=.12,left=.20,right=.96,wspace=.30)
 fig.suptitle('All 30 methods at 480 Hz: complete runs and failures',x=.20,y=.97,ha='left',fontsize=19,fontweight='bold')
 for ax,metric,factor,title,lo,hi in [(axs[0],'max_energy_rel',100,'Maximum physical energy error (%)',1e-10,1e5),(axs[1],'rms_position',100,'Position RMS (cm)',1e-3,1e4)]:
  values=np.array([[factor*summary[c,m,480][metric] if summary[c,m,480]['status']=='completed' else np.nan for c in cases] for m in methods]);cmap=plt.colormaps['viridis'].copy();cmap.set_bad('#ddd')
  im=ax.imshow(np.maximum(lo,values),aspect='auto',norm=LogNorm(lo,hi),cmap=cmap)
  ax.set_xticks(range(len(cases)),[NAMES[c].replace('Double spring: ','Double: ').replace(' with tangential velocity',' + velocity') for c in cases],rotation=40,ha='right',fontsize=9)
  ax.set_yticks(range(len(methods)),[LABEL.get(m,m.replace('_',' ')) for m in methods],fontsize=9);ax.set_title(title,loc='left',fontweight='bold')
  for i in range(len(methods)):
   for j in range(len(cases)):
    v=values[i,j];text='FAIL' if np.isnan(v) else f'{v:.2g}';color='#111' if np.isnan(v) else 'white' if math.log10(max(lo,v))<.5*(math.log10(lo)+math.log10(hi)) else '#111';ax.text(j,i,text,ha='center',va='center',fontsize=7.5,color=color)
  fig.colorbar(im,ax=ax,fraction=.025,pad=.015)
 fig.text(.20,.025,'Every numerical cell represents a complete 10-second run. FAIL means the run hit the safety limit; accepted-prefix errors are not presented as full-run accuracy.\nBoth panels use logarithmic colors. All summaries, attempted energies on failure, and 240/960 Hz runs are preserved in the CSV.',fontsize=9)
 fig.savefig(args.output/'all_methods.png',dpi=150,facecolor='white');plt.close(fig)
 info['conditioning']={'rows':len(cond),'singular_rows':sum(r['status']!='valid' for r in cond),'max_scaling_position_difference':max(float(r['scaling_position_difference']) for r in valid),'max_dual_gradient_relative_error':max(float(r['dual_gradient_relative_error']) for r in valid)}
 if args.long_results:
  long_summary=summaries(args.long_results);long_groups,long_col=traces(args.long_results);info['long']=audit(long_groups,long_col,long_summary,models,False);info['long']['runs']=len(long_summary);info['long']['uncertified_cases']=sorted({r['case'] for r in long_summary.values() if not r['reference_valid']})
  fig,axs=plt.subplots(3,3,figsize=(17,12));methods=['xpbd','EL_bounded','LE_rx_bounded','reservoir_fast','midpoint_4'];decorate(fig,'60-second check: energy and motion over many oscillations',methods,'480 Hz · traces shown every 8th step; summary metrics use every step. The small double-spring reference passes refinement.\nLarge double-spring long-run RMS is not certified; this figure uses the validated small-swing case.')
  for i,c in enumerate(['spring_release','pendulum_rest','double_small']):
   for j,metric in enumerate(['energy','motion','error']):graph(axs[i,j],c,methods,metric,long_groups,long_col);axs[i,j].set_title(NAMES[c],loc='left',fontweight='bold')
  fig.savefig(args.output/'long_horizon.png',dpi=150,facecolor='white');plt.close(fig)
 # Independent ledger-only trajectory equality, not just an energy identity.
 equality=0.
 for c in NAMES:
  a=groups[c,'xpbd',480];b=groups[c,'xpbd_ledger',480]
  for k in ['x1','y1','vx1','vy1','x2','y2','vx2','vy2']:equality=max(equality,float(np.max(np.abs(a[:,col[k]]-b[:,col[k]]))))
 assert equality==0;info['ledger_only_state_difference']=equality
 (args.output/'audit.json').write_text(json.dumps(info,indent=2)+'\n')
 (args.output/'selected_results.json').write_text(json.dumps([r for r in summary.values() if r['hz']==480 and r['method'] in MAIN],indent=2)+'\n')
 print('Saved',len(list(args.output.glob('*.png'))),'figures; audit:',args.output/'audit.json',flush=True)
if __name__=='__main__':main()
