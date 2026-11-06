#!/usr/bin/env python3
"""Audit backward-Euler candidates and render physical energy/motion comparisons."""
import argparse,csv,json,math,os
from array import array
from collections import defaultdict
from pathlib import Path
os.environ.setdefault('MPLCONFIGDIR','/tmp/xpbd-implicit-matplotlib')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.colors import LogNorm
import numpy as np
NAMES={'spring_release':'Spring release','pendulum_rest':'Rigid pendulum','double_small':'Double spring: small swing','double_large':'Double spring: large swing','spring_orbit':'Spring orbit','pendulum_kick':'Kicked pendulum','double_mixed':'Double spring: mixed motion'}
LABEL={'xpbd':'XPBD','xpbd_4':'XPBD, 4 passes','EL_soft_10':'Compliant EL, energy only','LE_soft_10':'Compliant LE, energy only','EL_angular_10':'Compliant EL + angular balance','LE_angular_10':'Compliant LE + angular balance','work_soft_10':'Compliant work, no angular row','work_angular':'Compliant work + angular, no repayment','work_repay':'Work + repayment, no angular row','work_repay_angular':'Work + angular, τ=0.25 s','workL_s1_tau005':'Work + angular + repayment (strength 1)','workL_s10_tau005':'Work + angular + repayment (strength 10)','workL_soft_01':'Work + repayment, soft angular (0.1)','workL_soft_1':'Work + repayment, soft angular (1)','workL_soft_10':'Work + repayment, soft angular (10)'}
COLORS={'xpbd':'#89959e','xpbd_4':'#89959e','EL_angular_10':'#d19024','work_soft_10':'#a56aab','work_angular':'#7f68a0','work_repay':'#cc6573','work_repay_angular':'#458896','workL_s1_tau005':'#176d44','workL_s10_tau005':'#1c8c9e','workL_soft_01':'#d16a29','workL_soft_1':'#ca9a4a','workL_soft_10':'#b5a561'}
def read(path):
 with path.open() as f:return list(csv.DictReader(f))
def summary(path):return {(r['case'],r['method'],int(r['hz'])):{k:v if k in ('case','method','status') else float(v) for k,v in r.items()} for r in read(path/'summary.csv')}
def traces(path):
 groups=defaultdict(lambda:array('d'))
 with (path/'traces.csv').open() as f:
  reader=csv.reader(f);head=next(reader);col={k:i for i,k in enumerate(head[3:])}
  for r in reader:groups[r[0],r[1],int(r[2])].extend(map(float,r[3:]))
 return {k:np.asarray(v).reshape(-1,len(col)) for k,v in groups.items()},col
def audit(groups,col,stats,models,full):
 result={'groups':len(groups),'rows':0,'max_energy_reconstruction':0.,'max_rms_reconstruction':0.,'max_energy_gain':0.,'max_ledger_error':0.,'max_kinematic_error':0.}
 for key,a in groups.items():
  s=stats[key];m=models[key[0]];n=int(m['count']);at=lambda k:a[:,col[k]];E=np.zeros(len(a));squared=np.zeros(len(a));kin=0.
  assert at('step')[0]==0 and at('step')[-1]==s['completed'],key
  for i in range(1,n+1):
   x,y,vx,vy=[at(k+str(i)) for k in ['x','y','vx','vy']];mass=m['m'+str(i)];E+=mass*(.5*(vx*vx+vy*vy)-m['gx']*x-m['gy']*y)
   dx=x-(at('x'+str(i-1)) if i>1 else 0);dy=y-(at('y'+str(i-1)) if i>1 else 0);alpha=m['alpha'+str(i)]
   if alpha>0:E+=(np.hypot(dx,dy)-m['l'+str(i)])**2/(2*alpha)
   squared+=(x-at('ref_x'+str(i)))**2+(y-at('ref_y'+str(i)))**2
   if full:kin=max(kin,float(np.max(np.abs(np.diff(x)*key[2]-vx[1:]))),float(np.max(np.abs(np.diff(y)*key[2]-vy[1:]))))
  err=float(np.max(np.abs(E-at('energy'))))/s['energy_scale'];assert err<1e-8,(key,err)
  gain=float(np.max((E-s['initial_energy'])/s['energy_scale']));assert gain<=1.05e-10,(key,gain)
  if full:
   assert np.all(np.diff(at('step'))==1)
   rms=math.sqrt(float(np.mean(squared[1:]/n))) if len(a)>1 else 0;delta=abs(rms-s['rms_position']);assert delta<1e-10,(key,delta)
   maximum=float(np.max(np.abs(E[1:]-s['initial_energy'])))/s['energy_scale'];assert math.isclose(maximum,s['max_energy_rel'],abs_tol=1e-10,rel_tol=1e-8)
   result['max_rms_reconstruction']=max(result['max_rms_reconstruction'],delta)
  if key[1].startswith(('work','xpbd')):
   ledger=float(np.max(np.abs(E+at('reservoir')-s['initial_energy'])))/s['energy_scale'];assert ledger<1e-8,(key,ledger);result['max_ledger_error']=max(result['max_ledger_error'],ledger)
  result['max_energy_reconstruction']=max(result['max_energy_reconstruction'],err);result['max_energy_gain']=max(result['max_energy_gain'],gain);result['max_kinematic_error']=max(result['max_kinematic_error'],kin);result['rows']+=len(a)
 return result
def main():
 p=argparse.ArgumentParser();p.add_argument('results',type=Path);p.add_argument('--long-results',type=Path,required=True);p.add_argument('--stress-long',type=Path,required=True);p.add_argument('--output',type=Path,required=True);args=p.parse_args();args.output.mkdir(parents=True,exist_ok=True)
 models={r['case']:{k:float(v) for k,v in r.items() if k!='case'} for r in read(args.results/'cases.csv')};stats=summary(args.results);groups,col=traces(args.results);long_stats=summary(args.long_results);long_groups,long_col=traces(args.long_results)
 info={'short':audit(groups,col,stats,models,True),'long':audit(long_groups,long_col,long_stats,models,False),'short_runs':len(stats),'long_runs':len(long_stats)}
 stress=read(args.results/'stability.csv');stress_long=read(args.stress_long/'stability.csv');methods=list(dict.fromkeys(r['method'] for r in stress));info['stability']={}
 for method in methods:
  a=[r for r in stress if r['method']==method];b=[r for r in stress_long if r['method']==method]
  info['stability'][method]={'one_second_completed':sum(r['status']=='completed' for r in a),'ten_second_completed':sum(r['status']=='completed' for r in b),'per_duration_runs':len(a),'max_accepted_energy_gain':max(float(r['max_energy_gain']) for r in a+b)}
 info['short_failures']=[{k:r[k] for k in ['case','method','hz','completed','steps']} for r in stats.values() if r['status']!='completed'];info['long_failures']=[{k:r[k] for k in ['case','method','completed','steps']} for r in long_stats.values() if r['status']!='completed']
 plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False})
 def style(method):return {'color':COLORS[method],'ls':'--' if method in ('workL_s10_tau005','xpbd_4') else '-','lw':1.7 if method=='workL_s1_tau005' else 1.1}
 def decorate(fig,title,methods,foot,reference=False):
  fig.subplots_adjust(top=.84,bottom=.10,left=.075,right=.98,wspace=.29,hspace=.43);fig.suptitle(title,x=.075,y=.98,ha='left',fontweight='bold',fontsize=19)
  handles=[Line2D([],[],label=LABEL[m],**style(m)) for m in methods]
  if reference:handles+=[Line2D([],[],label='Physical reference',color='#111',ls=':')]
  fig.legend(handles=handles,loc='upper left',bbox_to_anchor=(.07,.945),ncol=3,frameon=False,fontsize=9);fig.text(.075,.025,foot,fontsize=9,color='#333')
 def draw(ax,case,methods,metric,source=groups,index=col,angular=False):
  for method in methods:
   a=source[case,method,480];s=stats[case,method,480];at=lambda k:a[:,index[k]]
   if metric=='energy':v=100*np.abs(at('energy')-s['initial_energy'])/s['energy_scale']
   elif metric=='error':v=np.sqrt(sum((at('x'+str(i))-at('ref_x'+str(i)))**2+(at('y'+str(i))-at('ref_y'+str(i)))**2 for i in range(1,int(s['count'])+1))/s['count'])
   elif metric=='angular':v=np.abs(at('angular')-at('angular')[0]-at('torque_impulse'))
   elif case.startswith('pendulum'):v=np.arctan2(at('x1'),-at('y1'))
   else:v=at('x1')-(models[case]['l1'] if case=='spring_release' else 0)
   ax.plot(at('time'),v,**style(method))
  ax.grid(alpha=.2);ax.set_xlabel('Time (s)')
  if metric=='energy':ax.set_yscale('symlog',linthresh=1e-3);ax.set_ylim(bottom=0);ax.set_ylabel('Physical energy error (%) · symlog')
  elif metric=='error':ax.set_yscale('symlog',linthresh=1e-5);ax.set_ylim(bottom=0);ax.set_ylabel('Position error vs reference (m)')
  elif metric=='angular':ax.set_yscale('symlog',linthresh=1e-8);ax.set_ylim(bottom=0);ax.set_ylabel('Angular balance defect (kg m²/s)')
  else:
   a=source[case,'xpbd',480];at=lambda k:a[:,index[k]];v=np.arctan2(at('ref_x1'),-at('ref_y1')) if case.startswith('pendulum') else at('ref_x1')-(models[case]['l1'] if case=='spring_release' else 0)
   ax.plot(at('time'),v,':',color='#111',lw=1.1);ax.set_ylabel('Angle from vertical (rad)' if case.startswith('pendulum') else 'x − rest length (m)' if case=='spring_release' else 'First mass x (m)')
 def panel(name,title,cases,methods,metrics,foot,source=groups,index=col):
  fig,axs=plt.subplots(len(cases),len(metrics),figsize=(17,12),squeeze=False);decorate(fig,title,methods,foot,'motion' in metrics)
  for i,c in enumerate(cases):
   for j,metric in enumerate(metrics):draw(axs[i,j],c,methods,metric,source,index);axs[i,j].set_title(NAMES[c],loc='left',fontweight='bold')
  fig.savefig(args.output/name,dpi=150,facecolor='white');plt.close(fig)
 candidates=['xpbd','EL_angular_10','work_angular','workL_s1_tau005','workL_s10_tau005']
 for method in candidates:assert info['stability'][method]['one_second_completed']==info['stability'][method]['ten_second_completed']==96
 panel('implicit_candidates.png','Backward-Euler candidates: compliant work, angular balance and repayment',['spring_release','pendulum_rest','double_large'],candidates,['energy','motion','error'],'480 Hz · 10 seconds · endpoint velocity reconstruction · physical energy excludes reservoirs. All shown methods pass both stiffness/perturbation screens.\nAccepted steps obey a physical-energy ceiling; this is not a proof of unconditional solve completion or exact physical energy conservation.')
 panel('recommended_candidate.png','Implicit local work balance: the strongest common setting in this sweep',['spring_release','pendulum_rest','double_large'],['xpbd','workL_s1_tau005'],['energy','motion','error'],'480 Hz · 10 seconds · work-balance strength 1 · reservoir repayment τ=0.05 s · hard angular-balance row.\nBoth methods retain endpoint velocity reconstruction. The correction reduces physical damping; it does not eliminate it. Accepted states satisfy an energy ceiling.')
 panel('moving_and_coupled.png','Orbit and coupled motion: why angular balance matters',['spring_orbit','double_small','double_mixed'],['xpbd','work_soft_10','work_repay','work_angular','workL_s1_tau005','workL_s10_tau005'],['energy','motion','error'],'480 Hz · 10 seconds · every method uses the backward-Euler objective, compliant rows, linearization and scalar Gauss–Seidel.\nEnergy repair without angular control can strongly distort motion. Black dotted curves are independent refined references.')
 panel('angular_compliance.png','Angular compliance trades torque balance for energy correction',['pendulum_rest','spring_orbit','double_small'],['xpbd','workL_s1_tau005','workL_soft_01','workL_soft_1','workL_soft_10'],['energy','error','angular'],'All repayment variants here use work strength 1 and τ=0.05 s. Softer angular balance helps the rigid pendulum but worsens orbit dynamics.\nHard length plus hard angular balance restricts the rigid pendulum\'s next position; energy cannot be adjusted independently with endpoint kinematics.')
 panel('long_horizon.png','60-second check of the surviving implicit work candidates',['spring_release','pendulum_rest','double_small'],['xpbd','workL_s1_tau005','workL_s10_tau005','workL_soft_01'],['energy','motion','error'],'480 Hz · summary metrics use every step; traces shown every 8th step. These three long-run references pass refinement.\nLong large-swing double-pendulum RMS remains uncertified. Conserved augmented energy does not eliminate damping or phase error.',long_groups,long_col)
 # Frequency refinement and consistent per-step row-visit accounting.
 fig,axs=plt.subplots(3,3,figsize=(17,12));ms=['xpbd','work_angular','workL_s1_tau005','workL_s10_tau005'];decorate(fig,'Small-step refinement of implicit work corrections',ms,'Each curve uses the same settings at 240/480/960 Hz; no per-case fit. Angular balance includes gravity torque.\nCompliant energy/work corrections still dissipate physical energy; reservoirs record and attempt repayment of accepted work defects.')
 for i,c in enumerate(['spring_release','pendulum_rest','double_small']):
  for j,(metric,factor,label) in enumerate([('max_energy_rel',100,'Maximum physical energy error (%)'),('rms_position',100,'Position RMS (cm)'),('max_angular_balance',1,'Angular balance defect (kg m²/s)')]):
   for method in ms:axs[i,j].loglog([240,480,960],[max(1e-14,factor*stats[c,method,h][metric]) for h in [240,480,960]],marker='o',**style(method))
   axs[i,j].set_xticks([240,480,960],['240','480','960']);axs[i,j].set_xlabel('Small steps per second');axs[i,j].set_ylabel(label);axs[i,j].set_title(NAMES[c],loc='left',fontweight='bold');axs[i,j].grid(alpha=.2)
 fig.savefig(args.output/'refinement.png',dpi=150,facecolor='white');plt.close(fig)
 # All configurations, with failed attempts excluded from full-run metrics.
 methods_all=list(dict.fromkeys(k[1] for k in stats));cases=list(NAMES);fig,axs=plt.subplots(1,2,figsize=(20,16));fig.subplots_adjust(top=.92,bottom=.13,left=.20,right=.96,wspace=.30);fig.suptitle('All 38 implicit configurations at 480 Hz',x=.20,y=.97,ha='left',fontsize=19,fontweight='bold')
 for ax,key,factor,label,lo,hi in [(axs[0],'max_energy_rel',100,'Physical energy error (%)',1e-3,100),(axs[1],'rms_position',100,'Position RMS (cm)',.01,100)]:
  values=np.array([[factor*stats[c,m,480][key] if stats[c,m,480]['status']=='completed' else np.nan for c in cases] for m in methods_all]);cmap=plt.colormaps['viridis'].copy();cmap.set_bad('#ddd');im=ax.imshow(np.maximum(lo,values),aspect='auto',norm=LogNorm(lo,hi),cmap=cmap)
  ax.set_xticks(range(7),[NAMES[c] for c in cases],rotation=40,ha='right',fontsize=8);ax.set_yticks(range(len(methods_all)),[LABEL.get(m,m.replace('_',' ')) for m in methods_all],fontsize=8);ax.set_title(label,loc='left',fontweight='bold')
  for i in range(len(methods_all)):
   for j in range(7):v=values[i,j];ax.text(j,i,'FAIL' if np.isnan(v) else f'{v:.2g}',ha='center',va='center',fontsize=7,color='#111' if np.isnan(v) or v>math.sqrt(lo*hi) else 'white')
  fig.colorbar(im,ax=ax,fraction=.025,pad=.015)
 fig.text(.20,.025,'FAIL means no admissible next state was found; no violating state was committed. Failed accepted prefixes are not treated as full-run accuracy.\nCompletion at nominal stiffness does not imply passing the separate stiffness screen. These tables include excluded controls for auditability.',fontsize=9)
 fig.savefig(args.output/'all_methods.png',dpi=150,facecolor='white');plt.close(fig)
 # Stability/completion and actual solver visits, not an I/O-contaminated timing benchmark.
 fig,axs=plt.subplots(1,2,figsize=(17,11));fig.subplots_adjust(top=.88,bottom=.11,left=.30,right=.96,wspace=.30);fig.suptitle('Implicit candidates: safe completion and solve budget',x=.08,y=.98,ha='left',fontsize=19,fontweight='bold');y=np.arange(len(methods))
 a=[info['stability'][m]['one_second_completed'] for m in methods];b=[info['stability'][m]['ten_second_completed'] for m in methods];axs[0].barh(y-.18,a,height=.35,color='#8db6a7',label='1-second stress');axs[0].barh(y+.18,b,height=.35,color='#176d44',label='10-second stress');axs[0].set_yticks(y,[LABEL.get(m,m.replace('_',' ')) for m in methods],fontsize=8);axs[0].invert_yaxis();axs[0].set_xlim(0,100);axs[0].set_xlabel('Safely completed scenarios / 96');axs[0].legend(fontsize=9,frameon=False);axs[0].grid(axis='x',alpha=.2)
 visits=[stats['double_small',m,480]['scalar_visits']/stats['double_small',m,480]['completed'] for m in methods];axs[1].barh(y,visits,color='#587b9a');axs[1].set_yticks(y,[]);axs[1].invert_yaxis();axs[1].set_xlabel('Scalar visits / step (small double spring)');axs[1].grid(axis='x',alpha=.2)
 fig.text(.08,.025,'Stress cases vary stiffness ×1/100/10,000, light-mass ratio ×1/100, 240/960 Hz and an angle perturbation of 10⁻⁸.\nNo adaptive substeps are used. Backtracking blends corrections toward XPBD; inability to advance is recorded explicitly. Visits include the fallback baseline calculation.',fontsize=9)
 fig.savefig(args.output/'stability_and_cost.png',dpi=150,facecolor='white');plt.close(fig)
 # Baseline reproduction against the previous independent study.
 old={(r['case'],int(r['hz'])):r for r in read(args.results.parent.parent/'traditional_cases/results/summary.csv') if r['method']=='xpbd'};baseline_error=0.
 for (c,m,h),r in stats.items():
  if m=='xpbd':
   for field in ['max_energy_rel','rms_position','max_angular_balance']:baseline_error=max(baseline_error,abs(r[field]-float(old[c,h][field])))
 assert baseline_error==0;info['baseline_reproduction_error']=baseline_error
 chosen='workL_s1_tau005';ratios=[]
 for (c,m,h),r in stats.items():
  if m==chosen:ratios.append({'case':c,'hz':h,'rms_ratio':r['rms_position']/stats[c,'xpbd',h]['rms_position']})
 info['default_candidate']=chosen;info['default_candidate_short_improvements']=sum(r['rms_ratio']<1 for r in ratios);info['default_candidate_rms_ratios']=ratios
 info['long_uncertified_cases']=sorted({r['case'] for r in long_stats.values() if not r['reference_valid']})
 (args.output/'audit.json').write_text(json.dumps(info,indent=2)+'\n');print(json.dumps({k:v for k,v in info.items() if k not in ('stability','short_failures','default_candidate_rms_ratios')},indent=2));print('Figures:',len(list(args.output.glob('*.png'))))
if __name__=='__main__':main()
