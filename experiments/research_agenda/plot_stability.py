#!/usr/bin/env python3
import csv,json,os
from pathlib import Path
os.environ.setdefault('MPLCONFIGDIR','/tmp/xpbd-agenda-matplotlib')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
root=Path(__file__).resolve().parent
rows=list(csv.DictReader((root/'results/stability.csv').open()))
eq=list(csv.DictReader((root/'results/acceptance_equivalence.csv').open()))
methods=['xpbd','midpoint_4','guarded_midpoint'];info={}
for m in methods:
 rr=[r for r in rows if r['method']==m];info[m]={'runs':len(rr),'completed':sum(r['status']=='completed' for r in rr),'growth_or_failure':sum(float(r['max_energy_ratio'])>1.01 or r['status']!='completed' for r in rr),'maximum_absolute_energy_error':max(float(r['max_absolute_energy_error']) for r in rr),'maximum_substeps_per_step':max(int(r['substeps'])/int(r['steps']) for r in rr),'deepest':max(int(r['deepest']) for r in rr)}
assert info['guarded_midpoint']['completed']==96 and info['guarded_midpoint']['maximum_absolute_energy_error']<=1e-6
assert all(float(r['max_state_difference'])==0 and int(r['rejections'])==0 for r in eq)
info['nominal_equivalence']=eq
fig,axs=plt.subplots(2,2,figsize=(15,10));fig.subplots_adjust(top=.85,bottom=.12,wspace=.28,hspace=.4)
fig.suptitle('Stability screen: fixed four-pass midpoint is excluded',x=.08,y=.97,ha='left',fontsize=20,fontweight='bold')
fig.text(.08,.91,'96 scenarios per method · spring / orbit / pendulum / double spring · stiffness ×1 to ×10,000 · mass ratio to 100',fontsize=10)
cases=['spring_release','spring_orbit','pendulum_rest','double_large'];names=['Spring release','Spring orbit','Rigid pendulum','Double spring'];x=np.arange(4)
for k,(method,label,color) in enumerate([('xpbd','XPBD','#8b969e'),('midpoint_4','Fixed four-pass: rejected candidate','#bf5f42'),('guarded_midpoint','Guarded midpoint','#176d44')]):
 count=[sum(r['case']==c and (float(r['max_energy_ratio'])>1.01 or r['status']!='completed') for r in rows if r['method']==method) for c in cases]
 axs[0,0].bar(x+(k-1)*.23,count,width=.23,color=color,label=label)
axs[0,0].set_xticks(x,names);axs[0,0].set_ylabel('Runs with >1% energy growth or failure');axs[0,0].set_ylim(0,10);axs[0,0].legend(fontsize=8,frameon=False);axs[0,0].grid(axis='y',alpha=.2)
for hz,ls in [(240,'-'),(960,'--')]:
 for ratio,color in [(1,'#2585af'),(100,'#176d44')]:
  select=[r for r in rows if r['method']=='guarded_midpoint' and int(r['hz'])==hz and float(r['light_mass_ratio'])==ratio]
  xx=[1,100,10000];err=[max(float(r['max_absolute_energy_error']) for r in select if float(r['stiffness_factor'])==s) for s in xx];cost=[max(int(r['substeps'])/int(r['steps']) for r in select if float(r['stiffness_factor'])==s) for s in xx]
  axs[0,1].loglog(xx,err,ls=ls,color=color,marker='o',label=f'{hz} Hz, mass ratio {ratio}')
  axs[1,0].loglog(xx,cost,ls=ls,color=color,marker='o')
axs[0,1].axhline(1e-6,color='#bf5f42',ls=':',label='Accepted-state bound');axs[0,1].set_xlabel('Stiffness multiplier');axs[0,1].set_ylabel('Guarded: maximum |E−E_initial| / E_scale');axs[0,1].legend(fontsize=8,frameon=False);axs[0,1].grid(alpha=.2)
axs[1,0].set_xlabel('Stiffness multiplier');axs[1,0].set_ylabel('Guarded: accepted substeps / requested step');axs[1,0].grid(alpha=.2)
axs[1,1].axis('off');axs[1,1].text(0,1,'Acceptance checks',fontsize=15,fontweight='bold',va='top')
axs[1,1].text(0,.86,'Global physical energy error ≤ 10⁻⁶ of energy scale\nLocal kinetic/potential work defect ≤ 10⁻⁹\nPrimal and constitutive residual checks\nRod length error ≤ 10⁻⁹ × rest length\nFinite positions and velocities\n\nReject → retry as two half steps; depth limit 12.\nExhausted retries return failure without committing a state.\n\n96/96 guarded stress runs completed.\nNominal 60-second trajectories: guard changes no states.\nThis does not guarantee convergence for every problem.',fontsize=11,va='top',linespacing=1.6)
fig.text(.08,.025,'One-second stress runs screen boundedness, not high-stiffness trajectory accuracy. Maximum guarded subdivision depth: 7.\nExact implicit midpoint is stable for a linear oscillator; that statement does not establish stability of a finite GS solve of nonlinear equations.',fontsize=9)
out=root/'plots';out.mkdir(exist_ok=True);fig.savefig(out/'stability.png',dpi=160,facecolor='white');(out/'stability_audit.json').write_text(json.dumps(info,indent=2)+'\n');print(json.dumps(info,indent=2))
