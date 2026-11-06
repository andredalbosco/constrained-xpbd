from pathlib import Path
import csv, math, json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parent
OUT=ROOT/'figures';OUT.mkdir(exist_ok=True)
def read(path):
    with open(ROOT/path) as f:return list(csv.DictReader(f))
def num(rows,key):return np.array([float(r[key]) for r in rows])
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
colors={'obj_midpotential_b1':'#d1495b','obj_quarter_b1':'#00798c','canonical_dg':'#edae49','sym_energy':'#6a4c93','obj_ls_dg_E':'#30638e'}
labels={'obj_midpotential_b1':'Midpoint potential, no old reaction','obj_quarter_b1':'Quarter objective, old reaction','canonical_dg':'Canonical discrete gradient','sym_energy':'Symmetric energy projection','obj_ls_dg_E':'Mechanical residual + energy'}
fig,ax=plt.subplots(1,3,figsize=(15,4.3))
for panel,c in enumerate(['double_moderate','triple_moderate']):
    rows=read('objectives_order_results/summary.csv')
    for m in ['obj_midpotential_b1','obj_quarter_b1','obj_ls_dg_E']:
        r=sorted([x for x in rows if x['case']==c and x['method']==m],key=lambda x:float(x['h']))
        ax[panel].loglog(num(r,'h'),num(r,'position_rms'),'o-',label=labels[m],color=colors[m])
    ax[panel].set(title=c.replace('_',' ').title()+' — 1 second',xlabel='Step size (s)',ylabel='Position RMS (m)');ax[panel].grid(alpha=.2,which='both')
for m in ['obj_midpotential_b1','obj_quarter_b1','canonical_dg','sym_energy']:
    r=read('modal_transfer/'+m+'.csv')
    ax[2].semilogy(num(r,'t'),np.maximum(num(r,'mode0_energy')+num(r,'mode1_energy'),1e-12),color=colors[m],label=labels[m])
r=read('modal_transfer/reference.csv');ax[2].semilogy(num(r,'t'),np.maximum(num(r,'mode0_energy')+num(r,'mode1_energy'),1e-12),'--',color='black',lw=1,label='Refined reference')
ax[2].set(title='Triple pendulum — highest initial mode',xlabel='Time (s)',ylabel='Fraction of initial modal energy in other modes',ylim=(1e-10,2));ax[2].grid(alpha=.2,which='both')
unique={}
for a in ax:
    handles,ls=a.get_legend_handles_labels()
    for handle,label in zip(handles,ls):unique.setdefault(label,handle)
fig.legend(list(unique.values()),list(unique),loc='lower center',ncol=3,bbox_to_anchor=(.5,-.04));fig.tight_layout(rect=(0,.10,1,1));fig.savefig(OUT/'reaction_and_mode_transfer.png',dpi=180,bbox_inches='tight');fig.savefig(OUT/'reaction_and_mode_transfer.pdf',bbox_inches='tight');plt.close(fig)

fig,ax=plt.subplots(1,3,figsize=(15,4))
r=read('momentum_results/summary.csv');chosen=['rattle','cart_dg','canonical_midpoint','canonical_dg','canonical_dg_momentum','canonical_avf_momentum','sym_energy','sym_energy_momentum']
for m in chosen:
    rr=[x for x in r if x['method']==m and x['case']=='double_free' and abs(float(x['h'])-1/30)<1e-8][0]
    ax[0].scatter(max(float(rr['max_energy_relative']),1e-15),max(float(rr['max_angular_change']),1e-15),label=m,s=45)
ax[0].set(xscale='log',yscale='log',xlabel='Max relative energy error',ylabel='Max scaled angular-momentum change',title='Free double chain, 30 Hz, 10 s');ax[0].grid(alpha=.2)
r=read('hidden_modes.csv')
for c,col in [('double_moderate','#00798c'),('triple_moderate','#d1495b')]:
    rr=[x for x in r if x['case']==c and x['source']=='genuine_previous_step']
    ax[1].loglog(num(rr,'h'),num(rr,'doubling_gap'),'o-',color=col,label=c)
    ax[1].axhline(float(rr[0]['predicted_doubling_limit']),color=col,ls='--',lw=1)
ax[1].set(xlabel='Trial step h (s)',ylabel='Scaled one-step / two-half-step discrepancy',title='Inherited Cartesian DG normal mode');ax[1].grid(alpha=.2,which='both');ax[1].legend(fontsize=8)
r=read('accuracy_results/summary.csv')
for c,col in [('double_moderate','#00798c'),('triple_moderate','#d1495b'),('triple_contrast','#6a4c93')]:
    for m,mark in [('rattle_accuracy','o'),('sym_energy_accuracy','s'),('canonical_dg_accuracy','^')]:
        rr=[x for x in r if x['case']==c and x['method']==m][0]
        ax[2].scatter(float(rr['residual_evaluations']),float(rr['position_rms']),color=col,marker=mark,s=55,label=f'{c}: {m}')
ax[2].set(xscale='log',yscale='log',xlabel='Residual evaluations (whole trajectory)',ylabel='Position RMS (m)',title='Accuracy subdivision, 30-Hz output, 10 s');ax[2].grid(alpha=.2)
ax[0].legend(fontsize=7,loc='best');ax[2].text(.03,.04,'Circles: RATTLE; squares: symmetric projection\nTriangles: canonical DG\nColours: double / triple / contrast triple',transform=ax[2].transAxes,fontsize=8)
fig.tight_layout();fig.savefig(OUT/'invariants_hidden_modes_and_cost.png',dpi=180);fig.savefig(OUT/'invariants_hidden_modes_and_cost.pdf');plt.close(fig)

fig,ax=plt.subplots(1,2,figsize=(11,4))
r=read('modal_transfer/summary.csv');weights=[0,.5,1,1.5,2]
for beta,col in [(.1,'#00798c'),(1,'#d1495b'),(10,'#edae49')]:
    rr=[next(x for x in r if x['method']==f'obj_reaction{c:.6f}_b{beta:.6f}') for c in weights]
    ax[0].semilogy(weights,np.maximum(num(rr,'max_other_mode_fraction'),1e-12),'o-',color=col,label=f'beta = {beta:g}')
    ax[1].semilogy(weights,num(rr,'max_energy_relative'),'o-',color=col,label=f'beta = {beta:g}')
ax[0].set(xlabel='Fraction of the old reaction included',ylabel='Max energy fraction in other modes',title='Changing compliance does not repair omitted reaction')
ax[1].set(xlabel='Fraction of the old reaction included',ylabel='Max relative total-energy error',title='All variants conserve total energy closely')
for a in ax:a.grid(alpha=.2,which='both');a.legend()
fig.tight_layout();fig.savefig(OUT/'reaction_weight_beta.png',dpi=180);fig.savefig(OUT/'reaction_weight_beta.pdf');plt.close(fig)

folders=['results','order_results','accuracy_results','stress_results','long_results','momentum_results','metric_results','objectives_results','objectives_order_results','reaction_order_results','equilibria_results','equilibria_stable_results','momentum_stable_results']
summary=[]
for folder in folders:
    rows=read(folder+'/summary.csv');summary.append({'suite':folder,'runs':len(rows),'complete':sum(r['status']=='complete' for r in rows),'uncertified':sum(r['reference_certified']!='1' for r in rows)})
(ROOT/'counts.json').write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps(summary,indent=2))
