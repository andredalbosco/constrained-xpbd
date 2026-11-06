#!/usr/bin/env python3
"""Audit cold-start XPBD, energy schedules, and two-mass spring-pendulum traces."""
import argparse
import collections
import csv
import hashlib
import json
import math
import os
from pathlib import Path
os.environ.setdefault('MPLCONFIGDIR','/tmp/xpbd-traditional-matplotlib')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.colors import ListedColormap,BoundaryNorm
from matplotlib.patches import Rectangle
import numpy as np

CASES=['spring_release','spring_orbit','pendulum_rest','pendulum_kick','double_small','double_large','double_mixed']
NAMES={'spring_release':'Released spring','spring_orbit':'Orbiting spring','pendulum_rest':'Rigid pendulum from rest',
       'pendulum_kick':'Rigid pendulum with velocity','double_small':'Double spring pendulum: small swing',
       'double_large':'Double spring pendulum: large swing','double_mixed':'Double spring pendulum: mixed motion'}
METHODS=['xpbd','frozen_xpbd','EL_raw','EL_direction','EL_repair_cap','LE_raw','LE_direction','LE_repair_cap',
         'jacobi_raw','jacobi_direction','jacobi_repair_cap','joint_raw','joint_direction','joint_repair_cap',
         'EL_raw_gate','EL_raw_gate_cap','LEL_repair_cap','ELE_repair_cap']
SELECTED=['xpbd','frozen_xpbd','EL_raw_gate','EL_raw_gate_cap','EL_repair_cap','LE_repair_cap','joint_repair_cap']
COLORS={'xpbd':'#8a949e','frozen_xpbd':'#4d6074','EL_raw_gate':'#c87524','EL_raw_gate_cap':'#b13c60',
        'EL_repair_cap':'#387daf','LE_repair_cap':'#138674','joint_repair_cap':'#7d57a5'}
LABELS={'xpbd':'Ordinary XPBD: fresh spring rows','frozen_xpbd':'Material GS: frozen rows',
        'EL_raw_gate':'EL: raw energy-direction gate','EL_raw_gate_cap':'EL: energy gate + cap',
        'EL_repair_cap':'EL: repair + cap','LE_repair_cap':'LE: repair + cap',
        'joint_repair_cap':'Joint: repair + cap'}
PALETTE=['#eaf3e8','#cde3c7','#f0efbf','#f3d9a3','#eeb692','#dd8b88','#bd5e79']
BOUNDS=[0,.1,1,5,20,50,100,1e300]


def read(path):
    with path.open(newline='') as f:
        for r in csv.DictReader(f):
            for k in r:
                if k not in {'case','method','status'}:
                    r[k]=float(r[k]);assert math.isfinite(r[k]),(path,k,r)
            yield r


def key(r):return r['case'],r['method'],int(r['hz'])
def vals(tr,f):return np.asarray([r[f] for r in tr])
def complete(r):return r['completed']==r['steps']
def metrics(r):return 100*r['max_energy_rel'],100*r['rms_position']/r['reference_excursion']
def label(m):
    if m in LABELS:return LABELS[m]
    for suffix,title in [('repair_cap','repair + cap'),('direction','combined direction gate'),('raw','raw')]:
        if m.endswith('_'+suffix):return m[:-len(suffix)-1].replace('jacobi','Jacobi').replace('joint','Joint')+': '+title
    return m


def energy(c,r,reference=False):
    p='ref_' if reference else '';v=0.;last=(0.,0.)
    for i in range(1,int(c['count'])+1):
        x,y=r[p+f'x{i}'],r[p+f'y{i}'];vx,vy=r[p+f'vx{i}'],r[p+f'vy{i}']
        v+=.5*c[f'm{i}']*(vx*vx+vy*vy)-c[f'm{i}']*(c['gx']*x+c['gy']*y)
        if c[f'alpha{i}']>0:v+=(math.hypot(x-last[0],y-last[1])-c[f'l{i}'])**2/(2*c[f'alpha{i}'])
        last=(x,y)
    return v


def minimum(c):
    value=0.;down=0.;g=math.hypot(c['gx'],c['gy'])
    for i in range(int(c['count']),0,-1):
        down+=c[f'm{i}'];t=down*g;value-=t*c[f'l{i}']+.5*c[f'alpha{i}']*t*t
    return value


class Data:
    def __init__(self,path):
        self.path=path;self.rows=list(read(path/'summary.csv'));self.summary={key(r):r for r in self.rows}
        assert len(self.rows)==len(self.summary)
        self.cases={r['case']:r for r in read(path/'cases.csv')}
        self.hz=sorted({int(r['hz']) for r in self.rows})
        self.methods=[m for m in METHODS if any(r['method']==m for r in self.rows)]
        assert len(self.methods)==18
        assert set(self.summary)=={(c,m,hz) for c in CASES for m in self.methods for hz in self.hz}
        assert all(math.isclose(r['h'],1/r['hz']) and r['steps']==round(r['duration']*r['hz']) for r in self.rows)
        self.audit={'summary_rows':len(self.rows),'summary_sha256':hashlib.sha256((path/'summary.csv').read_bytes()).hexdigest(),
                    'guarded_stops':sum(not complete(r) for r in self.rows),'max_scaled_energy_reconstruction_error':0.,
                    'max_scaled_reference_energy_reconstruction_error':0.,'max_position_rms_reconstruction_error':0.,
                    'max_angular_balance_reconstruction_error':0.}
        for c in self.cases.values():
            assert c['count'] in [1,2]
            assert math.isclose(energy(c,c),c['initial_energy'],rel_tol=2e-12,abs_tol=2e-12)
            assert math.isclose(c['initial_energy']-minimum(c),c['energy_scale'],rel_tol=2e-12,abs_tol=2e-12)
        self.traces=collections.defaultdict(list)
        for r in read(path/'traces.csv'):
            c=self.cases[r['case']]
            for ref,field in [(False,'energy'),(True,'reference_energy')]:
                e=energy(c,r,ref);err=abs(e-r[field])/max(1,abs(e))
                af='max_scaled_reference_energy_reconstruction_error' if ref else 'max_scaled_energy_reconstruction_error'
                self.audit[af]=max(self.audit[af],err);assert err<1e-10,(key(r),r['step'],field,err)
            ell=sum(c[f'm{i}']*(r[f'x{i}']*r[f'vy{i}']-r[f'y{i}']*r[f'vx{i}']) for i in range(1,int(c['count'])+1))
            assert math.isclose(ell,r['angular'],rel_tol=2e-12,abs_tol=2e-10)
            self.traces[key(r)].append(r)
        for k,tr in self.traces.items():
            s=self.summary[k];c=self.cases[k[0]];n=int(c['count'])
            assert len(tr)==int(s['completed'])+1,(k,len(tr),s['completed'])
            assert all(t['step']==i and math.isclose(t['time'],i/s['hz'],abs_tol=1e-11) for i,t in enumerate(tr))
            ep=max(abs(t['energy']-c['initial_energy']) for t in tr)/c['energy_scale']
            assert math.isclose(ep,s['max_energy_rel'],rel_tol=1e-10,abs_tol=1e-11)
            dist=sum((vals(tr,f'x{i}')-vals(tr,f'ref_x{i}'))**2+(vals(tr,f'y{i}')-vals(tr,f'ref_y{i}'))**2 for i in range(1,n+1))/n
            rms=float(np.sqrt(np.mean(dist[1:])))
            self.audit['max_position_rms_reconstruction_error']=max(self.audit['max_position_rms_reconstruction_error'],abs(rms-s['rms_position']))
            assert math.isclose(rms,s['rms_position'],rel_tol=2e-11,abs_tol=1e-11)
            excursion=sum((vals(tr,f'ref_x{i}')[1:]-c[f'x{i}'])**2+(vals(tr,f'ref_y{i}')[1:]-c[f'y{i}'])**2 for i in range(1,n+1))/n
            exc=float(np.sqrt(np.mean(excursion)))
            assert math.isclose(exc,s['reference_excursion'],rel_tol=2e-11,abs_tol=1e-11)
            angular=np.max(np.abs(vals(tr,'angular')-tr[0]['angular']-vals(tr,'torque_impulse')))
            self.audit['max_angular_balance_reconstruction_error']=max(self.audit['max_angular_balance_reconstruction_error'],abs(angular-s['max_angular_balance']))
            assert math.isclose(angular,s['max_angular_balance'],rel_tol=2e-10,abs_tol=1e-9)
        assert set(self.traces)=={(c,m,480) for c in CASES for m in SELECTED}
        self.audit['trace_groups']=len(self.traces)
        self.audit['reference_checks']=list(read(path/'reference_checks.csv'))
        assert len(self.audit['reference_checks'])==len(CASES)*len(self.hz)
        assert all(r['max_position_difference']<1e-7 and r['max_velocity_difference']<1e-6 and r['max_energy_rel']<1e-8 for r in self.audit['reference_checks'])
        self.audit['baseline_position_refinement']={c:[self.get(c,'xpbd',hz)[0]['rms_position'] for hz in self.hz] for c in CASES}
        assert all(all(b<a for a,b in zip(v,v[1:])) for v in self.audit['baseline_position_refinement'].values())
        self.audit['baseline_position_error_decreases_with_refinement']=True

    def get(self,c,m,hz=480):return self.summary[(c,m,hz)],self.traces.get((c,m,hz),[])


def save(fig,out,name):fig.savefig(out/name,dpi=160,facecolor='white');plt.close(fig)
def footer(fig,t):fig.text(.045,.014,t,fontsize=9,color='#445563',va='bottom')
def legend(fig,methods,y=.90,reference=True):
    hh=[Line2D([],[],color=COLORS[m],lw=2,ls='--' if m=='frozen_xpbd' else '-',label=label(m)) for m in methods]
    if reference:hh.append(Line2D([],[],color='#172b38',ls=(0,(4,3)),lw=1.3,label='Refined RK4 reference'))
    fig.legend(handles=hh,loc='upper left',bbox_to_anchor=(.04,y),ncol=3,frameon=False,fontsize=9.5,columnspacing=2)


def scale(ax,arrays,ref,positive=False):
    v=np.concatenate(arrays+[ref]);lo=float(v.min());hi=float(v.max());typical=max(1.,float(np.max(np.abs(ref))))
    if max(abs(lo),abs(hi))>8*typical:
        ax.set_yscale('symlog',linthresh=typical)
        ticks=[0.]+[s*10.**p for s in [-1,1] for p in range(-3,20) if 10.**p>=typical and lo<=s*10.**p<=hi]
        ax.set_yticks(sorted(ticks));ax.text(.015,.95,'symlog · full range',transform=ax.transAxes,va='top',fontsize=8,color='#995522')
    else:
        pad=max(.065*(hi-lo),.01);ax.set_ylim(max(0,lo-pad) if positive else lo-pad,hi+pad)


def draw(ax,data,case,field,methods=SELECTED):
    c=data.cases[case];b,bt=data.get(case,'xpbd');es=c['energy_scale'];vmin=minimum(c)
    def values(tr,ref=False):
        if field=='energy':return (vals(tr,'reference_energy' if ref else 'energy')-vmin)/es
        if field=='position_error':
            if ref:return np.zeros(len(tr))
            return np.sqrt(sum((vals(tr,f'x{i}')-vals(tr,f'ref_x{i}'))**2+(vals(tr,f'y{i}')-vals(tr,f'ref_y{i}'))**2 for i in range(1,int(c['count'])+1))/c['count'])
        if field=='angle':return np.unwrap(np.arctan2(vals(tr,('ref_' if ref else '')+'x1'),-vals(tr,('ref_' if ref else '')+'y1')))
        return vals(tr,('ref_' if ref else '')+field)
    arrays=[]
    for m in methods:
        s,tr=data.get(case,m);t=vals(tr,'time');v=values(tr);arrays.append(v)
        ax.plot(t,v,color=COLORS[m],lw=1.15,ls='--' if m=='frozen_xpbd' else '-',alpha=.93)
        if not complete(s):ax.plot(t[-1],v[-1],marker='x',color=COLORS[m],ms=8)
    ref=values(bt,True)
    if field!='position_error':ax.plot(vals(bt,'time'),ref,color='#172b38',ls=(0,(4,3)),lw=1.2)
    scale(ax,arrays,ref,positive=field in ['energy','position_error'])
    ax.set_xlim(0,b['duration']);ax.set_xlabel('Time (s)',fontsize=9);ax.grid(alpha=.22)
    title={'energy':'(E − Vmin) / (E₀ − Vmin)','position_error':'Position error across bodies (m)','angle':'Pendulum angle (rad)'}.get(field)
    if title is None:title=f'Mass {field[-1]}: {field[0]} position (m)'
    ax.set_title(title,loc='left',fontsize=11)


def primary(data,out,case):
    if case=='double_small':fields=['energy','position_error','x1','y1','x2','y2'];shape=(3,2);size=(17,16.5)
    elif case=='pendulum_rest':fields=['energy','angle','x1','y1'];shape=(2,2);size=(17,12.5)
    else:fields=['energy','x1','y1','position_error'];shape=(2,2);size=(17,12.5)
    fig,axs=plt.subplots(*shape,figsize=size)
    fig.subplots_adjust(left=.07,right=.98,top=.79,bottom=.10,hspace=.4,wspace=.22)
    fig.suptitle(NAMES[case],x=.05,y=.975,ha='left',fontsize=24,weight='bold')
    fig.text(.05,.943,'480 Hz · one pass per substep · all multipliers reset each substep · same physical energy target',fontsize=11)
    legend(fig,SELECTED,.925)
    for ax,f in zip(axs.flat,fields):draw(ax,data,case,f)
    notes=[]
    for m in ['xpbd','EL_raw_gate_cap','LE_repair_cap']:
        r,_=data.get(case,m);e,p=metrics(r);notes.append(f'{label(m)}: max E error {e:.3g}%, position RMS {p:.3g}% of reference excursion')
    footer(fig,'\n'.join(notes)+'\nOrdinary XPBD refreshes each spring row after its preceding position correction. Other material rows remain frozen for the pass. × marks a guarded stop; full ranges are retained.')
    save(fig,out,f'{case}_energy_positions_480hz.png')


def number(v):
    if v>=10000:return f'{v:.0e}'
    if v>=100:return f'{v:.0f}'
    if v>=10:return f'{v:.1f}'
    if v>=.01:return f'{v:.2f}'
    return '0' if v==0 else f'{v:.0e}'


def scorecard(data,out,hz):
    fig,axs=plt.subplots(1,2,figsize=(24,15.5))
    fig.subplots_adjust(left=.20,right=.99,top=.865,bottom=.19,wspace=.52)
    fig.suptitle(f'Cold-start methods across all cases · {hz} Hz',x=.04,y=.972,ha='left',fontsize=25,weight='bold')
    fig.text(.04,.94,'One pass per small substep · same energy target · ordinary XPBD refreshes spring rows; the other schedules freeze them',fontsize=12)
    xlabels=['Spring\nrelease','Spring\norbit','Rigid\nrest','Rigid\nkick','Double\nsmall','Double\nlarge','Double\nmixed']
    for metric,ax in enumerate(axs):
        values=np.zeros((len(data.methods),len(CASES)));stops=np.zeros_like(values,dtype=bool)
        for i,m in enumerate(data.methods):
            for j,c in enumerate(CASES):
                s,_=data.get(c,m,hz);values[i,j]=metrics(s)[metric];stops[i,j]=not complete(s)
        ax.imshow(values,cmap=ListedColormap(PALETTE),norm=BoundaryNorm(BOUNDS,len(PALETTE)),aspect='auto',interpolation='nearest')
        ax.set_yticks(range(len(data.methods)),[label(m) for m in data.methods],fontsize=9)
        ax.set_xticks(range(len(CASES)),xlabels,fontsize=10);ax.tick_params(length=0,pad=6)
        ax.set_xticks(np.arange(-.5,len(CASES),1),minor=True);ax.set_yticks(np.arange(-.5,len(data.methods),1),minor=True)
        ax.grid(which='minor',color='white',lw=1.2);ax.tick_params(which='minor',length=0)
        for i in range(len(data.methods)):
            for j in range(len(CASES)):
                if stops[i,j]:
                    ax.add_patch(Rectangle((j-.5,i-.5),1,1,fc='#d9e0e5',hatch='///',ec='#9ba6af',lw=0));text='STOP'
                else:text=number(values[i,j])
                ax.text(j,i,text,ha='center',va='center',fontsize=9,color='#152c3b')
        for spine in ax.spines.values():spine.set_visible(False)
        ax.set_title('Maximum energy error / (E₀−Vmin), %' if metric==0 else 'Position RMS / reference excursion, %',loc='left',fontsize=14,pad=15)
    fig.legend(handles=[Rectangle((0,0),1,1,fc=c,label=l) for c,l in zip(PALETTE,['<0.1%','0.1–1%','1–5%','5–20%','20–50%','50–100%','≥100%'])],
               loc='lower center',bbox_to_anchor=(.5,.10),ncol=7,frameon=False,fontsize=11)
    footer(fig,'STOP means the run did not complete; no partial-prefix error is ranked. Small energy error is not evidence of correct motion.\nCombined-correction repair limits predicted repair to 0.5 of the material-only residual and caps the extra displacement at 0.1 × motion scale. Energy gate acts before the material response.')
    save(fig,out,f'all_methods_energy_motion_{hz}hz.png')


def holdouts(data,out):
    cases=['spring_orbit','pendulum_kick','double_large','double_mixed']
    methods=['xpbd','EL_raw_gate_cap','LE_repair_cap','joint_repair_cap']
    fig=plt.figure(figsize=(20,19.5))
    fig.text(.055,.976,'Holdouts: energy and both-body trajectory fidelity',fontsize=24,weight='bold')
    fig.text(.055,.95,'480 Hz · one pass per substep · fixed settings for every case · joint and LE use the same acceptance filter',fontsize=11)
    legend(fig,methods,.932)
    grid=fig.add_gridspec(4,3,left=.065,right=.98,top=.85,bottom=.09,hspace=.66,wspace=.24)
    for i,c in enumerate(cases):
        fields=['energy','x1','x2' if data.cases[c]['count']==2 else 'y1']
        for j,f in enumerate(fields):
            ax=fig.add_subplot(grid[i,j]);draw(ax,data,c,f,methods)
            if j==0:ax.text(0,1.2,NAMES[c],transform=ax.transAxes,fontsize=13,weight='bold')
        r,_=data.get(c,'LE_repair_cap');e,p=metrics(r)
        ax=fig.axes[-3];ax.text(0,-.35,f'LE repair + cap: energy {e:.3g}%, position RMS {p:.3g}% of reference excursion',transform=ax.transAxes,fontsize=9,color='#285c54')
    footer(fig,'For the double pendulum, x₁ and x₂ show both moving masses; the full reference-based RMS in the scorecards includes x and y for each body.\nFull value ranges are shown; a symlog label marks any required range compression. Ordinary XPBD remains the motion-control comparison.')
    save(fig,out,'holdout_energy_positions_480hz.png')


def convergence(data,out):
    methods=['xpbd','frozen_xpbd','EL_raw_gate_cap','EL_repair_cap','LE_repair_cap','joint_repair_cap']
    fig,axs=plt.subplots(2,4,figsize=(22,12.5));fig.subplots_adjust(left=.055,right=.985,top=.80,bottom=.11,hspace=.40,wspace=.28)
    fig.suptitle('Position convergence under smaller cold-start steps',x=.045,y=.973,ha='left',fontsize=23,weight='bold')
    fig.text(.045,.942,'Each point is a separate run; one pass per substep throughout',fontsize=11)
    legend(fig,methods,.925,reference=False)
    for ax,c in zip(axs.flat,CASES):
        for m in methods:
            rr=[data.get(c,m,h)[0] for h in data.hz]
            yy=[metrics(r)[1] if complete(r) else float('nan') for r in rr]
            ax.plot(data.hz,yy,color=COLORS[m],lw=1.7,marker='o',ms=4,ls='--' if m=='frozen_xpbd' else '-')
        ax.set_xscale('log',base=2);ax.set_xticks(data.hz,[str(h) for h in data.hz]);ax.set_yscale('log')
        ax.set_title(NAMES[c],loc='left',fontsize=11);ax.set_xlabel('Substep frequency (Hz)');ax.set_ylabel('Normalized position RMS (%)');ax.grid(alpha=.2,which='both')
    axs.flat[-1].axis('off')
    footer(fig,'RMS is averaged over all moving bodies and all completed time steps, then divided by the reference RMS excursion from its initial configuration.\nFailed runs are omitted from connecting curves and marked STOP in the complete scorecards; baseline validation is reported separately.')
    save(fig,out,'position_convergence.png')


def overview(data,out):
    methods=['xpbd','EL_repair_cap','LE_repair_cap','EL_raw_gate_cap']
    fig,axs=plt.subplots(3,2,figsize=(16,14.5))
    fig.subplots_adjust(left=.075,right=.98,top=.82,bottom=.11,hspace=.61,wspace=.23)
    fig.suptitle('Cold-start XPBD: energy correction versus motion accuracy',x=.05,y=.974,ha='left',fontsize=23,weight='bold')
    fig.text(.05,.944,'480 Hz · one pass per small substep · multipliers reset every substep · 10 seconds',fontsize=11)
    legend(fig,methods,.925,reference=False)
    for i,c in enumerate(['spring_release','pendulum_rest','double_small']):
        meta=data.cases[c];es=meta['energy_scale'];e0=meta['initial_energy'];n=int(meta['count'])
        largest=0.
        for m in methods:
            _,tr=data.get(c,m);t=vals(tr,'time');energy_error=100*np.abs(vals(tr,'energy')-e0)/es
            position_error=np.sqrt(sum((vals(tr,f'x{j}')-vals(tr,f'ref_x{j}'))**2+(vals(tr,f'y{j}')-vals(tr,f'ref_y{j}'))**2 for j in range(1,n+1))/n)
            axs[i,0].plot(t,energy_error,color=COLORS[m],lw=1.45)
            axs[i,1].plot(t,position_error,color=COLORS[m],lw=1.45)
            largest=max(largest,float(energy_error.max()))
        axs[i,0].set_yscale('symlog',linthresh=.001)
        axs[i,0].set_yticks([0]+[10.**p for p in range(-3,4) if 10.**p<=largest*1.3])
        axs[i,0].set_ylim(0,largest*1.12)
        axs[i,0].set_title(NAMES[c],loc='left',fontsize=14,weight='bold')
        axs[i,0].set_ylabel('100 |E − E₀| / (E₀ − Vmin) · symlog',fontsize=10)
        axs[i,1].set_title('Position error versus refined RK4 (m)',loc='left',fontsize=12)
        axs[i,1].set_ylim(bottom=0)
        for ax in axs[i]:ax.set_xlim(0,10);ax.set_xlabel('Time (s)');ax.grid(alpha=.2)
        b,_=data.get(c,'xpbd');r,_=data.get(c,'LE_repair_cap')
        note=f'XPBD → LE repair + cap: max energy error {100*b["max_energy_rel"]:.3g}% → {100*r["max_energy_rel"]:.3g}%; position RMS {b["rms_position"]:.4f} m → {r["rms_position"]:.4f} m'
        axs[i,0].text(0,-.35,note,transform=axs[i,0].transAxes,fontsize=9,color='#285c54')
    spike=100*data.get('spring_release','EL_raw_gate')[0]['max_energy_rel']
    footer(fig,f'Position error uses both Cartesian coordinates and, for the double pendulum, the RMS across both masses. Lower is better on both columns.\nThe uncapped EL energy-direction gate is excluded from this overview: its released-spring energy spike reaches {spike:.3g}%; full ranges appear in the individual figure.')
    save(fig,out,'overview.png')


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('results',type=Path)
    p.add_argument('--output',type=Path,default=Path('results/traditional_cases/plots'))
    args=p.parse_args();args.output.mkdir(parents=True,exist_ok=True)
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    d=Data(args.results)
    for c in ['spring_release','pendulum_rest','double_small']:primary(d,args.output,c)
    for hz in d.hz:scorecard(d,args.output,hz)
    holdouts(d,args.output);convergence(d,args.output);overview(d,args.output)
    (args.output/'plot_audit.json').write_text(json.dumps(d.audit,indent=2)+'\n')
    print(json.dumps({k:v for k,v in d.audit.items() if k!='reference_checks'},indent=2))
    print('Figures: '+', '.join(f.name for f in sorted(args.output.glob('*.png'))))


if __name__=='__main__':main()
