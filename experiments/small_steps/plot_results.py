#!/usr/bin/env python3
"""Audit and plot EL small-step experiments; one sweep per substep throughout."""
import argparse
import collections
import csv
import hashlib
import json
import math
import os
from pathlib import Path

os.environ.setdefault('MPLCONFIGDIR','/tmp/xpbd-small-steps-matplotlib')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
import numpy as np

CASES=['vertical_oscillation','small_swing','swing_rest','swing_tangent','coupled_motion']
NAMES=['Vertical oscillation','Small swing','Large swing from rest','Swing with tangential velocity','Coupled radial and swinging motion']
FAMILIES=['baseline','unguarded','directional','relax','cap','soft','effective','effective_cap']
FAMILY_COLORS={'baseline':'#939da7','unguarded':'#cf7b22','directional':'#156cb0','relax':'#b13b59',
               'cap':'#28856b','soft':'#9053a1','effective':'#43a2c2','effective_cap':'#a3a136'}
TEXT={'case','method','family','gate','status'}
VMIN=-10.2911805
PALETTE=['#f2f7f0','#d7e7cd','#f5edbf','#f1d1a0','#eaa998','#d87a84']
MAIN=['baseline','cap_0.001000','soft_10.000000','soft_1000.000000']
DIFFICULT=['baseline','unguarded','directional','relax_0.100000']
CONVERGENCE=['baseline','unguarded','directional','relax_0.100000','cap_0.001000','cap_0.100000','soft_10.000000','soft_1000.000000']
METHOD_COLORS={'baseline':'#939da7','unguarded':'#cf7b22','directional':'#156cb0','relax_0.100000':'#b13b59',
               'cap_0.001000':'#28856b','cap_0.100000':'#7e9d32','soft_10.000000':'#9053a1','soft_1000.000000':'#503593'}


def read(path):
    with path.open(newline='') as f: rows=list(csv.DictReader(f))
    for r in rows:
        for k in r:
            if k not in TEXT:
                r[k]=float(r[k])
                assert math.isfinite(r[k]),(path,k)
    return rows


def key(r): return r['case'],r['method'],int(r['substeps'])
def vals(tr,k): return np.array([r[k] for r in tr])
def complete(r): return r['completed_steps']==r['requested_steps']
def metrics(r): return 100*r['max_abs_energy_error']/r['Escale'],100*r['rms_position_error']/r['reference_rms_excursion']
def label(r):
    f=r['family']
    if f=='baseline': return 'Ordinary XPBD'
    if f=='unguarded': return 'Unguarded EL'
    if f=='directional': return 'Directional EL'
    if f=='effective': return 'Effective-direction EL'
    if f=='relax': return f"Directional + relaxation ω={r['omega']:g}"
    if f=='cap': return f"Directional + cap γ={r['cap']:g}"
    if f=='effective_cap': return f"Effective direction + cap γ={r['cap']:g}"
    return f"Directional + compliance β/Escale={r['beta_ratio']:g}"


class Data:
    def __init__(self,sweep,selected):
        self.rows=read(sweep/'summary.csv');self.summary={key(r):r for r in self.rows}
        assert len(self.rows)==len(self.summary)
        self.methods=list(dict.fromkeys(r['method'] for r in self.rows))
        expected={(c,m,n) for c in CASES for m in self.methods for n in [4,8,16]}
        assert set(self.summary)==expected
        assert all(r['sweeps']==1 and math.isclose(r['h'],1/(60*r['substeps'])) for r in self.rows)
        assert all(r['requested_steps']==600*r['substeps'] for r in self.rows)
        assert all(r['scalar_solves']==r['requested_steps']*(1 if r['family']=='baseline' else 2) for r in self.rows)
        assert all(r['evaluations']==r['requested_steps'] for r in self.rows)
        self.meta={m:next(r for r in self.rows if r['method']==m) for m in self.methods}
        self.traces=collections.defaultdict(list)
        self.audit={'source_runs':len(self.rows),'unique_source_keys':len(self.summary),
                    'source_sha256':hashlib.sha256((sweep/'summary.csv').read_bytes()).hexdigest(),
                    'max_scaled_energy_reconstruction_error':0.,'max_position_rms_reconstruction_error':0.,
                    'all_sweeps_equal_one':True,'guarded_stops':sum(not complete(r) for r in self.rows)}
        selected_rows=read(selected/'summary.csv')
        for r in selected_rows: assert r==self.summary[key(r)],('Selected rerun differs',key(r))
        for r in read(selected/'traces.csv'):
            assert r['sweeps']==1 and r['substeps']==8
            E=.5*(r['vx']**2+r['vy']**2)+50*(math.hypot(r['x'],r['y'])-1)**2+9.81*r['y']
            error=abs(E-r['energy'])/max(1,abs(E))
            self.audit['max_scaled_energy_reconstruction_error']=max(self.audit['max_scaled_energy_reconstruction_error'],error)
            assert error<5e-12
            assert math.isclose(r['angular_momentum'],r['x']*r['vy']-r['y']*r['vx'],abs_tol=1e-12,rel_tol=1e-12)
            self.traces[key(r)].append(r)
        for k,tr in self.traces.items():
            r=self.summary[k];assert len(tr)==int(r['completed_steps'])+1
            assert all(t['step']==i and math.isclose(t['time'],i*t['h'],abs_tol=1e-12) for i,t in enumerate(tr))
            ep=max(abs(t['energy']-t['target_energy']) for t in tr)
            assert math.isclose(ep,r['max_abs_energy_error'],rel_tol=1e-11,abs_tol=1e-11)
            dist=np.hypot(vals(tr,'x')-vals(tr,'ref_x'),vals(tr,'y')-vals(tr,'ref_y'))
            rms=float(np.sqrt(np.mean(dist[1:]**2)))
            self.audit['max_position_rms_reconstruction_error']=max(self.audit['max_position_rms_reconstruction_error'],abs(rms-r['rms_position_error']))
            assert math.isclose(rms,r['rms_position_error'],rel_tol=1e-12,abs_tol=1e-12)
            excursion=float(np.sqrt(np.mean((vals(tr,'ref_x')[1:]-tr[0]['x'])**2+(vals(tr,'ref_y')[1:]-tr[0]['y'])**2)))
            assert math.isclose(excursion,r['reference_rms_excursion'],rel_tol=1e-12,abs_tol=1e-12)
            for field in ['energy_attempts','energy_accepted','directional_skips','capped','singular_updates']:
                assert tr[-1][field]==r[field]
                assert all(b[field]>=a[field] for a,b in zip(tr,tr[1:]))
        # Prove the plotted-away duplicate, including every recorded state.
        state_fields=['x','y','vx','vy','energy']
        for case in CASES:
            for n in [4,8,16]:
                a=self.summary[(case,'directional',n)];b=self.summary[(case,'effective',n)]
                for f in ['max_abs_energy_error','rms_position_error','path_ratio','directional_skips']:
                    assert a[f]==b[f]
            a=self.traces[(case,'directional',8)];b=self.traces[(case,'effective',8)]
            assert all(all(x[f]==y[f] for f in state_fields) for x,y in zip(a,b)) and len(a)==len(b)
        self.audit['selected_summary_rows']=len(selected_rows)
        self.audit['audited_trace_groups']=len(self.traces)
        self.audit['raw_effective_state_identical_at_480hz']=True

    def get(self,case,m,n=8): return self.summary[(case,m,n)],self.traces.get((case,m,n),[])


def selection(data):
    configs=[]
    for m in data.methods:
        rows=[r for r in data.rows if r['method']==m]
        ep=max(metrics(r)[0] for r in rows);pos=max(metrics(r)[1] for r in rows)
        configs.append({'method':m,'family':rows[0]['family'],'complete_all':all(complete(r) for r in rows),
                        'worst_energy_error_percent':ep,'worst_normalized_position_rms_percent':pos,
                        'qualifies_five_percent_energy':all(complete(r) for r in rows) and ep<=5})
    chosen={}
    for fam in FAMILIES:
        candidates=[r for r in configs if r['family']==fam]
        qualified=[r for r in candidates if r['qualifies_five_percent_energy']]
        pool=qualified or [r for r in candidates if r['complete_all']]
        best=min(pool,key=lambda r:(r['worst_normalized_position_rms_percent'],r['worst_energy_error_percent'],r['method']))
        chosen[fam]=dict(best,selection_pool='energy-qualified' if qualified else 'fallback: no energy-qualified candidate')
    return {'rule':'One fixed parameter per family across all 5 cases × 3 substep counts. Require completion and max energy error/Escale <= 5% in every condition; among eligible candidates minimize the worst position RMS/reference RMS excursion. If none is eligible, minimize that same worst motion error among completed candidates. Ties use worst energy error then method id.',
            'energy_threshold_percent':5,'conditions_per_configuration':15,
            'position_normalization':'reference_rms_excursion = RMS norm(x_ref(t)-x_initial), excluding t=0',
            'qualifying_configurations':[r['method'] for r in configs if r['qualifies_five_percent_energy']],
            'representatives':chosen,'all_configurations':configs,
            'additional_illustrative_methods':['cap_0.100000','soft_10.000000'],
            'additional_method_reason':'Show the tradeoff between preserving energy and following the reference; these are not replacements for the fixed fallback representatives.',
            'effective_duplicate':'Effective and raw directional gates have identical states in all selected 480 Hz traces; convergence summaries also match. Only raw directional curves are drawn.'}


def method_legend(fig,data,methods,y,ncol=3,size=10):
    handles=[Line2D([],[],color=METHOD_COLORS[m],lw=2,label=label(data.meta[m])) for m in methods]
    handles.append(Line2D([],[],color='#182b3b',ls=(0,(4,3)),lw=1.3,label='Refined RK4 reference'))
    fig.legend(handles=handles,loc='upper left',bbox_to_anchor=(.047,y),ncol=ncol,frameon=False,fontsize=size,columnspacing=2)


def set_scale(ax,arrays,reference,energy=False):
    allv=np.concatenate(arrays+[reference]);lo=float(np.min(allv));hi=float(np.max(allv))
    typical=max(1.,float(np.max(np.abs(reference))))
    if max(abs(lo),abs(hi))>8*typical:
        ax.set_yscale('symlog',linthresh=typical)
        ticks=[0.]+[s*10.**p for s in [-1,1] for p in range(0,16) if 10.**p>=typical and lo<=s*10.**p<=hi]
        ax.set_yticks(sorted(ticks))
        ax.text(.02,.95,'symlog · full range',transform=ax.transAxes,va='top',fontsize=9,color='#8b4d20')
    else:
        pad=max(.07*(hi-lo),.02 if not energy else .01)
        ax.set_ylim((max(0,lo-pad) if energy else lo-pad),hi+pad)
    return lo,hi


def trajectory_figure(data,out,methods,difficult=False):
    fig=plt.figure(figsize=(20,23),facecolor='white')
    title='Raw, relaxed, and directional EL: full failure ranges' if difficult else 'Small steps: energy and motion tradeoffs'
    fig.text(.055,.971,title,fontsize=24,weight='bold')
    fig.text(.055,.947,'8 substeps per 60 Hz frame = 480 Hz · exactly one linearized EL sweep per substep · 10 seconds · fixed parameters across all five cases',fontsize=11)
    method_legend(fig,data,methods,.932,ncol=3)
    grid=fig.add_gridspec(5,3,left=.06,right=.98,top=.848,bottom=.064,hspace=.85,wspace=.25)
    for ri,case in enumerate(CASES):
        axs=[fig.add_subplot(grid[ri,c]) for c in range(3)]
        _,base=data.get(case,'baseline');time=vals(base,'time');es=base[0]['Escale']
        reference=[(vals(base,'ref_energy')-VMIN)/es,vals(base,'ref_x'),vals(base,'ref_y')]
        drawn=[[] for _ in range(3)]
        for m in methods:
            s,tr=data.get(case,m);t=vals(tr,'time')
            vv=[(vals(tr,'energy')-VMIN)/es,vals(tr,'x'),vals(tr,'y')]
            for c,ax in enumerate(axs):
                ax.plot(t,vv[c],color=METHOD_COLORS[m],lw=1.15,alpha=.94)
                drawn[c].append(vv[c])
                if not complete(s):ax.plot(t[-1],vv[c][-1],marker='X',color=METHOD_COLORS[m],ms=7)
        for c,ax in enumerate(axs):
            ax.plot(time,reference[c],color='#182b3b',ls=(0,(4,3)),lw=1.1,zorder=6)
            set_scale(ax,drawn[c],reference[c],c==0)
            ax.set(xlim=(0,10),xlabel='Time (s)',ylabel=['(E − Vmin) / (E₀ − Vmin)','x (m)','y (m)'][c],title=['Energy above minimum / initial','Horizontal position','Vertical position'][c])
            ax.grid(alpha=.17)
        p=axs[0].get_position();fig.text(.06,p.y1+.034,NAMES[ri],fontsize=14,weight='bold')
        notes=[]
        for m in methods:
            s,_=data.get(case,m);e,_=metrics(s)
            short={'baseline':'XPBD','cap_0.001000':'cap .001','soft_10.000000':'soft 10','soft_1000.000000':'soft 1000','unguarded':'unguarded','directional':'directional','relax_0.100000':'relax .1'}[m]
            notes.append(f"{short}: {e:.3g}% / {s['rms_position_error']:.3g} m")
        fig.text(.06,p.y0-.032,'Max energy error / position RMS:  '+'      |      '.join(notes),fontsize=9.4,color='#425667')
    fig.text(.06,.017,'All plotted curves use one sweep. Energy scale E₀−Vmin, Vmin = −10.2911805 J. Effective-direction EL is identical to directional EL here and is drawn once.',fontsize=9.8)
    name='difficult_variants_480hz.png' if difficult else 'small_steps_energy_and_motion_480hz.png'
    fig.savefig(out/name,dpi=155,facecolor='white');plt.close(fig)
    return name


def scorecard(data,out,n):
    count=len(data.methods);fig=plt.figure(figsize=(21,14.5),facecolor='white')
    fig.text(.026,.968,f'Every parameter · {n} substeps/frame · {60*n} Hz · one sweep/substep',fontsize=22,weight='bold')
    fig.text(.026,.936,f'{count} configurations × 5 cases, each simulated for 10 s. Completed runs remain visible even when their errors are large.',fontsize=11)
    left=.334;bottom=.18;height=.66
    labels=fig.add_axes([.026,bottom,.30,height]);ax=fig.add_axes([left,bottom,.592,height]);cost=fig.add_axes([.938,bottom,.046,height])
    for a in [labels,ax,cost]:a.set_ylim(count-.5,-.5);a.set_xticks([]);a.set_yticks([]);a.set_frame_on(False)
    labels.set_xlim(0,1);ax.set_xlim(-.5,9.5);cost.set_xlim(0,1)
    labels.text(0,-1.1,'Configuration / fixed parameter',fontsize=11,weight='bold',clip_on=False)
    cost.text(.5,-1.1,'Solves /\nframe',ha='center',fontsize=10,weight='bold',clip_on=False)
    short=['Vertical','Small swing','Large swing','Tangential','Coupled']
    for c,name in enumerate(short):
        ax.text(2*c+.5,-2.8,name,ha='center',fontsize=11,weight='bold',clip_on=False)
        ax.text(2*c,-1.1,'E error\n(%)',ha='center',fontsize=10,weight='bold',clip_on=False)
        ax.text(2*c+1,-1.1,'Position RMS\n(% excursion)',ha='center',fontsize=9.5,weight='bold',clip_on=False)
    displayed=[]
    for i,m in enumerate(data.methods):
        if i%2==0:labels.add_patch(Rectangle((0,i-.5),1,1,fc='#f2f5f7',ec='none'))
        labels.text(.01,i,label(data.meta[m]),fontsize=10,va='center')
        costs=[]
        for c,case in enumerate(CASES):
            r,_=data.get(case,m,n);displayed.append(key(r));costs.append(r['scalar_solves']/600)
            if not complete(r):
                ax.add_patch(Rectangle((2*c-.5,i-.5),2,1,fc='#5e6670',ec='white',lw=.4));ax.text(2*c+.5,i,'STOP',ha='center',va='center',color='white',weight='bold');continue
            for j,(value,bounds) in enumerate(zip(metrics(r),[[.1,1,5,20,100],[1,5,10,25,100]])):
                color=PALETTE[sum(value>b for b in bounds)]
                ax.add_patch(Rectangle((2*c+j-.5,i-.5),1,1,fc=color,ec='white',lw=.4))
                ax.text(2*c+j,i,f'{value:.3g}',ha='center',va='center',fontsize=9.2)
        cost.text(.5,i,f'{min(costs):g}' if min(costs)==max(costs) else f'{min(costs):g}–{max(costs):g}',ha='center',va='center',fontsize=10)
        if i==count-1 or data.meta[data.methods[i+1]]['family']!=data.meta[m]['family']:
            for a in [labels,ax,cost]:a.axhline(i+.5,color='#7b8791',lw=.7)
    for c in range(1,5):ax.axvline(2*c-.5,color='#7b8791',lw=.7)
    for y,title,bins in [(.14,'Energy bins (%)',['≤0.1','≤1','≤5','≤20','≤100','>100']),(.118,'Position bins (% excursion)',['≤1','≤5','≤10','≤25','≤100','>100'])]:
        fig.text(.026,y,title,fontsize=10,weight='bold')
        for i,(color,b) in enumerate(zip(PALETTE,bins)):
            x=.19+i*.066;fig.add_artist(Rectangle((x,y-.004),.062,.018,transform=fig.transFigure,fc=color,ec='white'))
            fig.text(x+.031,y+.004,b,fontsize=9.7,ha='center',va='center')
    fig.text(.026,.084,'Energy error = 100 max|E−E₀|/(E₀−Vmin). Position = 100 RMS|x−x_ref| / RMS|x_ref−x_initial|. Lower is better within each metric.',fontsize=10)
    fig.text(.026,.061,'Solves/frame counts scalar-solve attempts per 60 Hz frame; it is not a runtime benchmark. The 5% energy criterion is a selection aid, not a correctness guarantee.',fontsize=10)
    fig.text(.026,.038,'Every setting uses one EL sweep per substep. Raw/effective gates coincide in this study. No candidate keeps max energy error ≤5% in all 15 conditions.',fontsize=10)
    name=f'all_parameters_substeps{n}.png';fig.savefig(out/name,dpi=170,facecolor='white');plt.close(fig)
    return name,displayed


def pareto(data,out,selected):
    fig,axs=plt.subplots(1,2,figsize=(19,8.5));fig.subplots_adjust(left=.075,right=.98,bottom=.19,top=.80,wspace=.23)
    fig.suptitle('Parameter tradeoffs: energy preservation and motion fidelity',fontsize=22,weight='bold',y=.96)
    fig.text(.075,.905,'Every tested parameter is included. One fixed configuration is evaluated across cases; no per-case tuning. Smaller values are better.',fontsize=11)
    representatives={r['method'] for r in selected['representatives'].values()}
    for col,substeps in enumerate([None,8]):
        ax=axs[col]
        for m in data.methods:
            rows=[r for r in data.rows if r['method']==m and (substeps is None or r['substeps']==substeps)]
            if not all(complete(r) for r in rows):continue
            e=max(metrics(r)[0] for r in rows);p=max(metrics(r)[1] for r in rows);family=data.meta[m]['family']
            ax.scatter(e,p,s=65 if m in representatives else 30,color=FAMILY_COLORS[family],marker='s' if m in representatives else 'o',edgecolors='#263746' if m in representatives else 'none',zorder=4,alpha=.85)
            if m in ['baseline','unguarded','directional','relax_0.100000','cap_0.001000','soft_10.000000','soft_1000.000000']:
                short={'baseline':'XPBD','unguarded':'unguarded','directional':'directional / effective','relax_0.100000':'relax .1','cap_0.001000':'cap .001','soft_10.000000':'soft 10','soft_1000.000000':'soft 1000'}[m]
                offsets={'baseline':(7,8),'cap_0.001000':(12,-14),'soft_1000.000000':(7,-15),'soft_10.000000':(-18,14),'unguarded':(7,7),'directional':(7,7),'relax_0.100000':(7,7)}
                ax.annotate(short,(e,p),xytext=offsets[m],textcoords='offset points',fontsize=9)
        ax.axvline(5,color='#8b4d20',ls=':',lw=1.3,label='5% energy reference')
        ax.set_xscale('log');ax.set_yscale('log');ax.grid(alpha=.17,which='both')
        ax.set(xlabel='Worst max energy error (%)',ylabel='Worst position RMS / reference excursion (%)',title='All 15 conditions: five cases × 4/8/16 substeps' if substeps is None else 'Five cases at 8 substeps (480 Hz)')
    handles=[Line2D([],[],marker='o',ls='',color=FAMILY_COLORS[f],label=f.replace('_',' ')) for f in FAMILIES]
    fig.legend(handles=handles,loc='lower center',bbox_to_anchor=(.5,.08),ncol=4,frameon=False,fontsize=10)
    fig.text(.075,.025,'Squares: family representatives from the disclosed rule. The fallback minimizes worst normalized motion error when no candidate meets the energy criterion; it does not establish success.',fontsize=10)
    name='parameter_tradeoffs.png';fig.savefig(out/name,dpi=170,facecolor='white');plt.close(fig);return name


def convergence(data,out):
    fig=plt.figure(figsize=(16,22),facecolor='white')
    fig.text(.07,.971,'Convergence with one sweep per substep',fontsize=24,weight='bold')
    fig.text(.07,.946,'4 / 8 / 16 substeps per 60 Hz frame = 240 / 480 / 960 Hz. Same parameter values at every resolution and in every case.',fontsize=10.5)
    handles=[Line2D([],[],color=METHOD_COLORS[m],marker='o',lw=1.5,label=label(data.meta[m])) for m in CONVERGENCE]
    fig.legend(handles=handles,loc='upper left',bbox_to_anchor=(.06,.932),ncol=2,frameon=False,fontsize=9.5)
    grid=fig.add_gridspec(5,2,left=.09,right=.97,top=.815,bottom=.065,hspace=.65,wspace=.28)
    for ri,case in enumerate(CASES):
        axs=[fig.add_subplot(grid[ri,c]) for c in range(2)]
        for m in CONVERGENCE:
            rows=[data.get(case,m,n)[0] for n in [4,8,16]]
            for c,ax in enumerate(axs):
                ys=[metrics(r)[c] if complete(r) else np.nan for r in rows]
                ax.plot([4,8,16],ys,color=METHOD_COLORS[m],marker='o',lw=1.15,ms=4)
        for c,ax in enumerate(axs):
            ax.set_xscale('log',base=2);ax.set_yscale('log');ax.set_xticks([4,8,16],['4','8','16']);ax.grid(alpha=.18,which='both')
            ax.set(xlabel='Substeps per 60 Hz frame (one sweep each)',ylabel=['Max energy error (%)','Position RMS / ref. excursion (%)'][c])
            if c==0:ax.axhline(5,color='#8b4d20',ls=':',lw=.9)
        p=axs[0].get_position();fig.text(.09,p.y1+.022,NAMES[ri],fontsize=13,weight='bold')
    fig.text(.09,.028,'Counted cost per frame: XPBD 4 / 8 / 16 scalar solves; EL variants 8 / 16 / 32. All use 4 / 8 / 16 row evaluations. These are operation counts, not measured runtime.',fontsize=9)
    name='convergence_and_cost.png';fig.savefig(out/name,dpi=160,facecolor='white');plt.close(fig);return name


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('sweep',type=Path);p.add_argument('selected',type=Path);p.add_argument('--output',type=Path)
    args=p.parse_args();out=args.output or args.selected/'plots';out.mkdir(parents=True,exist_ok=True)
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False,'axes.edgecolor':'#aab5bd','text.color':'#193044','axes.labelcolor':'#304757','axes.titlecolor':'#193044'})
    data=Data(args.sweep,args.selected);chosen=selection(data)
    (out/'selection.json').write_text(json.dumps(chosen,indent=2)+'\n')
    images=[];shown=[]
    for n in [4,8,16]:name,keys=scorecard(data,out,n);images.append(name);shown.extend(keys)
    assert len(shown)==len(set(shown))==len(data.rows) and set(shown)==set(data.summary)
    images.append(pareto(data,out,chosen))
    images.append(trajectory_figure(data,out,MAIN))
    images.append(trajectory_figure(data,out,DIFFICULT,True))
    images.append(convergence(data,out))
    data.audit.update({'images':images,'scorecard_run_entries':len(shown),'no_missing_or_repeated_scorecard_keys':True,'configurations':len(data.methods),'qualifying_configurations':chosen['qualifying_configurations']})
    (out/'plot_audit.json').write_text(json.dumps(data.audit,indent=2)+'\n')
    print(json.dumps(data.audit,indent=2))


if __name__=='__main__':main()
