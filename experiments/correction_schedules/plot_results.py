#!/usr/bin/env python3
"""Audit and plot frozen XPBD/energy schedules crossed with acceptance rules."""
import argparse
import collections
import csv
import hashlib
import json
import math
import os
from pathlib import Path

os.environ.setdefault('MPLCONFIGDIR', '/tmp/xpbd-correction-matplotlib')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm, ListedColormap
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
import numpy as np

TEXT={'case','method','schedule','policy','status'}
CASES=['vertical_oscillation','small_swing','swing_rest','swing_tangent','coupled_motion',
       'spring_release','spring_orbit','pendulum_rest','pendulum_kick']
NAMES=['Vertical spring + gravity','Small elastic swing','Elastic swing from rest',
       'Elastic swing with velocity','Coupled elastic motion','Spring release',
       'Spring orbit','Rigid pendulum from rest','Rigid pendulum with velocity']
SCHEDULES=['jacobi','EL','LE','joint','guarded_joint','alternating','LEL','ELE',
           'spectral','spectral_p1','spectral_unbounded']
SCHED_LABELS=['Jacobi','E → L','L → E','Joint','Guarded joint','Alternate EL / LE',
              'L → E → L','E → L → E','Spectral p=½, gain≤4',
              'Spectral p=1, gain≤4','Spectral p=½, uncapped']
POLICIES=['raw','direction','cap','direction_cap','repair','repair_cap']
POLICY_LABELS=['Raw','Direction','Cap','Direction\n+ cap','Fractional\nrepair','Repair\n+ cap']
COLORS={'baseline':'#8c98a2','EL_raw':'#ce7326','LE_raw':'#a13f67',
        'LE_repair_cap':'#168270','spectral_repair_cap':'#7656a7',
        'joint_repair_cap':'#1c709f','joint_raw':'#c19a20',
        'alternating_repair_cap':'#c45676','LE_direction_cap':'#65a23e',
        'EL_repair_cap':'#649ec9'}
LABELS={'baseline':'Ordinary XPBD','EL_raw':'EL, raw','LE_raw':'LE, raw',
        'LE_repair_cap':'LE, fractional repair + cap',
        'spectral_repair_cap':'Spectral ½, fractional repair + cap',
        'joint_raw':'Joint, raw','joint_repair_cap':'Joint, fractional repair + cap',
        'alternating_repair_cap':'Alternating, fractional repair + cap',
        'LE_direction_cap':'LE, direction + cap','EL_repair_cap':'EL, fractional repair + cap'}
MAIN=['baseline','EL_raw','LE_raw','LE_repair_cap','spectral_repair_cap']
BOUNDARIES=[0,.1,1,5,20,50,100,1e300]
PALETTE=['#eaf3e8','#cde3c7','#f0efbf','#f3d9a3','#eeb692','#dd8b88','#bd5e79']


def read(path):
    with path.open(newline='') as f:
        for row in csv.DictReader(f):
            for k in row:
                if k not in TEXT:
                    row[k]=float(row[k]);assert math.isfinite(row[k]),(path,k,row)
            yield row


def key(r): return r['case'],r['method'],int(r['substeps']),r['duration']
def complete(r): return r['completed_steps']==r['requested_steps']
def vals(tr,f): return np.asarray([r[f] for r in tr])
def metrics(r):
    return 100*r['max_abs_energy_error']/r['Escale'],100*r['rms_position_error']/r['reference_rms_excursion']
def minimum(r):
    return -r['mass']*r['gravity']*r['length']-(0 if r['pendulum'] else .5*r['mass']**2*r['gravity']**2*r['alpha'])
def physical_energy(r):
    spring=0 if r['pendulum'] else (math.hypot(r['x'],r['y'])-r['length'])**2/(2*r['alpha'])
    return .5*r['mass']*(r['vx']**2+r['vy']**2)+spring+r['mass']*r['gravity']*r['y']


class Data:
    def __init__(self,sweep,selected,duration):
        self.rows=list(read(sweep/'summary.csv'))
        self.summary={key(r):r for r in self.rows};assert len(self.summary)==len(self.rows)
        self.duration=duration
        self.methods=list(dict.fromkeys(r['method'] for r in self.rows if r['duration']==duration))
        self.cases=[c for c in CASES if any(r['case']==c for r in self.rows)]
        self.substeps=sorted({int(r['substeps']) for r in self.rows})
        self.meta={m:next(r for r in self.rows if r['method']==m) for m in self.methods}
        expected={(c,m,n,duration) for c in self.cases for m in self.methods for n in self.substeps}
        assert {k for k in self.summary if k[-1]==duration}==expected
        assert all(r['linearizations']==1 and math.isclose(r['h'],1/(60*r['substeps'])) for r in self.rows)
        assert all(r['requested_steps']==round(r['duration']/r['h']) for r in self.rows)
        selected_rows=list(read(selected/'summary.csv'))
        for r in selected_rows: assert r==self.summary[key(r)],('Selected rerun differs',key(r))
        self.traces=collections.defaultdict(list)
        self.audit={'source_runs':len(self.rows),'unique_keys':len(self.summary),
                    'source_sha256':hashlib.sha256((sweep/'summary.csv').read_bytes()).hexdigest(),
                    'guarded_stops':sum(not complete(r) for r in self.rows),
                    'max_scaled_energy_reconstruction_error':0.,'max_rms_reconstruction_error':0.}
        for r in read(selected/'traces.csv'):
            e=physical_energy(r);err=abs(e-r['energy'])/max(1,abs(e))
            self.audit['max_scaled_energy_reconstruction_error']=max(self.audit['max_scaled_energy_reconstruction_error'],err)
            assert err<5e-12
            ell=r['mass']*(r['x']*r['vy']-r['y']*r['vx'])
            assert math.isclose(ell,r['angular_momentum'],abs_tol=2e-11,rel_tol=2e-12)
            assert math.isclose(r['target_energy']-minimum(r),r['Escale'],rel_tol=1e-12,abs_tol=1e-12)
            self.traces[key(r)].append(r)
        for k,tr in self.traces.items():
            s=self.summary[k];assert len(tr)==int(s['completed_steps'])+1
            assert all(r['step']==i and math.isclose(r['time'],i*r['h'],abs_tol=1e-12) for i,r in enumerate(tr))
            ep=max(abs(r['energy']-r['target_energy']) for r in tr)
            assert math.isclose(ep,s['max_abs_energy_error'],rel_tol=2e-11,abs_tol=1e-11)
            rms=float(np.sqrt(np.mean((vals(tr,'x')[1:]-vals(tr,'ref_x')[1:])**2+(vals(tr,'y')[1:]-vals(tr,'ref_y')[1:])**2)))
            self.audit['max_rms_reconstruction_error']=max(self.audit['max_rms_reconstruction_error'],abs(rms-s['rms_position_error']))
            assert math.isclose(rms,s['rms_position_error'],rel_tol=2e-12,abs_tol=1e-12)
            if complete(s):
                for f in ['energy_attempts','accepted','direction_skips','repair_skips','repair_scaled','capped','singular_updates']:
                    assert tr[-1][f]==s[f],(k,f)
                    assert all(b[f]>=a[f] for a,b in zip(tr,tr[1:]))
        self.audit['selected_runs']=len(selected_rows)
        self.audit['audited_trace_groups']=len(self.traces)
        self.audit['all_linearizations_equal_one']=True
        refs=list(read(sweep/'reference_checks.csv'))
        self.audit['reference_max_position_difference']=max(r['max_position_difference'] for r in refs)
        self.audit['reference_max_energy_error']=max(r['max_energy_error'] for r in refs)

    def get(self,c,m,n=8):
        k=(c,m,n,self.duration);return self.summary[k],self.traces.get(k,[])


def stamp(fig,text): fig.text(.04,.013,text,fontsize=9,color='#465561',va='bottom')
def save(fig,out,name):
    fig.savefig(out/name,dpi=160,facecolor='white');plt.close(fig)


def number(x):
    if x>=10000:return f'{x:.0e}'
    if x>=100:return f'{x:.0f}'
    if x>=10:return f'{x:.1f}'
    if x>=.01:return f'{x:.2f}'
    if x==0:return '0'
    return f'{x:.0e}'


def heatmap(ax,values,stops,ylabels=True):
    cmap=ListedColormap(PALETTE);norm=BoundaryNorm(BOUNDARIES,len(PALETTE))
    ax.imshow(values,cmap=cmap,norm=norm,aspect='auto',interpolation='nearest')
    ax.set_xticks(range(len(POLICIES)),POLICY_LABELS,fontsize=9)
    ax.set_yticks(range(len(SCHEDULES)),SCHED_LABELS if ylabels else [],fontsize=9)
    ax.tick_params(length=0,pad=5)
    ax.set_xticks(np.arange(-.5,len(POLICIES),1),minor=True)
    ax.set_yticks(np.arange(-.5,len(SCHEDULES),1),minor=True)
    ax.grid(which='minor',color='white',lw=1.3);ax.tick_params(which='minor',length=0)
    for i in range(values.shape[0]):
        for j in range(values.shape[1]):
            if stops[i,j]:
                ax.add_patch(Rectangle((j-.5,i-.5),1,1,facecolor='#dae0e5',hatch='///',edgecolor='#a4acb4',lw=0))
                t='STOP'
            else:t=number(values[i,j])
            ax.text(j,i,t,ha='center',va='center',fontsize=8.7,color='#142a39')
    for s in ax.spines.values():s.set_visible(False)


def legend_bins(fig,y=.05):
    labels=['<0.1%','0.1–1%','1–5%','5–20%','20–50%','50–100%','≥100%']
    fig.legend(handles=[Rectangle((0,0),1,1,fc=c,label=l) for c,l in zip(PALETTE,labels)],
               loc='lower center',bbox_to_anchor=(.5,y),ncol=7,frameon=False,fontsize=10)


def scorecards(data,out,n):
    for metric in [0,1]:
        fig,axs=plt.subplots(3,3,figsize=(23,20))
        fig.subplots_adjust(left=.135,right=.985,top=.895,bottom=.115,hspace=.48,wspace=.37)
        kind='Maximum physical energy error' if metric==0 else 'Position error relative to reference motion'
        fig.suptitle(f'{kind} · {n*60} Hz',x=.04,y=.975,ha='left',fontsize=25,weight='bold')
        fig.text(.04,.948,f'{n} substeps per 60 Hz frame · one frozen linearization per substep · {data.duration:g} s · identical equations across schedules',fontsize=12)
        fig.text(.04,.926,'Cells: 100 max|E−E₀| / (E₀−Vmin)' if metric==0 else 'Cells: 100 position RMS / reference RMS excursion from its initial position',fontsize=12)
        for ax,c in zip(axs.flat,data.cases):
            values=np.zeros((len(SCHEDULES),len(POLICIES)));stops=np.zeros_like(values,dtype=bool)
            for i,s in enumerate(SCHEDULES):
                for j,p in enumerate(POLICIES):
                    r,_=data.get(c,s+'_'+p,n);values[i,j]=metrics(r)[metric];stops[i,j]=not complete(r)
            heatmap(ax,values,stops)
            b,_=data.get(c,'baseline',n)
            ax.set_title(f'{NAMES[CASES.index(c)]}\nOrdinary XPBD: {number(metrics(b)[metric])}%',loc='left',fontsize=12,pad=13)
        legend_bins(fig,.065)
        stamp(fig,'Fractional repair limits predicted repair to ½ of the material-only linearized residual; cap limits the extra pair displacement to 0.1 × motion scale.\nSTOP = incomplete run; no prefix is ranked. Completion and low energy error do not establish correct motion. Spectral is a local dual-block preconditioner.')
        save(fig,out,f'schedule_policy_{"energy" if metric==0 else "position"}_{n*60}hz.png')


def overall(data,out):
    fig,axs=plt.subplots(1,2,figsize=(18,10.5))
    fig.subplots_adjust(left=.15,right=.985,top=.81,bottom=.19,wspace=.43)
    fig.suptitle('Schedule × acceptance: worst result over all 27 conditions',x=.04,y=.97,ha='left',fontsize=22,weight='bold')
    fig.text(.04,.925,f'9 cases × 4 / 8 / 16 substeps per 60 Hz frame · {data.duration:g} seconds · same original objective and hard endpoint energy',fontsize=11)
    for metric,ax in enumerate(axs):
        v=np.zeros((len(SCHEDULES),len(POLICIES)));stops=np.zeros_like(v,dtype=bool)
        for i,s in enumerate(SCHEDULES):
            for j,p in enumerate(POLICIES):
                rows=[r for r in data.rows if r['method']==s+'_'+p and r['duration']==data.duration]
                v[i,j]=max(metrics(r)[metric] for r in rows);stops[i,j]=any(not complete(r) for r in rows)
        heatmap(ax,v,stops)
        baseline=max(metrics(r)[metric] for r in data.rows if r['method']=='baseline' and r['duration']==data.duration)
        ax.set_title(('Maximum energy error / (E₀−Vmin)' if metric==0 else 'Position RMS / reference excursion')+f'\nOrdinary XPBD worst: {number(baseline)}%',loc='left',fontsize=13,pad=14)
    legend_bins(fig,.105)
    stamp(fig,'All cells are percentages; all 1,809 runs contribute. STOP = at least one incomplete condition. These are separate worst cases, not one trajectory.\nEnergy and position use different normalizations. Fractional repair = 0.5; displacement cap = 0.1. One frozen linearization per substep.')
    save(fig,out,'schedule_policy_worst_cases.png')


def axis_scale(ax,arrays,reference,positive=False):
    v=np.concatenate(arrays+[reference]);lo=float(np.min(v));hi=float(np.max(v));typical=max(1,float(np.max(np.abs(reference))))
    if max(abs(lo),abs(hi))>8*typical:
        ax.set_yscale('symlog',linthresh=typical)
        ticks=[0.]+[s*10.**p for s in [-1,1] for p in range(-3,20) if 10.**p>=typical and lo<=s*10.**p<=hi]
        ax.set_yticks(sorted(ticks));ax.text(.015,.96,'symlog · full range',transform=ax.transAxes,va='top',fontsize=8,color='#995522')
    else:
        pad=max(.06*(hi-lo),.01);ax.set_ylim(max(0,lo-pad) if positive else lo-pad,hi+pad)


def method_legend(fig,methods,y,reference=True):
    hh=[Line2D([],[],color=COLORS[m],lw=2,label=LABELS[m]) for m in methods]
    if reference:hh.append(Line2D([],[],color='#172c3a',ls=(0,(4,3)),lw=1.2,label='Refined RK4 reference'))
    fig.legend(handles=hh,loc='upper left',bbox_to_anchor=(.045,y),ncol=3,frameon=False,fontsize=10,columnspacing=2.5)


def trajectories(data,out,cases,suffix):
    fig=plt.figure(figsize=(20,4*len(cases)+3.0))
    fig.text(.055,.974,'Changing the frozen solve: energy and actual motion',fontsize=23,weight='bold')
    fig.text(.055,.95,f'8 substeps per 60 Hz frame = 480 Hz · one linearization per substep · {data.duration:g} seconds · fractional repair ½ + cap 0.1',fontsize=11)
    method_legend(fig,MAIN,.936)
    grid=fig.add_gridspec(len(cases),3,left=.06,right=.98,top=.858,bottom=.105,hspace=.83,wspace=.25)
    for ri,c in enumerate(cases):
        axs=[fig.add_subplot(grid[ri,j]) for j in range(3)]
        b,bt=data.get(c,'baseline');es=b['Escale'];vmin=minimum(b)
        refs=[(vals(bt,'ref_energy')-vmin)/es,vals(bt,'ref_x'),vals(bt,'ref_y')]
        arrays=[[] for _ in range(3)]
        for m in MAIN:
            s,tr=data.get(c,m);tt=vals(tr,'time')
            vv=[(vals(tr,'energy')-vmin)/es,vals(tr,'x'),vals(tr,'y')]
            for j,ax in enumerate(axs):
                ax.plot(tt,vv[j],color=COLORS[m],lw=1.1,alpha=.94)
                arrays[j].append(vv[j])
                if not complete(s):ax.plot(tt[-1],vv[j][-1],marker='x',color=COLORS[m],ms=8)
        for j,ax in enumerate(axs):
            ax.plot(vals(bt,'time'),refs[j],color='#172c3a',ls=(0,(4,3)),lw=1.2,alpha=.9)
            axis_scale(ax,arrays[j],refs[j],positive=j==0)
            ax.set_xlim(0,data.duration);ax.grid(alpha=.2);ax.set_xlabel('Time (s)',fontsize=9)
            ax.set_title(['(E − Vmin) / (E₀ − Vmin)','x position (m)','y position (m)'][j],loc='left',fontsize=11)
        axs[0].text(0,1.24,NAMES[CASES.index(c)],transform=axs[0].transAxes,fontsize=14,weight='bold')
        r,_=data.get(c,'LE_repair_cap');e,p=metrics(r)
        axs[0].text(0,-.4,f'LE repair + cap: max energy {e:.3g}% · normalized position RMS {p:.3g}% · path ratio {r["path_ratio"]:.3f}',transform=axs[0].transAxes,fontsize=9,color='#285c54')
    stamp(fig,'Energy target = 1 on the normalized energy axis. Every sampled state is drawn, including the raw EL release spike; symlog panels show their full range.\nA flat position trace with conserved energy is a frozen solution. Curves are Cartesian positions, not position errors. × marks a guarded stop.')
    save(fig,out,f'schedules_energy_motion_{suffix}_480hz.png')


def diagnostics(data,out):
    methods=['EL_raw','LE_raw','joint_raw','LE_repair_cap']
    fig,axs=plt.subplots(len(methods),2,figsize=(18,19.5))
    fig.subplots_adjust(left=.075,right=.98,top=.88,bottom=.09,hspace=.55,wspace=.25)
    fig.suptitle('Spring release: what the frozen energy prediction misses',x=.055,y=.974,ha='left',fontsize=23,weight='bold')
    fig.text(.055,.944,'First 24 small substeps at 480 Hz · one frozen system per substep · signed residuals relative to the same E₀',fontsize=11)
    ehandles=[Line2D([],[],color=co,ls=ls,lw=1.5,label=lab) for co,ls,lab in [('#737f89','--','Material-only: linear prediction'),('#7f6d9e',':','Material-only: actual energy'),('#cb7528','--','Accepted pair: linear prediction'),('#146b61','-','Accepted pair: actual energy')]]
    fig.legend(handles=ehandles,loc='upper left',bbox_to_anchor=(.055,.922),ncol=2,frameon=False,fontsize=10)
    for i,m in enumerate(methods):
        _,tr=data.get('spring_release',m);tr=tr[1:25];step=vals(tr,'step');es=tr[0]['Escale']
        lines=[vals(tr,'H0')/es,vals(tr,'actual_H0')/es,vals(tr,'predicted_H')/es,(vals(tr,'energy')-tr[0]['target_energy'])/es]
        for y,co,ls in zip(lines,['#737f89','#7f6d9e','#cb7528','#146b61'],['--',':','--','-']):
            axs[i,0].plot(step,y,color=co,ls=ls,lw=1.6,marker='o',ms=2.5)
        ax=axs[i,0];ax.axhline(0,color='#374b59',lw=.6);ax.set_yscale('symlog',linthresh=1e-6)
        ax.set_ylabel('(E − E₀) / Escale · symlog');ax.set_title(LABELS[m],loc='left',fontsize=14,weight='bold')
        ax=axs[i,1]
        for field,co,ls,lab in [('raw_energy_displacement','#ce7326','-','Energy contribution |g Δη / m|'),('material_response_displacement','#a13f67','--','Material response |n Δλextra / m|'),('candidate_extra_displacement','#788796',':','Combined candidate extra'),('extra_displacement','#168270','-','Accepted combined extra')]:
            ax.plot(step,vals(tr,field),color=co,ls=ls,lw=1.5,marker='o',ms=2.2,label=lab)
        ax.set_yscale('symlog',linthresh=1e-12);ax.set_ylabel('Displacement magnitude (m) · symlog')
        if i==0:ax.legend(fontsize=8.5,frameon=False,loc='upper right')
        for ax in axs[i]:ax.set_xlabel('Substep');ax.grid(alpha=.2);ax.set_xlim(1,24)
        accepted=vals(tr,'blend')
        axs[i,1].text(0,-.3,f'Accepted pair fraction: min {accepted.min():.3g}, max {accepted.max():.3g}',transform=axs[i,1].transAxes,fontsize=10)
    stamp(fig,'The energy and material components can cancel. Acceptance acts on the combined pair difference from ordinary XPBD, not on the energy component alone.\nThe frozen linear energy residual is not the nonlinear physical energy change. Zero curves remain visible in the linear center of each symlog axis.')
    save(fig,out,'spring_release_prediction_and_cancellation.png')


def compact(data,out):
    methods=['baseline','LE_raw','LE_repair_cap','spectral_repair_cap']
    fig,axs=plt.subplots(3,2,figsize=(15,13.8))
    fig.subplots_adjust(left=.075,right=.98,top=.82,bottom=.11,hspace=.60,wspace=.21)
    fig.suptitle('Energy repair restores release motion, but distorts the orbit',x=.055,y=.974,ha='left',fontsize=21,weight='bold')
    fig.text(.055,.945,'480 Hz · 8 small substeps per 60 Hz frame · one frozen solve per substep · fixed repair fraction ½ and cap 0.1',fontsize=10.5)
    method_legend(fig,methods,.928)
    for i,c in enumerate(['vertical_oscillation','spring_release','spring_orbit']):
        b,bt=data.get(c,'baseline');vmin=minimum(b);es=b['Escale'];coordinate='y' if c=='vertical_oscillation' else 'x'
        for m in methods:
            _,tr=data.get(c,m);tt=vals(tr,'time')
            axs[i,0].plot(tt,(vals(tr,'energy')-vmin)/es,color=COLORS[m],lw=1.4)
            axs[i,1].plot(tt,vals(tr,coordinate),color=COLORS[m],lw=1.25)
        axs[i,0].plot(vals(bt,'time'),(vals(bt,'ref_energy')-vmin)/es,color='#172c3a',ls=(0,(4,3)),lw=1.2)
        axs[i,1].plot(vals(bt,'time'),vals(bt,'ref_'+coordinate),color='#172c3a',ls=(0,(4,3)),lw=1.2)
        axs[i,0].set_title(NAMES[CASES.index(c)],loc='left',fontsize=14,weight='bold')
        axs[i,0].set_ylabel('(E − Vmin) / (E₀ − Vmin)')
        axs[i,1].set_title(coordinate+' position (m)',loc='left',fontsize=12)
        r,_=data.get(c,'LE_repair_cap');e,p=metrics(r)
        axs[i,0].text(0,-.34,f'LE repair + cap: energy error {e:.3g}% · position RMS {p:.3g}% of reference excursion',transform=axs[i,0].transAxes,fontsize=9,color='#285c54')
        for ax in axs[i]:ax.set_xlim(0,data.duration);ax.set_xlabel('Time (s)');ax.grid(alpha=.2)
    stamp(fig,'The same accepted correction works well for released scalar motion but accelerates the orbit. Nearly constant total energy is insufficient.\nRaw EL spikes are shown separately in the full-suite and first-substep figures; this comparison contains only the four explicitly listed methods.')
    save(fig,out,'core_energy_and_motion_480hz.png')


def orbit_partition(data,out):
    methods=['baseline','LE_repair_cap','spectral_repair_cap']
    fig,axs=plt.subplots(2,3,figsize=(19,11.5))
    fig.subplots_adjust(left=.065,right=.98,top=.79,bottom=.10,hspace=.40,wspace=.25)
    fig.suptitle('The spring orbit keeps total energy by changing its distribution',x=.05,y=.97,ha='left',fontsize=21,weight='bold')
    fig.text(.05,.935,f'No gravity: angular momentum should stay constant · 480 Hz · one frozen solve per substep · {data.duration:g} seconds',fontsize=11)
    method_legend(fig,methods,.914)
    b,bt=data.get('spring_orbit','baseline');es=b['Escale'];ell0=bt[0]['angular_momentum']
    def fields(tr,reference=False):
        prefix='ref_' if reference else ''
        x=vals(tr,prefix+'x');y=vals(tr,prefix+'y');radius=np.hypot(x,y)
        ell=vals(tr,prefix+'angular_momentum');energy=vals(tr,prefix+'energy')
        potential=(radius-b['length'])**2/(2*b['alpha'])
        tangential=ell**2/(2*b['mass']*radius**2)
        radial=energy-potential-tangential
        if not reference:
            direct=.5*b['mass']*((x*vals(tr,'vx')+y*vals(tr,'vy'))/radius)**2
            assert np.max(np.abs(direct-radial))<5e-12
        assert np.min(radial)>-5e-11
        return [energy/es,ell/ell0,np.unwrap(np.arctan2(y,x)),radial/es,tangential/es,potential/es]
    ref=fields(bt,True)
    for m in methods:
        _,tr=data.get('spring_orbit',m)
        for ax,v in zip(axs.flat,fields(tr)):
            ax.plot(vals(tr,'time'),v,color=COLORS[m],lw=1.35)
    titles=['Total E / E₀','Angular momentum L / L₀','Unwrapped polar angle (rad)',
            'Radial kinetic energy / E₀','Tangential kinetic energy / E₀','Spring potential energy / E₀']
    for ax,title,v in zip(axs.flat,titles,ref):
        ax.plot(vals(bt,'time'),v,color='#172c3a',ls=(0,(4,3)),lw=1.2)
        ax.set_title(title,loc='left',fontsize=12);ax.set_xlabel('Time (s)');ax.set_xlim(0,data.duration);ax.grid(alpha=.2)
    final=data.get('spring_orbit','LE_repair_cap')[1][-1]
    drift=100*(final['angular_momentum']/ell0-1)
    delta_angle=fields(data.get('spring_orbit','LE_repair_cap')[1])[2][-1]-ref[2][-1]
    stamp(fig,f'LE repair + cap: final angular-momentum drift = {drift:.2f}%; phase lead = {delta_angle:.3f} rad. Total energy error remains below {metrics(data.get("spring_orbit","LE_repair_cap")[0])[0]:.4f}%.\nEnergy partitions use the recorded mass, radius, velocity, and angular momentum; the reference radial part is Eref − Uref − Lref²/(2mrref²), without numerical differentiation.')
    save(fig,out,'spring_orbit_energy_partition_and_momentum.png')


def convergence(data,out):
    methods=['baseline','EL_raw','LE_raw','LE_repair_cap','spectral_repair_cap','joint_repair_cap']
    fig,axs=plt.subplots(3,3,figsize=(18,16))
    fig.subplots_adjust(left=.07,right=.98,top=.84,bottom=.09,hspace=.49,wspace=.25)
    fig.suptitle('Does taking smaller substeps fix the motion?',x=.05,y=.972,ha='left',fontsize=23,weight='bold')
    fig.text(.05,.944,f'Position RMS / reference excursion · fixed parameters · {data.duration:g} seconds · one frozen system per substep',fontsize=11)
    method_legend(fig,methods,.926,reference=False)
    for ax,c in zip(axs.flat,data.cases):
        for m in methods:
            rr=[data.get(c,m,n)[0] for n in data.substeps]
            ax.plot([60*r['substeps'] for r in rr],[metrics(r)[1] for r in rr],color=COLORS[m],lw=1.7,marker='o',ms=4)
        ax.set_xscale('log',base=2);ax.set_xticks([240,480,960],['240','480','960'])
        ax.set_yscale('log');ax.set_xlabel('Substep frequency (Hz)');ax.set_ylabel('Normalized position RMS (%)')
        ax.set_title(NAMES[CASES.index(c)],loc='left',fontsize=12);ax.grid(alpha=.22,which='both')
    stamp(fig,'Each point uses a separate run and its own refined RK4 reference. A flat error curve as the substep shrinks indicates unresolved dynamics, even if energy is nearly conserved.\nNo parameter is tuned separately for a case or a substep count. Raw LE remains frozen in the released scalar spring.')
    save(fig,out,'schedule_motion_convergence.png')


def selection(data):
    rows=[]
    for m in data.methods:
        rr=[r for r in data.rows if r['method']==m and r['duration']==data.duration]
        rows.append({'method':m,'complete_all':all(complete(r) for r in rr),
                     'worst_energy_error_percent':max(metrics(r)[0] for r in rr),
                     'worst_normalized_position_rms_percent':max(metrics(r)[1] for r in rr),
                     'minimum_path_ratio':min(r['path_ratio'] for r in rr)})
    eligible=[r for r in rows if r['complete_all'] and r['worst_energy_error_percent']<=5]
    order=sorted(eligible,key=lambda r:(r['worst_normalized_position_rms_percent'],r['worst_energy_error_percent']))
    return {'conditions_per_method':len(data.cases)*len(data.substeps),'rule':'For descriptive ranking only: require every condition to finish and worst max energy error/Escale <= 5%; then sort by worst position RMS/reference excursion. This is not a claim of physical correctness.',
            'energy_qualified_ranking':order,'all_methods':rows,'displayed_methods':MAIN}


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('sweep',type=Path);p.add_argument('selected',type=Path)
    p.add_argument('--output',type=Path,default=Path('results/correction_schedules/plots'))
    p.add_argument('--duration',type=float,default=10.)
    p.add_argument('--all-substeps',action='store_true')
    args=p.parse_args();args.output.mkdir(parents=True,exist_ok=True)
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    data=Data(args.sweep,args.selected,args.duration)
    for n in data.substeps if args.all_substeps else [8]:scorecards(data,args.output,n)
    overall(data,args.output)
    trajectories(data,args.output,data.cases[:5],'elastic')
    trajectories(data,args.output,data.cases[5:],'spring_and_rod')
    diagnostics(data,args.output);convergence(data,args.output)
    compact(data,args.output);orbit_partition(data,args.output)
    (args.output/'selection.json').write_text(json.dumps(selection(data),indent=2)+'\n')
    (args.output/'plot_audit.json').write_text(json.dumps(data.audit,indent=2)+'\n')
    print(json.dumps(data.audit,indent=2))
    print('Figures: '+', '.join(p.name for p in sorted(args.output.glob('*.png'))))


if __name__=='__main__':main()
