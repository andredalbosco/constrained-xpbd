#!/usr/bin/env python3
"""Plot directional clipping against unguarded EL and ordinary XPBD."""
import argparse
import collections
import csv
import hashlib
import json
import math
import os
from pathlib import Path

os.environ.setdefault('MPLCONFIGDIR', '/tmp/xpbd-directional-matplotlib')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D

CASES = ['vertical_oscillation', 'small_swing', 'swing_rest', 'swing_tangent', 'coupled_motion']
NAMES = ['Vertical oscillation', 'Small swing', 'Large swing released from rest',
         'Swing with tangential velocity', 'Coupled radial and swinging motion']
INITIAL = ['r₀ = 1.2981 m, θ₀ = 0, v₀ = 0',
           'r₀ = 1.0981 m, θ₀ = 0.15 rad, v₀ = 0',
           'r₀ = 1.25 m, θ₀ = 0.65 rad, v₀ = 0',
           'r₀ = 1.0981 m, θ₀ = 0.65 rad, vₜ = 1.2 m/s',
           'r₀ = 1.35 m, θ₀ = 0.65 rad, vᵣ = −0.7 m/s, vₜ = 0.9 m/s']
METHODS = [('baseline', -1.), ('unguarded', -1.), ('directional', -1.)]
COLORS = ['#939da7', '#cf7b22', '#1469ab']
STYLES = ['-', '--', '-']
LABELS = ['Ordinary XPBD', 'EL · unguarded energy correction', 'EL · directional clipping']
VMIN = -10.2911805
TEXT = {'case', 'order', 'gradient', 'energy_policy', 'status'}


def read(path):
    with path.open(newline='') as f:
        rows = list(csv.DictReader(f))
    for r in rows:
        for k in r:
            if k not in TEXT:
                r[k] = float(r[k])
                if not math.isfinite(r[k]):
                    raise ValueError(f'Nonfinite {k} in {path}')
    return rows


def key(row):
    return row['case'], row['energy_policy'], row['gamma'], round(1/row['h']), int(row['sweeps'])


def values(trace, field):
    return np.array([r[field] for r in trace])


class Data:
    def __init__(self, directional_root, unguarded_root):
        roots=[directional_root,unguarded_root]
        self.rows = [r for root in roots for r in read(root/'summary.csv')]
        self.summary = {}
        for r in self.rows:
            if key(r) in self.summary:
                assert r==self.summary[key(r)], 'Baseline summary differs across studies'
            assert r['gamma']==-1 and r['oversized_skips']==0
            if r['energy_policy']=='directional':
                assert r['max_accepted_objective_work']<=0
            self.summary[key(r)]=r
        assert len(self.rows)==180 and len(self.summary)==135
        self.traces = collections.defaultdict(list)
        self.audit = {'summary_runs': len(self.rows), 'max_scaled_energy_reconstruction_error': 0.,
                      'max_scaled_angular_reconstruction_error': 0.,
                      'max_position_rms_reconstruction_error': 0.,
                      'source_summary_sha256': {str(root):hashlib.sha256((root/'summary.csv').read_bytes()).hexdigest() for root in roots}}
        raw_groups=set(); duplicate_baseline=[]
        for ri,root in enumerate(roots):
          for r in read(root/'traces.csv'):
            raw_groups.add((ri,*key(r)))
            assert r['sweeps'] == 5 and round(1/r['h']) in [60, 240]
            assert r['retain_rx'] == 0
            assert r['gradient'] == ('none' if r['order'] == 'baseline' else 'fresh')
            energy = .5*(r['vx']**2+r['vy']**2)+50*(math.hypot(r['x'],r['y'])-1)**2+9.81*r['y']
            discrepancy = abs(energy-r['energy'])/max(1, abs(energy))
            self.audit['max_scaled_energy_reconstruction_error'] = max(self.audit['max_scaled_energy_reconstruction_error'], discrepancy)
            assert discrepancy < 2e-12
            angular = r['x']*r['vy']-r['y']*r['vx']
            discrepancy = abs(angular-r['angular_momentum'])/max(1, abs(angular))
            self.audit['max_scaled_angular_reconstruction_error'] = max(self.audit['max_scaled_angular_reconstruction_error'], discrepancy)
            assert discrepancy < 2e-12
            if ri==1 and r['energy_policy']=='baseline':
                duplicate_baseline.append(r)
            else:
                self.traces[key(r)].append(r)
        for r in duplicate_baseline:
            assert r==self.traces[key(r)][int(r['step'])], 'Baseline trace differs across studies'
        for k, tr in self.traces.items():
            s = self.summary[k]
            assert len(tr) == int(s['completed_steps'])+1
            for i, r in enumerate(tr):
                assert r['step'] == i and math.isclose(r['time'], i*r['h'], abs_tol=1e-12)
                assert math.isclose(r['Escale'], r['target_energy']-VMIN, rel_tol=1e-12)
            error = max(abs(r['energy']-r['target_energy']) for r in tr)
            assert math.isclose(error, s['max_abs_energy_error'], abs_tol=2e-12, rel_tol=2e-11)
            p = np.hypot(values(tr,'x')-values(tr,'ref_x'), values(tr,'y')-values(tr,'ref_y'))
            rms = float(np.sqrt(np.mean(p[1:]**2)))
            self.audit['max_position_rms_reconstruction_error'] = max(self.audit['max_position_rms_reconstruction_error'], abs(rms-s['rms_position_error']))
            assert math.isclose(rms,s['rms_position_error'],abs_tol=1e-12,rel_tol=1e-12)
            assert tr[-1]['directional_skips']==s['directional_skips']
            assert tr[-1]['energy_attempts']==s['energy_attempts']
            assert tr[-1]['energy_accepted']==s['energy_accepted']
            for a,b in zip(tr,tr[1:]):
                assert b['directional_skips']>=a['directional_skips']
        self.audit['audited_trace_groups_including_duplicate_baselines'] = len(raw_groups)
        self.audit['unique_trace_groups'] = len(self.traces)
        self.audit['identical_baselines'] = True
        assert len(self.traces) == 30 and len(raw_groups)==40

    def get(self, case, method, hz):
        k = (case, *method, hz, 5)
        return self.summary[k], self.traces[k]

    def arrays(self, case, method, hz, normalize=False):
        s, tr = self.get(case, method, hz)
        e = values(tr, 'energy')
        if normalize:
            e = (e-VMIN)/s['Escale']
        return values(tr,'time'), [e, values(tr,'x'), values(tr,'y'),
                                  np.hypot(values(tr,'x')-values(tr,'ref_x'),values(tr,'y')-values(tr,'ref_y'))]


def legend(fig, y):
    handles = [Line2D([], [], color=c, ls=ls, lw=2, label=label) for c,ls,label in zip(COLORS,STYLES,LABELS)]
    handles.append(Line2D([], [], color='#1a2933', lw=1.4, ls=(0,(4,3)), label='Refined RK4 reference'))
    fig.legend(handles=handles, loc='upper left', bbox_to_anchor=(.045,y), ncol=4, frameon=False,
               handlelength=2.7, columnspacing=2.5, fontsize=11)


def shared_limits(data, case, col, normalize=False, energy_hz=None):
    allvalues = []
    for hz in [60,240]:
        if col==0 and energy_hz is not None and hz!=energy_hz: continue
        for method in METHODS:
            _, v = data.arrays(case,method,hz,normalize)
            allvalues.extend(v[col])
        s,tr = data.get(case,METHODS[0],hz)
        if col < 3:
            r = values(tr,['ref_energy','ref_x','ref_y'][col])
            if col == 0 and normalize:
                r=(r-VMIN)/s['Escale']
            allvalues.extend(r)
    lo,hi=min(allvalues),max(allvalues)
    pad=max(.07*(hi-lo),.015 if col else .01)
    if col == 3:
        return 0,hi+pad
    return lo-pad,hi+pad


def row_note(data,case,hz):
    parts=[]
    for m,label in zip(METHODS,['XPBD','unguarded EL','directional EL']):
        s,_=data.get(case,m,hz)
        if s['completed_steps'] != s['requested_steps']:
            parts.append(f'{label}: STOP')
        else:
            parts.append(f"{label}: {100*s['max_abs_energy_error']/s['Escale']:.3g}% / {s['rms_position_error']:.3g} m / {s['path_ratio']:.3f}×")
    s,_=data.get(case,('directional',-1.),hz)
    return 'Max energy error / position RMS / path ratio:  '+ '     |     '.join(parts)+f"     |     directional skips: {int(s['directional_skips'])}"


def energy_scale(ax,data,case,hz,normalize=False):
    s,_=data.get(case,('directional',-1.),hz)
    if s['max_abs_energy_error']>10*s['Escale']:
        threshold=1. if normalize else s['Escale']
        ax.set_yscale('symlog',linthresh=threshold)
        lo,hi=ax.get_ylim()
        ticks=[0.]+[sign*10.**p for sign in [-1,1] for p in range(-3,9)
                    if 10.**p>=threshold and lo<=sign*10.**p<=hi]
        ax.set_yticks(sorted(ticks))
        ax.text(.03,.94,'symlog axis · full spike shown',transform=ax.transAxes,
                va='top',fontsize=9,color='#8b4d20')


def full_figure(data,hz,out):
    fig=plt.figure(figsize=(22,23.5),facecolor='white')
    fig.text(.05,.973,f'Directional clipping · elastic pendulum · {hz} Hz',fontsize=25,weight='bold')
    fig.text(.05,.949,'Gate: discard the raw energy correction when ∇F · Δx_E > 0. F is the full implicit objective; no magnitude threshold.',fontsize=12)
    fig.text(.05,.929,'Fresh EL · dropped rₓ · five Gauss–Seidel sweeps · m = 1 kg, rest length = 1 m, k = 100 N/m, gravity = 9.81 m/s².',fontsize=11)
    legend(fig,.918)
    grid=fig.add_gridspec(5,4,left=.055,right=.985,top=.85,bottom=.067,hspace=.80,wspace=.26)
    for ri,case in enumerate(CASES):
        axs=[fig.add_subplot(grid[ri,c]) for c in range(4)]
        for mi,(method,color,style) in enumerate(zip(METHODS,COLORS,STYLES)):
            s,tr=data.get(case,method,hz)
            t,v=data.arrays(case,method,hz)
            for col,ax in enumerate(axs):
                ax.plot(t,v[col],color=color,ls=style,lw=1.25,alpha=.96)
                if s['completed_steps'] != s['requested_steps']:
                    ax.plot(t[-1],v[col][-1],marker='X',color=color,ms=7)
        s,tr=data.get(case,METHODS[0],hz)
        t=values(tr,'time')
        reference=[values(tr,'ref_energy'),values(tr,'ref_x'),values(tr,'ref_y'),np.zeros_like(t)]
        for col,(ax,r,title,ylabel) in enumerate(zip(axs,reference,['Mechanical energy','Horizontal position','Vertical position','Distance from reference'],['E (J)','x (m)','y (m)','|x − x_ref| (m)'])):
            ax.plot(t,r,color='#1a2933',ls=(0,(4,3)),lw=1.1,zorder=6)
            ax.set_title(title,fontsize=11,pad=9)
            ax.set(xlim=(0,10),ylim=shared_limits(data,case,col,energy_hz=hz),xlabel='Time (s)',ylabel=ylabel)
            ax.grid(alpha=.19)
        energy_scale(axs[0],data,case,hz)
        p=axs[0].get_position()
        fig.text(.055,p.y1+.033,NAMES[ri],fontsize=14,weight='bold')
        fig.text(.32,p.y1+.034,INITIAL[ri],fontsize=10.5,color='#465461')
        fig.text(.055,p.y0-.031,row_note(data,case,hz),fontsize=9.5,color='#344b5b')
    fig.text(.055,.012,'Energy error is normalized by E₀ − Vmin, Vmin = −10.2911805 J. Position axes match across timesteps; energy axes fit each timestep. No energy spikes are omitted.',fontsize=10.5)
    path=out/f'directional_vs_unguarded_h{hz}_sweeps5.png'
    fig.savefig(path,dpi=150,facecolor='white');plt.close(fig)
    return path


def failure_figure(data,out):
    fig=plt.figure(figsize=(18,15),facecolor='white')
    fig.text(.06,.969,'Directional clipping: improvement and remaining failures',fontsize=23,weight='bold')
    fig.text(.06,.945,'Discard Δx_E when ∇F · Δx_E > 0, using the full implicit-objective gradient. Fresh EL · dropped rₓ · five sweeps.',fontsize=11.5)
    legend(fig,.931)
    grid=fig.add_gridspec(3,3,left=.065,right=.975,top=.84,bottom=.08,hspace=.78,wspace=.28)
    rows=[('vertical_oscillation',60),('vertical_oscillation',240),('coupled_motion',240)]
    for ri,(case,hz) in enumerate(rows):
        axs=[fig.add_subplot(grid[ri,c]) for c in range(3)]
        coordinate='y' if case=='vertical_oscillation' else 'x'
        for method,color,style in zip(METHODS,COLORS,STYLES):
            s,tr=data.get(case,method,hz);t,v=data.arrays(case,method,hz,True)
            axs[0].plot(t,v[0],color=color,ls=style,lw=1.35)
            axs[1].plot(t,values(tr,coordinate),color=color,ls=style,lw=1.2)
            axs[2].plot(t,v[3],color=color,ls=style,lw=1.2)
        s,tr=data.get(case,METHODS[0],hz);t=values(tr,'time')
        axs[0].plot(t,(values(tr,'ref_energy')-VMIN)/s['Escale'],color='#1a2933',ls=(0,(4,3)),lw=1.2)
        axs[1].plot(t,values(tr,'ref_'+coordinate),color='#1a2933',ls=(0,(4,3)),lw=1.2)
        axs[2].plot(t,np.zeros_like(t),color='#1a2933',ls=(0,(4,3)),lw=1.2)
        lim=shared_limits(data,case,0,True,energy_hz=hz)
        axs[0].set(xlim=(0,10),ylim=(max(0,lim[0]),lim[1]),xlabel='Time (s)',ylabel='(E − Vmin) / (E₀ − Vmin)',title='Energy above minimum / initial')
        energy_scale(axs[0],data,case,hz,True)
        axs[1].set(xlim=(0,10),ylim=shared_limits(data,case,2 if coordinate=='y' else 1),xlabel='Time (s)',ylabel=coordinate+' (m)',title='Vertical position' if coordinate=='y' else 'Horizontal position')
        upper=max(float(np.max(data.arrays(case,m,hz)[1][3])) for m in METHODS)
        axs[2].set(xlim=(0,10),ylim=(0,upper*1.08),xlabel='Time (s)',ylabel='|x − x_ref| (m)',title='Distance from reference')
        for ax in axs: ax.grid(alpha=.19)
        p=axs[0].get_position()
        name=NAMES[CASES.index(case)]
        fig.text(.065,p.y1+.037,f'{name} · {hz} Hz',fontsize=14,weight='bold')
        sg,_=data.get(case,('directional',-1.),hz)
        note=f"Directional EL: {int(sg['directional_skips'])} discarded corrections; max energy error {100*sg['max_abs_energy_error']/sg['Escale']:.3g}%; position RMS {sg['rms_position_error']:.4f} m; path {sg['path_ratio']:.3f}× reference."
        fig.text(.065,p.y0-.042,note,fontsize=10.5,color='#465461')
    fig.text(.065,.013,'Energy uses endpoint velocity. Vmin = −10.2911805 J. The directional derivative is a first-order test; accepting a correction does not guarantee finite-step descent.',fontsize=9.8)
    path=out/'directional_energy_and_motion.png'
    fig.savefig(path,dpi=160,facecolor='white');plt.close(fig)
    return path


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('directional_results',type=Path)
    p.add_argument('unguarded_results',type=Path)
    p.add_argument('--output',type=Path)
    args=p.parse_args();out=args.output or args.directional_results/'plots';out.mkdir(parents=True,exist_ok=True)
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,
                         'axes.spines.right':False,'axes.edgecolor':'#aab5bd','text.color':'#193044',
                         'axes.labelcolor':'#304757','axes.titlecolor':'#193044'})
    data=Data(args.directional_results,args.unguarded_results)
    images=[full_figure(data,hz,out) for hz in [60,240]]
    images.append(failure_figure(data,out))
    data.audit['images']=[str(p) for p in images]
    data.audit['plotted_method_groups']=30
    data.audit['selected_methods']=[list(x) for x in METHODS]
    for case in CASES:
        for normalize in [False,True]:
            for col in range(4):
                lo,hi=shared_limits(data,case,col,normalize)
                for hz in [60,240]:
                    for method in METHODS:
                        _,v=data.arrays(case,method,hz,normalize)
                        assert np.min(v[col])>=lo and np.max(v[col])<=hi
    data.audit['axes_clipped']=False
    (out/'plot_audit.json').write_text(json.dumps(data.audit,indent=2)+'\n')
    print(json.dumps(data.audit,indent=2))


if __name__=='__main__':
    main()
