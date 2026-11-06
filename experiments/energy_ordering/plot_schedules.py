#!/usr/bin/env python3
"""Plot physical energy and Cartesian positions for frozen XPBD energy schedules."""
import argparse
import collections
import csv
import json
import math
import os
from pathlib import Path

os.environ.setdefault('MPLCONFIGDIR', '/tmp/xpbd-matplotlib')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np

CASES = ['spring_rest', 'spring_radial', 'spring_tangent', 'pendulum_rest', 'pendulum_tangent']
CASE_LABELS = {
    'spring_rest': 'Spring: released from rest',
    'spring_radial': 'Spring: initial radial velocity',
    'spring_tangent': 'Spring: initial tangential velocity',
    'pendulum_rest': 'Pendulum: released from rest',
    'pendulum_tangent': 'Pendulum: initial tangential velocity',
}
METHODS = ['baseline', 'jacobi', 'joint', 'gs_le', 'gs_el', 'gs_lel', 'gs_ele']
LABELS = {
    'baseline': 'Ordinary XPBD',
    'jacobi': 'Jacobi: independent multiplier updates',
    'joint': 'Joint multiplier solve',
    'gs_le': 'Gauss–Seidel: material → energy (LE)',
    'gs_el': 'Gauss–Seidel: energy → material (EL)',
    'gs_lel': 'Gauss–Seidel: material → energy → material (LEL)',
    'gs_ele': 'Gauss–Seidel: energy → material → energy (ELE)',
}
SHORT = {'baseline': 'XPBD', 'jacobi': 'Jacobi', 'joint': 'Joint', 'gs_le': 'GS LE', 'gs_el': 'GS EL', 'gs_lel': 'GS LEL', 'gs_ele': 'GS ELE'}
COLORS = {'baseline':'#60666e', 'jacobi':'#6953bc', 'joint':'#1769aa', 'gs_le':'#e07a14', 'gs_el':'#138765', 'gs_lel':'#a04d96', 'gs_ele':'#c54051'}
STYLES = {'baseline':'-', 'jacobi':'-', 'joint':'-', 'gs_le':'-', 'gs_el':'--', 'gs_lel':'-.', 'gs_ele':':'}
KEY = ['case', 'method', 'application', 'retain_rx', 'h', 'sweeps', 'softness']
TEXT = {'case', 'method', 'application', 'status'}


def read_csv(path):
    with path.open() as stream:
        rows = list(csv.DictReader(stream))
    for row in rows:
        row.setdefault('application','deferred')
        for field in row:
            if field not in TEXT:
                row[field] = float(row[field])
                if not math.isfinite(row[field]):
                    raise ValueError(f'Nonfinite field {field} in {path}: {row}')
    return rows


def key(row):
    return tuple(row[field] for field in KEY)


def close(a, b):
    return math.isclose(a, b, rel_tol=1e-11, abs_tol=1e-13)


def energy(row, reference=False):
    field = 'ref_energy' if reference else 'energy'
    if field in row:
        return row[field]
    x, y = (row['ref_x'], row['ref_y']) if reference else (row['x'], row['y'])
    kinetic = row['ref_kinetic'] if reference else row.get('kinetic', .5*(row['vx']**2+row['vy']**2))
    potential = 9.81*y if row['case'].startswith('pendulum') else 50*(math.hypot(x,y)-1)**2
    return kinetic+potential


class Data:
    def __init__(self, root, h, sweeps, inner):
        self.h, self.sweeps, self.inner = h, sweeps, inner
        all_rows = read_csv(root/'summary.csv')
        self.summaries = [r for r in all_rows if close(r['h'],h) and r['sweeps']==sweeps
                          and r.get('inner',inner)==inner and r['method'] in METHODS]
        self.methods=[m for m in METHODS if any(r['method']==m for r in self.summaries)]
        self.lookup = {key(r):r for r in self.summaries}
        if len(self.lookup) != len(self.summaries):
            raise ValueError('Duplicate summary keys')
        self.traces = collections.defaultdict(list)
        for row in read_csv(root/'traces.csv'):
            if key(row) in self.lookup:
                self.traces[key(row)].append(row)
        self.energy_audit = []
        for k, rows in self.traces.items():
            summary = self.lookup[k]
            if len(rows) != summary['completed_steps']+1:
                raise ValueError(f'Trace length mismatch {k}')
            if any(r['step'] != i or not close(r['time'],i*h) for i,r in enumerate(rows)):
                raise ValueError(f'Trace time/step mismatch {k}')
            E0 = rows[0].get('target_energy',energy(rows[0]))
            def physical_energy(r):
                potential = 9.81*r['y'] if r['case'].startswith('pendulum') else 50*(math.hypot(r['x'],r['y'])-1)**2
                return .5*(r['vx']**2+r['vy']**2)+potential
            error = max(abs(energy(r)-physical_energy(r)) for r in rows)
            scale = max(1,abs(E0),max(abs(energy(r)) for r in rows))
            if error > 2e-12*scale:
                raise ValueError(f'Energy is inconsistent with energy_error for {k}: {error}')
            self.energy_audit.append(error/scale)
            measured_error = max(abs(energy(r)-E0) for r in rows)
            if 'max_abs_energy_error' in summary and not math.isclose(measured_error,summary['max_abs_energy_error'],rel_tol=1e-10,abs_tol=1e-10):
                raise ValueError(f'Summary and trajectory energy disagree {k}')
            summary['max_abs_energy_error'] = measured_error

    def get(self, case, method, soft, retain, application='deferred'):
        selected = [r for r in self.summaries if r['case']==case and r['method']==method
                    and r['retain_rx']==(0 if method=='baseline' else retain)
                    and close(r['softness'],0 if method=='baseline' else soft)
                    and r['application']==application]
        if len(selected)!=1:
            raise ValueError(f'Expected one summary for {case,method,soft,retain}; got {len(selected)}')
        row = selected[0]
        return row,self.traces[key(row)]

    def reference(self, case):
        return self.get(case,'baseline',0,0)[1]


def coordinates(rows, reference=False):
    prefix = 'ref_' if reference else ''
    return [np.array([energy(r,reference) for r in rows]),
            np.array([r[prefix+'x'] for r in rows]),
            np.array([r[prefix+'y'] for r in rows])]


def scales(data, soft, retain):
    """Common bounded-run ranges; isolated runaway panels get marked symlog."""
    result = {}
    for case in CASES:
        ref = coordinates(data.reference(case), True)
        all_values = [v.tolist() for v in ref]
        for method in data.methods:
            applications=['deferred','immediate'] if method.startswith('gs_') else ['deferred']
            for application in applications:
                _, rows = data.get(case,method,soft,retain,application)
                for i, values in enumerate(coordinates(rows)):
                    extent=max(1.0,float(np.max(np.abs(ref[i]))))
                    allowed=(np.max(np.abs(values-rows[0]['target_energy'])) <= 10*rows[0]['Escale']) if i==0 else (np.max(np.abs(values)) <= 3*extent)
                    if allowed:
                        all_values[i].extend(values)
        for column in range(3):
            ref_extent = max(1.0, float(np.max(np.abs(ref[column]))))
            vals = all_values[column]+ref[column].tolist()
            lo, hi = min(vals), max(vals)
            spread = hi-lo
            pad = max(spread*.07, .05*ref_extent)
            ylim = (lo-pad,hi+pad)
            result[case,column] = {'symlog':False,'linthresh':ref_extent,
                                   'ylim':ylim}
    return result


def panel_scale(data, limits, case, column, trace_groups):
    scale=limits[case,column]
    ref=coordinates(data.reference(case),True)[column]
    values=ref.tolist()
    for rows in trace_groups:values.extend(coordinates(rows)[column])
    lo,hi=min(values),max(values)
    extent=scale['linthresh']
    sample=trace_groups[0][0]
    runaway=(max(abs(v-sample['target_energy']) for v in values) > 10*sample['Escale']) if column==0 else (max(abs(lo),abs(hi)) > 3*extent)
    if not runaway:return scale
    ylim=(lo*1.15 if lo < -extent else lo-.05*extent,
          hi*1.15 if hi > extent else hi+.05*extent)
    return {'symlog':True,'linthresh':extent,'ylim':ylim}


def setup_axis(ax, scale, duration, label):
    if scale['symlog']:
        ax.set_yscale('symlog',linthresh=scale['linthresh'])
        locator=ticker.SymmetricalLogLocator(base=10,linthresh=scale['linthresh'])
        locator.set_params(numticks=5)
        ticks=[v for v in locator.tick_values(*scale['ylim']) if v==0 or abs(v)>=scale['linthresh']]
        lo,hi=scale['ylim']
        visible=[v for v in ticks if lo<=v<=hi]
        if len(visible)<2:
            visible=sorted(set(visible+[scale['linthresh']*v for v in [-4,-2,-1,0,1,2,4]
                                        if lo<=scale['linthresh']*v<=hi]))
            ticks=visible if len(visible)>=2 else list(np.linspace(lo,hi,3))
            ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda value,_:f'{value:g}'))
        ax.yaxis.set_major_locator(ticker.FixedLocator(ticks))
        ax.yaxis.set_minor_locator(ticker.NullLocator())
        ax.text(.985,.96,'symlog',transform=ax.transAxes,ha='right',va='top',fontsize=7,color='#9c332a',
                bbox={'facecolor':'white','edgecolor':'none','alpha':.85,'pad':1})
    ax.set_ylim(*scale['ylim'])
    ax.set_xlim(0,duration)
    ax.grid(True,alpha=.18,linewidth=.65)
    ax.set_ylabel(label,fontsize=8)
    ax.tick_params(labelsize=7)


def draw_reference(axes, rows):
    times = [r['time'] for r in rows]
    for ax,values in zip(axes,coordinates(rows,True)):
        ax.plot(times,values,color='#747a83',ls='--',lw=1.15,alpha=.8,label='RK4 physical reference',zorder=1)
    E0 = rows[0].get('target_energy',energy(rows[0]))
    axes[0].axhline(E0,color='#1d2228',ls=':',lw=.85,label='Initial target energy',zorder=1)


def draw_trace(axes, summary, rows, method, muted=False, application_label=False):
    times = [r['time'] for r in rows]
    color = '#b9bdc4' if muted else COLORS[method]
    style=STYLES[method]
    label=SHORT[method]
    if application_label:
        immediate=summary['application']=='immediate'
        style='--' if immediate else '-'
        color='#b42483' if immediate else color
        label='GS: immediate positions' if immediate else 'GS: summed positions'
    for ax,values in zip(axes,coordinates(rows)):
        ax.plot(times,values,color=color,ls=style,lw=1.0 if muted else 1.4,
                alpha=.8 if muted else .96,label=label,zorder=2 if muted else 3)
        if summary['completed_steps'] < summary['requested_steps']:
            ax.plot(times[-1],values[-1],marker='x',color=color,ms=7,mew=1.5,zorder=6)
            if not muted:
                ax.axvline(times[-1],color=color,ls=':',lw=.65,alpha=.6)
    return times[-1]


def condition_label(soft, retain):
    energy_label = 'Hard energy constraint' if soft==0 else f'Compliant energy: β/Escale = {soft:g}'
    return f"{energy_label}  ·  {'retain' if retain else 'drop'} rₓ"


def condition_slug(soft, retain):
    return f"{'hard' if soft==0 else 'soft_001'}_{'retain' if retain else 'drop'}_rx"


def make_method(data, method, soft, retain, out, limits):
    fig, axes = plt.subplots(5,3,figsize=(13.5,13.0),squeeze=False)
    labels = ['Total mechanical E [J]','Position x [m]','Position y [m]']
    for row,case in enumerate(CASES):
        summary, trace = data.get(case,method,soft,retain)
        baseline, base_trace = data.get(case,'baseline',0,0)
        draw_reference(axes[row],base_trace)
        if method!='baseline':
            draw_trace(axes[row],baseline,base_trace,'baseline',muted=True)
        final = draw_trace(axes[row],summary,trace,method,application_label=method.startswith('gs_'))
        trace_groups=[trace]
        if method.startswith('gs_'):
            immediate_summary,immediate_trace=data.get(case,method,soft,retain,'immediate')
            draw_trace(axes[row],immediate_summary,immediate_trace,method,application_label=True)
            trace_groups.append(immediate_trace)
        duration = summary['h']*summary['requested_steps']
        for col,ax in enumerate(axes[row]):
            setup_axis(ax,panel_scale(data,limits,case,col,trace_groups),duration,labels[col])
            if row==4:ax.set_xlabel('Time [s]',fontsize=9)
        status = f"stopped at {final:.3g} s" if summary['completed_steps']<summary['requested_steps'] else '10 s complete'
        max_error = summary['max_abs_energy_error']/summary['Escale']
        axes[row,0].set_title(CASE_LABELS[case],loc='left',fontsize=10,pad=7,weight='bold')
        axes[row,1].set_title(f"max |E − E*|/Escale = {max_error:.3g}",loc='left',fontsize=8,pad=7)
        if method.startswith('gs_'):
            immediate_error=immediate_summary['max_abs_energy_error']/immediate_summary['Escale']
            axes[row,1].set_title(f"max |E − E*|/Escale: summed {max_error:.3g}, immediate {immediate_error:.3g}",loc='left',fontsize=7,pad=7)
        issues = int(summary['incompatible_updates'])
        text = status + (f' · {issues} unsatisfied singular updates' if issues else '')
        if method.startswith('gs_'):
            immediate_time=immediate_summary['completed_steps']*data.h
            immediate_issues=int(immediate_summary['incompatible_updates'])
            text=f'summed {final:g} s; immediate {immediate_time:g} s'
            if issues or immediate_issues:text+=f' · singular failures {issues}/{immediate_issues}'
        axes[row,2].set_title(text,loc='left',fontsize=8,pad=7,color='#963a2b' if issues or final<duration else '#333333')
    handles, labels = axes[0,0].get_legend_handles_labels()
    fig.legend(handles,labels,loc='upper center',bbox_to_anchor=(.5,.958),ncol=5,frameon=False,fontsize=8)
    detail = 'Material λ only' if method=='baseline' else condition_label(soft,retain)
    fig.suptitle(LABELS[method]+'\n'+detail,fontsize=15,y=.995)
    fig.text(.5,.014,
        f'h = 1/60 s; {data.sweeps} outer sweeps; {data.inner} frozen pass per sweep. '
        'Curvature dropped; gradients refreshed only between outer sweeps.\n'
        'Pendulum gravity potential is 9.81 y, zero at the anchor. Crosses mark the last accepted state before a guard.\n'
        'Bounded panels share axis ranges; runaway panels use marked symlog ranges. E is actual mechanical energy, not its residual.',
        ha='center',fontsize=7.5,color='#454a50')
    fig.tight_layout(rect=(.005,.065,.995,.935),h_pad=1.8,w_pad=1.5)
    slug = 'baseline' if method=='baseline' else method+'_'+condition_slug(soft,retain)
    fig.savefig(out/(slug+'.png'),dpi=170,facecolor='white')
    plt.close(fig)
    return slug


def make_overview(data, soft, retain, out, limits, cases=CASES, suffix=''):
    nrows = len(cases)
    fig,axes=plt.subplots(nrows,3,figsize=(13.5,2.45*nrows+1.35),squeeze=False)
    for row,case in enumerate(cases):
        base,base_trace=data.get(case,'baseline',0,0)
        draw_reference(axes[row],base_trace)
        for method in data.methods:
            summary,trace=data.get(case,method,soft,retain)
            draw_trace(axes[row],summary,trace,method)
        for col,ax in enumerate(axes[row]):
            setup_axis(ax,limits[case,col],base['requested_steps']*data.h,['Total mechanical E [J]','Position x [m]','Position y [m]'][col])
            if row==nrows-1:ax.set_xlabel('Time [s]',fontsize=9)
        axes[row,0].set_title(CASE_LABELS[case],loc='left',fontsize=10,pad=7,weight='bold')
    handles,labels=axes[0,0].get_legend_handles_labels()
    fig.legend(handles,labels,loc='upper center',bbox_to_anchor=(.5,.955),ncol=8,frameon=False,fontsize=8)
    fig.suptitle('Frozen linearization: all solve schedules\n'+condition_label(soft,retain),fontsize=14,y=.995)
    fig.text(.5,.018,f'h = 1/60 s; {data.sweeps} outer sweeps; one frozen pass/sweep. Overlapping curves can hide identical solutions.\n'
             'Gray dashed: RK4. Dotted horizontal: initial energy. Crosses: last accepted state before a guard. Symlog axes marked.',
             ha='center',fontsize=7.5,color='#454a50')
    fig.tight_layout(rect=(.005,.065,.995,.91 if nrows<4 else .935),h_pad=1.7,w_pad=1.5)
    slug='overview_'+condition_slug(soft,retain)+suffix
    fig.savefig(out/(slug+'.png'),dpi=170,facecolor='white')
    plt.close(fig)
    return slug


def make_case(data, case, soft, retain, out, limits):
    """One row per method, including independent Jacobi when present."""
    nrows=len(data.methods)
    fig,axes=plt.subplots(nrows,3,figsize=(12.8,1.85*nrows+1.7),squeeze=False)
    base,base_trace=data.get(case,'baseline',0,0)
    for row,method in enumerate(data.methods):
        summary,trace=data.get(case,method,soft,retain)
        draw_reference(axes[row],base_trace)
        if method!='baseline':draw_trace(axes[row],base,base_trace,'baseline',muted=True)
        final=draw_trace(axes[row],summary,trace,method)
        for col,ax in enumerate(axes[row]):
            setup_axis(ax,panel_scale(data,limits,case,col,[trace]),base['requested_steps']*data.h,['E [J]','x [m]','y [m]'][col])
            if row==nrows-1:ax.set_xlabel('Time [s]',fontsize=9)
        axes[row,0].set_title(LABELS[method],loc='left',fontsize=9,weight='bold',color=COLORS[method])
        if summary['completed_steps']<summary['requested_steps']:
            axes[row,1].set_title(f'Guard: stopped at {final:.3g} s',loc='left',fontsize=8,color='#963a2b')
    fig.suptitle(CASE_LABELS[case]+'\n'+condition_label(soft,retain)+' · summed positions',fontsize=14,y=.995)
    fig.text(.5,.952,'Each row is a different solve schedule. Gray dashed: RK4. Light gray: ordinary XPBD. Dotted horizontal: initial energy.',ha='center',fontsize=8)
    fig.text(.5,.017,f'h = 1/60 s; {data.sweeps} outer sweeps; {data.inner} frozen pass per sweep. Actual total mechanical energy and both Cartesian coordinates.\n'
             'Bounded panels share axis ranges; runaway panels use marked symlog ranges. Crosses mark the last accepted state before a guard.',ha='center',fontsize=8)
    fig.tight_layout(rect=(.005,.055,.995,.94),h_pad=1.4,w_pad=1.4)
    slug=case+'_'+condition_slug(soft,retain)
    fig.savefig(out/(slug+'.png'),dpi=170,facecolor='white')
    plt.close(fig)
    return slug


def make_pair(data, cases, soft, retain, out, limits, suffix, methods=None):
    """One row per schedule, energy and horizontal position for two scenes."""
    methods=methods or data.methods
    nrows=len(methods)
    comparison=('jacobi' in methods and 'joint' not in methods)
    fig,axes=plt.subplots(nrows,4,figsize=(15.5,1.9*nrows+1.5),squeeze=False)
    for row,method in enumerate(methods):
        for case_index,case in enumerate(cases):
            summary,trace=data.get(case,method,soft,retain)
            base,base_trace=data.get(case,'baseline',0,0)
            times=[r['time'] for r in trace]
            reference_times=[r['time'] for r in base_trace]
            for local_col,value_index in enumerate([0,1]):
                col=case_index*2+local_col
                ax=axes[row,col]
                ax.plot(reference_times,coordinates(base_trace,True)[value_index],color='#747a83',ls='--',lw=1.05,alpha=.8)
                if method!='baseline':
                    ax.plot(reference_times,coordinates(base_trace)[value_index],color='#c3c6cc',lw=.9,alpha=.9)
                vals=coordinates(trace)[value_index]
                ax.plot(times,vals,color=COLORS[method],lw=1.4)
                if value_index==0:ax.axhline(trace[0]['target_energy'],color='#1d2228',ls=':',lw=.75)
                setup_axis(ax,panel_scale(data,limits,case,value_index,[trace]),base['requested_steps']*data.h,'E [J]' if value_index==0 else 'x [m]')
                if summary['completed_steps']<summary['requested_steps']:
                    ax.plot(times[-1],vals[-1],marker='x',color=COLORS[method],ms=7,mew=1.4)
                    ax.axvline(times[-1],color=COLORS[method],ls=':',lw=.65,alpha=.6)
                    if local_col==0:
                        ax.text(.97,.09,f'guard at {times[-1]:.2g} s',ha='right',transform=ax.transAxes,fontsize=7,color='#963a2b')
                if row==0:
                    ax.set_title(CASE_LABELS[case]+('\nEnergy' if value_index==0 else '\nHorizontal position'),fontsize=10,pad=8)
                if row==nrows-1:ax.set_xlabel('Time [s]',fontsize=8)
        axes[row,0].text(-.23,.5,SHORT[method],transform=axes[row,0].transAxes,rotation=90,
                         va='center',ha='center',weight='bold',fontsize=10,color=COLORS[method])
    title='Jacobi versus Gauss–Seidel' if comparison else 'Frozen linearization: all solve choices'
    fig.suptitle(title+'\n'+condition_label(soft,retain)+' · summed positions',fontsize=15,y=.997)
    note='Jacobi uses the old residual for both multiplier solves. GS updates the next residual with the preceding correction.' if 'jacobi' in methods else 'GS uses the preceding correction when solving the next multiplier. L = material λ; E = energy η.'
    fig.text(.5,.92 if nrows<=4 else .947,note,ha='center',fontsize=9)
    fig.text(.5,.018,'Gray dashed: RK4 physical reference. Light gray: ordinary XPBD. Dotted horizontal: initial energy.\n'
             f'h = 1/60 s; {data.sweeps} outer sweeps, one frozen pass per sweep. All curves show accepted states; crosses mark guards.\n'
             'Actual signed mechanical energy. Position is Cartesian x(t); individual figures also show y(t). Bounded panels share ranges; runaway panels use marked symlog.',
             ha='center',fontsize=8)
    fig.tight_layout(rect=(.025,.065,.995,.935),h_pad=1.15,w_pad=1.25)
    slug=('jacobi_vs_gs_' if comparison else 'overview_')+condition_slug(soft,retain)+'_'+suffix
    fig.savefig(out/(slug+'.png'),dpi=170,facecolor='white')
    plt.close(fig)
    return slug


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('results',type=Path)
    parser.add_argument('--output',type=Path,help='Default: RESULTS/plots')
    parser.add_argument('--h',type=float,default=1/60)
    parser.add_argument('--sweeps',type=int,default=5)
    parser.add_argument('--inner',type=int,default=1)
    args=parser.parse_args()
    if args.output is None:args.output=args.results/'plots'
    args.output.mkdir(parents=True,exist_ok=True)
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.spines.top':False,'axes.spines.right':False,
                         'axes.labelcolor':'#353a40','axes.edgecolor':'#c0c3c8','text.color':'#242b32','savefig.dpi':170})
    data=Data(args.results,args.h,args.sweeps,args.inner)
    figures=[]
    for soft in [0,.01]:
        for retain in [0,1]:
            limits=scales(data,soft,retain)
            for method in data.methods:
                if method=='baseline' and (soft or retain):continue
                figures.append(make_method(data,method,soft,retain,args.output,limits))
            if soft==0 and retain==1:
                figures.append(make_pair(data,['spring_rest','pendulum_rest'],soft,retain,args.output,limits,'rest'))
                figures.append(make_pair(data,['spring_tangent','pendulum_tangent'],soft,retain,args.output,limits,'moving'))
                for case in CASES:
                    figures.append(make_case(data,case,soft,retain,args.output,limits))
            if 'jacobi' in data.methods:
                comparison=['baseline','jacobi','gs_le','gs_el']
                figures.append(make_pair(data,['spring_rest','pendulum_rest'],soft,retain,args.output,limits,'rest',comparison))
                figures.append(make_pair(data,['spring_tangent','pendulum_tangent'],soft,retain,args.output,limits,'moving',comparison))
    audit={'trace_groups':len(data.traces),'summary_groups':len(data.summaries),
           'max_scaled_energy_reconstruction_difference':max(data.energy_audit),'figures':figures,
           'formula':'E = 0.5 |v|^2 + 50 (|x|-1)^2 for spring; E = 0.5 |v|^2 + 9.81 y for pendulum',
           'all_panels':'Actual energy, Cartesian x, Cartesian y; time is seconds',
           'scales':'Shared across bounded schedules per case and condition; marked symlog when energy error exceeds 10 Escale or a coordinate exceeds 3 times reference length scale',
           'methods':data.methods,
           'schedule_definitions':{'jacobi':'Each multiplier correction uses the same old frozen residual; sum the two corrections once.',
                                   'gauss_seidel':'Each following multiplier solve uses the preceding correction in its residual, even when positions are summed at the end.'},
           'condition_values':{'softness':[0,.01],'retain_rx':[0,1],'h':args.h,'sweeps':args.sweeps,'inner':args.inner}}
    (args.output/'plot_audit.json').write_text(json.dumps(audit,indent=2)+'\n')
    print(json.dumps({k:v for k,v in audit.items() if k!='figures'},indent=2))
    print(f'Wrote {len(figures)} figures (PNG) to {args.output}')


if __name__=='__main__':main()
