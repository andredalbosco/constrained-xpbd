#!/usr/bin/env python3
"""Compare spring work-midpoint with ordinary midpoint and baseline XPBD."""
import argparse
import collections
import csv
import json
import math
import os
from pathlib import Path

os.environ.setdefault('MPLCONFIGDIR','/tmp/xpbd-work-midpoint-matplotlib')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

CASES=['spring_rest','spring_radial','spring_tangent']
TITLES={'spring_rest':'Spring released from rest','spring_radial':'Spring with radial velocity',
        'spring_tangent':'Spring with tangential velocity'}
METHODS=['standard_xpbd','midpoint_xpbd','work_midpoint']
LABELS={'standard_xpbd':'Standard XPBD · 1 sweep','midpoint_xpbd':'Midpoint · 1 sweep',
        'work_midpoint':'Work midpoint · 5 sweeps'}
COLORS={'standard_xpbd':'#747b83','midpoint_xpbd':'#1769aa','work_midpoint':'#138765'}
TEXT={'case','method','phase_kind','status'}
FLOOR=1e-15


def csv_rows(path):
    with path.open() as stream:
        for row in csv.DictReader(stream):
            for key in row:
                if key not in TEXT:
                    row[key]=float(row[key])
                    if not math.isfinite(row[key]):raise ValueError(f'Nonfinite field {key} in {path}')
            row.setdefault('duration',10.)
            yield row


def key(row):return tuple(row[k] for k in ['case','method','h','sweeps','duration'])


def load(root,h,selected,summary_lookup):
    summaries={key(row):row for row in csv_rows(root/'summary.csv')}
    counts=collections.Counter()
    max_error=collections.defaultdict(float)
    energy_difference=0.
    angular_difference=0.
    for row in csv_rows(root/'traces.csv'):
        settings=key(row)
        summary=summaries[settings]
        E0=summary['E0']
        energy=.5*(row['vx']**2+row['vy']**2)+50*(math.hypot(row['x'],row['y'])-1)**2
        e_difference=abs(energy-row['energy'])/max(1,abs(row['energy']))
        l_difference=abs(row['x']*row['vy']-row['y']*row['vx']-row['angular_momentum'])
        energy_difference=max(energy_difference,e_difference)
        angular_difference=max(angular_difference,l_difference)
        if e_difference>2e-12 or l_difference>2e-12:raise ValueError('Physical energy/angular reconstruction failed')
        if not math.isclose(row['time'],counts[settings]*row['h'],rel_tol=1e-12,abs_tol=1e-12):
            raise ValueError('Nonsequential trajectory time')
        counts[settings]+=1
        max_error[settings]=max(max_error[settings],abs(energy-E0)/E0)
        sweeps=5 if row['method']=='work_midpoint' else 1
        if row['method'] in METHODS and row['case'] in CASES and math.isclose(row['h'],h,rel_tol=1e-12) \
                and row['sweeps']==sweeps and row['duration']==10:
            selected[row['case'],row['method']].append(row)
            summary_lookup[row['case'],row['method']]=summary
    for settings,count in counts.items():
        summary=summaries[settings]
        expected=summary.get('completed_steps',summary.get('steps'))
        if count!=expected+1:raise ValueError('Trajectory length does not match summary')
        if not math.isclose(max_error[settings],summary['max_relative_energy_error'],rel_tol=1e-8,abs_tol=1e-12):
            raise ValueError('Summary energy error does not match trajectory')
    return {'root':str(root),'validated_trajectory_groups':len(counts),
            'max_scaled_energy_reconstruction_difference':energy_difference,
            'max_angular_reconstruction_difference':angular_difference}


def decorate(ax,ylabel,log=False):
    ax.set_ylabel(ylabel,fontsize=9)
    ax.set_xlabel('Time [s]',fontsize=9)
    ax.set_xlim(0,10)
    if log:ax.set_yscale('log');ax.set_ylim(bottom=FLOOR/2)
    ax.grid(True,alpha=.18,linewidth=.7)
    ax.tick_params(labelsize=8)


def series(rows,summary,kind):
    if kind=='energy':return [max(abs(r['energy']-summary['E0'])/summary['E0'],FLOOR) for r in rows]
    if kind=='angular':return [max(abs(r['angular_momentum']-summary['L0']),FLOOR) for r in rows]
    return [r['x'] for r in rows]


def trace(ax,rows,summary,kind,method):
    ax.plot([r['time'] for r in rows],series(rows,summary,kind),color=COLORS[method],
            lw=1.3,label=LABELS[method])


def footer(fig,h):
    fig.text(.5,.016,f'h = {h:g} s, full ten-second trajectories. Energy uses stored endpoint velocity. Different sweep counts are explicit; this is not an equal-cost comparison.\n'
             'Energy and angular error axes are logarithmic; values below 10⁻¹⁵ are drawn at 10⁻¹⁵. No trajectory is clipped or rescaled.',
             ha='center',fontsize=8)


def overview(selected,summaries,out,h):
    fig,axes=plt.subplots(3,3,figsize=(14,9),squeeze=False)
    for index,case in enumerate(CASES):
        ref=selected[case,'midpoint_xpbd']
        axes[index,1].plot([r['time'] for r in ref],[r['ref_x'] for r in ref],color='#252c33',ls='--',lw=1.0,label='RK4 reference')
        for method in METHODS:
            rows=selected[case,method]
            summary=summaries[case,method]
            for col,kind in enumerate(['energy','x','angular']):trace(axes[index,col],rows,summary,kind,method)
        for col,label in enumerate(['|E − E₀| / E₀','Horizontal position x [m]','|L − L₀| [kg m²/s]']):
            decorate(axes[index,col],label,log=col!=1)
        axes[index,0].set_title(TITLES[case],loc='left',weight='bold',fontsize=11,pad=9)
        work=summaries[case,'work_midpoint']
        axes[index,1].set_title(f'Work midpoint RMS position error: {work["rms_position_error"]:.4g} m',loc='left',fontsize=9,pad=9)
        axes[index,2].set_title('Angular momentum about the fixed anchor',loc='left',fontsize=9,pad=9)
    handles,labels=axes[0,1].get_legend_handles_labels()
    fig.legend(handles,labels,loc='upper center',bbox_to_anchor=(.5,.944),ncol=4,frameon=False,fontsize=9)
    fig.suptitle('Spring: energy, motion, and angular momentum\nStandard XPBD, midpoint, and work midpoint',fontsize=15,y=.995)
    footer(fig,h)
    fig.tight_layout(rect=(.005,.07,.995,.91),h_pad=2,w_pad=1.6)
    name=f'work_midpoint_spring_h{round(1/h)}.png'
    fig.savefig(out/name,dpi=170,facecolor='white');plt.close(fig)
    return name


def tangent(selected,summaries,out,h):
    case='spring_tangent'
    fig,axes=plt.subplots(2,2,figsize=(12,9),squeeze=False)
    ref=selected[case,'midpoint_xpbd']
    axes[0,1].plot([r['time'] for r in ref],[r['ref_x'] for r in ref],color='#252c33',ls='--',lw=1.0,label='RK4 reference')
    axes[1,0].plot([r['ref_x'] for r in ref],[r['ref_y'] for r in ref],color='#252c33',ls='--',lw=1.0,label='RK4 reference')
    for method in METHODS:
        rows=selected[case,method];summary=summaries[case,method]
        trace(axes[0,0],rows,summary,'energy',method)
        trace(axes[0,1],rows,summary,'x',method)
        trace(axes[1,1],rows,summary,'angular',method)
        axes[1,0].plot([r['x'] for r in rows],[r['y'] for r in rows],color=COLORS[method],lw=1.2,label=LABELS[method],alpha=.9)
    decorate(axes[0,0],'|E − E₀| / E₀',True)
    decorate(axes[0,1],'Horizontal position x [m]')
    decorate(axes[1,1],'|L − L₀| [kg m²/s]',True)
    axes[1,0].plot(ref[0]['x'],ref[0]['y'],'o',color='#252c33',ms=4)
    axes[1,0].set_xlabel('x [m]');axes[1,0].set_ylabel('y [m]')
    axes[1,0].set_aspect('equal',adjustable='datalim');axes[1,0].grid(True,alpha=.18)
    axes[1,0].tick_params(labelsize=8)
    for ax,title in zip(axes.flat,['Mechanical energy error','Endpoint position','2D path; black dot marks the start','Angular momentum drift']):
        ax.set_title(title,loc='left',fontsize=11,pad=10,weight='bold')
    handles,labels=axes[0,1].get_legend_handles_labels()
    fig.legend(handles,labels,loc='upper center',bbox_to_anchor=(.5,.944),ncol=4,frameon=False,fontsize=9)
    fig.suptitle('2D spring with initial tangential velocity\nWork balance correction versus ordinary midpoint',fontsize=15,y=.995)
    footer(fig,h)
    fig.tight_layout(rect=(.005,.07,.995,.91),h_pad=2,w_pad=1.6)
    name=f'work_midpoint_tangent_h{round(1/h)}.png'
    fig.savefig(out/name,dpi=170,facecolor='white');plt.close(fig)
    return name


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('midpoint_results',type=Path)
    parser.add_argument('work_results',type=Path)
    parser.add_argument('--output',type=Path,help='Default: WORK_RESULTS/plots')
    parser.add_argument('--h',type=float,default=1/60)
    args=parser.parse_args()
    if args.output is None:args.output=args.work_results/'plots'
    args.output.mkdir(parents=True,exist_ok=True)
    selected=collections.defaultdict(list);summaries={}
    audits=[load(args.midpoint_results,args.h,selected,summaries),load(args.work_results,args.h,selected,summaries)]
    if len(selected)!=9:raise ValueError(f'Expected nine representative traces, found {len(selected)}')
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.spines.top':False,'axes.spines.right':False,
                         'axes.edgecolor':'#c0c4c8','text.color':'#252c33','axes.labelcolor':'#333d46'})
    figures=[overview(selected,summaries,args.output,args.h),tangent(selected,summaries,args.output,args.h)]
    metrics=[{k:s[k] for k in ['case','method','h','sweeps','max_relative_energy_error','max_angular_error',
                            'rms_position_error','final_phase_error']} for s in summaries.values()]
    audit={'source_validation':audits,'selected_duration':10,'figures':figures,'selected_metrics':metrics,
           'counts':'Baseline and midpoint use one sweep; work midpoint uses five sweeps.',
           'error_display_floor':FLOOR,'velocity':'Stored endpoint velocity, not chord velocity.'}
    (args.output/'work_plot_audit.json').write_text(json.dumps(audit,indent=2)+'\n')
    print(json.dumps(audit,indent=2))


if __name__=='__main__':main()
