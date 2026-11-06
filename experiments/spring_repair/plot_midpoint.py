#!/usr/bin/env python3
"""Plot independently verified spring midpoint versus standard XPBD trajectories."""
import argparse
import collections
import csv
import json
import math
import os
from pathlib import Path

os.environ.setdefault('MPLCONFIGDIR','/tmp/xpbd-midpoint-matplotlib')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

CASES=['spring_rest','spring_radial','spring_tangent']
LABELS={'spring_rest':'Spring released from rest','spring_radial':'Spring with initial radial velocity',
        'spring_tangent':'Spring with initial tangential velocity'}
METHODS=['standard_xpbd','midpoint_xpbd']
COLORS={'standard_xpbd':'#747b83','midpoint_xpbd':'#1769aa'}
NAMES={'standard_xpbd':'Standard XPBD','midpoint_xpbd':'Midpoint XPBD'}
KEY=['case','method','h','sweeps']


def read(path):
    with path.open() as stream:rows=list(csv.DictReader(stream))
    for row in rows:
        for key in row:
            if key not in ['case','method','phase_kind']:
                row[key]=float(row[key])
                if not math.isfinite(row[key]):raise ValueError(f'Nonfinite {key}: {path}')
    return rows


def key(row):return tuple(row[k] for k in KEY)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('results',type=Path)
    parser.add_argument('--output',type=Path,help='Default: RESULTS/plots')
    parser.add_argument('--h',type=float,default=1/60)
    parser.add_argument('--sweeps',type=int,default=1)
    args=parser.parse_args()
    if args.output is None:args.output=args.results/'plots'
    args.output.mkdir(parents=True,exist_ok=True)
    summaries={key(r):r for r in read(args.results/'summary.csv')}
    trajectories=collections.defaultdict(list)
    worst_reconstruction=0.0
    for row in read(args.results/'traces.csv'):
        potential=50*(math.hypot(row['x'],row['y'])-1)**2
        reconstructed=.5*(row['vx']**2+row['vy']**2)+potential
        difference=abs(reconstructed-row['energy'])/max(1,abs(row['energy']))
        worst_reconstruction=max(worst_reconstruction,difference)
        if difference>2e-12:raise ValueError('Endpoint energy reconstruction failed')
        angular=row['x']*row['vy']-row['y']*row['vx']
        if abs(angular-row['angular_momentum'])>1e-12:raise ValueError('Angular momentum reconstruction failed')
        trajectories[key(row)].append(row)
    for settings,rows in trajectories.items():
        summary=summaries[settings]
        if len(rows)!=int(summary['steps'])+1:raise ValueError('Trace length mismatch')
        energy_error=max(abs(r['energy']-summary['E0']) for r in rows)/summary['E0']
        if not math.isclose(energy_error,summary['max_relative_energy_error'],rel_tol=1e-11,abs_tol=1e-13):
            raise ValueError('Trace and summary energy errors disagree')
    selected={settings:rows for settings,rows in trajectories.items()
              if math.isclose(settings[2],args.h,rel_tol=1e-12) and settings[3]==args.sweeps}
    if len(selected)!=6:raise ValueError(f'Expected six representative trajectories, got {len(selected)}')
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.spines.top':False,'axes.spines.right':False,
                         'axes.edgecolor':'#c0c4c8','text.color':'#252c33','axes.labelcolor':'#333d46'})
    fig,axes=plt.subplots(3,3,figsize=(14,9),squeeze=False)
    audit_rows=[]
    for row_index,case in enumerate(CASES):
        midpoint=selected[(case,'midpoint_xpbd',args.h,float(args.sweeps))]
        times=[r['time'] for r in midpoint]
        axes[row_index,0].plot(times,[r['ref_energy'] for r in midpoint],color='#252c33',ls='--',lw=1.1,label='RK4 reference',zorder=1)
        axes[row_index,1].plot(times,[r['ref_x'] for r in midpoint],color='#252c33',ls='--',lw=1.1,label='RK4 reference',zorder=1)
        axes[row_index,2].axhline(0,color='#252c33',ls='--',lw=1.1,label='RK4 reference',zorder=1)
        for method in METHODS:
            rows=selected[(case,method,args.h,float(args.sweeps))]
            summary=summaries[key(rows[0])]
            for col,field in enumerate(['energy','x','phase_error']):
                axes[row_index,col].plot(times,[r[field] for r in rows],color=COLORS[method],
                                         lw=1.35,label=NAMES[method],zorder=2)
            audit_rows.append({k:summary[k] for k in ['case','method','h','sweeps','max_relative_energy_error',
                               'max_angular_error','rms_position_error','final_phase_error','phase_kind']})
        axes[row_index,0].set_title(LABELS[case],loc='left',fontsize=11,weight='bold',pad=9)
        axes[row_index,0].set_ylabel('Total mechanical energy [J]')
        axes[row_index,1].set_ylabel('Horizontal position x [m]')
        axes[row_index,1].set_title('Endpoint position',loc='left',fontsize=9,pad=9)
        kind='Orbital angle' if case=='spring_tangent' else 'Harmonic phase'
        axes[row_index,2].set_ylabel('Phase difference versus RK4 [rad]')
        axes[row_index,2].set_title(f'{kind}; midpoint final drift {midpoint[-1]["phase_error"]:.4g} rad',loc='left',fontsize=9,pad=9)
        for col,ax in enumerate(axes[row_index]):
            ax.grid(True,alpha=.18,linewidth=.7)
            ax.set_xlim(0,10)
            ax.tick_params(labelsize=8)
            if row_index==2:ax.set_xlabel('Time [s]')
    handles,labels=axes[0,0].get_legend_handles_labels()
    fig.legend(handles,labels,loc='upper center',bbox_to_anchor=(.5,.943),ncol=3,frameon=False)
    fig.suptitle(f'Standard XPBD versus midpoint XPBD — spring\nStep h = {args.h:g} s; {args.sweeps} constraint sweep per timestep',fontsize=15,y=.995)
    fig.text(.5,.015,'Energy uses the stored endpoint velocity; midpoint reconstructs vₙ₊₁ = 4(z − xₙ)/h − vₙ. No energy multiplier or velocity rescaling.\n'
             'Harmonic phase: unwrapped atan2(−vᵣ/10, r−1). Orbital phase: unwrapped atan2(y, x). All runs complete the full ten seconds.',
             ha='center',fontsize=8)
    fig.tight_layout(rect=(.005,.07,.995,.91),h_pad=2,w_pad=1.6)
    filename=f'midpoint_spring_h{round(1/args.h)}_sweeps{args.sweeps}.png'
    fig.savefig(args.output/filename,dpi=170,facecolor='white')
    plt.close(fig)
    audit={'all_trajectory_groups':len(trajectories),'max_scaled_endpoint_energy_reconstruction_difference':worst_reconstruction,
           'plotted':audit_rows,'figure':filename,
           'velocity_warning':'Midpoint endpoint velocity differs from the chord (x_new-x_old)/h; energy uses endpoint velocity.'}
    (args.output/'plot_audit.json').write_text(json.dumps(audit,indent=2)+'\n')
    print(json.dumps(audit,indent=2))
    print(args.output/filename)


if __name__=='__main__':main()
