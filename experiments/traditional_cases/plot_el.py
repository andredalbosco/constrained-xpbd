#!/usr/bin/env python3
"""Focused EL comparisons from C++ study --el-focus output."""
import argparse
import csv
import json
import math
import os
from collections import defaultdict
from pathlib import Path
os.environ.setdefault('MPLCONFIGDIR', '/tmp/xpbd-traditional-matplotlib')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import SymmetricalLogLocator
import numpy as np
from plot_results import energy as physical_energy

METHODS = ['xpbd', 'EL_raw', 'EL_raw_gate', 'EL_raw_gate_cap', 'EL_repair_cap']
LABELS = {'xpbd': 'Ordinary XPBD', 'EL_raw': 'EL: no clipping',
          'EL_raw_gate': 'EL: your direction clip',
          'EL_raw_gate_cap': 'EL: direction clip + cap',
          'EL_repair_cap': 'EL: bounded energy repair + cap'}
COLORS = {'xpbd': '#7a8795', 'EL_raw': '#c33d45', 'EL_raw_gate': '#dd8a20',
          'EL_raw_gate_cap': '#ab478a', 'EL_repair_cap': '#2478ac'}
CASES = {'spring_release': 'Spring released from rest', 'pendulum_rest': 'Rigid pendulum released from rest',
         'spring_orbit': 'Spring with tangential velocity', 'pendulum_kick': 'Rigid pendulum with tangential velocity'}

def read(path):
    with path.open(newline='') as f:
        for row in csv.DictReader(f):
            yield {k: v if k in {'case', 'method', 'status'} else float(v) for k, v in row.items()}

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('results', type=Path)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    rows = list(read(args.results / 'summary.csv'))
    summaries = {(r['case'], r['method'], int(r['hz'])): r for r in rows}
    models = {r['case']: r for r in read(args.results / 'cases.csv')}
    traces = defaultdict(list)
    for row in read(args.results / 'traces.csv'):
        traces[row['case'], row['method']].append(row)
    audit = {'summary_runs': len(rows), 'trace_groups': len(traces),
             'max_scaled_energy_reconstruction': 0., 'max_rms_reconstruction': 0.}
    for key, tr in traces.items():
        c, m = key
        summary = summaries[c, m, 480]
        assert summary['completed'] == summary['steps'] == len(tr) - 1
        assert all(r['step'] == i and abs(r['time'] - i/480) < 1e-14 for i,r in enumerate(tr))
        for r in tr:
            delta = abs(physical_energy(models[c], r) - r['energy']) / summary['energy_scale']
            audit['max_scaled_energy_reconstruction'] = max(audit['max_scaled_energy_reconstruction'], delta)
            assert delta < 1e-10
        rms = math.sqrt(sum((r['x1']-r['ref_x1'])**2+(r['y1']-r['ref_y1'])**2 for r in tr[1:])/(len(tr)-1))
        audit['max_rms_reconstruction'] = max(audit['max_rms_reconstruction'], abs(rms-summary['rms_position']))
        assert math.isclose(rms, summary['rms_position'], rel_tol=1e-12, abs_tol=1e-12)
        em = max(abs(r['energy']-summary['initial_energy'])/summary['energy_scale'] for r in tr[1:])
        assert math.isclose(em, summary['max_energy_rel'], rel_tol=1e-12, abs_tol=1e-12)
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,
                         'axes.spines.top':False,'axes.spines.right':False})
    for c in CASES:
        fig, axes = plt.subplots(len(METHODS), 2, figsize=(14, 13.7))
        fig.subplots_adjust(top=.875, bottom=.11, hspace=.64, wspace=.22, left=.08, right=.97)
        fig.suptitle(CASES[c] + ': EL variants', x=.08, y=.977, ha='left', fontsize=22, fontweight='bold')
        fig.text(.08, .946, '480 Hz · 10 seconds · one frozen EL pass per substep · multipliers reset every substep', fontsize=11)
        fig.legend(handles=[Line2D([], [], color='#555', lw=2, label='Simulated motion'),
                            Line2D([], [], color='black', ls='--', lw=1.3, label='Reference motion / exact energy target')],
                   loc='upper left', bbox_to_anchor=(.075,.932), frameon=False, ncol=2)
        for row, method in enumerate(METHODS):
            tr = traces[c, method]; summary = summaries[c, method, 480]
            t = np.array([r['time'] for r in tr]); e = np.array([100*(r['energy']-summary['initial_energy'])/summary['energy_scale'] for r in tr])
            left, right = axes[row]
            left.plot(t,e,color=COLORS[method],lw=1.15)
            left.axhline(0,color='black',ls='--',lw=1)
            if np.max(np.abs(e)) > 1000:
                left.set_yscale('symlog',linthresh=1)
                locator=SymmetricalLogLocator(base=10,linthresh=1); locator.set_params(numticks=7)
                left.yaxis.set_major_locator(locator)
                left.set_ylim(float(e.min())-.05*max(1,abs(float(e.min()))),float(e.max())*1.1)
                left.set_ylabel('Energy error (%) · symlog')
            else:left.set_ylabel('Energy error (%)')
            left.set_title(LABELS[method]+f"\nMaximum energy error: {100*summary['max_energy_rel']:.4g}%", loc='left', fontsize=11)
            if c.startswith('spring'):
                offset=models[c]['l1'] if c=='spring_release' else 0
                p = np.array([r['x1'] - offset for r in tr])
                ref = np.array([r['ref_x1'] - offset for r in tr])
                right.set_ylabel('x − rest length (m)' if offset else 'x coordinate (m)')
                if np.max(np.abs(p)) > 3:
                    right.set_yscale('symlog',linthresh=.25)
                    right.set_ylabel('x − rest length (m) · symlog')
            else:
                p = np.array([math.atan2(r['x1'],-r['y1']) for r in tr])
                ref = np.array([math.atan2(r['ref_x1'],-r['ref_y1']) for r in tr])
                right.set_ylabel('Angle from vertical (rad)')
            right.plot(t,p,color=COLORS[method],lw=1.2)
            right.plot(t,ref,color='black',ls='--',lw=.95,alpha=.8)
            right.set_title(f"Position RMS: {100*summary['rms_position']:.4g} cm",loc='left',fontsize=11)
            for ax in (left,right):
                ax.set_xlim(0,10); ax.grid(alpha=.2)
                if row==len(METHODS)-1:ax.set_xlabel('Time (s)')
        fig.text(.08,.026,'Energy error is normalized by initial energy above minimum potential. Each row has its own vertical range; large spikes remain visible.\n'
                 'Direction clip rejects an energy contribution opposed to objective descent. Cap = 10% of estimated motion. Repair limit = 50% of predicted deficit.\n'
                 'Position RMS uses both Cartesian coordinates. Dashed black is the reference; energy target is zero error.',fontsize=9,color='#444')
        fig.savefig(args.output/(c+'_EL_comparison.png'),dpi=160,facecolor='white')
        plt.close(fig)
    (args.output/'audit.json').write_text(json.dumps(audit,indent=2)+'\n')
    print(json.dumps(audit,indent=2))
    print('Figures: '+', '.join(str(p) for p in sorted(args.output.glob('*.png'))))
if __name__ == '__main__':main()
