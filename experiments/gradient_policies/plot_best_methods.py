#!/usr/bin/env python3
"""Compare useful endpoint-energy solver candidates against motion references."""
import argparse
import csv
import json
import math
import os
from collections import defaultdict
from pathlib import Path

os.environ.setdefault('MPLCONFIGDIR', '/tmp/xpbd-full-results-mpl')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

CASES = ['spring_rest', 'spring_radial', 'spring_tangent', 'pendulum_rest', 'pendulum_tangent']
TITLES = ['Spring released from rest', 'Spring with radial velocity', 'Spring with tangential velocity',
          'Pendulum released from rest', 'Pendulum with initial velocity']
BASE = ('baseline', 'none', -1., 0)
EL4 = ('EL', 'fresh', 4., 0)
LE025 = ('LE', 'fresh', .25, 0)
LE1 = ('LE', 'fresh', 1., 0)
EL1 = ('EL', 'fresh', 1., 0)
ELFULL = ('EL', 'fresh', -1., 1)
ELFROZEN = ('EL', 'frozen_step', 1., 0)
LEFROZEN = ('LE', 'frozen_step', 1., 0)
OLDSTEP = ('EL', 'previous_step', 1., 1)
SELECTION = {
    'spring_rest': [EL4, LE025, ELFROZEN, BASE],
    'spring_radial': [EL4, LE025, ELFROZEN, BASE],
    'spring_tangent': [EL4, LE1, OLDSTEP, BASE],
    'pendulum_rest': [LE1, ELFULL, ELFROZEN, BASE],
    'pendulum_tangent': [LE1, EL1, LEFROZEN, BASE],
}
NOTES = [
    'EL / current / gamma 4 has the lowest position RMS among all 45 tested configurations at these settings.',
    'EL / current / gamma 4 also leads here; LE / current / gamma 0.25 is a useful alternative.',
    'No satisfactory orbit: tiny energy error coexists with large position and angular-momentum errors.',
    'LE / current / gamma 1 balances energy and motion; unbounded EL with retained r_x tracks position most closely.',
    'LE / current / gamma 1 has much smaller energy error; frozen-gradient LE tracks position more closely.',
]
COLORS = ['#0077b6', '#e67c1d', '#984ea3', '#9ba4ae']


def method(row):
    return row['order'], row['gradient'], float(row['gamma']), int(row['retain_rx'])


def label(m):
    order, gradient, gamma, retain = m
    if order == 'baseline':
        return 'Ordinary XPBD'
    g = {'fresh': 'current', 'frozen_step': 'predictor', 'lag_sweep': 'previous sweep',
         'previous_step': 'previous step'}[gradient]
    bound = 'no cutoff' if gamma < 0 else f'gamma {gamma:g}'
    return f'{order} / {g} / {bound} / {"retain" if retain else "drop"} r_x'


def load(root):
    summary = list(csv.DictReader((root / 'summary.csv').open()))
    lookup = {}
    for r in summary:
        if math.isclose(float(r['h']), 1/60) and r['sweeps'] == '5':
            lookup[(r['case'], method(r))] = r
    traces = defaultdict(list)
    for r in csv.DictReader((root / 'traces.csv').open()):
        assert math.isclose(float(r['h']), 1/60) and r['sweeps'] == '5'
        traces[(r['case'], method(r))].append(r)
    return summary, lookup, traces


def numbers(trace, name):
    return np.array([float(row[name]) for row in trace])


def metrics(r):
    return 100 * float(r['max_abs_energy_error']) / float(r['Escale']), float(r['rms_position_error'])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('results', type=Path)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    out = args.output or args.results / 'full_results'
    out.mkdir(parents=True, exist_ok=True)
    rows, lookup, traces = load(args.results)
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 9,
                         'axes.spines.top': False, 'axes.spines.right': False,
                         'axes.edgecolor': '#bdc7cf', 'text.color': '#182e40',
                         'axes.labelcolor': '#263f50', 'axes.titlecolor': '#182e40'})
    fig = plt.figure(figsize=(23, 23), facecolor='white')
    grid = fig.add_gridspec(5, 4, left=.055, right=.985, top=.83, bottom=.075,
                           hspace=.82, wspace=.25)
    fig.text(.055, .971, 'Which constrained-XPBD methods work best?', fontsize=25, weight='bold')
    fig.text(.055, .945, 'Energy and motion tell different stories', fontsize=17, color='#526c7c')
    fig.text(.055, .918, 'Original endpoint-energy formulation  |  60 Hz  |  5 Gauss-Seidel sweeps  |  10 seconds  |  Reference: refined RK4', fontsize=12)
    fig.text(.055, .895, 'Representative candidates below; the companion scorecards include all 45 configurations at all six timestep / sweep settings.', fontsize=11)
    selected_rows = []
    for row_index, case in enumerate(CASES):
        axes = [fig.add_subplot(grid[row_index, col]) for col in range(4)]
        baseline = traces[(case, BASE)]
        t = numbers(baseline, 'time')
        es = float(baseline[0]['Escale'])
        target = float(baseline[0]['target_energy'])
        shift = es - target
        reference = [(numbers(baseline, 'ref_energy') + shift) / es,
                     numbers(baseline, 'ref_x'), numbers(baseline, 'ref_y'), np.zeros_like(t)]
        for ax, values in zip(axes, reference):
            ax.plot(t, values, color='#1a2933', ls=(0, (3, 2)), lw=1.35, zorder=5)
        handles, legends = [], []
        for index, m in enumerate(SELECTION[case]):
            r, trace = lookup[(case, m)], traces[(case, m)]
            assert r['completed_steps'] == r['requested_steps'], 'Selected a stopped run'
            assert len(trace) == 601
            ti = numbers(trace, 'time')
            energy = (numbers(trace, 'energy') + shift) / es
            x, y = numbers(trace, 'x'), numbers(trace, 'y')
            distance = np.hypot(x - numbers(trace, 'ref_x'), y - numbers(trace, 'ref_y'))
            e, rms = metrics(r)
            assert math.isclose(float(np.sqrt(np.mean(distance[1:]**2))), rms, rel_tol=2e-12)
            assert math.isclose(float(np.max(np.abs(energy-1)))*100, e, rel_tol=2e-7, abs_tol=2e-12)
            for col, values in enumerate([energy, x, y, distance]):
                line, = axes[col].plot(ti, values, color=COLORS[index], lw=1.45,
                                      ls='--' if index == 2 else '-', alpha=.93,
                                      zorder=2 if m == BASE else 3)
                if col == 0: handles.append(line)
            legends.append(f'{label(m)}\nmax energy error {e:.3g}% | position RMS {rms:.4f} m')
            selected_rows.append({'case': case, 'method': label(m), 'max_energy_error_percent': e,
                                  'position_rms_m': rms, 'path_ratio': float(r['path_ratio'])})
        pos = axes[0].get_position()
        fig.text(.055, pos.y1 + .047, TITLES[row_index], fontsize=16, weight='bold')
        fig.text(.055, pos.y1 + .032, NOTES[row_index], fontsize=10.2, color='#526c7c')
        fig.legend(handles, legends, loc='lower left', bbox_to_anchor=(.052, pos.y1 + .004),
                   bbox_transform=fig.transFigure, frameon=False, ncol=4, fontsize=8.4,
                   handlelength=2.8, columnspacing=2.1, borderaxespad=0)
        for col, ax in enumerate(axes):
            ax.set_xlim(0, 10)
            ax.grid(alpha=.18, color='#7390a4')
            if row_index == 4:
                ax.set_xlabel('Time [s]')
            ax.set_ylabel(['Energy above minimum / initial', 'Position x [m]', 'Position y [m]',
                           'Distance to reference [m]'][col])
        axes[0].set_ylim(-.04, 1.34 if case == 'spring_radial' else 1.10)
        if case in ('spring_rest', 'spring_radial'):
            axes[2].set_ylim(-.025, .025)
        axes[3].set_ylim(bottom=0)
    fig.text(.055, .04, 'Dashed black: continuous reference. Energy ratio 1 is the target. For the pendulum, energy is measured above the lowest gravitational potential.', fontsize=10)
    fig.text(.055, .025, 'EL = energy then material; LE = material then energy, within the same frozen linearization. gamma bounds the raw energy displacement relative to h|v| + h²|f|/m.', fontsize=10)
    fig.text(.055, .010, 'Candidates illustrate tradeoffs, not one universal winner. All shown trajectories complete 10 s; failed and frozen alternatives remain visible in the exhaustive scorecards.', fontsize=10)
    fig.savefig(out / 'best_methods_trajectories.png', dpi=160, facecolor='white')
    fig.savefig(out / 'best_methods_trajectories.pdf', facecolor='white')
    plt.close(fig)
    with (out / 'selected_metrics.csv').open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=selected_rows[0].keys())
        writer.writeheader(); writer.writerows(selected_rows)
    audit = {'input_runs': len(rows), 'selected_trajectory_entries': len(selected_rows),
             'trace_samples_per_entry': 601, 'energy_and_position_metrics_reconstructed': True,
             'h': 1/60, 'sweeps': 5, 'duration_s': 10,
             'figures': ['best_methods_trajectories.png', 'best_methods_trajectories.pdf']}
    (out / 'trajectory_plot_audit.json').write_text(json.dumps(audit, indent=2)+'\n')
    print(json.dumps(audit, indent=2))


if __name__ == '__main__':
    main()
