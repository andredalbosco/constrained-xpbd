#!/usr/bin/env python3
"""Validate and summarize the direct-energy experiment; plot accepted trajectories."""
import argparse
import collections
import csv
import json
import math
import os
from pathlib import Path

CASES = ['spring_rest', 'spring_radial', 'spring_tangent', 'pendulum_rest', 'pendulum_tangent']
METHODS = ['baseline', 'joint', 'immediate', 'frozen', 'reverse', 'symmetric']
COLORS = dict(zip(METHODS, ['#68737b', '#2775b6', '#e48422', '#36966c', '#9656a2', '#cb4b75']))
KEY = ['case', 'method', 'retain_rx', 'inner', 'h', 'sweeps', 'softness']
TEXT = {'case', 'method', 'status'}


def read_csv(path):
    with path.open() as stream:
        rows = list(csv.DictReader(stream))
    if not rows:
        raise ValueError(f'No rows in {path}')
    for row in rows:
        for field in row:
            if field not in TEXT:
                row[field] = float(row[field])
                if not math.isfinite(row[field]):
                    raise ValueError(f'Nonfinite {path.name} field {field}: {row}')
    return rows


def key(row):
    return tuple(row[k] for k in KEY)


def complete(row):
    return row['completed_steps'] == row['requested_steps']


def write_csv(path, rows):
    with path.open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def close(a, b):
    return abs(a-b) < 1e-12


def select(rows, **values):
    return [r for r in rows if all(close(r[k], v) if isinstance(v, (float, int)) else r[k] == v for k, v in values.items())]


def one(rows, **values):
    selected = select(rows, **values)
    if len(selected) != 1:
        raise ValueError(f'Expected one row for {values}, found {len(selected)}')
    return selected[0]


def label(row):
    result = row['method']
    if row['method'] in ['joint', 'frozen']:
        result += ' + R' if row['retain_rx'] else ' - R'
    if row['method'] == 'frozen':
        result += f" (inner {int(row['inner'])})"
    return result


def plot_traces(rows, traces, directory, softness, filename, cases):
    import matplotlib.pyplot as plt
    figure, axes = plt.subplots(3, len(cases), figsize=(5*len(cases), 9.5), squeeze=False)
    for column, case in enumerate(cases):
        reference_drawn = False
        for method in METHODS:
            row = one(rows, case=case, method=method, retain_rx=0, inner=1,
                      h=1/60, sweeps=5, softness=0 if method == 'baseline' else softness)
            trajectory = traces[key(row)]
            times = [r['time'] for r in trajectory]
            energy = [max(abs(r['energy_error'])/r['Escale'], 1e-16)
                      if r['step'] or softness==0 else math.nan for r in trajectory]
            kinetic = [r['kinetic']/r['Escale'] for r in trajectory]
            position = [max(math.hypot(r['x']-r['ref_x'], r['y']-r['ref_y']), 1e-16)
                        if r['step'] else math.nan for r in trajectory]
            legend = method if complete(row) else f"{method}: guard at {row['completed_time']:.3g}s"
            for axis, values in zip(axes[:, column], [energy, kinetic, position]):
                axis.plot(times, values, color=COLORS[method], label=legend, linewidth=1.4)
                if not complete(row):
                    axis.plot(times[-1], values[-1], 'x', color=COLORS[method], markersize=7)
            if not reference_drawn:
                axes[1, column].plot(times, [r['ref_kinetic']/r['Escale'] for r in trajectory],
                                     color='black', linestyle=':', linewidth=1.2, label='RK4 reference')
                reference_drawn = True
        axes[0, column].set_title(case.replace('_', ' '))
        axes[0, column].set_yscale('log')
        axes[1, column].set_yscale('symlog', linthresh=.01)
        axes[2, column].set_yscale('log')
        axes[2, column].set_xlabel('Time [s]')
        for axis in axes[:, column]:
            axis.set_xlim(0, 10)
            axis.grid(True, alpha=.22)
        axes[0, column].legend(fontsize=7, loc='best')
    axes[0, 0].set_ylabel('|Actual energy error| / energy scale')
    axes[1, 0].set_ylabel('Kinetic energy / energy scale')
    axes[2, 0].set_ylabel('Position error versus RK4 [length]')
    figure.suptitle(f'Energy softness beta/Escale = {softness:g}; h=1/60, five sweeps, drop R, frozen inner=1', fontsize=12)
    figure.text(.5, .012, 'Accepted states only; crosses mark guards and lines stop there. Black dotted line: reference kinetic energy. Log values below 1e-16 are floored.',
                ha='center', fontsize=8)
    figure.tight_layout(rect=(0, .035, 1, .96))
    figure.savefig(directory/filename, dpi=170)
    plt.close(figure)


def plot_controls(rows, directory):
    import matplotlib.pyplot as plt
    figure, axes = plt.subplots(2, 3, figsize=(14, 8.2))
    cases = ['spring_rest', 'spring_tangent', 'pendulum_tangent']
    for column, case in enumerate(cases):
        for retain in [0, 1]:
            for softness, marker in [(0, 'o'), (.01, 's'), (1, '^')]:
                series = sorted(select(rows, case=case, method='frozen', h=1/60, sweeps=5,
                                       retain_rx=retain, softness=softness), key=lambda r:r['inner'])
                linestyle = '-' if retain else '--'
                color = {0:'#2775b6', .01:'#e48422', 1:'#36966c'}[softness]
                legend = f"soft={softness:g}, {'retain' if retain else 'drop'} R"
                for index, field in enumerate(['max_relative_energy_error', 'max_rx']):
                    # Exclude partial trajectories from connecting lines. A guard is marked
                    # separately at its accepted-prefix metric, so its time budget is visible.
                    values = [max(r[field], 1e-16) if complete(r) else math.nan for r in series]
                    axes[index, column].plot([r['inner'] for r in series], values,
                                            marker=marker, linestyle=linestyle, color=color, label=legend)
                    for r in series:
                        if not complete(r):
                            axes[index, column].plot(r['inner'], max(r[field], 1e-16), 'x', color=color, markersize=9)
                            axes[index, column].annotate(f"{r['completed_time']:.2g}s", (r['inner'], max(r[field],1e-16)), fontsize=6)
        axes[0, column].set_title(case.replace('_', ' '))
        for axis in axes[:, column]:
            axis.set_xscale('log');axis.set_yscale('log');axis.set_xticks([1,5,20], ['1','5','20'])
            axis.grid(True, alpha=.22)
        axes[1, column].set_xlabel('Inner passes per frozen block')
    axes[0,0].set_ylabel('Max |actual energy error| / energy scale')
    axes[1,0].set_ylabel('Max full position residual |R|/h²')
    axes[0,0].legend(fontsize=7)
    figure.suptitle('Frozen linear-system passes and retained position residual; h=1/60, five outer sweeps', fontsize=12)
    figure.text(.5,.013,'More inner passes increase scalar solve counts; they are not equal-work comparisons. Crosses show guarded accepted prefixes, with completed duration.',ha='center',fontsize=8)
    figure.tight_layout(rect=(0,.035,1,.96));figure.savefig(directory/'residual_controls.png',dpi=170);plt.close(figure)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('results', nargs='?', type=Path, default=Path('results/energy_ordering'))
    parser.add_argument('--output', type=Path, help='Default: RESULTS/analysis')
    parser.add_argument('--no-plots', action='store_true')
    args=parser.parse_args()
    if args.output is None:args.output=args.results/'analysis'
    args.output.mkdir(parents=True,exist_ok=True)
    rows=read_csv(args.results/'summary.csv')
    references=read_csv(args.results/'reference_checks.csv')
    if len({key(r) for r in rows}) != len(rows):
        raise ValueError('Duplicate summary settings')
    required=['max_abs_energy_error','Escale','max_rx','path_length','requested_steps','completed_steps',
              'rank_deficient_blocks','late_samples']
    for row in rows:
        if any(field not in row for field in required):raise ValueError('Missing required summary fields')
        if row['status'].startswith('completed') != complete(row):raise ValueError(f'Status/duration mismatch: {row}')
        row['max_relative_energy_error']=row['max_abs_energy_error']/row['Escale']
        row['completed_time']=row['h']*row['completed_steps']
        row['full_duration']=complete(row)
        row['late_mean_kinetic_available']=row['late_samples']>0
        row['exact_frozen_motion']=complete(row) and row['path_length']<1e-12 and row['max_speed']<1e-12
        row['closed_row_fraction_of_requested']=row['rows_closed_steps']/row['requested_steps']
        row['closed_kkt_fraction_of_requested']=row['kkt_closed_steps']/row['requested_steps']
    traces=collections.defaultdict(list)
    for row in read_csv(args.results/'traces.csv'):
        traces[key(row)].append(row)
    lookup={key(r):r for r in rows}
    angular=[]
    for settings, trajectory in traces.items():
        row=lookup[settings]
        expected=int(row['completed_steps'])+1
        if len(trajectory)!=expected:raise ValueError(f'Trace length mismatch: {settings}')
        if any(int(r['step'])!=i for i,r in enumerate(trajectory)):raise ValueError('Nonsequential trace')
        if not close(trajectory[-1]['time'],row['completed_time']):raise ValueError('Final trace time mismatch')
        max_energy=max(abs(r['energy_error']) for r in trajectory)
        if not math.isclose(max_energy,row['max_abs_energy_error'],rel_tol=1e-12,abs_tol=1e-12):
            raise ValueError('Trace/summary energy mismatch')
        if row['case'].startswith('spring_'):
            # Every harness case has unit mass. The anchored spring force is
            # central, so angular momentum about the fixed anchor is physical.
            momenta=[r['x']*r['vy']-r['y']*r['vx'] for r in trajectory]
            angular.append({**{k:row[k] for k in KEY}, 'completed_time':row['completed_time'],
                            'status':row['status'], 'initial_angular_momentum':momenta[0],
                            'final_angular_momentum':momenta[-1],
                            'max_angular_momentum_change':max(abs(v-momenta[0]) for v in momenta)})
    representative=select(rows,h=1/60,sweeps=5,retain_rx=0,inner=1)
    controls=[r for r in rows if close(r['h'],1/60) and r['sweeps'] in [5,20]
              and r['method'] in ['joint','frozen']]
    write_csv(args.output/'representative.csv',representative)
    write_csv(args.output/'controls.csv',controls)
    write_csv(args.output/'spring_angular_momentum.csv',angular)
    write_csv(args.output/'all_derived.csv',rows)
    summary={
        'runs':len(rows),'statuses':dict(collections.Counter(r['status'] for r in rows)),
        'full_duration_runs':sum(complete(r) for r in rows),'guarded_prefixes':sum(not complete(r) for r in rows),
        'exact_frozen_full_runs':sum(r['exact_frozen_motion'] for r in rows),
        'hard_energy_runs':sum(r['softness']==0 and r['method']!='baseline' for r in rows),
        'runs_with_rank_deficient_blocks':sum(r['rank_deficient_blocks']>0 for r in rows),
        'runs_with_incompatible_local_updates':sum(r['incompatible_updates']>0 for r in rows),
        'hard_pendulum_runs_with_infeasible_endpoint_steps':sum(r['infeasible_steps']>0 for r in rows),
        'runs_without_late_kinetic_samples':sum(not r['late_mean_kinetic_available'] for r in rows),
        'max_reference_refinement_position_error':max(r['max_refinement_position_error'] for r in references),
        'max_reference_energy_error':max(r['max_reference_energy_error'] for r in references),
        'trace_groups':len(traces),'trace_rows':sum(len(v) for v in traces.values()),
        'all_numeric_csv_cells_finite':True,
    }
    (args.output/'audit.json').write_text(json.dumps(summary,indent=2)+'\n')
    if not args.no_plots:
        os.environ.setdefault('MPLCONFIGDIR','/tmp')
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        plt.rcParams.update({'font.size':9,'axes.spines.top':False,'axes.spines.right':False})
        plot_traces(rows,traces,args.output,0,'hard_energy_ordering.png',
                    ['spring_rest','spring_tangent','pendulum_tangent'])
        plot_traces(rows,traces,args.output,.01,'soft_energy_ordering.png',
                    ['spring_rest','pendulum_rest','pendulum_tangent'])
        plot_controls(rows,args.output)
    print(json.dumps(summary,indent=2))
    print('Analysis files:',args.output)


if __name__=='__main__':main()
