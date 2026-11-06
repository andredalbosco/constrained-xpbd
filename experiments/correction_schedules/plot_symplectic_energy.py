"""Plot and independently audit the symplectic-predictor/energy study."""
import json
from pathlib import Path
import numpy as np
import csv
from collections import defaultdict
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent
OUT = ROOT / 'symplectic_energy_plots'
OUT.mkdir(exist_ok=True)
def read_summary(path):
    with open(path) as f:
        return [{k: v if k in {'case','method','status'} else float(v) for k,v in r.items()} for r in csv.DictReader(f)]
short = read_summary(ROOT / 'symplectic_energy_results/summary.csv')
long = read_summary(ROOT / 'symplectic_energy_long/summary.csv')
traces = defaultdict(lambda: defaultdict(list))
with open(ROOT / 'symplectic_energy_results/traces.csv') as f:
    for row in csv.DictReader(f):
        target = traces[row['case'], row['method']]
        for key,value in row.items():
            if key not in {'case','method'}: target[key].append(float(value))
traces = {key: {col:np.array(v) for col,v in row.items()} for key,row in traces.items()}
def select(frame, case, method):
    return next(r for r in frame if r['case']==case and r['method']==method and r['substeps']==8)

cases = ['spring_release', 'vertical_oscillation', 'small_swing', 'swing_rest',
         'swing_tangent', 'coupled_motion', 'spring_orbit', 'circular_orbit']
labels = ['Released\nspring', 'Vertical\nspring', 'Small\nswing', 'Large\nswing',
          'Tangential\nswing', 'Coupled\nmotion', 'Eccentric\norbit', 'Circular\norbit']
methods = ['symplectic', 'E_half_cap', 'paired_LE_repair_cap']
names = ['Symplectic predictor alone', 'Predictor + bounded energy correction',
         'Warm LE + bounded energy correction']
colors = ['#3567a3', '#19947a', '#df8b33']

# Reconstruct all saved endpoint energies and RMS errors independently.
max_energy_discrepancy = 0.
max_rms_discrepancy = 0.
for (case, method), t in traces.items():
    row = select(short,case,method)
    assert len(t['step']) == int(row['completed']) + 1
    gravity = 0. if case in ['spring_release', 'spring_orbit', 'circular_orbit'] else 9.81
    E = .5 * (t['vx']**2 + t['vy']**2) + 50 * (np.hypot(t['x'], t['y']) - 1)**2 + gravity*t['y']
    discrepancy = float(np.max(np.abs(E-t['energy'])) / row['Escale'])
    max_energy_discrepancy = max(max_energy_discrepancy, discrepancy)
    rms = float(np.sqrt(np.mean(((t['x']-t['ref_x'])**2+(t['y']-t['ref_y'])**2)[1:])))
    max_rms_discrepancy = max(max_rms_discrepancy, abs(rms-row['rms_position']))
assert max_energy_discrepancy < 1e-7
assert max_rms_discrepancy < 1e-9
assert len(short) == 384 and len(long) == 192
for frame in [short, long]:
    assert all(r['status']=='completed' for r in frame if r['method'] in ['symplectic','xpbd'])
    for method, cap in [('E_half_cap',.1),('E_half_cap_small',.01),('E_skip',.1),('E_direction_cap',.1)]:
        assert all(r['max_correction_ratio'] <= cap+1e-10 for r in frame if r['method']==method)

plt.rcParams.update({'font.size': 10, 'axes.spines.top': False, 'axes.spines.right': False,
                     'axes.titleweight': 'bold', 'axes.grid': True, 'grid.alpha': .18})
fig, axs = plt.subplots(2, 2, figsize=(15, 10))
for row, (frame, duration) in enumerate([(short,10),(long,100)]):
    for col, (metric, title, ylabel) in enumerate([
        ('max_energy_percent','Maximum energy error','% of initial energy above minimum potential'),
        ('rms_position','Position error against reference','RMS distance (m)')]):
        ax=axs[row,col]
        for j,(method,color) in enumerate(zip(methods,colors)):
            data=[select(frame,case,method) for case in cases]
            assert all(r['status']=='completed' for r in data)
            ax.bar(np.arange(8)+(j-1)*.24, [r[metric] for r in data], .23, color=color, label=names[j])
        ax.set_yscale('log');ax.set_xticks(np.arange(8),labels)
        ax.set_title(f'{duration} seconds: {title}',loc='left');ax.set_ylabel(ylabel)
fig.suptitle('Bounded energy correction improves energy, but usually worsens motion', fontsize=19, x=.06, ha='left')
handles, legend_labels=axs[0,0].get_legend_handles_labels()
fig.legend(handles,legend_labels,loc='upper left',bbox_to_anchor=(.055,.955),ncol=3,frameon=False)
fig.text(.06,.025,'480 Hz · one linearized correction pass per substep · energy repair fraction ½ · displacement cap 0.1 × predicted motion\nLower bars are better. All displayed runs completed; raw failures and other methods are retained in the CSV files.',fontsize=10)
fig.subplots_adjust(left=.075,right=.985,top=.87,bottom=.13,hspace=.43,wspace=.25)
fig.savefig(OUT/'energy_and_motion.png',dpi=180);plt.close(fig)

fig,axs=plt.subplots(3,2,figsize=(14,12))
for row,case in enumerate(['spring_release','spring_orbit','coupled_motion']):
    summary=select(short,case,'symplectic')
    for method,color,name in zip(methods,colors,names):
        t=traces[case,method]
        E0=t['energy'][0]
        axs[row,0].plot(t['time'],100*(t['energy']-E0)/summary['Escale'],color=color,lw=1.3,label=name)
        axs[row,1].plot(t['time'],t['x'],color=color,lw=1.3,label=name)
    t=traces[case,'symplectic']
    axs[row,0].axhline(0,color='#333333',ls='--',lw=1,label='Reference energy')
    axs[row,1].plot(t['time'],t['ref_x'],color='#222222',ls='--',lw=1.3,label='Refined RK4 reference')
    axs[row,0].set_title(case.replace('_',' ').capitalize()+': energy',loc='left')
    axs[row,1].set_title('Position',loc='left')
    axs[row,0].set_ylabel('Energy error (%)');axs[row,1].set_ylabel('x (m)')
    for ax in axs[row]:ax.set_xlabel('Time (s)')
fig.suptitle('Accurate total energy does not preserve the trajectory',fontsize=19,x=.06,ha='left')
handles,legend_labels=axs[0,1].get_legend_handles_labels()
fig.legend(handles,legend_labels,loc='upper left',bbox_to_anchor=(.055,.958),ncol=2,frameon=False)
fig.subplots_adjust(left=.075,right=.985,top=.86,bottom=.065,hspace=.4,wspace=.22)
fig.savefig(OUT/'energy_and_position_traces.png',dpi=180);plt.close(fig)
audit=dict(short_runs=len(short),long_runs=len(long),trace_groups=len(traces),
           max_scaled_energy_reconstruction_error=max_energy_discrepancy,max_position_RMS_discrepancy=max_rms_discrepancy,
           short_explicit_stops=sum(r['status']!='completed' for r in short),long_explicit_stops=sum(r['status']!='completed' for r in long),
           all_baselines_completed=True,all_recorded_correction_caps_respected=True)
(OUT/'audit.json').write_text(json.dumps(audit,indent=2)+'\n')
print(json.dumps(audit,indent=2))
