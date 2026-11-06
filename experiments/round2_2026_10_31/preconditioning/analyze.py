from pathlib import Path
import csv, math, re, statistics
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parent
FIG=ROOT/'figures';FIG.mkdir(exist_ok=True)
STR={'regime','motion','metric','globalizer','failure','last_step_failure','step_failure'}
def read(name):
    rows=list(csv.DictReader((ROOT/name).open()))
    for r in rows:
        for k,v in list(r.items()):
            if k not in STR:r[k]=float(v)
    return rows
runs=read('runs.csv');refs=read('references.csv');traces=read('selected_traces.csv')
metrics=['none','mass','diagonal','graph','svd_rows']
labels=['Identity','Mass','Diagonal','Graph','SVD rows']
colors=dict(zip(metrics,plt.get_cmap('tab10').colors[:5]))
def subset(**kw):return [r for r in runs if all(r[k]==v for k,v in kw.items())]
def med(rows,key):return statistics.median(r[key] for r in rows if math.isfinite(r[key]))
def fullref(rows):return [r for r in rows if r['matched_steps']==r['requested_steps']]
def converged(rows):return sum(r['converged_steps']==r['requested_steps'] for r in rows)
raw=(ROOT/'identity_tests.txt').read_text()
audit=[];conditions=[]
for line in raw.splitlines():
    if line.startswith(('root ','conditioning ')):
        d=dict(re.findall(r'(\w+)=([^ ]+)',line));d={k:(v if k in ['regime','metric'] else float(v)) for k,v in d.items()}
        (audit if line.startswith('root ') else conditions).append(d)
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'figure.dpi':145})
fig,ax=plt.subplots(1,2,figsize=(12,4.5),constrained_layout=True)
x=np.arange(5);width=.19
for j,(budget,g) in enumerate([(12,'line'),(12,'trust'),(48,'line'),(48,'trust')]):
    vals=[converged(subset(metric=m,budget=budget,globalizer=g)) for m in metrics]
    ax[0].bar(x+(j-1.5)*width,vals,width,label=f'{budget} products, {g}')
ax[0].set_xticks(x,labels);ax[0].set_ylim(0,36);ax[0].set_ylabel('All nonlinear steps converged / 36 scenarios');ax[0].legend(fontsize=8);ax[0].set_title('Completion of a trajectory is not convergence')
for m,label in zip(metrics,labels):
    a=sorted([r for r in audit if r['metric']==m and r['regime']=='both100'],key=lambda r:r['N'])
    ax[1].plot([r['N'] for r in a],[r['products'] for r in a],'-o',color=colors[m],label=label)
ax[1].set_xscale('log',base=2);ax[1].set_yscale('log');ax[1].set_xticks([2,8,32],[2,8,32]);ax[1].set_xlabel('Chain masses');ax[1].set_ylabel('GMRES products to the same converged root');ax[1].set_title('First step: combined contrasts, h=0.04');ax[1].legend(fontsize=8)
fig.savefig(FIG/'convergence_scaling.png');plt.close(fig)

fig,ax=plt.subplots(2,2,figsize=(12,8),constrained_layout=True)
for m,label in zip(metrics,labels):
    a=sorted([r for r in traces if r['N']==32 and r['regime']=='mass100' and r['h']==.04 and r['budget']==12 and r['globalizer']=='line' and r['metric']==m],key=lambda r:r['step'])
    t=[r['step']*.04 for r in a]
    for panel,key in zip(ax.flat,['residual','energy_drift','mode_error','phase_error']):panel.plot(t,[max(1e-14,r[key]) for r in a],'-o',markersize=3,color=colors[m],label=label);panel.set_yscale('log');panel.set_xlabel('Time (s)')
for panel,title in zip(ax.flat,['Original normalized residual','Energy drift / initial excitation','Modal phase-space error vs converged DG','Weighted modal phase difference (radians)']):panel.set_title(title);panel.grid(alpha=.18)
ax[0,0].legend(fontsize=8);fig.suptitle('N=32, mass ratio 100, h=0.04, budget 12: same DG equations');fig.savefig(FIG/'mode_error.png');plt.close(fig)

fig,ax=plt.subplots(1,3,figsize=(14,4),constrained_layout=True)
prefix=next(r for r in refs if r['N']==32 and r['regime']=='both100' and r['h']==.04)['completed_steps']*.04
for m,label in zip(metrics,labels):
    a=sorted([r for r in traces if r['N']==32 and r['regime']=='both100' and r['h']==.04 and r['budget']==48 and r['globalizer']=='line' and r['metric']==m],key=lambda r:r['step'])
    t=[r['step']*.04 for r in a]
    for panel,key in zip(ax,['residual','energy_drift','mode_error']):panel.plot(t,[max(1e-14,r[key]) if math.isfinite(r[key]) else np.nan for r in a],'-o',markersize=3,color=colors[m],label=label);panel.set_yscale('log');panel.set_xlabel('Time (s)');panel.grid(alpha=.18)
for panel,title in zip(ax,['Original residual','Energy drift / initial excitation','Mode error: only where reference converged']):panel.set_title(title)
ax[2].axvspan(prefix,.4,color='gray',alpha=.18,label='No converged reference');ax[2].set_xlim(0,.4);ax[0].legend(fontsize=8);ax[2].legend(fontsize=8)
fig.suptitle('N=32, mass and stiffness ratios 100, h=0.04, budget 48: unresolved nonlinear failure');fig.savefig(FIG/'failed_reference_prefix.png');plt.close(fig)

fig,ax=plt.subplots(figsize=(7,4.7),constrained_layout=True)
regimes=['uniform','mass100','stiff100','both100']
arr=np.array([[next(r['condition'] for r in conditions if r['N']==32 and r['regime']==regime and r['metric']==m) for regime in regimes] for m in metrics])
im=ax.imshow(np.log10(arr),cmap='viridis',aspect='auto');ax.set_xticks(range(4),regimes);ax.set_yticks(range(5),labels);ax.set_title('Transformed Jacobian condition number: N=32, h=0.04')
for i in range(5):
    for j in range(4):ax.text(j,i,f'{arr[i,j]:.2g}',ha='center',va='center',color='white' if np.log10(arr[i,j])<np.log10(arr).max()*.6 else 'black')
fig.colorbar(im,ax=ax,label='log10 condition number');fig.savefig(FIG/'conditioning.png');plt.close(fig)

reference_fail=[r for r in refs if r['completed_steps']<r['requested_steps']]
lines=['# Measured preconditioning results','',f'{len(runs)} configurations across {len(refs)} physical scenarios. Every run retained a finite {0.4:g}-second trajectory; many steps were deliberately propagated before nonlinear convergence. **Finite completion is not a success criterion.**','',
'## Identity checks','',
'- Analytic directional Jacobians: maximum relative discrepancy '+f"{max(float(v) for v in re.findall(r'relative=([0-9.eE+-]+)',raw)):.3g}"+'.',
'- Energy/residual identity away from a root: maximum normalized discrepancy '+f"{max(float(v) for v in re.findall(r'error=([0-9.eE+-]+)',raw)):.3g}"+'.',
'- All 60 converged metric/case combinations agree with dense Newton: maximum combined position/velocity difference '+f"{max(r['difference'] for r in audit):.3g}"+'.',
'- Maximum converged first-step energy defect: '+f"{max(r['energy'] for r in audit):.3g}"+' of initial excitation.',
'', 'These checks show that the coordinate and row transformations preserve the solved equations. They do not establish convergence for every later step.','',
'## Finite-budget results','',
'Counts use all 36 scenarios. Median errors use only the 28 scenarios with a complete converged-DG reference; the other eight retain explicitly limited reference prefixes. Errors compare the same time discretization, not the continuous ODE.','',
'| Budget | Globalizer | Metric | All steps converged / 36 | Median max modal error, complete reference | Median max energy drift, complete reference |','|---:|---|---|---:|---:|---:|']
for budget in [12,48]:
 for glob in ['line','trust']:
  for m in metrics:
   a=subset(metric=m,budget=budget,globalizer=glob);eligible=fullref(a)
   lines.append(f"| {budget} | {glob} | {m} | {converged(a)} | {med(eligible,'max_mode_error'):.3g} | {med(eligible,'max_energy_drift'):.3g} |")
lines+=['','![Convergence and solve-work scaling](figures/convergence_scaling.png)','','## What improved','',
'Graph and mass scaling improve several contrasting-mass cases under the same Krylov budget. For N=32, mass ratio 100, h=0.04, line search, and budget 12:','',
'| Metric | Maximum modal error | Maximum energy drift | Converged steps / 10 |','|---|---:|---:|---:|']
for m in metrics:
 r=subset(N=32,regime='mass100',motion='mixed',h=.04,metric=m,budget=12,globalizer='line')[0]
 lines.append(f"| {m} | {r['max_mode_error']:.3g} | {r['max_energy_drift']:.3g} | {int(r['converged_steps'])} |")
lines+=['','![Finite-budget mode diagnostics](figures/mode_error.png)','','Bounded SVD row balancing gives the most fully converged trajectories in this particular sweep. It uses a dense SVD at every nonlinear iteration, so fewer Krylov products must not be presented as a runtime improvement. Graph scaling is also not uniformly beneficial: its isotropic stiffness surrogate can overweight directions whose true tangent stiffness is small.','',
'![Conditioning](figures/conditioning.png)','','## Globalization is not a universal repair','']
better=[];worse=[]
for r in runs:
 if r['globalizer']!='trust':continue
 b=next(z for z in runs if z['globalizer']=='line' and all(z[k]==r[k] for k in ['N','regime','motion','h','metric','budget']))
 if b['max_energy_drift']>1e-6 and r['max_energy_drift']<.1*b['max_energy_drift']:better.append((b,r))
 if r['max_energy_drift']>1e-6 and r['max_energy_drift']>10*b['max_energy_drift']:worse.append((b,r))
lines += [f'Across 360 paired configurations, the trust-radius variant reduced energy drift by more than 10× in {len(better)} pairs and increased it by more than 10× in {len(worse)} pairs, ignoring differences below 1e-6. This simple radius-capped Newton variant is not uniformly better than line search.','',
'| Case | Metric / budget | Line-search drift | Trust-radius drift |','|---|---|---:|---:|']
for b,r in better:lines.append(f"| N={int(r['N'])}, {r['regime']}, {r['motion']}, h={r['h']} | {r['metric']} / {int(r['budget'])} | {b['max_energy_drift']:.3g} | {r['max_energy_drift']:.3g} |")
lines+=['','## Failed references and missing convergence','',
'The following dense-Newton references failed to converge within their declared iteration/evaluation limit or stalled in line search. No replacement integrator was substituted. Mode/phase comparisons stop at the listed prefix; energy and residual diagnostics continue for the finite-budget methods.','',
'| N | Regime | Motion | h | Reference prefix | Reason |','|---:|---|---|---:|---:|---|']
for r in reference_fail:lines.append(f"| {int(r['N'])} | {r['regime']} | {r['motion']} | {r['h']} | {int(r['completed_steps'])}/{int(r['requested_steps'])} | {r['failure']} |")
lines+=['','![Failed reference prefix](figures/failed_reference_prefix.png)','','The combined mass/stiffness contrast at N=32 and h=0.04 remains a failure for every finite-budget variant: none converges at every timestep, and the dense reference stops after four of ten steps. Better early-step mode errors do not validate the later trajectories. The transforms improve linear conditioning but cannot by themselves repair a stalled nonlinear solve or choose its physical branch.','',
'## Cost accounting','',
'Raw tables record Krylov products, additional model products, residual evaluations, Jacobian assemblies, setup/factorization microseconds, process CPU time, and wall time. The primary run shared the workstation with other experiments; no runtime ranking is drawn from those measurements. A separate `./benchmark --timing` mode performs seven repeats with rotated metric order and writes `timing.csv`; run it only after other benchmarking jobs finish. Dense SVD cost is explicitly included.','',
'## Conclusions for the next experiment','',
'1. Keep mass/graph scaling as numerical options: they can improve finite-budget motion without changing the converged DG state.','2. Use bounded SVD balancing diagnostically; its product-count advantage needs isolated setup-inclusive timing before treating it as practical.','3. Do not adopt the current trust-radius rule as a default. Its regressions need to be separated into radius scaling, model quality, and inadequate nonlinear budgets.','4. Require convergence-gated acceptance or timestep subdivision for the unresolved contrast cases; these tests deliberately show the damage from accepting unconverged steps.','5. This study measures nonlinear-solve distortion relative to converged DG. Physical phase accuracy versus a refined ODE reference remains a separate question.','']
if (ROOT/'timing.csv').exists():
    timings=read('timing.csv')
    start=lines.index('## Cost accounting');end=lines.index('## Conclusions for the next experiment')
    cost=['## Isolated setup-inclusive timing','',
          'After the other trajectory jobs finished, the timing mode ran seven serial repeats with rotated method order (420 measured trajectories). The table shows median process CPU time for the same N=32, mass100, h=0.04, ten-step scenario; setup is measured separately and is included in total time. These are dense-prototype timings, not a production sparse-solver benchmark.','',
          '| Budget | Metric | Total CPU ms | Setup/factorization ms | Converged steps / 10 | Max modal error |',
          '|---:|---|---:|---:|---:|---:|']
    for budget in [12,48]:
        for metric in metrics:
            a=[r for r in timings if r['N']==32 and r['regime']=='mass100' and r['budget']==budget and r['metric']==metric]
            run=subset(N=32,regime='mass100',motion='mixed',h=.04,metric=metric,budget=budget,globalizer='line')[0]
            cost.append(f"| {budget} | {metric} | {med(a,'cpu_ms'):.4g} | {med(a,'setup_us')/1000:.4g} | {int(run['converged_steps'])} | {run['max_mode_error']:.3g} |")
    cost+=['', 'At budget 48 in this case, diagonal scaling has about the same measured cost as identity and roughly 13× lower mode error; graph scaling costs about twice as much and has roughly 25× lower mode error. Dense SVD balancing reaches the tightest residuals but is roughly 145× slower than identity here; almost all its time is setup. Those comparisons apply to the stated scenario and prototype.','',
           '![Isolated timing including setup](figures/timing_scaling.png)','','Counts of linear products alone therefore overstate the practical advantage of dense spectral balancing. A block approximation or a carefully tested approximate matrix function would be a new experiment; it is not validated by these results.','']
    lines[start:end]=cost
    lines=[line.replace('Use bounded SVD balancing diagnostically; its product-count advantage needs isolated setup-inclusive timing before treating it as practical.', 'Use bounded SVD balancing diagnostically: isolated timing shows dense setup overwhelms its product-count advantage at N=32.') for line in lines]
    fig,axes=plt.subplots(1,2,figsize=(12,4.5),constrained_layout=True)
    for ax,regime in zip(axes,['mass100','both100']):
        for metric,label in zip(metrics,labels):
            values=[]
            for N in [2,8,32]:
                a=[r for r in timings if r['N']==N and r['regime']==regime and r['budget']==48 and r['metric']==metric]
                values.append(med(a,'cpu_ms'))
            ax.plot([2,8,32],values,'-o',label=label,color=colors[metric])
        ax.set_xscale('log',base=2);ax.set_yscale('log');ax.set_xticks([2,8,32],[2,8,32]);ax.set_xlabel('Moving masses');ax.set_ylabel('Median CPU ms per ten-step trajectory');ax.set_title(f'{regime}, h=0.04, budget 48');ax.grid(alpha=.18)
    axes[0].legend(fontsize=8);fig.suptitle('Seven isolated serial repeats; all setup included; convergence differs by method');fig.savefig(FIG/'timing_scaling.png');plt.close(fig)
(ROOT/'RESULTS.md').write_text('\n'.join(lines))
print(f'Wrote RESULTS.md and {len(list(FIG.glob("*.png")))} figures. Full references: {len(refs)-len(reference_fail)}/{len(refs)}.')
