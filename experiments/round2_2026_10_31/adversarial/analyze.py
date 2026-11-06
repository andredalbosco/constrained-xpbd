from pathlib import Path
import csv,json,collections,math
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
P=Path(__file__).resolve().parent
D=[]
for folder in ['results','long_results','angles_results','angles_long_results']:
 for r in csv.DictReader((P/folder/'summary.csv').open()):
  r['source']=folder
  for k in list(r):
   if k not in ['case','method','status','reason','source']:r[k]=float(r[k])
  D.append(r)
(P/'figures').mkdir(exist_ok=True)
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'figure.dpi':150})
def trace(case,method,folder='results'):
 f=P/folder/f'trace_{case}_{method}_0.033333.csv'
 return np.genfromtxt(f,delimiter=',',names=True)
fig,axs=plt.subplots(2,2,figsize=(12,8),layout='constrained')
for name,label in [('rattle','RATTLE'),('dg64','Cartesian DG'),('sym_energy','Symmetric energy projection')]:
 t=trace('rod_separatrix_below',name);axs[0,0].plot(t['t'],np.unwrap(np.arctan2(t['x0'],-t['x1'])),label=label)
t=trace('rod_separatrix_below','angle_dg','angles_results');axs[0,0].plot(t['t'],np.unwrap(np.arctan2(t['x0'],-t['x1'])),label='Angular DG',ls='--');axs[0,0].plot(t['t'],np.unwrap(np.arctan2(t['reference_x0'],-t['reference_x1'])),label='Refined reference',color='black',ls=':',lw=2)
axs[0,0].set(title='Just below pendulum separatrix · 30 Hz',xlabel='Time (s)',ylabel='Unwrapped angle (rad)');axs[0,0].legend(fontsize=8)
rows=list(csv.DictReader((P/'radial_mode.csv').open()))
for name,label in [('tangent','Initially tangent'),('after_dg_0.01','After one DG step'),('injected_radial_0.01','Injected radial speed 0.01')]:
 a=[r for r in rows if r['initial']==name];axs[0,1].loglog([float(r['h']) for r in a],[float(r['step_doubling_error']) for r in a],'.-',label=label)
axs[0,1].loglog([float(r['h']) for r in a],[float(r['tolerance']) for r in a],'k--',label='Acceptance tolerance')
axs[0,1].set(title='Radial velocity defeats step-doubling',xlabel='Trial step h (s)',ylabel='Scaled state discrepancy');axs[0,1].legend(fontsize=8)
for method,label in [('dg64','Raw DG / 64 iterations'),('dg_adaptive','DG / residual-gated subdivision'),('dg_accuracy','DG / accuracy subdivision')]:
 t=trace('radial_stiff',method);e0=504.5;axs[1,0].semilogy(t['t'],np.maximum(abs(t['energy']-e0)/e0,1e-15),label=label)
axs[1,0].set(title='Stiff radial spring · 30 Hz · k=100,000',xlabel='Time (s)',ylabel='Relative energy error');axs[1,0].legend(fontsize=8)
ms=['dg64','dg_tangent_energy','rattle_energy','post_energy','sym_energy','angle_dg'];vals=[]
for m in ms:
 r=next(r for r in D if r['case']=='rod_release' and r['method']==m and r['duration']==10 and abs(r['h']-1/30)<1e-12)
 vals.append(r['max_local_reversal'])
# A failed backward solve is a category, not a numerical error of 1e9.
for i,v in enumerate(vals):
 if v==1e9:axs[1,1].text(i,1e-1,'backward\nfailed',ha='center',va='center',fontsize=8)
 else:axs[1,1].bar(i,max(v,1e-16))
axs[1,1].set_yscale('log');axs[1,1].set_ylim(1e-16,1);axs[1,1].set_xticks(range(len(ms)),['Cartesian\nDG','DG + tangent\n+ energy','RATTLE\n+ energy','Post\nprojection','Symmetric\nprojection','Angular\nDG']);axs[1,1].set(title='Pendulum release · reversal over first 32 steps',ylabel='Scaled forward/backward discrepancy')
fig.savefig(P/'figures'/'new_failures_and_controls.png');plt.close(fig)

def fmt(v):return f'{v:.3g}'
lines=['# Adversarial mechanics and pendulum compositions','',f'{len(D):,} trajectory configurations: '+', '.join(f'{v} {k}' for k,v in collections.Counter(r['status'] for r in D).items())+'. A completed trajectory may contain unconverged steps; these remain explicitly counted.','',f"{sum(r['reference_certified']==0 for r in D)} rows have uncertified references and are excluded from trajectory ranking. All pendulum references pass the two-resolution check.",'','## Pendulum just above the separatrix','', 'Initial angle 0, angular speed 1.001 × 2√9.81, unit mass and length. Ten seconds at 30 Hz. This is deliberately sensitive to small energy changes.','', '| Method | Position RMS (m) | Max relative energy error | Max scaled tangency | Status |','|---|---:|---:|---:|---|']
for m in ['rattle','dg64','dg_tangent','dg_tangent_energy','sym_geom','sym_energy','post_energy','angle_dg','rattle_accuracy']:
 r=next(r for r in D if r['case']=='rod_separatrix_above' and r['method']==m and r['duration']==10 and abs(r['h']-1/30)<1e-12)
 lines.append(f"| {m} | {fmt(r['position_rms'])} | {fmt(r['max_energy_relative'])} | {fmt(r['max_tangency_scaled'])} | {r['status']} |")
lines+=['','RATTLE is a useful structural baseline, but at30 Hz the just-below-separatrix case changes from reference libration into numerical rotation (the candidate angle reaches22.37 rad while the reference stays between−3.053 and3.053 rad). Symmetric projection and the angular discrete-gradient control substantially reduce that specific error. This does not establish superiority at equal cost.','', '## Finite convergence and extreme steps','', '| Case / h | Method | Max relative energy error | Unconverged accepted steps | Rejected attempts |','|---|---|---:|---:|---:|']
for case,h in [('radial_stiff',1/30),('elastic_stiff',.25),('rod_fast_spin',.25)]:
 for m in ['dg64','dg_adaptive']:
  r=next(r for r in D if r['case']==case and r['method']==m and r['duration']==10 and abs(r['h']-h)<1e-12)
  lines.append(f"| {case} / {h:.5g} | {m} | {fmt(r['max_energy_relative'])} | {int(r['unconverged_steps'])} | {int(r['rejected_attempts'])} |")
lines+=['','Residual gating rescues energy when an incomplete nonlinear solve is the cause. It does not certify temporal accuracy: the fully converged fast-spin Cartesian DG trajectory still has 1.41 m position RMS at h=0.25. Increasing Newton iterations cannot repair the wrong temporal resolution.','', '## The inherited radial mode','',r'For a rigid rod, exact old/new length and Cartesian midpoint kinematics imply $n^T(v_{n+1}+v_n)=O(h)$. If the old velocity contains a finite normal component, the small-step limit is a reflection $v_{n+1}\to v_n-2n(n^Tv_n)$. This preserves kinetic energy. Two half steps reflect twice, while one full step reflects once. The resulting velocity discrepancy tends to $2|n^Tv_n|$, rather than zero.','', 'The probe starts from a genuine DG output, not only an artificial off-manifold input. After h=0.01 from pendulum release, radial velocity is 8.8082e-6. At the next trial h≈2.03e-6 it becomes −8.8098e-6; the scaled step-doubling discrepancy is 1.2459e-5 against tolerance 2.03e-10, while the solver residual is 3.03e-11.','', 'The accuracy controller fails on 24/28 ten-second gravity-pendulum configurations; both free-spin cases at each step size work because they do not generate this radial defect. Two additional stiff-spring configurations hit the controller depth limit. Those latter failures are not attributed to the rigid hidden mode.','', '## Projection and reversibility','', 'The symmetric projection experiment follows Hairer’s same-multiplier before/after construction around an implicit midpoint step of the acceleration-level pendulum ODE. Position, tangency, and optionally energy are enforced together in fixed dimensionless coordinates. Ordinary postprojection is a separate control. [Primary source](https://www.unige.ch/~hairer/preprints/symproj.pdf).','', 'For release at 30 Hz, symmetric energy projection has local reversal discrepancy 1.37e-15; postprojection 2.11e-6; DG followed by tangent/energy repair 2.30e-4. RATTLE plus energy-only velocity rescaling encounters a backward step whose position requires negative kinetic energy. A small energy error or exact tangency does not by itself imply reversibility.','', 'Angular DG uses the same potential in angle coordinates and midpoint angular velocity. Rod length and tangency hold by construction; the exact potential secant closes energy. It intentionally changes Cartesian midpoint kinematics. All 42 angle-control trajectories complete, including h=0.25 fast rotation; that case has 0.161 m RMS versus 1.41 m for Cartesian DG. This is a single-pendulum coordinate control, not a general chain solver.','', '## Numerical limits','', '- Relative energy uses initial excitation energy above the hanging minimum; it is not divided by signed total energy.','- Reversal compares one forward/backward step from each of the first 32 states. The CSV sentinel 1e9 means a backward solve failed, not a measured discrepancy of that magnitude.','- Midpoint and DG impulse defects are macro-step diagnostics; adaptive compositions obey those equations at each accepted substep, not necessarily at the macro-step.','- Accuracy subdivision accepts two half steps when their scaled discrepancy from one full step is below 1e-4 |h|. It has a depth cap14 and keeps the two-half-step state without extrapolation.','- The residual controller subdivides only on solver/geometry failure, to depth10. State-dependent subdivision is not claimed time reversible.','- Timings here were collected concurrently and include the first32 backward probes; use them only as diagnostic workload, not performance rankings.','- This experiment does not cover collisions, friction, damping, external time-dependent forcing, or rigid chains.','', '![New failure mechanisms and controls](figures/new_failures_and_controls.png)','']
(P/'RESULTS.md').write_text('\n'.join(lines))
(P/'counts.json').write_text(json.dumps({'trajectories':len(D),'status':dict(collections.Counter(r['status'] for r in D)),'uncertified_reference_rows':sum(r['reference_certified']==0 for r in D)},indent=2))
print(P/'RESULTS.md')
