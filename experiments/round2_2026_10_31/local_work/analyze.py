#!/usr/bin/env python3
"""Summarize the local-work experiment without dropping inexact solves."""
import csv
import statistics as st
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parent
def rows(name):
    p = ROOT/'results'/name
    return list(csv.DictReader(p.open())) if p.exists() else []
def n(r, key): return float(r[key])
def fmt(x): return f'{x:.4g}'
def table(headers, data):
    return '\n'.join(['| '+' | '.join(headers)+' |', '| '+' | '.join(['---']*len(headers))+' |'] + ['| '+' | '.join(map(str, r))+' |' for r in data])

P=rows('prediction_accuracy.csv'); R=rows('prediction_reference_checks.csv'); A=rows('local_solver_comparison.csv'); D=rows('frozen_target_dynamics.csv'); L=rows('local_solver_dynamics.csv')
lines=['# Local-work results', '', 'Generated from the complete CSVs by `python3 analyze.py`. These are research prototypes. The original work-target comparison retains its position objective; a separate final cross experiment also tests the repaired quarter objective.', '']
lines += ['## Start-only prediction', '', f'{len(R)} physical state/timestep probes, {len(P)} per-incidence target comparisons. Maximum normalized 128→256 reference-work refinement: {fmt(max(n(r,"max_normalized_work_refinement") for r in R))}.', '', 'The following is mean absolute per-incidence work error divided by the initial energy above the minimum. It pools the five compliant initial conditions and includes their gravitational work incidences; rigid-rod work is reported separately.', '']
laws=['frozen_power','acceleration_chord','taylor2','predicted_midpoint','predicted_dg','implicit_dg','implicit_midpoint_dg']
data=[]
for scenario,h in [('nominal',1/120),('nominal',1/30),('stiff',1/120),('stiff',1/30)]:
    rr=[r for r in P if r['scenario']==scenario and abs(n(r,'h')-h)<1e-12 and r['case']!='rod_kick']
    data.append([scenario,'1/'+str(round(1/h))]+[fmt(st.fmean(n(r,'normalized_error') for r in rr if r['target']==law)) for law in laws])
lines += [table(['Scenario','h']+laws,data),'','`stiff` means k×100. The start-only midpoint linear solve improves the useful-step-size prediction error and avoids the explicit predictor’s large stiff excursion. The backward-Euler linear predictor is more biased at small h but can be more robust at the largest stiff step. Neither gives an identity at the actual optimized endpoint.','']
rr=[r for r in P if r['case']=='rod_kick' and r['potential']=='0']
lines += [table(['Rigid-rod target','Max absolute error'],[[law,fmt(max(n(r,'absolute_error') for r in rr if r['target']==law))] for law in sorted({r['target'] for r in rr})]),'','The Taylor cancellation assumes a feasible rod and tangent old velocity. `projected_dg` uses a feasible predicted rod endpoint. The two linear-solve predictors alias this projected construction for the rod.','']
lines += ['## Local augmented-Lagrangian solves','','All work variables and all vertex multipliers are shared. The methods minimize the same explicit-W objective; only their solve schedules differ. A solved row is not counted as convergence unless the full feasibility and stationarity tolerances pass.','']
data=[]
for budget in ['32','128','512']:
    for method in ['joint_budget','local_AL']:
        rr=[r for r in A if r['method']==method and r['budget']==budget]
        data.append([method,budget,f'{sum(r["converged"]=="1" for r in rr)}/{len(rr)}',fmt(st.median(n(r,'primal_residual') for r in rr)),fmt(st.median(n(r,'stationarity') for r in rr))])
lines += [table(['Method','Evaluation budget','KKT converged','Median feasibility','Median stationarity'],data),'']
data=[]
for rho in ['0.1','1','10','100']:
    for adaptive in ['0','1']:
        rr=[r for r in A if r['method']=='local_AL' and r['budget']=='512' and r['rho_initial']==rho and r['adaptive']==adaptive]
        data.append([rho,'progress adaptive' if adaptive=='1' else 'fixed',f'{sum(r["converged"]=="1" for r in rr)}/{len(rr)}',fmt(st.median(n(r,'primal_residual') for r in rr))])
lines += [table(['Initial rho','Penalty schedule','KKT converged at512','Median feasibility'],data),'','ρ≈10 works well with this particular normalization. These are full-AD evaluation budgets, not a claim of equal FLOPs or distributed speedup. The dense joint reference is much cheaper on the small fixtures.','']
branches={}
for r in A:
    if r['method']=='local_AL' and r['converged']=='1' and n(r,'objective_gap')>1e-5:
        if r['fixture'] not in branches or n(r,'objective_gap')>n(branches[r['fixture']],'objective_gap'):branches[r['fixture']]=r
lines += ['Different converged branches remain possible:', '',table(['Fixture','Objective gap','Local reduced curvature','Joint reduced curvature'],[[k,fmt(n(r,'objective_gap')),fmt(n(r,'reduced_min_eigenvalue')),fmt(n(r,'joint_reduced_min_eigenvalue'))] for k,r in sorted(branches.items())]),'','Positive reduced curvature confirms that these are different local minima, not merely unfinished residual corrections. The joint reference is a checked local minimum, not a certified global optimum.','']
lines += ['## Frozen targets in dynamics','',f'{len(D)} one-second runs, statuses: {dict(Counter(r["status"] for r in D))}. {sum(int(r["nonconverged_steps"])>0 for r in D)} runs contain at least one finite but inexact optimization step.', '', 'For nominal stiffness, h=1/60, and soft endpoint kinematics:', '']
data=[]
for case in ['scalar_release','radial_orbit','rod_kick','elastic_mixed','chain_mixed']:
    for target in ['frozen_power','predicted_dg','implicit_midpoint_dg','live_dg']:
        rr=[r for r in D if r['case']==case and r['target']==target and r['kinematics']=='soft_endpoint' and r['kfactor']=='1' and abs(n(r,'h')-1/60)<1e-12]
        if rr:
            r=rr[0];data.append([case,target,fmt(n(r,'rms_position')),fmt(n(r,'first_displacement_ratio')),fmt(n(r,'max_energy_error')),r['nonconverged_steps']])
lines += [table(['Case','Target','RMS x error','First displacement / reference','Max energy error','Inexact steps'],data),'','More accurate predicted work does not consistently improve motion. For example, the scalar first-step displacement remains biased despite a converged objective, and the elastic/chain soft-midpoint trajectories retain substantial phase error. Work-target quality and mechanical/kinematic consistency are separate questions.','']
if L:
    lines += ['## Local-solver rollouts','',f'{len(L)} half-second rollouts at h=1/60. Statuses: {dict(Counter(r["status"] for r in L))}. Multipliers reset each time step; no trajectory is silently removed for inexact solves.','']
    data=[]
    for method,budget,rho in [('joint','128','0'),('local_AL','128','10'),('local_AL','512','10'),('local_AL','512','100')]:
        rr=[r for r in L if r['method']==method and r['budget']==budget and r['rho']==rho]
        if rr:data.append([method,budget,rho,f'{sum(r["status"]=="completed" for r in rr)}/{len(rr)}',fmt(st.median(n(r,'max_energy_error') for r in rr)),fmt(st.median(n(r,'rms_position') for r in rr))])
    lines += [table(['Method','Budget','rho','All steps converged','Median max energy error','Median RMS x error'],data),'','Every recorded prefix is available in `results/local_traces/`; summary rows retain unfinished solves and their residuals.','']
C=rows('objective_work_cross.csv'); S=rows('scalar_stalls.csv')
if C:
    lines += ['## Corrected objective and start-only work', '', '24 one-second runs at h=1/60, with unchanged explicit-W balances, beta=gamma=1, and soft midpoint velocity penalty. Joint solves get 128 evaluations; local AL gets 512 at fixed rho10. Six derivative checks pass.', '']
    data=[]
    for case in ['scalar_release','radial_orbit','chain_mixed']:
        for objective,target in [('original','frozen_power'),('quarter','frozen_power'),('quarter','implicit_midpoint_dg')]:
            r=next(x for x in C if x['case']==case and x['objective']==objective and x['target']==target and x['method']=='joint')
            data.append([case,objective,target,fmt(n(r,'rms_position')),r['nonconverged_steps'],fmt(n(r,'max_primal'))])
    lines += [table(['Case','Position objective','Start target','RMS x error','Inexact joint steps','Max balance'],data),'','The repaired objective makes the more accurate start-only work target useful for radial motion and the chain. Their reported joint trajectories converge at every step; local AL closely matches them, with 4/1 mildly inexact steps for corrected chain frozen-power/predicted-DG. Scalar joint results include two substantial line-search stalls and must not be read as solved-objective accuracy.','']
if S:
    lines += ['The scalar diagnosis replays all four stalled old states without changing their objective or constraints. Increasing the BE-initialized budget from128 to4096 does not help: the same line search stalls after69–112 evaluations. The constraint Jacobian retains minimum singular value0.1887–0.1892, so this is not loss of constraint rank. A start-only linear-midpoint initialization converges in1–3 evaluations; rho10 AL converges in395–453. An independent dense feasible-circle scan agrees with the successful local minimum. The stalled state approaches near-zero reduced curvature while still infeasible; it is not a verified stationary branch. See [CROSS_OBJECTIVE.md](CROSS_OBJECTIVE.md) and `results/scalar_stalls.csv`.','']
lines += ['## Limits and failure categories','','- Predictor error: explicit Taylor endpoints can overshoot stiff motion; the reference-work error exposes this before optimization.','- Geometric failure: unprojected rod predictions generally assign nonzero rod work; tangent power/Taylor cancellation and feasible predicted DG are distinct fixes.','- Inexact solve: finite trajectory steps with unresolved feasibility/stationarity are counted and retained.','- Local branch: positive-curvature local minima can differ by order and rho even after convergence.','- Physical-model error: converged energy/work balances can coexist with biased displacement, phase, or velocity because the objective/kinematics remain unchanged.','- Scalability is untested: the local algorithm restricts primal variables but uses a dense AD evaluator for this prototype.','']
(ROOT/'RESULTS.md').write_text('\n'.join(lines))
print('wrote RESULTS.md')
