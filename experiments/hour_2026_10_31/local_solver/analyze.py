#!/usr/bin/env python3
import csv
import statistics as st
from collections import Counter
from pathlib import Path
ROOT=Path(__file__).resolve().parent
def read(name):
    p=ROOT/'results'/name
    return list(csv.DictReader(p.open())) if p.exists() else []
def num(r,k):return float(r[k])
def fmt(v):return f'{v:.4g}'
def table(head,rows):return '\n'.join(['| '+' | '.join(head)+' |','| '+' | '.join(['---']*len(head))+' |']+['| '+' | '.join(map(str,r))+' |' for r in rows])
S=read('small.csv');C=read('chains.csv');B=read('branch_comparison.csv');P=read('block_comparison.csv');V=read('validation.csv');L=read('large_validation.csv');R=read('reference_checks.csv');ST=read('small_stable.csv');OLD=read('long_block_comparison.csv');NEW=read('long_block_stable.csv');SP=read('iteration_spectrum.csv');SV=read('spectrum_validation.csv');IV=read('incremental_validation.csv');W=read('warm_block_dynamics.csv');ROD=read('rod_reaction_control.csv')
lines=['# Numerical local-work results','','Generated from retained CSVs. All options solve the same final quarter objective with explicit W, beta=gamma=1, and a frozen start-only midpoint-predicted DG target. Finite unconverged steps are accepted deliberately, counted, and never classified as solved objectives.','']
if V:
    lines += [f'{sum(r["passed"]=="1" for r in V)}/{len(V)} five-environment analytic-vs-AD derivative checks pass.']
if L:lines += [f'{sum(r["passed"]=="1" for r in L)}/{len(L)} larger-chain checks pass (N8 against AD, N32 directional finite differences).']
if R:lines += [f'{len(R)} reference refinement checks: maximum position difference {fmt(max(num(r,"position_refinement") for r in R))}, velocity difference {fmt(max(num(r,"velocity_refinement") for r in R))}, normalized energy difference {fmt(max(num(r,"relative_energy_refinement") for r in R))}.']
if S:
    lines += ['','## Initialization and globalization','','81 one-second trajectories cover the five environments, with nominal and k×10 compliant cases. The rod has one nominal case. Each row below therefore covers nine trajectories.','']
    data=[]
    for method in sorted({r['method'] for r in S}):
        a=[r for r in S if r['method']==method]
        data.append([method,f'{sum(r["status"]=="completed" for r in a)}/{len(a)}',sum(int(r['inexact_steps']) for r in a),fmt(max(num(r,'max_primal') for r in a)),fmt(st.mean(num(r,'mean_evaluations') for r in a))])
    lines += [table(['Method','All steps converged','Inexact steps','Worst balance','Mean evaluations/step'],data),'','The midpoint initializer is the decisive repair in these tests. Continuation, residual damping, and the merit safeguard do not make the BE initializer uniformly reliable. More solver iterations cannot repair a line search that has already stopped. Warm multipliers add little on the small cases because the midpoint predictor already lands close to the relevant branch.','']
    data=[]
    for r in sorted([r for r in S if r['method']=='kkt_midpoint'],key=lambda r:(r['case'],r['scenario'])):data.append([r['case'],r['scenario'],fmt(num(r,'rms_position')),fmt(num(r,'max_energy_error')),fmt(num(r,'mean_evaluations'))])
    lines += [table(['Case','Regime','RMS x error','Max normalized energy error','Evaluations/step'],data),'','The scalar harmonic predictor is itself the constrained stationary solution, hence one evaluation. It still has the phase error of midpoint mechanics, especially when stiff; zero energy error is not exact dynamics.','']
if B:
    lines += ['## Same-state branch comparison','',f'{len(B)} method comparisons use common physical old states at three phases and two time steps. Each objective gap is against the midpoint-initialized joint reference for that exact old state. Warm-only duplicates are omitted.','']
    data=[]
    for method in sorted({r['method'] for r in B}):
        a=[r for r in B if r['method']==method];con=[r for r in a if r['converged']=='1'];other=[r for r in con if abs(num(r,'objective_gap'))>1e-6]
        data.append([method,f'{len(con)}/{len(a)}',len(other),sum(num(r,'reduced_curvature')<0 for r in con),fmt(max([abs(num(r,'objective_gap')) for r in con] or [0]))])
    lines += [table(['Method','KKT converged','Objective gap above1e-6','Negative reduced curvature','Largest converged objective gap'],data),'','A converged KKT point is classified further by reduced curvature. Negative curvature is a stationary saddle, and positive curvature is evidence of a local minimum only. Different old states along separately evolving trajectories are not used for objective-gap claims.','']
if C:
    lines += ['## Shared multipliers on longer chains','',f'{len(C)} eight-step trajectories cover N2/8/32, nominal/stiffness×10/alternating-mass×10, fixed evaluation budgets, and warm/cold shared multipliers.','']
    data=[]
    for N in ['2','8','32']:
        for method in ['kkt_midpoint_cold','kkt_midpoint_warm','AL_cold_128','AL_warm_128','AL_cold_512','AL_warm_512','AL_rhoN_cold','AL_rhoN_warm']:
            a=[r for r in C if r['N']==N and r['method']==method]
            if a:data.append([N,method,f'{sum(r["status"]=="completed" for r in a)}/{len(a)}',fmt(max(num(r,'max_primal') for r in a)),fmt(max(num(r,'max_stationarity') for r in a)),fmt(st.mean(num(r,'mean_evaluations') for r in a))])
    lines += [table(['N','Method','All steps converged','Worst balance','Worst stationarity','Mean evaluations'],data),'','At a fixed evaluation budget the original split partition performs fewer complete passes as N grows. With N32 and128 evaluations, it can update the dual only after a partial three-pass primal cycle. Warm duals and rho=10N do not remove the unresolved primal/dual modes. Rho tuning is numerical; it does not change the constrained objective.','']
if P:
    lines += ['## Grouping position and velocity in local blocks','',f'{len(P)} same-state AL solves compare vertex-full and edge-full blocks, rho1/10/100, and256/1024 evaluations. Objective gaps use a converged joint solve of exactly the same objective and old state.','']
    data=[]
    for N in ['2','8','32']:
        for partition in ['vertex_full','edge_full']:
            for budget in ['256','1024']:
                a=[r for r in P if r['N']==N and r['partition']==partition and r['budget']==budget]
                if a:data.append([N,partition,budget,f'{sum(r["converged"]=="1" for r in a)}/{len(a)}',fmt(st.median(num(r,'primal') for r in a)),fmt(st.median(num(r,'stationarity') for r in a))])
    lines += [table(['N','Block','Budget','KKT converged','Median balance','Median stationarity'],data),'']
if ST:
    lines += ['## Stable line-search differences','','The exact incremental objective/balance difference changes arithmetic only. It avoids cancellation of nearly equal total gravitational/objective values in the Armijo test.', '',f'{sum(r["passed"]=="1" for r in IV)}/{len(IV)} incremental-difference checks pass across all five small environments and an N8 mass-contrast chain. The 45 additional one-second trajectories retain all baseline CSVs.','']
    data=[]
    for method in sorted({r['method'] for r in ST}):
        a=[r for r in ST if r['method']==method]
        data.append([method,f'{sum(r["status"]=="completed" for r in a)}/{len(a)}',sum(int(r['inexact_steps']) for r in a),fmt(st.mean(num(r,'mean_evaluations') for r in a))])
    lines += [table(['Method','All steps converged','Inexact steps','Mean evaluations/step'],data),'','In one/two-mass systems, an edge-full block can contain the entire primal system. These small-case successes are not evidence that purely local solves scale.','']
if OLD and NEW:
    lines += [f'Longer same-state comparisons converge in {sum(r["converged"]=="1" for r in OLD)}/{len(OLD)} cases with ordinary subtraction and {sum(r["converged"]=="1" for r in NEW)}/{len(NEW)} with stable differences. Both remaining stable failures are N32 mass-contrast/rho10 cases; fixed budgets still matter.','']
    data=[]
    for N,scenario,partition,rho in [('8','nominal','edge_full','100'),('32','nominal','edge_full','100'),('32','mass10','edge_full','100'),('32','nominal','global','100'),('32','mass10','global','100')]:
        pick=lambda rows:next(r for r in rows if r['N']==N and r['scenario']==scenario and r['partition']==partition and r['rho']==rho)
        a,b=pick(OLD),pick(NEW);data.append([N,scenario,partition,a['evaluations']+(' (inexact)' if a['converged']=='0' else ''),b['evaluations'],fmt(num(b,'primal')),fmt(num(b,'stationarity'))])
    lines += [table(['N','Regime','Block','Ordinary evaluations','Stable evaluations','Stable balance','Stable stationarity'],data),'','The global block is a control that solves the full primal matrix on each visit. It is not a local method or a claim of equal factorization cost. The tiny joint KKT solve remains much cheaper on these fixtures.','']
if SP:
    lines += ['## Local iteration spectrum','',f'{sum(r["passed"]=="1" for r in SV)}/{len(SV)} finite-difference checks validate the derived local full-step AL iteration matrix, with maximum relative error {fmt(max(num(r,"directional_relative_error") for r in SV))}. These local asymptotic predictions assume the nearby full-step regime; they do not guarantee nonlinear convergence or predict an arbitrary initial error vector.','']
    data=[]
    for partition,rho in [('split','10'),('split','100'),('edge_full','100'),('edge_full','1000'),('global','100'),('global','1000')]:
        r=next(x for x in SP if x['N']=='32' and x['partition']==partition and x['rho']==rho)
        data.append([partition,rho,fmt(num(r,'spectral_radius')),fmt(num(r,'predicted_cycles_for_1e6_reduction')),fmt(num(r,'predicted_block_visits'))])
    lines += [table(['N32 block','rho','Spectral radius','Cycles for six decades','Block visits for six decades'],data),'','The slow modes remain even after the line-search arithmetic fix. Increasing rho can improve dual correction while worsening local primal coupling. See [TURNING_AND_LOCALITY.md](TURNING_AND_LOCALITY.md) for the derivation and the scalar turning-point trigger.','']
if W:
    lines += ['## Warm starts after a converged joint first step','',f'{len(W)} retained trajectories compare cold duals, ordinary warm duals, and one joint first step followed by warm local AL. All use stable differences, edge-full blocks, and rho100. The joint first step is counted and labeled as a hybrid initialization.','']
    data=[]
    for N in ['8','32']:
        for budget in ['512','2048']:
            for initialization in ['cold','warm','joint_bootstrap']:
                a=[r for r in W if r['N']==N and r['budget']==budget and r['initialization']==initialization]
                if a:data.append([N,budget,initialization,f'{sum(r["status"]=="completed" for r in a)}/{len(a)}',sum(int(r['inexact_steps']) for r in a),fmt(max(num(r,'max_primal') for r in a)),fmt(st.mean(num(r,'mean_evaluations') for r in a))])
    lines += [table(['N','Budget','Dual initialization','All steps converged','Inexact steps','Worst balance','Mean evaluations'],data),'']
if ROD:
    lines += ['## Single-rod reaction control','','All main hour rod runs include the old physical reaction in the quarter predictor. A separate omitted-reaction control also appears second order for this single-degree-of-freedom pendulum because energy, position, and tangency already restrict its motion strongly. It must not be extrapolated to multi-rod chains.','']
    data=[]
    for r in ROD:data.append(['included' if r['reaction']=='1' else 'omitted',str(round(1/num(r,'h'))),fmt(num(r,'rms_position')),fmt(num(r,'max_energy_error')),r['inexact_steps']])
    lines += [table(['Old reaction','Hz','Two-second RMS x error','Max energy error','Inexact steps'],data),'','The fixed per-step residual tolerance accumulates to energy error about 1.5e-6 at 240 Hz. Both the convergence order and that accumulation are retained, rather than treating a green residual as exact conservation.','']
lines += ['## Boundaries','','- Work-target accuracy, local energy feasibility, stationary-objective convergence, and trajectory accuracy are distinct measurements.','- Some inexact large-chain updates stay closer to the predictor and therefore to the high-resolution trajectory. That does not make them accurate solves of the fixed objective.','- The chain family changes physical lengths, stiffnesses, and total normalization with N. It does not establish replication invariance or a uniform continuum limit.','- Dense analytic derivative storage and dense KKT solves remain in this prototype. Recorded timings are not isolated scaling benchmarks.','- Positive reduced curvature supports local minimality; it never certifies a global optimum.','']
(ROOT/'RESULTS.md').write_text('\n'.join(lines))
print('wrote RESULTS.md')
