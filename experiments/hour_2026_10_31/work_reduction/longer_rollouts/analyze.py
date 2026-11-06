from pathlib import Path
import csv,sys,math,json
root=Path(__file__).resolve().parent
old=Path(sys.argv[1]) if len(sys.argv)>1 else root.parent/'rollout_results'
def read(p):
 return list(csv.DictReader(p.open())) if p.exists() else []
def key(r): return (r['N'],r['stiffness'],r['mass_ratio'],round(float(r['h']),12))
refs=read(root/'references/references.csv');cert={key(r):r['certified']=='1' for r in refs}
allrows=[];prefix=[]
for directory in ['gn_polish','full','exact_parallel']:
 for r in read(root/directory/'trajectories.csv'):
  allrows.append(dict(r,reference_certified=int(cert.get(key(r),False)),budget_exhausted=int(int(r['failed_evaluations'])>=128)))
  name=f"{r['N']}_{float(r['stiffness']):.6f}_{float(r['h']):.6f}_{r['method']}.csv"
  a=read(old/name);b=read(root/directory/name);b=[x for x in b if float(x['t'])<=1+1e-12]
  mx=0;evalmatch=True
  for x,y in zip(a,b):
   for c in x:
    if c=='evaluations':evalmatch &= int(x[c])==int(y[c])
    else:mx=max(mx,abs(float(x[c])-float(y[c]))/(1+abs(float(x[c]))))
  prefix.append(dict(N=r['N'],stiffness=r['stiffness'],h=r['h'],method=r['method'],old_frames=len(a),new_frames=len(b),max_scaled_difference=mx,evaluations_equal=int(evalmatch),passed=int(len(a)==len(b) and mx<1e-7 and evalmatch),has_accepted_prefix=int(bool(a))))
for name,rows in [('combined.csv',allrows),('prefix_checks.csv',prefix)]:
 if rows:
  with (root/name).open('w') as f:w=csv.DictWriter(f,fieldnames=rows[0].keys());w.writeheader();w.writerows(rows)
lines=['# Five-second work-reduction follow-up','',f'{len(allrows)}/24 requested trajectories have recorded terminal outcomes. All runs stop at the first rejected nonlinear step. Physical reference certification covers the full five-second horizon, independently of candidate completion.','', '| Method | Complete five seconds | Completed with certified reference | Recorded configurations |','|---|---:|---:|---:|']
for method in ['reduced_GN_polish','full','reduced_exact']:
 r=[r for r in allrows if r['method']==method];lines.append(f"| {method} | {sum(x['complete']=='1' for x in r)}/8 | {sum(x['complete']=='1' and x['reference_certified'] for x in r)}/8 | {len(r)}/8 |")
lines+=['',f"The original 128/256 reference pair certifies {sum(r['certified']=='1' for r in refs)}/{len(refs)} configurations. Both N32 stiffness×10 references fail over five seconds; their physical trajectory RMS remains **uncertified**, regardless of a candidate's solver convergence or energy error.",'','| Failed method | N / stiffness | h | Accepted time (s) | Reason | Rejected primal | Rejected stationarity |','|---|---|---:|---:|---|---:|---:|']
for r in allrows:
 if r['complete']!='1':lines.append(f"| {r['method']} | {r['N']} / {r['stiffness']} | {float(r['h']):.6g} | {int(r['steps'])*float(r['h']):.6g} | {r['failure_reason']}{' (128-evaluation budget exhausted)' if r['budget_exhausted'] else ''} | {float(r['failed_primal']):.3g} | {float(r['failed_stationarity']):.3g} |")
fine=read(root/'references_fine/references.csv')
if fine:
 lines+=['','Targeted 512/1024-substep reference follow-up (the candidate-error baseline is still the original256-substep reference):','','| N / stiffness | h | Position gap | Scaled velocity gap | Energy drift | Certified |','|---|---:|---:|---:|---:|---|']
 for r in fine:lines.append(f"| {r['N']} / {r['stiffness']} | {float(r['h']):.6g} | {float(r['position_gap']):.3g} | {float(r['scaled_velocity_gap']):.3g} | {float(r['energy_error']):.3g} | {r['certified']} |")
lines+=['',f"First-second regression: {sum(p['passed'] for p in prefix)}/{len(prefix)} saved-prefix checks pass; {sum(p['has_accepted_prefix'] for p in prefix)} contain accepted states. Maximum scaled numeric discrepancy is {max((p['max_scaled_difference'] for p in prefix),default=0):.3g}. Zero-frame comparisons only check matching immediate failure. Per-step evaluation counts are compared separately and retained.",'','These are bounded follow-up observations. They do not establish long-duration universal success, a general-mesh result, or a runtime ranking. Full residuals, failed prefixes, and raw accuracy columns remain in the CSVs. See [REPRODUCE.md](REPRODUCE.md) for exact commands, dependencies, and the unchanged formulation.']
(root/'RESULTS.md').write_text('\n'.join(lines)+'\n')
print('\n'.join(lines))
