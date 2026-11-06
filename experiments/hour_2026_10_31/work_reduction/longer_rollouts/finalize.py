from pathlib import Path
import csv,json
root=Path(__file__).resolve().parent
rows=list(csv.DictReader((root/'combined.csv').open()));by={(r['N'],r['stiffness'],round(float(r['h']),12),r['method']):r for r in rows}
status=[]
for method,directory in [('reduced_GN_polish','gn_polish'),('full','full'),('reduced_exact','exact_parallel')]:
 for N in [8,32]:
  for k in [1,10]:
   for h in [1/30,1/60]:
    key=(str(N),str(k),round(h,12),method);r=by.get(key);p=root/directory/f'{N}_{k:.6f}_{h:.6f}_{method}.csv';frames=0;seconds=0
    if r:state='complete' if r['complete']=='1' else 'solver_failure';frames=int(r['steps']);seconds=frames*h
    elif p.exists():
     state='interrupted_by_wall_clock_cutoff';raw=p.read_bytes();(p.with_suffix('.interrupted_raw.txt')).write_bytes(raw)
     if not raw.endswith(b'\n'):raw=raw[:raw.rfind(b'\n')+1];p.write_bytes(raw)
     trace=list(csv.DictReader(p.open()));frames=len(trace);seconds=float(trace[-1]['t']) if trace else 0
    else:state='not_started_due_wall_clock_cutoff'
    status.append(dict(N=N,stiffness=k,h=h,method=method,status=state,saved_accepted_steps=frames,saved_accepted_time=seconds,requested_time=5))
(root/'run_status.json').write_text(json.dumps(dict(cutoff_utc='2026-10-31 19:02:30 UTC',all_processes_stopped=True,runs=status),indent=2)+'\n')
inter=[r for r in status if r['status'].startswith('interrupted')];notstarted=[r for r in status if r['status'].startswith('not_started')]
p=root/'RESULTS.md';s=p.read_text();s+='\n## Wall-clock cutoff\n\nAll 8 GN-polish and all 8 exact-reduced configurations reached terminal outcomes. Five full controls completed (all four N8 cases and N32 stiffness1 at30Hz). Full N32 stiffness1 at60Hz was interrupted by the wall-clock cutoff before any accepted trace row was flushed; both full N32 stiffness10 controls were not started. These three missing outcomes are **not solver failures**. The original interrupted file is retained, but no unsaved state or accuracy is inferred. See `run_status.json` for every requested configuration.\n\nThe first full N32 terminal row arrived between an interim analysis and interruption; the final report was regenerated after all processes stopped. All21 terminal-outcome first-second comparisons pass exactly, including evaluation counts; two have no accepted states because of immediate solver failure. No prefix comparison is possible for the interrupted empty trace.\n';p.write_text(s)
print(s)
