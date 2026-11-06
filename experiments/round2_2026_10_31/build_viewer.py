#!/usr/bin/env python3
"""Build a dependency-free trajectory viewer from the recorded round-two CSVs.

Run only after the benchmark writers have finished. This reads existing data;
it does not integrate, extrapolate, resample, or repair any trajectory.
"""
import csv
import datetime
import json
import math
from pathlib import Path
import re

ROOT = Path(__file__).resolve().parent

CASES = {
    "scalar_release": ("Scalar spring · release", 1),
    "scalar_mixed": ("Scalar spring · moving", 1),
    "radial_release": ("2D central spring · release", 2),
    "radial_orbit": ("2D central spring · orbit", 2),
    "rod_release": ("Rigid pendulum · release", 2),
    "rod_kick": ("Rigid pendulum · kick", 2),
    "elastic_release": ("Elastic pendulum · release", 2),
    "elastic_mixed": ("Elastic pendulum · mixed motion", 2),
    "chain_mixed": ("Two-spring chain · mixed motion", 2),
    "chain_large": ("Two-spring chain · large motion", 2),
    "rod_separatrix_below": ("Rigid pendulum · just below separatrix", 2),
    "rod_separatrix_above": ("Rigid pendulum · just above separatrix", 2),
    "rod_fast_spin": ("Rigid pendulum · fast spin", 2),
}

LABELS = {
    "quarter_b1": "Quarter objective · β = 1",
    "quarter_h2": "Quarter objective · β ∝ h²",
    "midpotential_b1": "Midpoint-potential objective · β = 1",
    "midpotential_h2": "Midpoint-potential objective · β ∝ h²",
    "original_endpoint_b1": "Original endpoint objective · β = 1",
    "original_midpoint_b1": "Original midpoint objective · β = 1",
    "original_midpoint_b01": "Original midpoint objective · β = 0.1",
    "original_midpoint_b001": "Original midpoint objective · β = 0.01",
    "original_midpoint_h2": "Original midpoint objective · β ∝ h²",
    "original_midpoint_h2_c01": "Original midpoint objective · β ∝ 0.1 h²",
    "original_midpoint_h4": "Original midpoint objective · β ∝ h⁴",
    "half_half_bad": "Half/half objective · inconsistent control",
    "projection_xpbd_E": "Fixed XPBD target + energy projection",
    "projection_verlet_E": "Fixed Verlet/RATTLE target + energy",
    "projection_verlet_EL": "Fixed Verlet/RATTLE target + energy, L*",
    "projection_midpoint_E": "Fixed DG target + energy projection",
    "projection_midpoint_EL": "Fixed DG target + energy, L*",
    "ls_midpoint": "Midpoint mechanical residual objective",
    "ls_midpoint_E": "Midpoint residual objective + energy",
    "ls_dg": "DG mechanical residual objective",
    "ls_dg_E": "DG residual objective + energy",
    "ls_dg_EL": "DG residual objective + energy, L*",
    "dg4": "Cartesian midpoint DG · 4 iterations",
    "dg16": "Cartesian midpoint DG · 16 iterations",
    "dg64": "Cartesian midpoint DG · 64 iterations",
    "dg_adaptive": "DG · convergence-gated subdivision",
    "dg_accuracy": "DG · accuracy-controlled subdivision",
    "rattle": "RATTLE",
    "rattle_adaptive": "RATTLE · convergence-gated subdivision",
    "rattle_accuracy": "RATTLE · accuracy-controlled subdivision",
    "rattle_energy": "RATTLE + energy projection",
    "shake": "Position-projected Verlet · no final velocity projection",
    "ode_midpoint": "Cartesian ODE implicit midpoint",
    "implicit_midpoint": "Implicit midpoint",
    "post_geom": "Midpoint + endpoint geometry projection",
    "post_energy": "Midpoint + endpoint geometry/energy projection",
    "sym_geom": "Symmetric midpoint geometry projection",
    "sym_energy": "Symmetric midpoint geometry/energy projection",
    "sym_energy_adaptive": "Symmetric energy projection · gated subdivision",
    "dg_tangent": "Cartesian DG + endpoint tangency projection",
    "dg_tangent_energy": "Cartesian DG + tangency/energy projection",
    "angle_dg": "Angle-coordinate DG",
}


def number(value):
    try:
        n = float(value)
        return n if math.isfinite(n) else None
    except (TypeError, ValueError):
        return None


def compact(value):
    n = number(value)
    # Nine significant digits are ample for drawing. Scientific summary metrics
    # retain their CSV precision and must be used for numerical comparisons.
    return float(f"{n:.9g}") if n is not None else None


def metadata(value):
    try:
        float(value)
        return number(value)
    except (TypeError, ValueError):
        return value


def read_csv(path):
    with path.open(newline="") as f:
        return list(csv.DictReader(f))


def trace(path, kind):
    if not path.exists():
        return {"frames": [], "reference": [], "dof": 0,
                "issue": "Trace file is absent; no trajectory has been invented."}
    rows = read_csv(path)
    if not rows:
        return {"frames": [], "reference": [], "dof": 0,
                "issue": "No accepted output sample was recorded."}
    columns = sorted((k for k in rows[0] if re.fullmatch(r"x\d+", k)),
                     key=lambda k: int(k[1:]))
    t_key = "time" if kind == "objective" else "t"
    prefix = "ref_" if kind == "objective" else "reference_"
    frames, refs = [], []
    issue = ""
    for row in rows:
        t = number(row[t_key])
        pos = [compact(row[k]) for k in columns]
        if t is None or any(v is None for v in pos):
            issue = "A nonfinite output row occurred; only its finite prefix is displayed."
            break
        if frames and t <= frames[-1][0]:
            raise ValueError(f"Non-increasing trace times: {path}")
        frames.append([compact(t), compact(row.get("energy")), *pos])
        ref = [compact(row.get(prefix + k)) for k in columns]
        if all(v is not None for v in ref):
            refs.append([compact(t), *ref])
    return {"frames": frames, "reference": refs, "dof": len(columns), "issue": issue}


def make_scenario(case, kind, duration, hz):
    title, dim = CASES[case]
    return {"id": f"{kind}:{case}:{duration:g}", "case": case,
            "title": title, "kind": kind, "duration": duration,
            "hz": hz, "dim": dim, "rigid": case.startswith("rod_"),
            "methods": [], "reference": [], "referenceCertified": False,
            "group": "Nominal · objective comparisons" if kind == "objective" else
                     f"Rigid stress cases · {duration:g} seconds"}


def add_method(scenario, row, path, summary_path):
    objective = scenario["kind"] == "objective"
    name = row["variant"] if objective else row["method"]
    source = trace(path, scenario["kind"])
    frames = source.pop("frames")
    reference = source.pop("reference")
    n = lambda key: number(row.get(key))
    certified = n("reference_valid" if objective else "reference_certified") == 1
    method = {"id": name, "label": LABELS.get(name, name), "frames": frames,
              "dof": source["dof"], "traceIssue": source["issue"],
              "status": row.get("status", "unknown"),
              "reason": row.get("reason", ""),
              "completed": n("completed" if objective else "steps_done"),
              "requested": n("steps" if objective else "steps_requested"),
              "inexact": n("unconverged_steps"), "referenceCertified": certified,
              "energyError": n("max_energy_error" if objective else "max_energy_relative"),
              "positionError": n("rms_position" if objective else "position_rms"),
              "residual": n("max_stationarity" if objective else "max_solver_residual"),
              "balance": n("max_balance") if objective else None,
              "lengthError": n("max_rod" if objective else "max_length_relative"),
              "tangencyError": n("max_tangency" if objective else "max_tangency_scaled"),
              "summary": str(summary_path.relative_to(ROOT)),
              "trace": str(path.relative_to(ROOT)),
              "raw": {k: metadata(v) for k, v in row.items()}}
    if frames and frames[-1][0] > scenario["duration"] + 1e-5:
        raise ValueError(f"Trace extends beyond matching summary duration: {path}")
    if method["status"] == "complete" and method["completed"] == method["requested"]:
        if not frames or abs(frames[-1][0] - scenario["duration"]) > 1e-5:
            raise ValueError(f"Complete summary has an incomplete trace (writer active?): {path}")
    scenario["methods"].append(method)
    if reference and (certified, len(reference)) > (scenario["referenceCertified"], len(scenario["reference"])):
        scenario["reference"] = reference
        scenario["referenceCertified"] = certified
    return method


def build_data():
    scenarios = []
    summary = ROOT / "objectives/summary.csv"
    objective_rows = read_csv(summary)
    for case in list(CASES)[:10]:
        scenario = make_scenario(case, "objective", 10, 120)
        rows = [r for r in objective_rows if r["case"] == case and
                number(r["duration"]) == 10 and number(r["hz"]) == 120]
        if len(rows) != 22:
            raise ValueError(f"Expected 22 matched 10-second variants for {case}; got {len(rows)}")
        for row in rows:
            add_method(scenario, row, ROOT / "objectives/traces" /
                       f"{case}__{row['variant']}.csv", summary)
        scenario["outputHz"] = 60
        scenarios.append(scenario)
    for duration, cases, folders in [
        (10, ("rod_release", "rod_separatrix_below", "rod_separatrix_above", "rod_fast_spin"),
         ("results", "angles_results")),
        (60, ("rod_release", "rod_fast_spin"), ("long_results", "angles_long_results")),
    ]:
        for case in cases:
            scenario = make_scenario(case, "rigid", duration, 30)
            for folder in folders:
                summary = ROOT / "adversarial" / folder / "summary.csv"
                if not summary.exists():
                    continue
                for row in read_csv(summary):
                    if row["case"] != case or number(row["duration"]) != duration:
                        continue
                    if not math.isclose(number(row["h"]) or 0, 1 / 30, rel_tol=1e-7):
                        continue
                    path = summary.parent / f"trace_{case}_{row['method']}_{float(row['h']):.6f}.csv"
                    add_method(scenario, row, path, summary)
            if not scenario["methods"]:
                raise ValueError(f"No recorded rigid methods for {case}, {duration}s")
            scenario["outputHz"] = 30
            scenarios.append(scenario)
    for scenario in scenarios:
        ids = [m["id"] for m in scenario["methods"]]
        if len(ids) != len(set(ids)):
            raise ValueError(f"Duplicate methods in {scenario['id']}")
        if scenario["reference"]:
            scenario["start"] = scenario["reference"][0][0]
        else:
            scenario["start"] = min((m["frames"][0][0] for m in scenario["methods"] if m["frames"]), default=0)
    return {"built": datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds"),
            "scenarios": scenarios}


HTML = r'''<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>Constrained XPBD · recorded trajectory comparisons</title>
<style>
:root{color-scheme:dark;--bg:#10151c;--panel:#18212c;--text:#e8edf3;--muted:#a4b1c1;--line:#334254;--accent:#8ec5ff;--ok:#8dd5ac;--warn:#f4ca72;--bad:#f6a5a0}
*{box-sizing:border-box}body{margin:0;background:var(--bg);color:var(--text);font:15px/1.5 system-ui,sans-serif}main{max-width:1700px;margin:auto;padding:28px}h1{font-size:25px;font-weight:650;margin:0 0 8px}p{margin:8px 0;color:var(--muted)}a{color:var(--accent)}button,select,input{font:inherit}button,select{background:#202d3d;color:var(--text);border:1px solid #4a5d73;border-radius:6px;padding:9px 11px;min-height:42px}button{cursor:pointer}button:hover{background:#2b3c50}select:focus,button:focus,input:focus{outline:2px solid var(--accent);outline-offset:2px}label{color:var(--muted)}.toolbar{display:flex;flex-wrap:wrap;gap:16px;align-items:end;margin:22px 0 12px}.field{display:flex;flex-direction:column;gap:5px}.field:first-child{flex:1;min-width:280px}.field select{max-width:100%}.line{display:flex;gap:14px;align-items:center;flex-wrap:wrap}.transport{padding:15px;background:var(--panel);border:1px solid var(--line);border-radius:10px}.scrub{flex:1;min-width:180px;accent-color:var(--accent)}.time{font-variant-numeric:tabular-nums;min-width:116px}#scenario-note{min-height:48px;font-size:14px;margin:12px 0 20px}.panels{display:grid;grid-template-columns:repeat(3,minmax(0,1fr));gap:16px}.panel{border:1px solid var(--line);border-radius:10px;overflow:hidden;background:var(--panel)}.panel-top{padding:14px;min-height:129px}.panel-top select{width:100%;font-size:13px}.method-id{font:12px ui-monospace,monospace;overflow-wrap:anywhere;color:var(--muted);margin-top:8px}.status{display:inline-block;font-size:12px;border-radius:4px;padding:3px 7px;margin-top:8px;background:#263b33;color:var(--ok)}.status.inexact{background:#443a26;color:var(--warn)}.status.failed{background:#452e32;color:var(--bad)}canvas{display:block;width:100%;height:340px;background:#111a24}.sample-note{min-height:66px;padding:10px 14px;font-size:13px;color:var(--muted);border-bottom:1px solid var(--line)}.sample-note.warning{color:var(--warn)}.metrics{display:grid;grid-template-columns:1fr 1fr;gap:12px;padding:14px}.metric span{display:block;color:var(--muted);font-size:12px}.metric strong{font:15px ui-monospace,monospace}.details{font-size:13px;padding:0 14px 14px;color:var(--muted)}.details p{margin:6px 0}summary{cursor:pointer;color:var(--accent)}pre{white-space:pre-wrap;word-break:break-word;font:11px/1.5 ui-monospace,monospace;max-height:260px;overflow:auto}.foot{margin-top:22px;font-size:13px}.legend{display:inline-flex;gap:16px;flex-wrap:wrap;font-size:13px}.dot{display:inline-block;width:9px;height:9px;border-radius:50%;margin-right:6px;background:#8894a4}.dot.candidate{background:#77b5ff}.checkbox{display:flex;gap:7px;align-items:center;font-size:13px}.checkbox input{accent-color:var(--accent)}@media(max-width:1050px){main{padding:18px}.panels{grid-template-columns:1fr}.panel-top{min-height:0}canvas{height:360px}.sample-note{min-height:0}.metrics{grid-template-columns:repeat(4,1fr)}}@media(max-width:560px){.metrics{grid-template-columns:1fr 1fr}h1{font-size:21px}.field:first-child{min-width:100%}canvas{height:300px}}
</style></head><body><main>
<h1>Constrained XPBD · recorded trajectory comparisons</h1>
<p>Compare motion, numerical convergence, and conservation separately. Every animated point comes from the saved CSV traces.</p>
<div class="toolbar"><label class="field">Environment / initial condition<select id="scenario"></select></label><label class="field">Playback speed<select id="speed"><option value="0.1">0.1×</option><option value="0.25">0.25×</option><option value="0.5">0.5×</option><option value="1" selected>1×</option><option value="2">2×</option><option value="4">4×</option></select></label><label class="field">Camera<select id="camera"><option value="reference">Shared reference bounds</option><option value="selected">Fit selected recorded trajectories</option></select></label></div>
<div class="transport line"><button id="play" type="button">Play</button><button id="reset" type="button">Restart</button><input id="scrub" class="scrub" aria-label="Simulation time" type="range" min="0" max="10" step="0.001" value="0"><output id="time" class="time">0.000 / 10 s</output><label class="checkbox"><input id="ghost" type="checkbox" checked>Reference ghost</label><label class="checkbox"><input id="trails" type="checkbox" checked>2 s trails</label></div>
<p id="scenario-note"></p><div class="panels" id="panels"></div>
<div class="foot"><div class="legend"><span><i class="dot candidate"></i>Selected candidate</span><span><i class="dot"></i>Recorded reference (certification shown above)</span></div>
<p>Playback holds each recorded sample until the next sample. It never extrapolates or fills a failed prefix. Nominal objectives were solved at 120 Hz and recorded at 60 Hz; the rigid stress cases shown here use 30 Hz macro steps and output. Adaptive methods can take additional internal substeps. The first saved frame is after t = 0.</p>
<p>Metric cards report full-run CSV summaries, including any recorded failed prefix. “Completed, inexact” means the horizon was reached but at least one step missed the solver threshold. Objective stationarity and mechanical solver residuals have different definitions. A small energy error does not establish accurate motion. Drawing coordinates use nine significant digits; CSV summary metrics retain their original precision.</p>
<p>* L denotes angular momentum and is imposed only where the objective experiment enables that constraint (the 2D central spring without gravity). All links below are optional; playback is fully self-contained.</p>
<p><a href="objectives/README.md">Objective definitions and results</a> · <a href="adversarial/RESULTS.md">Rigid stress results</a> · <a href="preconditioning/RESULTS.md">Preconditioning results</a> · <a href="preconditioning/figures/convergence_scaling.png">Preconditioning figure</a>. The preconditioning study contains diagnostics, not animation traces.</p><p id="provenance"></p></div>
</main><script id="data" type="application/json">__DATA__</script><script>
'use strict';
const DATA=JSON.parse(document.getElementById('data').textContent), $=id=>document.getElementById(id);
const colors=['#77b5ff','#efb66a','#b5a0fb'];
let scenario=DATA.scenarios[0], time=0, playing=false, lastStamp=null, bounds=null;
const panels=[];
function format(x){if(x===null||x===undefined||!Number.isFinite(x))return 'unavailable';if(x===0)return '0';return Math.abs(x)<.001||Math.abs(x)>=10000?x.toExponential(3):x.toPrecision(4)}
function status(m){if(m.status!=='complete'||m.completed!==m.requested)return ['failed','Stopped / failed'];if(m.inexact>0)return ['inexact','Completed, inexact'];return ['','Completed']}
function methodById(id){return scenario.methods.find(m=>m.id===id)||scenario.methods[0]}
function sampleIndex(frames,t){let lo=0,hi=frames.length;while(lo<hi){const mid=(lo+hi)>>1;if(frames[mid][0]<=t+1e-8)lo=mid+1;else hi=mid}return lo-1}
function points(row,offset){if(!row)return [];const result=[];for(let k=offset;k<row.length;k+=scenario.dim)result.push([row[k],scenario.dim===1?0:row[k+1]]);return result}
function chooseBounds(){let pts=[[0,0]];for(const r of scenario.reference)pts.push(...points(r,1));if($('camera').value==='selected')for(const p of panels)for(const r of p.method.frames)pts.push(...points(r,2));if(pts.length===1)for(const p of panels)if(p.method.frames.length)pts.push(...points(p.method.frames[0],2));let xmin=Infinity,xmax=-Infinity,ymin=Infinity,ymax=-Infinity;for(const [x,y]of pts){if(!Number.isFinite(x)||!Number.isFinite(y))continue;xmin=Math.min(xmin,x);xmax=Math.max(xmax,x);ymin=Math.min(ymin,y);ymax=Math.max(ymax,y)}let width=Math.max(xmax-xmin,.5),height=Math.max(ymax-ymin,.5);if(scenario.dim===1)height=Math.max(height,.8);bounds={cx:(xmin+xmax)/2,cy:(ymin+ymax)/2,width:width*1.25,height:height*1.25};}
function draw(p){const canvas=p.canvas,rect=canvas.getBoundingClientRect(),ratio=window.devicePixelRatio||1;if(canvas.width!==Math.round(rect.width*ratio)||canvas.height!==Math.round(rect.height*ratio)){canvas.width=Math.round(rect.width*ratio);canvas.height=Math.round(rect.height*ratio)}const c=canvas.getContext('2d');c.setTransform(ratio,0,0,ratio,0,0);const w=rect.width,h=rect.height;c.clearRect(0,0,w,h);const scale=Math.min((w-44)/bounds.width,(h-44)/bounds.height);const screen=([x,y])=>[w/2+(x-bounds.cx)*scale,h/2-(y-bounds.cy)*scale];
  c.strokeStyle='#243241';c.lineWidth=1;c.beginPath();const [ox,oy]=screen([0,0]);c.moveTo(0,oy);c.lineTo(w,oy);c.moveTo(ox,0);c.lineTo(ox,h);c.stroke();
  function chain(pts,color,ghost=false){c.strokeStyle=color;c.fillStyle=color;c.lineWidth=ghost?1.5:2.4;c.setLineDash(ghost?[5,5]:[]);let prev=[0,0];for(const pt of pts){const a=screen(prev),b=screen(pt);c.beginPath();c.moveTo(...a);if(!scenario.rigid&&!ghost){const dx=b[0]-a[0],dy=b[1]-a[1],len=Math.hypot(dx,dy);if(len>20){const px=-dy/len*4,py=dx/len*4;for(let j=1;j<12;j++){const s=j/12,amp=j%2?1:-1;c.lineTo(a[0]+dx*s+px*amp,a[1]+dy*s+py*amp)}}}c.lineTo(...b);c.stroke();c.beginPath();c.arc(b[0],b[1],ghost?5:7,0,Math.PI*2);ghost?c.stroke():c.fill();prev=pt}c.setLineDash([])}
  function trail(frames,offset,color){const last=sampleIndex(frames,time);if(last<0)return;const first=Math.max(0,sampleIndex(frames,Math.max(0,time-2)));const count=(frames[last].length-offset)/scenario.dim;c.strokeStyle=color;c.lineWidth=1.3;for(let body=0;body<count;body++){c.beginPath();let started=false;for(let j=first;j<=last;j++){const pt=screen(points(frames[j],offset)[body]);if(!started){c.moveTo(...pt);started=true}else c.lineTo(...pt)}c.stroke()}}
  if($('trails').checked){if($('ghost').checked)trail(scenario.reference,1,'#64728580');trail(p.method.frames,2,p.color+'70')}
  const ri=sampleIndex(scenario.reference,time),i=sampleIndex(p.method.frames,time),frame=i>=0?p.method.frames[i]:null;
  if($('ghost').checked&&ri>=0)chain(points(scenario.reference[ri],1),'#8c99aa',true);
  let outside=false;if(frame){const pts=points(frame,2);outside=pts.some(pt=>{const [x,y]=screen(pt);return x<8||x>w-8||y<8||y>h-8});chain(pts,p.color)}
  c.fillStyle='#aab8c9';c.beginPath();c.arc(ox,oy,4,0,Math.PI*2);c.fill();c.font='12px system-ui';c.fillText('fixed anchor',Math.max(8,Math.min(w-80,ox+9)),Math.max(15,Math.min(h-10,oy-9)));
  const bar=Math.pow(10,Math.floor(Math.log10(Math.max(bounds.width/5,1e-9))));c.strokeStyle='#a5b1c2';c.beginPath();c.moveTo(16,h-19);c.lineTo(16+bar*scale,h-19);c.stroke();c.fillText(format(bar)+' m',16,h-27);
  const fs=p.method.frames,stopped=fs.length&&time>fs[fs.length-1][0]+1/scenario.outputHz*.6&&fs[fs.length-1][0]<scenario.duration-1e-5;let note='';
  if(!fs.length)note='No recorded trajectory: '+(p.method.reason||p.method.traceIssue||'the solve stopped before its first output.');else if(!frame)note='No saved sample yet. First recorded state: '+fs[0][0].toFixed(4)+' s.';else if(stopped)note='Stopped at '+fs[fs.length-1][0].toFixed(4)+' s. Showing the last recorded state; subsequent motion is unavailable.';else note='Recorded state: '+frame[0].toFixed(4)+' s · energy '+format(frame[1])+' J.';
  if(outside)note+=' Candidate leaves the reference view; select “Fit selected recorded trajectories” to inspect it.';p.note.textContent=note;p.note.className='sample-note'+(stopped||outside||!fs.length?' warning':'');
}
function metric(label,value){const box=document.createElement('div');box.className='metric';const l=document.createElement('span');l.textContent=label;const v=document.createElement('strong');v.textContent=value;box.append(l,v);return box}
function updatePanel(p,id){p.method=methodById(id);p.select.value=p.method.id;const m=p.method,[cls,label]=status(m);p.id.textContent=m.id;p.badge.className='status '+cls;p.badge.textContent=label;p.metrics.replaceChildren(metric('Max relative energy error',format(m.energyError)),metric('Position RMS (m)',m.referenceCertified?format(m.positionError):'uncertified ref.'),metric(scenario.kind==='objective'?'Max KKT stationarity':'Max solver residual',format(m.residual)),metric('Inexact solver steps',m.inexact===null?'unavailable':String(m.inexact)));p.details.replaceChildren();const info=document.createElement('p');info.textContent=(m.completed??'?')+' / '+(m.requested??'?')+' steps recorded in summary. Reference: '+(m.referenceCertified?'certified.':'uncertified; position errors should not be ranked.')+(m.reason?' Stop reason: '+m.reason+'.':'')+(m.traceIssue?' '+m.traceIssue:'');p.details.append(info);if(scenario.kind==='objective'){const feasibility=document.createElement('p');feasibility.textContent='Max constraint balance: '+format(m.balance)+'.';p.details.append(feasibility)}const detail=document.createElement('details'),summary=document.createElement('summary'),pre=document.createElement('pre');summary.textContent='Exact summary row and source files';pre.textContent=JSON.stringify(m.raw,null,2);detail.append(summary,pre);const links=document.createElement('p');for(const [label,path]of [['Summary CSV',m.summary],['Trajectory CSV',m.trace]]){const a=document.createElement('a');a.textContent=label;a.href=path;links.append(a,document.createTextNode('  '))}detail.append(links);p.details.append(detail)}
for(let n=0;n<3;n++){const panel=document.createElement('section');panel.className='panel';const top=document.createElement('div');top.className='panel-top';const select=document.createElement('select');select.setAttribute('aria-label','Candidate '+(n+1)+' method');const id=document.createElement('div');id.className='method-id';const badge=document.createElement('span');const canvas=document.createElement('canvas');canvas.setAttribute('aria-label','Candidate '+(n+1)+' recorded trajectory');const note=document.createElement('div'),metrics=document.createElement('div'),details=document.createElement('div');metrics.className='metrics';details.className='details';top.append(select,id,badge);panel.append(top,canvas,note,metrics,details);$('panels').append(panel);const p={select,id,badge,canvas,note,metrics,details,color:colors[n]};panels.push(p);select.addEventListener('change',()=>{updatePanel(p,select.value);chooseBounds();render()})}
let group=null;for(const s of DATA.scenarios){if(!group||group.label!==s.group){group=document.createElement('optgroup');group.label=s.group;$('scenario').append(group)}const option=document.createElement('option');option.value=s.id;option.textContent=s.title+' · '+s.duration+' s';group.append(option)}
function setScenario(id){scenario=DATA.scenarios.find(s=>s.id===id)||DATA.scenarios[0];$('scenario').value=scenario.id;time=0;playing=false;$('play').textContent='Play';$('scrub').max=scenario.duration;const defaults=scenario.kind==='objective'?['original_midpoint_h2','quarter_b1','ls_dg_E']:['dg64','rattle','angle_dg'];panels.forEach((p,i)=>{p.select.replaceChildren();for(const m of scenario.methods){const opt=document.createElement('option');opt.value=m.id;opt.textContent=m.label+(status(m)[0]?' ['+status(m)[1]+']':'');p.select.append(opt)}updatePanel(p,defaults[i])});const refEnd=scenario.reference.length?scenario.reference[scenario.reference.length-1][0]:0;$('scenario-note').textContent=scenario.methods.length+' methods · '+scenario.hz+' Hz '+(scenario.kind==='objective'?'solver':'macro step')+' · '+scenario.outputHz+' Hz saved frames · '+scenario.duration+' s requested. Reference '+(scenario.referenceCertified?'certified':'uncertified')+', recorded through '+refEnd.toFixed(3)+' s. All panels share one camera and simulation clock.';chooseBounds();render()}
function render(){$('scrub').value=time;$('time').textContent=time.toFixed(3)+' / '+scenario.duration+' s';for(const p of panels)draw(p)}
$('scenario').addEventListener('change',()=>setScenario($('scenario').value));$('camera').addEventListener('change',()=>{chooseBounds();render()});$('ghost').addEventListener('change',render);$('trails').addEventListener('change',render);$('play').addEventListener('click',()=>{if(time>=scenario.duration)time=0;playing=!playing;lastStamp=null;$('play').textContent=playing?'Pause':'Play';render()});$('reset').addEventListener('click',()=>{time=0;lastStamp=null;render()});$('scrub').addEventListener('input',()=>{time=Number($('scrub').value);lastStamp=null;render()});window.addEventListener('resize',render);
function tick(stamp){if(playing&&lastStamp!==null){time=Math.min(scenario.duration,time+Math.min((stamp-lastStamp)/1000,.2)*Number($('speed').value));if(time>=scenario.duration){playing=false;$('play').textContent='Play'}render()}lastStamp=stamp;requestAnimationFrame(tick)}
$('provenance').textContent='Generated from local recorded data at '+DATA.built+'. '+DATA.scenarios.length+' scenarios; '+DATA.scenarios.reduce((n,s)=>n+s.methods.length,0)+' recorded method/scenario entries.';
setScenario(DATA.scenarios.find(s=>s.case==='radial_orbit'&&s.kind==='objective').id);requestAnimationFrame(tick);
</script></body></html>'''


if __name__ == "__main__":
    data = build_data()
    encoded = json.dumps(data, ensure_ascii=False, separators=(",", ":"), allow_nan=False).replace("</", "<\\/")
    target = ROOT / "viewer.html"
    target.write_text(HTML.replace("__DATA__", encoded))
    counts = [(s["id"], len(s["methods"]), sum(len(m["frames"]) for m in s["methods"])) for s in data["scenarios"]]
    print(json.dumps({"output": str(target), "bytes": target.stat().st_size, "scenarios": counts}, indent=2))
