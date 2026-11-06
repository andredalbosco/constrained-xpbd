#!/usr/bin/env python3
"""Embed finalized trace CSVs in an offline motion-first comparison viewer.

No web dependencies. By default, include trace/long tiers when they exist and
otherwise fall back to smoke traces. Regenerate after more runs finish.
"""
import argparse
import csv
import datetime
import json
import math
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parent


def number(value, default=0.0):
    try:
        x = float(value)
        return x if math.isfinite(x) else default
    except (ValueError, TypeError):
        return default


def metadata(row):
    strings = {'id', 'tier', 'case', 'family', 'work', 'treatment', 'kinematics',
               'schedule', 'curvature', 'status', 'reason'}
    return {k: (v if k in strings else number(v)) for k, v in row.items() if k}


def read_trace(path, limit):
    with path.open() as f:
        reader = csv.DictReader(f)
        fields = reader.fieldnames or []
        dims = len([k for k in fields if k.startswith('x') and k[1:].isdigit()])
        if not dims:
            return [], 0, 0
        cols = ['time', 'energy', 'reference_energy']
        cols += [f'x{i}' for i in range(dims)] + [f'ref_x{i}' for i in range(dims)]
        samples = []
        for row in reader:
            try:
                values = [float(row[k]) for k in cols]
            except (KeyError, ValueError, TypeError):
                continue  # A partially written last CSV row is not a frame.
            if all(math.isfinite(v) for v in values):
                # Preserve energy precision; coordinates need much less than
                # the source's 16 digits for a pixel-based animation.
                samples.append(values[:3] + [float(format(x, '.10g')) for x in values[3:]])
    original = len(samples)
    if original > limit:
        indices = {round(i * (original - 1) / (limit - 1)) for i in range(limit)}
        # Preserve extremes so a short energy spike or runaway is not hidden.
        for col in range(1, 3 + dims):
            indices.add(min(range(original), key=lambda i: samples[i][col]))
            indices.add(max(range(original), key=lambda i: samples[i][col]))
        samples = [samples[i] for i in sorted(indices)]
    return samples, dims, original


def select_representatives(entries, limit):
    groups = defaultdict(list)
    for m, path in entries:
        group = tuple(m.get(k) for k in ('case', 'tier', 'hz', 'duration', 'stiffness_factor', 'mass_ratio', 'perturbation'))
        groups[group].append((m, path))
    chosen = []
    for entries in groups.values():
        keep = []
        def score(item):
            m = item[0]
            return (m['status'] != 'completed', m.get('rms_position', 1e6), abs(m.get('iterations', 0)-16))
        def take(predicate, key=score):
            available = [p for p in entries if p not in keep and predicate(p[0])]
            if available and len(keep) < limit:
                keep.append(min(available, key=key))
        take(lambda m: m['family'] == 'xpbd')
        take(lambda m: m['family'] == 'force_momentum' and m['work'] == 'dg_average')
        take(lambda m: m['family'] == 'independent' and m['treatment'] == 'global' and m['kinematics'] == 'soft_endpoint')
        take(lambda m: m['family'] == 'independent' and m['treatment'] == 'global' and m['kinematics'] == 'soft_midpoint')
        # Preserve stopped traces and a frozen/excess-motion completed example.
        take(lambda m: m['status'] != 'completed', key=lambda p: (p[0].get('completed', 0) == 0, -p[0].get('max_energy_rel', 0)))
        take(lambda m: m['status'] == 'completed' and (m.get('path_ratio', 1) < .05 or m.get('path_ratio', 1) > 3), key=lambda p: p[0].get('path_ratio', 1))
        take(lambda m: m['family'] == 'rattle')
        take(lambda m: m['family'] == 'force_momentum' and m['work'] == 'force_average')
        take(lambda m: m['family'] == 'force_momentum' and m['work'] == 'separate_average')
        take(lambda m: m['family'] == 'independent' and m['work'] == 'start')
        take(lambda m: m['family'] == 'force_momentum' and m['work'] == 'dg_end')
        take(lambda m: m['family'] == 'EL')
        for item in sorted(entries, key=score):
            if len(keep) >= limit:
                break
            if item not in keep:
                keep.append(item)
        chosen.extend(keep)
    return chosen


def collect(results, max_frames, per_setup):
    summaries = sorted(results.rglob('summary.csv'))
    preferred = [p for p in summaries if p.parent.name in {'traces', 'long'}]
    if preferred:
        summaries = preferred
    runs, skipped, entries = [], 0, []
    for summary in summaries:
        with summary.open() as f:
            for row in csv.DictReader(f):
                if not row.get('id') or not row.get('status'):
                    continue
                trace = summary.parent / 'traces' / (row['id'] + '.csv')
                if not trace.exists():
                    continue
                meta = metadata(row)
                meta['source'] = str(trace.relative_to(ROOT)) if trace.is_relative_to(ROOT) else str(trace)
                meta['uid'] = str(summary.parent.relative_to(results)) + '/' + row['id']
                entries.append((meta, trace))
    for meta, trace in select_representatives(entries, per_setup):
        try:
            samples, dims, original = read_trace(trace, max_frames)
        except OSError:
            skipped += 1
            continue
        if not samples:
            skipped += 1
            continue
        meta['recorded_end'] = samples[-1][0]
        meta['original_frames'] = original
        runs.append({'m': meta, 'd': dims, 's': samples})
    return {'built': datetime.datetime.now(datetime.timezone.utc).isoformat(),
            'runs': runs, 'catalog': [m for m, _ in entries], 'skipped': skipped,
            'note': 'Finalized summary rows only; trajectories end at their last recorded frame.'}


HTML = r'''<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>Constrained XPBD — motion comparisons</title>
<link rel="icon" href="data:image/svg+xml,%3Csvg xmlns=%27http://www.w3.org/2000/svg%27 viewBox=%270 0 32 32%27%3E%3Crect width=%2732%27 height=%2732%27 rx=%276%27 fill=%27%23142c38%27/%3E%3Cpath d=%27M16 5v7l7 12%27 fill=%27none%27 stroke=%27%23e2a17d%27 stroke-width=%272%27/%3E%3Ccircle cx=%2723%27 cy=%2724%27 r=%274%27 fill=%27%23e2a17d%27/%3E%3C/svg%3E">
<style>
:root{color-scheme:light;--ink:#172c38;--muted:#657783;--line:#dbe3e6;--paper:#f2f5f5;--accent:#bd553c;--ref:#788b97}
*{box-sizing:border-box}body{margin:0;background:var(--paper);font:14px/1.45 system-ui,-apple-system,sans-serif;color:var(--ink)}
header{padding:28px max(24px,calc((100vw - 1680px)/2));background:#142c38;color:#f6fafb;border-bottom:5px solid #d98b62}h1{font-size:27px;letter-spacing:-.6px;margin:0 0 6px;font-weight:650}header p{margin:0;color:#c2d0d6;font-size:14px}.eyebrow{font:11px/1.4 ui-monospace,monospace;text-transform:uppercase;letter-spacing:1.8px;color:#dfae91;margin-bottom:10px}
main{max-width:1730px;margin:auto;padding:22px 24px 30px}.controls{display:flex;gap:14px;flex-wrap:wrap;align-items:end;background:white;padding:17px;border:1px solid var(--line);border-radius:12px}.field{display:grid;gap:5px;min-width:160px}.field.wide{flex:1;min-width:250px}label,.small-label{font-size:11px;font-weight:650;letter-spacing:.6px;text-transform:uppercase;color:var(--muted)}select,button{font:inherit;border:1px solid #cdd7dd;border-radius:7px;background:#fff;color:var(--ink);padding:9px 10px;min-height:39px}select{max-width:100%;cursor:pointer}button{cursor:pointer;font-weight:600}button.primary{background:#193c4b;border-color:#193c4b;color:white;min-width:94px}button:hover{filter:brightness(.97)}button:disabled{opacity:.4;cursor:default}
.timeline{display:grid;grid-template-columns:auto auto 1fr auto auto;gap:13px;align-items:center;padding:18px 1px 13px}.timeline input{width:100%;accent-color:#bd553c}.time{font:13px ui-monospace,monospace;min-width:100px;font-variant-numeric:tabular-nums}.hint{color:var(--muted);font-size:12px}.info{display:flex;justify-content:space-between;gap:12px;margin:0 0 17px;flex-wrap:wrap}.panes{display:grid;grid-template-columns:repeat(var(--pane-count,4),minmax(0,1fr));gap:14px;align-items:start}.pane{background:white;border:1px solid var(--line);border-radius:12px;overflow:hidden}.pane-head{padding:15px 14px 12px;border-bottom:1px solid #e8edef;min-height:166px}.pane-head h2{margin:0 0 8px;font-size:15px}.pane-head select{width:100%;font-size:12px;padding:7px;min-height:36px}.method-description{font-size:11px;line-height:1.6;margin-top:9px;color:#50636e;min-height:36px;overflow-wrap:anywhere}.status{display:flex;align-items:center;gap:6px;font-size:11px;margin-top:7px;color:#566d78}.dot{width:6px;height:6px;border-radius:50%;background:#6c8898;display:inline-block}.dot.stopped{background:#bd553c}.view{position:relative;background:linear-gradient(#fbfcfc,#f2f6f7)}canvas{display:block;width:100%}.world{aspect-ratio:1/1.04}.frame-status{position:absolute;top:10px;left:11px;right:11px;font-size:11px;pointer-events:none;color:#6f7f88}.frame-status.ended{background:#fff0e6ed;color:#923d28;border:1px solid #ecc6af;border-radius:6px;padding:8px}.pane-body{padding:12px 14px 15px}.observation{font-size:12px;min-height:36px;margin-bottom:12px;color:#344f5d}.observation.alert{color:#9b442f}.metrics{display:grid;grid-template-columns:1fr 1fr;gap:9px;margin-bottom:14px}.metric{padding:8px 9px;background:#f3f6f7;border-radius:6px}.metric b{display:block;font:15px ui-monospace,monospace;font-variant-numeric:tabular-nums;margin-top:3px}.metric span{font-size:10px;color:var(--muted)}.chart-title{display:flex;justify-content:space-between;align-items:center;margin:11px 0 4px;font-size:11px;font-weight:600}.chart{height:105px}.path{height:130px}.legend{font-size:10px;color:var(--muted)}details{margin-top:12px;font-size:11px;color:#617583}summary{cursor:pointer}details p{margin:8px 0;overflow-wrap:anywhere}.source{color:#586f7d}.empty{padding:65px 25px;text-align:center;border:1px dashed #b5c4cc;border-radius:12px;background:white}.empty h2{font-size:22px}.footer{margin-top:21px;padding-top:17px;border-top:1px solid var(--line);font-size:12px;color:var(--muted)}.footer a{color:#365a70}.hidden{display:none!important}.badge{padding:3px 7px;border:1px solid #d7e1e6;border-radius:4px;font:10px ui-monospace,monospace}.reference .pane-head{background:#f5f8fa}.reference .pane-head h2{margin-top:2px}.reference .method-description{margin-top:12px}.reference .pane-head select{visibility:hidden}.section-title{font-size:12px;font-weight:600;color:#617583;margin-bottom:10px}
@media(max-width:1250px){.panes{grid-template-columns:repeat(2,minmax(0,1fr))}.world{aspect-ratio:1.4/1}.pane-head{min-height:145px}}@media(max-width:690px){main{padding:15px}.panes{grid-template-columns:1fr}.timeline{grid-template-columns:auto auto 1fr}.timeline input{grid-column:1/-1;grid-row:2}.timeline .speed-label{display:none}.timeline .time{text-align:right}.controls{gap:10px}.field{flex:1;min-width:130px}header{padding:22px 18px}.world{aspect-ratio:1.3/1}}
.pane-head{height:250px}.method-filter{width:100%;font:11px system-ui;border:1px solid #dce4e8;border-radius:5px;padding:6px 8px;margin:0 0 6px;color:#475f6d}.reference .method-filter{visibility:hidden}.catalog{background:white;border:1px solid var(--line);border-radius:9px;padding:14px;margin-top:18px}.table-wrap{overflow:auto}table{width:100%;border-collapse:collapse;margin-top:12px;font-size:11px;white-space:nowrap}th,td{text-align:left;padding:7px 9px;border-bottom:1px solid #e6edef}th{color:#6b7f8c}td a{color:#305f76}

/* Keep the animation itself above the fold; detailed composition stays below. */
header{padding:15px max(20px,calc((100vw - 1680px)/2));border-bottom-width:3px}h1{font-size:24px;margin-bottom:3px}.eyebrow{margin-bottom:5px;font-size:10px}header p{font-size:12px}
main{padding:12px 20px 24px}.controls{padding:10px 12px;gap:12px}.field{gap:3px}select,button{padding:6px 9px;min-height:34px}.timeline{padding:10px 1px 9px;gap:10px}.info{margin-bottom:10px;font-size:11px}.hint{font-size:11px}
.pane-head{height:140px;min-height:140px;padding:9px 11px}.pane-head h2{font-size:13px;margin-bottom:5px}.pane-head select{font-size:11px;padding:6px;min-height:31px}.method-filter{padding:4px 6px;margin-bottom:5px;font-size:10px}.pane-head .method-description{display:none}.status{font-size:10px;margin-top:7px;line-height:1.3}.reference .method-filter,.reference .pane-head select{display:none}.reference .pane-head .method-description{display:block;margin-top:10px;font-size:11px;min-height:0;line-height:1.5}.reference .pane-head h2{margin-top:0}
</style></head><body>
<header><div class="eyebrow">Constrained XPBD / composition experiments</div><h1>Watch the motion. Then read the energy.</h1><p>Candidate and reference trajectories share initial conditions, a camera, and a clock.</p></header>
<main>
<div id="empty" class="empty hidden"><h2>No finalized traces yet</h2><p>Run the trace experiments, then regenerate this file with <code>python3 build_viewer.py</code>.</p><p>Partial runs appear only after their summary row is written.</p></div>
<div id="app"><div class="controls">
<div class="field"><label for="environment">Environment</label><select id="environment"></select></div>
<div class="field"><label for="initial">Initial condition</label><select id="initial"></select></div>
<div class="field wide"><label for="setup">Recorded run setup</label><select id="setup"></select></div>
<div class="field"><label for="paneCount">Candidate panes</label><select id="paneCount"><option value="3">Three comparisons</option><option value="2">Two comparisons</option><option value="1">One comparison</option></select></div>
</div>
<div class="timeline"><button id="play" class="primary">Play</button><button id="reset" aria-label="Return to start">↺</button><input id="scrub" aria-label="Simulation time" type="range" min="0" max="1" step="0.001" value="0"><span id="time" class="time">0.00 / 0.00 s</span><select id="speed" aria-label="Playback speed"><option value="0.1">0.1×</option><option value="0.25">0.25×</option><option value="0.5">0.5×</option><option value="1" selected>1×</option><option value="2">2×</option><option value="4">4×</option></select></div>
<div class="info"><span id="recording" class="hint"></span><span class="hint">Solid: candidate · dashed gray: reference · no extrapolation past a recorded endpoint</span></div>
<div class="panes" id="panes"></div>
<details class="catalog"><summary id="catalogSummary">All recorded runs in this setup</summary><div class="table-wrap"><table><thead><tr><th>Family</th><th>Work</th><th>Treatment</th><th>Kinematics</th><th>Schedule</th><th>Tangent</th><th>Status</th><th>Completed</th><th>Path ratio</th><th>RMS x</th><th>Animation / source</th></tr></thead><tbody id="catalogBody"></tbody></table></div></details>
</div><div class="footer"><p>Energy drift is (E − E₀)/(E₀ − E_min). Exact bookkeeping alone does not establish correct motion. A stopped or frozen trajectory remains visible as an observation.</p><p id="buildInfo"></p><a href="audits_notes.md">Work-identity and EL audit notes</a> · <a href="results/">Raw result files</a></div>
</main>
<script id="embedded-data" type="application/json">__DATA__</script>
<script>
'use strict';
const DATA=JSON.parse(document.getElementById('embedded-data').textContent);
const $=id=>document.getElementById(id), ESC=s=>String(s??'').replace(/[&<>"']/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
const ENV=[['scalar','Scalar spring'],['radial','Planar spring'],['rod','Rigid pendulum'],['elastic','Elastic pendulum'],['chain','Two-spring pendulum']];
const COLORS=['#bc563c','#247d8e','#8871aa'];
const ALL=DATA.runs, byID=new Map(ALL.map(r=>[r.m.uid,r]));
let selected=[],pool=[],reference=null,time=0,duration=0,playing=false,lastWall=0,camera=null,panes=[];
const envOf=r=>r.m.case.split('_')[0];
const setupOf=r=>[r.m.case,r.m.hz,r.m.duration,r.m.stiffness_factor,r.m.mass_ratio,r.m.perturbation,r.m.tier].join('|');
const fnum=(v,d=3)=>!Number.isFinite(v)?'—':Math.abs(v)>0&&(Math.abs(v)<.001||Math.abs(v)>10000)?v.toExponential(2):v.toFixed(d);
const nice=s=>String(s).replaceAll('_',' ');
function options(el,arr,value){el.innerHTML=arr.map(([v,t])=>`<option value="${ESC(v)}">${ESC(t)}</option>`).join('');if(arr.some(a=>a[0]===value))el.value=value;}
function methodName(r){let m=r.m,fam={independent:'Independent',force_momentum:'Force',xpbd:'XPBD',fepr_projection:'Projection'}[m.family]||m.family,kin={soft_endpoint:'soft endpoint',soft_midpoint:'soft midpoint',hard_endpoint:'endpoint',hard_midpoint:'midpoint',fixed_target:'fixed target'}[m.kinematics]||nice(m.kinematics),treat={global:'global E',prescribed:'fixed W',penalized:'soft W',impulse_closed:'impulses',none:''}[m.treatment]??nice(m.treatment);return `${fam} · ${treat? treat+' · ':''}${kin} · ${nice(m.work)} · ${nice(m.schedule)} ${m.iterations}i${m.tangency?' · tangent':' · no tangent'}${m.status==='completed'?'':' · STOPPED'}`;}
function fullDescription(r){let m=r.m;return `Work: ${nice(m.work)} · ${nice(m.treatment)}<br>Kinematics: ${nice(m.kinematics)} · schedule: ${nice(m.schedule)}<br>Tangency: ${m.tangency?'enforced':'not enforced'} · ${m.iterations} iterations`;}
function makePane(index){let el=document.createElement('article');el.className='pane'+(index<0?' reference':'');el.innerHTML=`<div class="pane-head"><h2>${index<0?'Reference trajectory':'Candidate '+(index+1)}</h2><input class="method-filter" placeholder="Filter loaded methods…" aria-label="Filter loaded methods"><select class="method" aria-label="${index<0?'Reference':'Candidate '+(index+1)+' method'}"></select><div class="method-description"></div><div class="status"></div></div><div class="view"><canvas class="world"></canvas><div class="frame-status"></div></div><div class="pane-body"><div class="observation"></div><div class="metrics"><div class="metric"><span>Recorded motion / reference</span><b class="path-metric"></b></div><div class="metric"><span>RMS position error</span><b class="error-metric"></b></div></div><div class="chart-title"><span>Energy drift</span><span class="energy-now"></span></div><canvas class="chart energy"></canvas><div class="chart-title"><span class="path-title">Position trajectory</span><span class="legend">shared coordinates</span></div><canvas class="chart path"></canvas><details><summary>Run details and source</summary><p class="details"></p></details></div>`;
 $('panes').appendChild(el);let p={el,index,world:el.querySelector('.world'),energy:el.querySelector('.energy'),path:el.querySelector('.path'),selector:el.querySelector('.method'),run:null};p.selector.addEventListener('change',()=>{selected[index]=byID.get(p.selector.value);refreshPanes();});el.querySelector('.method-filter').addEventListener('input',e=>{let q=e.target.value.toLowerCase(),items=pool.filter(r=>(methodName(r)+' '+r.m.treatment+' '+r.m.kinematics).toLowerCase().includes(q));let current=selected[index];if(current&&!items.includes(current))items.unshift(current);options(p.selector,items.map(r=>[r.m.uid,methodName(r)]),selected[index]?.m.uid);});return p;}
function updateInitial(){let env=$('environment').value,cases=[...new Set(ALL.filter(r=>envOf(r)===env).map(r=>r.m.case))].sort();let prefer={scalar:'scalar_release',radial:'radial_orbit',rod:'rod_release',elastic:'elastic_mixed',chain:'chain_mixed'}[env];options($('initial'),cases.map(c=>[c,nice(c.replace(env+'_',''))]),prefer);updateSetups();}
function updateSetups(){let c=$('initial').value,groups=new Map();ALL.filter(r=>r.m.case===c).forEach(r=>groups.set(setupOf(r),r));let entries=[...groups].sort((a,b)=>{let priority=r=>(r.m.tier==='traces'?0:r.m.tier==='long'?1:2)*1000+Math.abs(r.m.hz-120);return priority(a[1])-priority(b[1]);});options($('setup'),entries.map(([key,r])=>[key,`${r.m.duration}s · ${r.m.hz} Hz · k ×${r.m.stiffness_factor} · mass ratio ${r.m.mass_ratio} · amplitude ×${r.m.perturbation} · ${r.m.tier}`]));updatePool();}
function rank(r){let m=r.m;return (m.status==='completed'?0:100000)+(m.path_ratio<.05?10000:0)+Math.min(1000,m.rms_position??1000)+Math.abs((m.iterations||0)-16)*.0001;}
function choose(predicate,exclude=[]){return pool.filter(r=>predicate(r.m)&&!exclude.includes(r)).sort((a,b)=>rank(a)-rank(b))[0];}
function updatePool(){playing=false;time=0;pool=ALL.filter(r=>setupOf(r)===$('setup').value).sort((a,b)=>methodName(a).localeCompare(methodName(b)));reference=[...pool].sort((a,b)=>b.m.recorded_end-a.m.recorded_end)[0];
 selected[0]=choose(m=>m.family==='xpbd')||pool[0];selected[1]=choose(m=>m.family==='force_momentum'&&m.work==='dg_average')||choose(m=>m.family==='force_momentum',[selected[0]])||pool[1]||pool[0];selected[2]=choose(m=>m.family==='independent'&&m.treatment==='global'&&m.kinematics==='soft_endpoint')||choose(m=>m.family==='independent')||pool[2]||pool[0];
 panes.filter(p=>p.index>=0).forEach(p=>{p.el.querySelector('.method-filter').value='';options(p.selector,pool.map(r=>[r.m.uid,methodName(r)]),selected[p.index]?.m.uid);});camera=makeCamera(reference);refreshPanes();}
function observe(r,isRef){if(isRef)return r.m.reference_valid?'Reference refinement check passed.':'Reference refinement check did not pass; compare cautiously.';let m=r.m,notes=[];if(m.status!=='completed')notes.push(`${nice(m.status)}: ${nice(m.reason)}.`);if(m.path_ratio<.05)notes.push(`Very little motion: ${fnum(100*m.path_ratio,1)}% of the reference path.`);else if(m.path_ratio>3)notes.push(`Excess motion: ${fnum(m.path_ratio,1)}× the reference path.`);else if(m.final_energy_rel>.5)notes.push('Substantial energy gain.');else if(m.final_energy_rel<-.5)notes.push('Substantial energy loss.');else notes.push(`Reference path ratio ${fnum(m.path_ratio,2)}×; inspect phase and shape.`);return notes.join(' ');}
function refreshPanes(){if(!reference)return;duration=reference.m.recorded_end;let count=+$('paneCount').value;$('panes').style.setProperty('--pane-count',count+1);panes.forEach(p=>{let isRef=p.index<0;p.run=isRef?reference:selected[p.index];p.el.classList.toggle('hidden',!isRef&&p.index>=count);let r=p.run,m=r.m;p.el.querySelector('.method-description').innerHTML=isRef?`Refined RK4 trajectory<br>Same initial state and physical parameters<br>Refinement: ${m.reference_valid?'passed':'not certified'}`:fullDescription(r);p.el.querySelector('.status').innerHTML=`<span class="dot ${!isRef&&m.status!=='completed'?'stopped':''}"></span>${isRef?'Reference frames':nice(m.status)} · recorded ${fnum(m.recorded_end,2)} / ${fnum(m.duration,2)} s`;
 let obs=p.el.querySelector('.observation');obs.textContent=observe(r,isRef);obs.classList.toggle('alert',!isRef&&(m.status!=='completed'||m.path_ratio<.05||m.path_ratio>3));p.el.querySelector('.path-metric').textContent=isRef?'1.000×':fnum(m.path_ratio,3)+'×';p.el.querySelector('.error-metric').textContent=isRef?'—':fnum(m.rms_position,4);p.el.querySelector('.path-title').textContent=r.d===1?'Position versus time':'Position trajectory';
 p.el.querySelector('.details').innerHTML=`<b>${ESC(m.family)}</b><br>${fullDescription(r)}<br>${ESC(m.id)}<br>β = ${m.beta}; γ = ${m.gamma}; curvature ${ESC(m.curvature)}<br>rₓ ${m.keep_rx?'retained':'dropped'}; ${m.relinearize?'refreshed':'frozen'} linearization; ${m.reverse?'reverse':'forward'} order<br>Solver completed ${m.completed}/${m.steps} steps; recorded through ${fnum(m.recorded_end,4)} s.<br>Max normalized energy error: ${fnum(m.max_energy_rel)}<br>Reference refinement Δx: ${fnum(m.reference_dx)}; Δv: ${fnum(m.reference_dv)}<br>${m.original_frames} source frames; ${r.s.length} embedded frames, with extremes preserved.<br><a class="source" href="${ESC(m.source)}">Raw trace CSV</a>`;});
 $('scrub').max=String(Math.max(duration,.001));time=Math.min(time,duration);$('recording').textContent=`${pool.length} embedded representative methods · shared camera from reference · playback ends at ${fnum(duration,2)} s${duration+1e-6<reference.m.duration?' (partial recording)':''}`;updateCatalog();drawAll();}
function updateCatalog(){let rows=(DATA.catalog||[]).filter(m=>setupOf({m})===$('setup').value);$('catalogSummary').textContent=`All ${rows.length} recorded runs in this setup (${pool.length} animations embedded)`;$('catalogBody').innerHTML=rows.map(m=>`<tr><td>${ESC(m.family)}</td><td>${ESC(m.work)}</td><td>${ESC(m.treatment)}</td><td>${ESC(m.kinematics)}</td><td>${ESC(m.schedule)}</td><td>${m.tangency?'yes':'no'}</td><td>${ESC(m.status)}</td><td>${fnum(m.completed/m.hz,2)} / ${m.duration}s</td><td>${fnum(m.path_ratio,2)}×</td><td>${fnum(m.rms_position,4)}</td><td>${byID.has(m.uid)?'embedded · ':''}<a href="${ESC(m.source)}">CSV</a></td></tr>`).join('');}
function sample(r,t,isRef=false){let s=r.s;if(t<=s[0][0])return s[0];if(t>=s[s.length-1][0])return s[s.length-1];let lo=0,hi=s.length-1;while(hi-lo>1){let mid=(lo+hi)>>1;if(s[mid][0]<=t)lo=mid;else hi=mid;}let a=s[lo],b=s[hi],w=(t-a[0])/(b[0]-a[0]);return a.map((v,i)=>v+(b[i]-v)*w);}
function fit(canvas){let ratio=Math.min(devicePixelRatio||1,2),w=canvas.clientWidth,h=canvas.clientHeight;let iw=Math.round(w*ratio),ih=Math.round(h*ratio);if(canvas.width!==iw||canvas.height!==ih){canvas.width=iw;canvas.height=ih;}let ctx=canvas.getContext('2d');ctx.setTransform(ratio,0,0,ratio,0,0);ctx.clearRect(0,0,w,h);return {ctx,w,h};}
function makeCamera(r){let xs=[0],ys=[0];r.s.forEach(a=>{for(let i=0;i<r.d;i+=r.d===1?1:2){xs.push(a[3+r.d+i]);ys.push(r.d===1?0:a[3+r.d+i+1]);}});let xmin=Math.min(...xs),xmax=Math.max(...xs),ymin=Math.min(...ys),ymax=Math.max(...ys);let span=Math.max(xmax-xmin,ymax-ymin, r.d===1?1.6:2);return {cx:(xmin+xmax)/2,cy:(ymin+ymax)/2,span:span*1.32};}
function transform(w,h){let scale=Math.min(w,h-35)/camera.span;return {scale,to:(x,y)=>[w/2+(x-camera.cx)*scale,h/2+12-(y-camera.cy)*scale]};}
function line(ctx,a,b,color,width=1,dash=[]){ctx.beginPath();ctx.strokeStyle=color;ctx.lineWidth=width;ctx.setLineDash(dash);ctx.moveTo(a[0],a[1]);ctx.lineTo(b[0],b[1]);ctx.stroke();ctx.setLineDash([]);}
function disk(ctx,p,r,color){ctx.beginPath();ctx.fillStyle=color;ctx.arc(p[0],p[1],r,0,Math.PI*2);ctx.fill();}
function spring(ctx,a,b,color,rigid,dash=[]){if(rigid){line(ctx,a,b,color,2.5,dash);return;}let dx=b[0]-a[0],dy=b[1]-a[1],len=Math.hypot(dx,dy),nx=-dy/(len||1),ny=dx/(len||1);ctx.beginPath();ctx.strokeStyle=color;ctx.lineWidth=1.7;ctx.setLineDash(dash);ctx.moveTo(...a);for(let k=1;k<18;k++){let t=k/18,amp=k<2||k>16?0:(k%2?4:-4);ctx.lineTo(a[0]+dx*t+nx*amp,a[1]+dy*t+ny*amp);}ctx.lineTo(...b);ctx.stroke();ctx.setLineDash([]);}
function world(p){let {ctx,w,h}=fit(p.world),r=p.run,isRef=p.index<0,a=sample(r,time),T=transform(w,h),origin=T.to(0,0);ctx.font='10px system-ui';ctx.fillStyle='#91a0a8';
 let grid=camera.span>4?1:.5;for(let x=Math.ceil((camera.cx-camera.span/2)/grid)*grid;x<camera.cx+camera.span/2;x+=grid){let xx=T.to(x,0)[0];line(ctx,[xx,34],[xx,h],'#e6ecef');}for(let y=Math.ceil((camera.cy-camera.span/2)/grid)*grid;y<camera.cy+camera.span/2;y+=grid){let yy=T.to(0,y)[1];line(ctx,[0,yy],[w,yy],'#e6ecef');}
 if(r.d===1){line(ctx,T.to(-.1,0),T.to(camera.cx+camera.span/2,0),'#b7c6ce');let eq=T.to(1,0);line(ctx,[eq[0],eq[1]-16],[eq[0],eq[1]+16],'#b5c2c9',1,[3,3]);ctx.fillText('rest length',eq[0]-24,eq[1]+31);}
 disk(ctx,origin,4,'#344f5d');ctx.fillStyle='#6e818d';ctx.fillText('fixed',origin[0]+7,origin[1]-7);
 const positions=(row,ref)=>{let offset=3+(ref?r.d:0),out=[];for(let i=0;i<r.d;i+=r.d===1?1:2)out.push(T.to(row[offset+i],r.d===1?0:row[offset+i+1]));return out;};
 if(!isRef){let prev=origin;for(let point of positions(a,true)){spring(ctx,prev,point,'#9aabb5',r.m.case.startsWith('rod'),[4,4]);disk(ctx,point,5,'#acbac2');prev=point;}}
 let color=isRef?'#6b8391':COLORS[p.index],prev=origin,outside=false;for(let point of positions(a,isRef)){spring(ctx,prev,point,color,r.m.case.startsWith('rod'));let clipped=[Math.max(9,Math.min(w-9,point[0])),Math.max(39,Math.min(h-10,point[1]))];if(clipped[0]!==point[0]||clipped[1]!==point[1]){outside=true;disk(ctx,clipped,6,color);}else disk(ctx,point,8,color);prev=point;}
 if(['rod','elastic','chain'].includes(envOf(r))){ctx.fillStyle='#778d99';ctx.fillText('gravity ↓',w-66,h-12);}ctx.fillStyle='#788d99';ctx.fillText('1 m',12,h-12);line(ctx,[37,h-16],[37+T.scale,h-16],'#788d99',2);
 let ended=time>r.m.recorded_end+1e-8,box=p.el.querySelector('.frame-status');box.classList.toggle('ended',ended);box.textContent=ended?`Trace ended at ${fnum(r.m.recorded_end,3)} s. ${isRef?'No later reference frame.':nice(r.m.status)+': '+nice(r.m.reason)}`:outside?'Body outside the shared reference view':`t = ${fnum(Math.min(time,r.m.recorded_end),3)} s`;
}
function plotLine(ctx,points,color,width=1.5,dash=[]){if(!points.length)return;ctx.beginPath();ctx.strokeStyle=color;ctx.lineWidth=width;ctx.setLineDash(dash);points.forEach((p,i)=>i?ctx.lineTo(p[0],p[1]):ctx.moveTo(p[0],p[1]));ctx.stroke();ctx.setLineDash([]);}
function energyChart(p){let {ctx,w,h}=fit(p.energy),r=p.run,isRef=p.index<0,color=isRef?'#6b8391':COLORS[p.index],scale=Math.max(r.m.energy_scale,1e-12),E0=r.s[0][isRef?2:1],vals=r.s.map(a=>(a[isRef?2:1]-E0)/scale),peak=Math.max(1e-10,...vals.map(Math.abs)),pad=25,plotW=w-pad-4,plotH=h-28;let x=t=>pad+plotW*t/(duration||1),y=v=>12+plotH/2-v/(2*peak)*plotH;ctx.font='9px ui-monospace,monospace';ctx.fillStyle='#7a8d97';ctx.fillText(fnum(peak,2),0,10);ctx.fillText('0',8,y(0)+3);line(ctx,[pad,y(0)],[w,y(0)],'#dfe7eb');plotLine(ctx,r.s.map((a,i)=>[x(a[0]),y(vals[i])]),color);line(ctx,[x(time),9],[x(time),h-13],'#607b89',1,[2,3]);ctx.fillText('0',pad,h-1);ctx.fillText(fnum(duration,1)+'s',w-38,h-1);let current=sample(r,time),v=(current[isRef?2:1]-E0)/scale;p.el.querySelector('.energy-now').textContent=fnum(v,3);}
function pathChart(p){let {ctx,w,h}=fit(p.path),r=p.run,isRef=p.index<0,color=isRef?'#6b8391':COLORS[p.index];ctx.save();ctx.beginPath();ctx.rect(0,0,w,h);ctx.clip();let offset=3+(isRef?r.d:0);
 if(r.d===1){let vals=reference.s.map(a=>a[4]),lo=Math.min(...vals)-.08,hi=Math.max(...vals)+.08;let xy=(t,x)=>[18+(w-24)*t/(duration||1),h-12-(h-24)*(x-lo)/(hi-lo||1)];plotLine(ctx,r.s.map(a=>xy(a[0],a[4])),'#a7b7c1',1,[3,3]);plotLine(ctx,r.s.map(a=>xy(a[0],a[offset])),color);let a=sample(r,time);disk(ctx,xy(Math.min(time,r.m.recorded_end),a[offset]),3,color);
 }else{let sc=Math.min(w,h)/camera.span,to=(x,y)=>[w/2+(x-camera.cx)*sc,h/2-(y-camera.cy)*sc];for(let j=0;j<r.d;j+=2){if(!isRef)plotLine(ctx,r.s.map(a=>to(a[3+r.d+j],a[3+r.d+j+1])),'#a7b7c1',1,[3,3]);plotLine(ctx,r.s.map(a=>to(a[offset+j],a[offset+j+1])),color,1.2);let a=sample(r,time);disk(ctx,to(a[offset+j],a[offset+j+1]),3,color);}disk(ctx,to(0,0),3,'#607683');}ctx.restore();}
function drawAll(){if(!reference)return;$('play').textContent=playing?'Pause':'Play';$('time').textContent=`${fnum(time,2)} / ${fnum(duration,2)} s`;$('scrub').value=String(time);for(let p of panes)if(!p.el.classList.contains('hidden')){world(p);energyChart(p);pathChart(p);}}
function tick(now){if(playing){if(lastWall)time=Math.min(duration,time+(now-lastWall)/1000*+$('speed').value);if(time>=duration)playing=false;drawAll();}lastWall=now;requestAnimationFrame(tick);}
function init(){ $('buildInfo').textContent=`${ALL.length} embedded run traces. Generated ${new Date(DATA.built).toLocaleString()}. Playback uses recorded frames with linear interpolation; source CSVs retain the original samples.`;if(!ALL.length){$('empty').classList.remove('hidden');$('app').classList.add('hidden');return;}
 panes=[makePane(-1),makePane(0),makePane(1),makePane(2)];options($('environment'),ENV.filter(([key])=>ALL.some(r=>envOf(r)===key)));$('environment').addEventListener('change',updateInitial);$('initial').addEventListener('change',updateSetups);$('setup').addEventListener('change',updatePool);$('paneCount').addEventListener('change',refreshPanes);$('play').addEventListener('click',()=>{if(time>=duration)time=0;playing=!playing;lastWall=0;drawAll();});$('reset').addEventListener('click',()=>{time=0;playing=false;drawAll();});$('scrub').addEventListener('input',()=>{playing=false;time=+$('scrub').value;drawAll();});window.addEventListener('resize',drawAll);updateInitial();requestAnimationFrame(tick);}
init();window.viewerState=()=>({runCount:ALL.length,case:reference?.m.case,time,duration,playing,selected:selected.map(r=>r?.m.id),visiblePanes:panes.filter(p=>!p.el.classList.contains('hidden')).length});
</script></body></html>'''


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--results', type=Path, default=ROOT/'results')
    ap.add_argument('--out', type=Path, default=ROOT/'viewer.html')
    ap.add_argument('--max-frames', type=int, default=600)
    ap.add_argument('--per-setup', type=int, default=10)
    args = ap.parse_args()
    if args.max_frames < 2:
        ap.error('--max-frames must be at least 2')
    if args.per_setup < 3:
        ap.error('--per-setup must be at least 3')
    payload = collect(args.results, args.max_frames, args.per_setup)
    encoded = json.dumps(payload, separators=(',', ':'), allow_nan=False).replace('</', r'<\/')
    args.out.write_text(HTML.replace('__DATA__', encoded))
    print(f'{args.out}: {len(payload["runs"])} traces; {args.out.stat().st_size:,} bytes')


if __name__ == '__main__':
    main()
