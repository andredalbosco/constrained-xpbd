from pathlib import Path
import csv,json,math
ROOT=Path(__file__).resolve().parent
def read(p):
    with open(ROOT/p) as f:return list(csv.DictReader(f))
def points(rows,prefix,n):return [[round(float(r[prefix+str(j)]),10) for j in range(2*n)] for r in rows]
datasets=[]
for case,n in [('double_moderate',2),('triple_moderate',3)]:
    series=[]
    for folder,method,label in [('objectives_results','obj_midpotential_b1','Midpoint potential: omitted old reaction'),('objectives_results','obj_quarter_b1','Quarter objective: old reaction included'),('results','canonical_dg','Canonical discrete gradient')]:
        rows=read(f'{folder}/trace_{case}_{method}_0.033333.csv')
        series.append({'label':label,'x':points(rows,'x',n)})
    series.append({'label':'Refined reference','x':points(rows,'ref_x',n)})
    datasets.append({'name':case.replace('_',' ').title()+' — 10 seconds','n':n,'t':[float(r['t']) for r in rows],'series':series,'magnification':1,'caption':'30 Hz. References pass the paired-resolution check. All displayed trajectories completed.'})
series=[]
for method,label in [('obj_midpotential_b1','Midpoint potential: omitted old reaction'),('obj_quarter_b1','Quarter objective: old reaction included'),('canonical_dg','Canonical discrete gradient'),('reference','Refined reference')]:
    rows=read(f'modal_transfer/{method}.csv')
    series.append({'label':label,'x':points(rows,'x',3),'modes':[[float(r[f'mode{k}_energy']) for k in range(3)] for r in rows],'energy':[float(r['relative_energy']) for r in rows]})
datasets.append({'name':'Triple highest mode — 60 seconds','n':3,'t':[float(r['t']) for r in rows],'series':series,'magnification':8,'caption':'30 Hz. Initial highest-mode amplitude 0.02, lower-mode seed 2e-7. Angle magnification is a display aid; plotted modal energy and simulation data are unchanged.'})
# One shared viewport per experiment and display magnification, covering all saved frames.
for d in datasets:
    d['bounds']={}
    for zoom in (1,4,8,16):
        xmin=xmax=ymin=ymax=0.0
        for series in d['series']:
            for raw in series['x']:
                oldx=oldy=x=y=0.0
                for k in range(d['n']):
                    nx,ny=raw[2*k:2*k+2]
                    angle=math.atan2(nx-oldx,-(ny-oldy))*zoom
                    x+=math.sin(angle);y-=math.cos(angle)
                    oldx,oldy=nx,ny
                    xmin=min(xmin,x);xmax=max(xmax,x);ymin=min(ymin,y);ymax=max(ymax,y)
        d['bounds'][str(zoom)]=[xmin,xmax,ymin,ymax]
template=r'''<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>Rigid-chain reaction experiment</title><style>
body{background:#f6f7f9;color:#162a34;font:15px system-ui,sans-serif;margin:24px;max-width:1600px}h1{font-size:26px;margin-bottom:8px}p{max-width:1100px;line-height:1.45}.controls{display:flex;flex-wrap:wrap;align-items:center;gap:16px;background:white;padding:16px;border-radius:10px}select,button,input{font:inherit}button{padding:6px 15px}label{display:flex;align-items:center;gap:7px}#time{width:min(44vw,600px)}#panels{display:grid;grid-template-columns:repeat(4,minmax(180px,1fr));gap:12px;margin-top:18px}.panel{background:white;border:1px solid #dce3e8;border-radius:12px;padding:12px}.panel h2{font-size:15px;min-height:42px;margin:0}.panel canvas{width:100%;height:400px}.stats{font:13px ui-monospace,monospace;min-height:70px;line-height:1.5}#clock{font-variant-numeric:tabular-nums;min-width:80px}small{color:#566a76}@media(max-width:900px){#panels{grid-template-columns:repeat(2,1fr)}}
</style><h1>Rigid chains: where does the energy go?</h1><p>Compare the missing old-reaction predictor with the corrected objective. The highest-mode experiment shows almost complete transfer to lower modes even though total energy stays nearly constant.</p><div class="controls"><label>Experiment <select id="experiment"></select></label><button id="play">Play</button><label>Time <input id="time" aria-label="Simulation time" type="range" min="0" step="1" value="0"></label><span id="clock"></span><label>Speed <select id="speed"><option value=".25">0.25×</option><option value=".5">0.5×</option><option selected value="1">1×</option><option value="2">2×</option><option value="4">4×</option></select></label><label>Angle magnification <select id="zoom"><option>1</option><option>4</option><option>8</option><option>16</option></select></label></div><p id="caption"></p><div id="panels"></div><p><small>Each rod has unit length. All panels share one viewport fitted to the full experiment. Modal fractions use fixed linear modes and initial highest-mode energy as their denominator; the energy residual uses exact nonlinear energy. <a href="README.md">Read the study and validation limits</a>.</small></p><script>
const datasets=DATA_PLACEHOLDER;
const $=id=>document.getElementById(id);let selected=0,index=0,playing=false,last=0,elapsed=0;
datasets.forEach((d,i)=>{const o=document.createElement('option');o.value=i;o.textContent=d.name;$('experiment').append(o)});
function setup(){const d=datasets[selected];index=0;elapsed=0;$('time').max=d.t.length-1;$('time').value=0;$('zoom').value=d.magnification;$('caption').textContent=d.caption;$('panels').replaceChildren();d.series.forEach(s=>{const box=document.createElement('section');box.className='panel';const h=document.createElement('h2');h.textContent=s.label;const c=document.createElement('canvas');c.width=600;c.height=800;const stat=document.createElement('div');stat.className='stats';box.append(h,c,stat);$('panels').append(box)});render()}
function render(){const d=datasets[selected],zoom=Number($('zoom').value);$('clock').textContent=d.t[index].toFixed(2)+' s';$('time').value=index;[...$('panels').children].forEach((p,j)=>{const s=d.series[j],raw=s.x[index],c=p.querySelector('canvas'),g=c.getContext('2d'),w=c.width,h=c.height;g.clearRect(0,0,w,h);const [xmin,xmax,ymin,ymax]=d.bounds[zoom],scale=Math.min((w-60)/Math.max(xmax-xmin,.1),(h-70)/Math.max(ymax-ymin,.1)),ox=w/2-scale*(xmin+xmax)/2,oy=h/2+scale*(ymin+ymax)/2;g.strokeStyle='#e0e6eb';g.lineWidth=2;g.beginPath();g.moveTo(ox,oy);g.lineTo(ox,h-25);g.stroke();let prev=[0,0],sum=[0,0];g.strokeStyle=['#d1495b','#00798c','#edae49','#273b49'][j];g.fillStyle=g.strokeStyle;g.lineWidth=5;g.beginPath();g.moveTo(ox,oy);const drawn=[];for(let k=0;k<d.n;k++){const x=raw[2*k],y=raw[2*k+1],dx=x-prev[0],dy=y-prev[1],a=Math.atan2(dx,-dy)*zoom;sum=[sum[0]+Math.sin(a),sum[1]-Math.cos(a)];drawn.push([ox+scale*sum[0],oy-scale*sum[1]]);g.lineTo(...drawn[k]);prev=[x,y]}g.stroke();drawn.forEach(p=>{g.beginPath();g.arc(p[0],p[1],10,0,2*Math.PI);g.fill()});g.fillStyle='#162a34';g.beginPath();g.arc(ox,oy,7,0,2*Math.PI);g.fill();const stats=p.querySelector('.stats');if(s.modes){const e=s.modes[index];stats.textContent='Mode fractions: '+e.map(v=>v.toExponential(2)).join(' / ')+'\nOther modes: '+(e[0]+e[1]).toExponential(3)+'\nEnergy residual: '+s.energy[index].toExponential(2);stats.style.whiteSpace='pre-wrap'}else stats.textContent='All rods and masses shown at '+zoom+'× angle magnification.';});}
$('experiment').onchange=()=>{selected=Number($('experiment').value);setup()};$('zoom').onchange=render;$('play').onclick=()=>{playing=!playing;$('play').textContent=playing?'Pause':'Play';last=0};$('time').oninput=()=>{index=Number($('time').value);elapsed=0;render()};
function frame(t){if(playing&&last){elapsed+=(t-last)/1000*Number($('speed').value);const d=datasets[selected],dt=d.t[1]-d.t[0];if(elapsed>=dt){const step=Math.floor(elapsed/dt);elapsed-=step*dt;index=(index+step)%d.t.length;render()}}last=t;requestAnimationFrame(frame)}setup();requestAnimationFrame(frame);window.viewerState=()=>({experiment:datasets[selected].name,index,time:datasets[selected].t[index],playing,panels:datasets[selected].series.length});
</script></html>'''
(ROOT/'viewer.html').write_text(template.replace('DATA_PLACEHOLDER',json.dumps(datasets,separators=(',',':'))))
print('Viewer:',ROOT/'viewer.html',len(datasets),'experiments')
