"""Merge explicit reruns, generate measured tables and motion/energy figures."""
from pathlib import Path
import csv,json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation,FFMpegWriter
ROOT=Path(__file__).resolve().parent

def read(p):
    if not p.exists():return []
    with p.open() as f:return list(csv.DictReader(f))

def key(r):return tuple(r[k] for k in ('case','method','sweeps','hz'))

def merge():
    summaries={};states={}
    for folder,group in [('planar_results','endpoint'),('raw_recheck','endpoint'),('planar_results','dg'),('planar_results','fepr'),('planar_results','cpd'),('stress_results','endpoint'),('stress_results','fepr'),('stress_results','cpd'),('substep_results','fepr')]:
        for r in read(ROOT/folder/(group+'_summary.csv')):summaries[key(r)]=r
        current={}
        for r in read(ROOT/folder/(group+'_states.csv')):current.setdefault(key(r),[]).append(r)
        states.update(current)
    out=ROOT/'results';out.mkdir(exist_ok=True)
    with (out/'summary.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(next(iter(summaries.values()))));w.writeheader();w.writerows(summaries.values())
    return out,summaries,states


def plots(out,summaries,states):
    methods=[('xpbd','1','XPBD','#777777'),('nonlinear_EL','1','Exact E → L, 1 sweep','#d95f02'),('pd_true','8','Strict constrained PD','#a6761d'),('fepr_step_target','1','FEPR (step target)','#1b9e77'),('dg','1','Discrete gradient control','#7570b3')]
    fig,axs=plt.subplots(2,3,figsize=(14,7),sharex='col',layout='constrained')
    for col,case in enumerate(['spring_release','spring_orbit','pendulum_release']):
        reference=False
        for method,sweep,label,color in methods:
            rr=states.get((case,method,sweep,'60'),[])
            if not rr:continue
            t=np.array([float(r['time']) for r in rr]);x=np.array([[float(r['x']),float(r['y'])] for r in rr]);ref=np.array([[float(r['ref_x']),float(r['ref_y'])] for r in rr]);e=np.array([float(r['rel_energy']) for r in rr])
            def coord(z):
                if case=='spring_release':return z[:,0]-1
                if case=='spring_orbit':return np.unwrap(np.arctan2(z[:,1],z[:,0]))
                return np.unwrap(np.arctan2(z[:,0],-z[:,1]))
            if not reference:
                axs[0,col].plot(t,coord(ref),'k--',lw=1.5,label='Continuous reference');reference=True
            axs[0,col].plot(t,coord(x),color=color,lw=1.4,label=label)
            axs[1,col].semilogy(t,np.maximum(np.abs(e),1e-16),color=color,lw=1.3)
        axs[0,col].set_title(case.replace('_',' ').capitalize()+' · 60 Hz')
        axs[0,col].set_ylabel('Extension (m)' if col==0 else 'Angle (rad)')
        axs[1,col].set_ylabel('Absolute normalized energy error');axs[1,col].set_xlabel('Time (s)')
        for ax in axs[:,col]:ax.grid(alpha=.2)
    handles,labels=axs[0,0].get_legend_handles_labels();fig.legend(handles,labels,loc='outside lower center',ncol=3,fontsize=9)
    fig.suptitle('Exact energy can coexist with stopped or incorrect motion',fontsize=15)
    fig.savefig(out/'motion_energy.png',dpi=160);plt.close(fig)


def movie(out,states):
    methods=[('ref','1','Reference'),('xpbd','1','XPBD'),('nonlinear_EL','8','Exact E → L (8)'),('pd_true','8','Strict constrained PD'),('fepr_step_target','1','FEPR'),('dg','1','Discrete gradient')]
    fig,axs=plt.subplots(2,len(methods),figsize=(14,5),layout='constrained');artists=[]
    cases=['spring_release','pendulum_release']
    for row,case in enumerate(cases):
        for col,(method,sweep,label) in enumerate(methods):
            source='dg' if method=='ref' else method;rr=states.get((case,source,sweep,'60'),[])
            if not rr:continue
            prefix='ref_' if method=='ref' else ''
            pos=np.array([[float(r[prefix+'x']),float(r[prefix+'y'])] for r in rr])
            ax=axs[row,col];ax.set_aspect('equal');ax.set_xticks([]);ax.set_yticks([]);ax.set_title(label if row==0 else ('+ hard rod extension' if method.startswith('fepr') else ''),fontsize=10)
            if row==0:ax.set_xlim(-.05,1.4);ax.set_ylim(-.4,.4)
            else:ax.set_xlim(-1.05,1.05);ax.set_ylim(-1.1,.1)
            if col==0:ax.set_ylabel('Released spring' if row==0 else 'Released pendulum')
            ax.plot(0,0,'ks',ms=5);line,=ax.plot([],[],'-',color='#246aab',lw=1.8);dot,=ax.plot([],[],'o',color='#c34a36',ms=8)
            artists.append((line,dot,pos,row))
    title=fig.suptitle('Identical initial states · 60 Hz · normal-speed playback')
    def update(frame):
        j=2*frame
        for line,dot,pos,row in artists:
            p=pos[min(j,len(pos)-1)];line.set_data([0,p[0]],[0,p[1]]);dot.set_data([p[0]],[p[1]])
        title.set_text(f'Identical initial states · 60 Hz · t = {frame/30:.2f} s')
        return [a for entry in artists for a in entry[:2]]
    update(45);fig.savefig(out/'playback_poster.png',dpi=140)
    ani=FuncAnimation(fig,update,frames=151,interval=1000/30,blit=False)
    ani.save(out/'comparison.mp4',writer=FFMpegWriter(fps=30,bitrate=1800),dpi=120);plt.close(fig)


def report(out,summaries):
    methods=[('xpbd','1'),('raw_EL','1'),('nonlinear_EL','1'),('nonlinear_EL','8'),('pd_true','8'),('pd_surrogate','1'),('pd_inertia','1'),('fepr_step_target','1'),('fepr_base','1'),('dg','1'),('cpd','1'),('fepr_sub4','1'),('fepr_sub8','1')]
    lines=['# Measured results','', 'Five-second runs. Energy is normalized by initial energy above the minimum potential. Incomplete runs remain explicitly failed; their prefix errors are not comparable to full-horizon accuracy. `raw_EL` refreshes the material row after E. FEPR uses the paper fast update and a separate hard-rod extension for rigid cases.','']
    for case in ['spring_release','spring_orbit','spring_gravity','pendulum_release','pendulum_kick','spring_stiff','pendulum_large']:
        lines += ['## '+case,'','| Method | Hz | Sweeps | Full 5 s? | Max relative energy error | RMS position (m) | Path length (m) | Failure |','|---|---:|---:|---|---:|---:|---:|---|']
        for hz in ['30','60','240']:
            for meth,sweep in methods:
                r=summaries.get((case,meth,sweep,hz))
                if not r:continue
                lines.append(f"| {meth} | {hz} | {sweep} | {'yes' if r['completed']=='1' else 'NO'} | {float(r['max_rel_energy']):.4g} | {float(r['rms_position']):.5g} | {float(r['path_length']):.5g} | {r['reason']} |")
        lines+=['']
    lines += ['## Scope','', 'The single radial spring contour search currently supports the energy component not enclosing the anchor. Unsupported topology, missing contour/intersection, and nonfinite states report failure. These small cases do not establish robustness for coupled networks, contacts, or indefinite multibody energy surfaces. The strict PD controls are not the published relaxed-energy CPD algorithm. FEPR constraint convergence is not fixed-objective KKT convergence. Discrete-gradient controls change kinematics and use nonlinear root solves. No cross-method runtime claim is made from these Python prototypes.','']
    (out/'TABLES.md').write_text('\n'.join(lines))

if __name__=='__main__':
    out,summaries,states=merge();report(out,summaries);plots(out,summaries,states);movie(out,states)
    print('wrote',len(summaries),'runs to',out)
