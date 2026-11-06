"""Unforced relative-coordinate comparison; all trajectories are precomputed.

Full CPD has P/L/E rows. The no-L column is explicitly an ablation with
unchanged Kstar. Independent-v LE is read from the previous canonical archive.
"""
from pathlib import Path
import csv
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation,FFMpegWriter

ROOT=Path(__file__).resolve().parent
INDEPENDENT=Path("/home/blacknoir/Projects/constrained-xpbd/experiments/independent_velocity_2026_11_01/results")
CASES=("release","orbit","stiff")
COLORS=("#d55e00","#a6761d","#0072b2")
LABELS=("Full CPD: P / L / E","Angular-row ablation: P / E","Independent-v LE · βᵥ=1")


def read(path):
    with path.open(newline="") as f:return list(csv.DictReader(f))


def select(source,case,independent=False):
    rows=read(source/"summary.csv")
    candidates=[r for r in rows if r["case"]==case and r["hz"]=="60" and float(r["epsilon_ratio"])==1e8
                and (r["sweeps"]=="16" and float(r["beta_v"])==1. and r["order"]=="LE" if independent
                     else r["policy"]=="fixed" and r["iterations_budget"]=="16")]
    if len(candidates)!=1:raise ValueError(f"Expected one selected run: {source} {case}, got{len(candidates)}")
    row=candidates[0]
    records=read(source/"trajectories"/(row["run_id"]+".csv"))
    accepted=[r for r in records if r["accepted"]=="1"]
    d={k:np.array([float(r[k]) for r in accepted]) for k in accepted[0]}
    if independent:
        p=np.c_[d["x"],d["y"]]
        v=np.c_[d["vx"],d["vy"]]
        angular=d["angular_momentum"]
        radius=d["radius"]
        mismatch=d["velocity_chord_mismatch"]
    else:
        p=np.c_[d["x_0"]-d["x_1"],d["y_0"]-d["y_1"]]
        v=np.c_[d["vx_0"]-d["vx_1"],d["vy_0"]-d["vy_1"]]
        angular=d["angular_global"]
        radius=d["mode_length"]
        mismatch=np.zeros(len(p))
    return dict(row=row,time=d["time"],pos=p,velocity=v,speed=np.linalg.norm(v,axis=1),
                radius=radius,angular=angular,energy=d["energy"],rel_energy=d["rel_energy"],
                mismatch=mismatch,source=str(source))


def load_all():
    return {case:[select(ROOT/"results","free_"+case),
                  select(ROOT/"angular_ablation","free_"+case),
                  select(INDEPENDENT,"spring_"+case,True)] for case in CASES}


def scientific_plot(datasets,out):
    fig,axes=plt.subplots(3,3,figsize=(14,9),layout="constrained")
    for i,(case,runs) in enumerate(datasets.items()):
        for j,run in enumerate(runs):
            axes[i,0].plot(run["time"],run["radius"],color=COLORS[j],lw=1.3,label=LABELS[j])
            axes[i,1].plot(run["time"],run["angular"],color=COLORS[j],lw=1.3)
            axes[i,2].plot(run["time"],100*run["rel_energy"],color=COLORS[j],lw=1.3)
        axes[i,0].set_title("Released spring" if case=="release" else ("Stiff orbiting spring" if case=="stiff" else "Orbiting spring"))
        axes[i,0].set_ylabel("Spring separation (m)")
        axes[i,1].set_ylabel("Angular momentum")
        axes[i,2].set_ylabel("Actual energy error (%)")
        for ax in axes[i]:
            ax.set_xlabel("Time (s)");ax.grid(alpha=.22)
            ax.ticklabel_format(axis="y",style="plain",useOffset=False)
    handles,labels=axes[0,0].get_legend_handles_labels()
    fig.legend(handles,labels,loc="outside lower center",ncol=3,fontsize=9)
    fig.suptitle("Unforced spring modes · 60 Hz · 16 iterations/sweeps · ε ratio=10⁸\n"
                 "CPD: free two-body relative coordinate, reduced mass1. Independent-v: equivalent anchored mass1. No phase scoring.",fontsize=12)
    fig.savefig(out/"unforced_relative_modes.png",dpi=150)
    plt.close(fig)


def spring_points(p):
    length=np.linalg.norm(p)
    along=np.r_[0.,.08,np.linspace(.1,.9,19),.92,1.]
    side=np.zeros(len(along));side[2:-2]=.035*np.where(np.arange(19)%2,-1.,1.)
    points=along[:,None]*p+side[:,None]*np.array([-p[1],p[0]])/max(length,1e-12)
    return points[:,0],points[:,1]


def movie(datasets,out):
    fig,axes=plt.subplots(3,3,figsize=(14,11),layout="constrained")
    artists=[]
    for i,(case,runs) in enumerate(datasets.items()):
        for j,run in enumerate(runs):
            ax=axes[i,j];ax.set_aspect("equal");ax.set_xlim(-1.4,1.4);ax.set_ylim(-1.4,1.4)
            ax.set_xticks([]);ax.set_yticks([]);ax.plot(0,0,"ks",ms=4)
            if i==0:ax.set_title(LABELS[j],fontsize=11)
            if j==0:ax.set_ylabel(case.replace("_"," ").capitalize()+" spring",fontsize=11)
            line,=ax.plot([],[],color=COLORS[j],lw=1.5)
            dot,=ax.plot([],[],"o",color=COLORS[j],ms=8)
            trail,=ax.plot([],[],color=COLORS[j],lw=.9,alpha=.3)
            text=ax.text(.02,.98,"",transform=ax.transAxes,va="top",fontsize=9,
                         bbox=dict(facecolor="white",alpha=.8,edgecolor="none",pad=1.))
            tail=run["radius"][run["time"]>=4]
            radius_frozen=np.ptp(tail)<1e-6 and abs(tail[-1]-1)>.02
            position_frozen=float(run["row"]["path_length"])<1e-4
            artists.append(dict(run=run,line=line,dot=dot,trail=trail,text=text,
                                radius_frozen=radius_frozen,position_frozen=position_frozen,case=case,column=j))
    title=fig.suptitle("")
    fps=30
    def update(frame):
        t=frame/fps
        title.set_text(f"Relative spring motion · 60 Hz · 16 global steps / GS sweeps · ε ratio=10⁸ · t={t:.2f}s\n"
                       "r = x₀ − x₁ for CPD; reduced mass1. Independent-v column uses equivalent anchored coordinates.")
        for a in artists:
            run=a["run"]
            k=min(max(0,np.searchsorted(run["time"],t+1e-10,side="right")-1),len(run["time"])-1)
            p=run["pos"][k]
            a["line"].set_data(*spring_points(p));a["dot"].set_data([p[0]],[p[1]])
            start=np.searchsorted(run["time"],max(0.,t-.7))
            a["trail"].set_data(run["pos"][start:k+1,0],run["pos"][start:k+1,1])
            message=(f"length {run['radius'][k]:.6f}m\nspeed {run['speed'][k]:.3g}m/s\n"
                     f"L {run['angular'][k]:.6g}\nactual E {run['energy'][k]:.7g}J\n"
                     f"ΔE {100*run['rel_energy'][k]:+.3g}%")
            if a["column"]==2:message+=f"\n|v−chord| {run['mismatch'][k]:.3g}m/s"
            if a["position_frozen"]:message+="\nNearly frozen trajectory"
            elif a["radius_frozen"]:message+="\nNearly constant stretched radius"
            a["text"].set_text(message)
            a["text"].set_color("#a40000" if a["position_frozen"] or a["radius_frozen"] else "#222222")
        return []
    update(fps)
    fig.savefig(out/"unforced_poster.png",dpi=120)
    fig.set_layout_engine("none")
    animation=FuncAnimation(fig,update,frames=151,interval=1000/fps,blit=False)
    animation.save(out/"comparison.mp4",writer=FFMpegWriter(fps=fps,bitrate=3000,
                    extra_args=["-pix_fmt","yuv420p"]),dpi=120)
    plt.close(fig)


def main():
    out=ROOT/"results"/"figures";out.mkdir(exist_ok=True,parents=True)
    datasets=load_all();scientific_plot(datasets,out)
    manifest=[dict(case=case,column=LABELS[j],source=run["source"],run_id=run["row"]["run_id"])
              for case,runs in datasets.items() for j,run in enumerate(runs)]
    with (out/"unforced_selected_runs.csv").open("w",newline="") as f:
        writer=csv.DictWriter(f,fieldnames=list(manifest[0]));writer.writeheader();writer.writerows(manifest)
    movie(datasets,out)
    print("Wrote unforced original/ablation/independent-v comparison",flush=True)


if __name__=="__main__":main()
