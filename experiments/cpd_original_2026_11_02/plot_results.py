"""Read CPD campaign CSVs; distinguish actual energy from paper energy budget.

No trajectory integration or continuous-reference comparison is performed.
Original full free-body cases and anchored/hard-rod adaptations are labeled.
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, FFMpegWriter


ROOT = Path(__file__).resolve().parent
COLORS = {"paper": "#d55e00", "fixed": "#0072b2", "stationary": "#7b2cbf"}
LABELS = {"paper": "Paper constraint stop", "fixed": "Fixed 16 global steps", "stationary": "Full KKT stop"}


def read(path):
    with path.open(newline="") as f:
        return list(csv.DictReader(f))


def find(rows, case, policy, args):
    found = [r for r in rows if r["case"] == case and r["policy"] == policy
             and r["method"] == "CPD" and float(r["hz"]) == args.hz
             and float(r["epsilon_ratio"]) == args.epsilon_ratio and r["angular_active"] == "1"]
    return found[0] if found else None


def load(results, row):
    if row is None:
        return None
    records = read(results/"trajectories"/(row["run_id"]+".csv"))
    data = {key: np.asarray([float(r[key]) for r in records]) for key in records[0]}
    accepted = data["accepted"] > .5
    return {key: value[accepted] for key,value in data.items()}


def positions(data, scene):
    return np.stack([np.c_[data[f"x_{i}"],data[f"y_{i}"]] for i in range(len(scene["masses"]))],axis=1)


def scope(scene):
    if scene["family"] == "free":
        return "Full original P / L / E rows"
    if scene["family"] == "hard_rod_extension":
        return "Exact hard-rod extension"
    return "Anchored adaptation: P omitted, L retained"


def by_case(rows, args, scenes):
    return {case:[(policy,find(rows,case,policy,args)) for policy in ("paper","fixed","stationary")
                  if find(rows,case,policy,args) is not None] for case in scenes}


def traces(rows, scenes, args, out):
    for family in ("free","anchored","hard_rod_extension"):
        cases=[name for name,s in scenes.items() if s["family"] == family]
        if not cases:
            continue
        fig,axes=plt.subplots(len(cases),3,figsize=(15,3*len(cases)),squeeze=False,layout="constrained")
        for i,case in enumerate(cases):
            for policy in ("paper","fixed","stationary"):
                row=find(rows,case,policy,args)
                data=load(args.results,row)
                if data is None:
                    continue
                color=COLORS[policy]
                axes[i,0].plot(data["time"],data["mode_length"],color=color,lw=1.2,label=LABELS[policy])
                axes[i,1].plot(data["time"],100*data["rel_energy"],color=color,lw=1.2)
                axes[i,2].plot(data["time"],data["angular_global"],color=color,lw=1.2)
                if row["completed"] != "1":
                    for ax in axes[i]:
                        ax.axvline(float(row["accepted_time"]),color=color,ls=":",lw=.9)
            axes[i,0].set_title(case.replace("_"," "))
            axes[i,0].set_ylabel("Spring length / triangle gyration (m)")
            axes[i,1].set_ylabel("Actual mechanical energy error (%)")
            axes[i,2].set_ylabel("Angular momentum about origin")
            angular_values=np.concatenate([line.get_ydata() for line in axes[i,2].lines])
            angular_margin=max(.05,.1*float(np.max(np.abs(angular_values))))
            axes[i,2].set_ylim(float(np.min(angular_values))-angular_margin,
                               float(np.max(angular_values))+angular_margin)
            for ax in axes[i]:
                ax.set_xlabel("Time (s)")
                ax.grid(alpha=.22)
                ax.ticklabel_format(axis="y",style="plain",useOffset=False)
        axes[0,0].legend(fontsize=8)
        fig.suptitle(f"{scope(scenes[cases[0]])} · {args.hz} Hz · ε / initial energy scale={args.epsilon_ratio:g}\n"
                     "Actual energy includes the conservative potential of external constant forces; no phase scoring",fontsize=12)
        fig.savefig(out/(family+"_motion_energy_angular.png"),dpi=130)
        plt.close(fig)


def budget_and_stopping(rows,scenes,args,out):
    fig,axes=plt.subplots(len(scenes),3,figsize=(15,2.7*len(scenes)),squeeze=False,layout="constrained")
    for i,(case,scene) in enumerate(scenes.items()):
        row=find(rows,case,"paper",args)
        data=load(args.results,row)
        if data is None:
            continue
        scale=float(row["energy_scale"])
        axes[i,0].plot(data["time"],100*data["rel_energy"],label="Actual ΔE",color="black",lw=1.2)
        axes[i,0].plot(data["time"],100*data["work_approximation_defect"]/scale,label="Paper power − exact work",color="#d55e00",lw=.9)
        axes[i,0].plot(data["time"],100*data["cumulative_relaxation"]/scale,label="Energy relaxation",color="#7b2cbf",lw=.9)
        axes[i,0].plot(data["time"],100*data["cumulative_constraint_defect"]/scale,label="Accumulated closure defect",color="#0072b2",lw=.9)
        axes[i,0].set_title(case.replace("_"," ")+" · paper stopping")
        axes[i,0].set_ylabel("Mechanical energy budget components (%)")
        for policy in ("paper","fixed","stationary"):
            rr=find(rows,case,policy,args)
            dd=load(args.results,rr)
            if dd is None:
                continue
            axes[i,1].semilogy(dd["time"],np.maximum(dd["raw_constraints"],1e-18),color=COLORS[policy],lw=1.,label=LABELS[policy])
            axes[i,2].semilogy(dd["time"],np.maximum(dd["kkt_scaled"],1e-18),color=COLORS[policy],lw=1.)
        axes[i,1].axhline(1e-4,color="black",lw=.7,ls=":")
        axes[i,2].axhline(1e-9,color="black",lw=.7,ls=":")
        axes[i,1].set_ylabel("Raw conservation-row residual")
        axes[i,2].set_ylabel("Full scaled KKT residual")
        for ax in axes[i]:
            ax.grid(alpha=.2)
            ax.set_xlabel("Time (s)")
    axes[0,0].legend(fontsize=7)
    axes[0,1].legend(fontsize=7)
    fig.suptitle("Actual energy budget = power approximation defect + relaxation + conservation closure defect\n"
                 f"{args.hz} Hz · ε ratio={args.epsilon_ratio:g} · families and source policies remain distinct",fontsize=12)
    fig.savefig(out/"energy_budget_stopping.png",dpi=120)
    plt.close(fig)


def spring_points(a,b):
    difference=a-b
    radius=np.linalg.norm(difference)
    if radius <= 1e-12:
        return [a[0],b[0]],[a[1],b[1]]
    along=np.r_[0.,.08,np.linspace(.1,.9,19),.92,1.]
    side=np.zeros(len(along))
    side[2:-2]=.03*np.where(np.arange(19)%2,-1.,1.)
    points=b+along[:,None]*difference+side[:,None]*np.array([-difference[1],difference[0]])/radius
    return points[:,0],points[:,1]


def draw_state(ax,scene,pos,color):
    for edge in scene["edges"]:
        a=pos[edge["a"]]
        b=pos[edge["b"]] if edge["b"]>=0 else np.zeros(2)
        xx,yy=([a[0],b[0]],[a[1],b[1]]) if edge["hard"] else spring_points(a,b)
        ax.plot(xx,yy,color=color,lw=1.3)
        if edge["b"]<0:
            ax.plot(0,0,"ks",ms=4)
    ax.scatter(pos[:,0],pos[:,1],s=28,color=color,zorder=5)


def contacts(rows,scenes,args,out):
    times=(0.,.25,.5,1.,2.,3.,4.,5.)
    for case,entries in by_case(rows,args,scenes).items():
        if not entries:
            continue
        scene=scenes[case]
        datasets=[(p,r,load(args.results,r)) for p,r in entries]
        radius=max(1.1,min(3.,1.1*max(np.max(np.linalg.norm(positions(d,scene),axis=2)) for _,_,d in datasets)))
        fig,axes=plt.subplots(len(datasets),8,figsize=(18,3*len(datasets)),squeeze=False,layout="constrained")
        for i,(policy,row,data) in enumerate(datasets):
            pos=positions(data,scene)
            for j,t in enumerate(times):
                index=min(max(0,np.searchsorted(data["time"],t+1e-10,side="right")-1),len(data["time"])-1)
                ax=axes[i,j]
                ax.set_aspect("equal");ax.set_xlim(-radius,radius);ax.set_ylim(-radius,radius)
                ax.set_xticks([]);ax.set_yticks([]);ax.set_title(f"t={t:g} s",fontsize=9)
                if j==0:ax.set_ylabel(LABELS[policy],fontsize=9)
                for particle in range(len(scene["masses"])):
                    ax.plot(pos[:index+1,particle,0],pos[:index+1,particle,1],color=COLORS[policy],alpha=.25,lw=.6)
                draw_state(ax,scene,pos[index],COLORS[policy])
                text=(f"length {data['mode_length'][index]:.4f}\nE {data['energy'][index]:.5g} J\n"
                      f"ΔE {100*data['rel_energy'][index]:+.3g}%\nρ {data['rho'][index]:+.2g}\n"
                      f"KKT {data['kkt_scaled'][index]:.2g}")
                failed=row["completed"]!="1" and t>float(row["accepted_time"])
                frozen=float(row["path_length"])<1e-4
                tail=float(row["tail_path_length"])<1e-4 and t>=4.
                if failed:text+="\nFAILED · prefix held"
                elif frozen:text+="\nNearly frozen trajectory"
                elif tail:text+="\nNearly stationary tail"
                ax.text(.02,.98,text,transform=ax.transAxes,va="top",fontsize=7,
                        color="#a40000" if failed or frozen or tail else "#222222",
                        bbox=dict(facecolor="white",alpha=.8,edgecolor="none",pad=1.))
        fig.suptitle(case.replace("_"," ")+" · "+scope(scene)+f" · {args.hz} Hz · ε ratio={args.epsilon_ratio:g}\n"
                     "Archived accepted-state snapshots; raw constraint stopping and full stationarity are separate",fontsize=12)
        fig.savefig(out/(case+"_contacts.png"),dpi=130)
        plt.close(fig)


def tables(rows,scenes,args,out):
    lines=["# CPD measured comparison","",f"Selected {args.hz} Hz, epsilon ratio {args.epsilon_ratio:g}.","",
           "Actual energy is K+elastic−sum(f·x), while the algorithm constrains H=K+elastic using incoming external power. Policies and extensions remain explicit.","",
           "| Case | Family | Policy | Complete? | Max actual E error | Relative travel | Tail travel | Raw residual | KKT | Stationary steps | Failure |",
           "|---|---|---|---|---:|---:|---:|---:|---:|---:|---|"]
    selected=[]
    for case,entries in by_case(rows,args,scenes).items():
        for policy,row in entries:
            travel=float(row["relative_mode_path_length"])
            if not np.isfinite(travel):travel=float(row["path_length"])
            lines.append(f"| {case} | {row['family']} | {policy} | {row['completed']} | "
                         f"{100*float(row['max_accepted_rel_energy']):.5g}% | {travel:.5g} | "
                         f"{float(row['tail_path_length']):.5g} | {float(row['max_raw_constraints']):.3g} | "
                         f"{float(row['max_kkt_scaled']):.3g} | {100*float(row['stationary_fraction']):.2f}% | {row['reason']} |")
            selected.append(row)
    lines += ["","A completed run has reached its requested horizon; that flag does not establish plausible motion or stationarity.",""]
    (out/"SELECTED_RESULTS.md").write_text("\n".join(lines))
    if selected:
        with (out/"selected_runs.csv").open("w",newline="") as f:
            w=csv.DictWriter(f,fieldnames=list(selected[0]));w.writeheader();w.writerows(selected)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results",type=Path,default=ROOT/"results")
    parser.add_argument("--output",type=Path)
    parser.add_argument("--hz",type=int,default=60)
    parser.add_argument("--epsilon-ratio",type=float,default=1e8)
    parser.add_argument("--cases",default="")
    args=parser.parse_args()
    out=args.output or args.results/"figures"
    out.mkdir(parents=True,exist_ok=True)
    rows=read(args.results/"summary.csv")
    metadata=json.loads((args.results/"metadata.json").read_text())
    scenes={s["name"]:s for s in metadata["scenes"] if not args.cases or s["name"] in args.cases.split(",")}
    tables(rows,scenes,args,out)
    traces(rows,scenes,args,out)
    budget_and_stopping(rows,scenes,args,out)
    contacts(rows,scenes,args,out)
    print(f"Wrote CPD motion/energy inspection figures to {out}")


if __name__=="__main__":main()
