#!/usr/bin/env python3
"""Analyze objective-revision experiments without ranking incomplete references."""
import csv
import json
import math
import os
from pathlib import Path
os.environ.setdefault("MPLCONFIGDIR", "/tmp/xpbd-objectives-matplotlib")
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "analysis"
OUT.mkdir(exist_ok=True)
def read(name):
    with (ROOT/name).open() as f:
        return list(csv.DictReader(f))
def write(name, rows):
    with (OUT/name).open("w", newline="") as f:
        w=csv.DictWriter(f, fieldnames=list(rows[0]) if rows else ["empty"])
        w.writeheader(); w.writerows(rows)
local, summary, landscape = read("local.csv"), read("summary.csv"), read("landscape.csv")
near, moving, rod = read("near_rest.csv"), read("moving_angles.csv"), read("rod_kinematics.csv")
orders=[]
for case in sorted({r["case"] for r in local}):
    for variant in sorted({r["variant"] for r in local}):
        rows=sorted([r for r in local if r["case"]==case and r["variant"]==variant], key=lambda r:int(r["hz"]))
        last=rows[-2:]
        out={"case":case,"variant":variant,"last_hz":last[-1]["hz"],"last_two_converged":all(r["converged"]=="1" for r in last)}
        for field in ["position_error","velocity_error"]:
            x,y=(float(r[field]) for r in last)
            out[field+"_order"]=math.log(x/y,2) if min(x,y)>1e-14 else "roundoff_limited"
        out.update(last_position_error=last[-1]["position_error"],last_velocity_error=last[-1]["velocity_error"],last_acceleration_ratio=last[-1]["acceleration_ratio"])
        orders.append(out)
write("local_orders.csv",orders)
failures=[]
for r in summary:
    r["eligible_motion_comparison"] = r["status"]=="complete" and r["completed"]==r["steps"] and r["reference_valid"]=="1"
    reasons=[]
    if r["completed"]!=r["steps"] or r["status"]!="complete":reasons.append("incomplete")
    if r["reference_valid"]!="1":reasons.append("uncertified_reference")
    if int(r["unconverged_steps"])>0:reasons.append("some_steps_missed_residual_threshold")
    if reasons:failures.append(dict(r,failure_categories=";".join(reasons)))
write("trajectory_comparison.csv",summary);write("limitations.csv",failures)
counts={"variants":len({r["variant"]for r in summary}),"local_cases":len(local),"dense_scalar_cases":len(landscape),"trajectories":len(summary),"full_trajectories":sum(r["completed"]==r["steps"]for r in summary),"certified_full_trajectories":sum(r["eligible_motion_comparison"]for r in summary),"trajectories_with_unconverged_steps":sum(int(r["unconverged_steps"])>0 for r in summary),"near_rest_trajectories":len(near),"moving_angle_probes":len(moving),"rod_kinematic_probes":len(rod)}
(OUT/"counts.json").write_text(json.dumps(counts,indent=2)+"\n")
plt.rcParams.update({"font.size":10,"axes.spines.top":False,"axes.spines.right":False,"figure.dpi":130,"savefig.dpi":180})
chosen={"original_midpoint_b1":("Old objective, β=1","#D55E00"),"original_midpoint_h2":("Old objective, β=(hω)²","#E69F00"),"quarter_b1":("Quarter predictor / potential","#0072B2"),"midpotential_b1":("Midpoint potential","#56B4E9"),"projection_verlet_E":("Verlet target + energy","#CC79A7"),"ls_dg_E":("DG mechanical residual + energy","#009E73")}
fig,axes=plt.subplots(2,2,figsize=(11,8),constrained_layout=True)
for ax,case in zip(axes.flat,["scalar_release","scalar_mixed","radial_orbit","chain_mixed"]):
    for v,(label,color)in chosen.items():
        rr=sorted([r for r in local if r["case"]==case and r["variant"]==v],key=lambda r:int(r["hz"]))
        ax.loglog([1/float(r["hz"]) for r in rr],[max(float(r["velocity_error"]),1e-14) for r in rr],"o-",label=label,color=color,markersize=3)
    ax.set(title=case.replace("_"," "),xlabel="Step size h (s)",ylabel="One-step velocity error (m/s)");ax.grid(alpha=.2,which="both")
    ax.set_xticks([1/320,1/160,1/80,1/40,1/20],["1/320","1/160","1/80","1/40","1/20"]);ax.xaxis.set_minor_formatter(NullFormatter())
handles,labels=axes[0,0].get_legend_handles_labels();fig.legend(handles,labels,loc="outside lower center",ncol=2,frameon=False)
fig.suptitle("Correcting the objective restores local consistency\nOrders are measured against a high-resolution physical reference, at turning and moving states.",fontsize=12)
fig.savefig(OUT/"local_consistency.png",bbox_inches="tight");fig.savefig(OUT/"local_consistency.pdf",bbox_inches="tight");plt.close(fig)

fig,axes=plt.subplots(1,2,figsize=(12,4.8),constrained_layout=True)
for variant in ["original_midpoint_b1","original_midpoint_h2","quarter_b1"]:
    rr=sorted([r for r in moving if r["variant"]==variant and int(r["hz"])==5120 and abs(float(r["speed"])-.8)<1e-8],key=lambda r:float(r["angle_degrees"]))
    axes[0].plot([float(r["angle_degrees"])for r in rr],[float(r["acceleration_error_to_physics"])for r in rr],"o-",label=chosen[variant][0],color=chosen[variant][1])
theta=np.linspace(0,180,200);axes[0].plot(theta,10*np.sin(theta*np.pi/180),"k--",label="Predicted old-objective limit: 10 sin θ")
axes[0].set(xlabel="Angle between velocity and force (degrees)",ylabel="Acceleration error (m/s²)",title="Speed 0.8 m/s; h=1/5120 s")
for speed in [.1,.8,2.]:
    rr=sorted([r for r in moving if r["variant"]=="original_midpoint_b1" and abs(float(r["speed"])-speed)<1e-8 and float(r["angle_degrees"])==90],key=lambda r:int(r["hz"]))
    axes[1].loglog([1/float(r["hz"])for r in rr],[float(r["acceleration_error_to_predicted_branch"])for r in rr],"o-",label=f"speed {speed:g} m/s")
axes[1].set(xlabel="Step size h (s)",ylabel="Error relative to predicted biased acceleration",title="The error tends to the derived biased branch")
for ax in axes:ax.grid(alpha=.2,which="both");ax.legend(fontsize=8,frameon=False)
fig.suptitle("The old midpoint objective doubles the force component perpendicular to velocity",fontsize=12)
fig.savefig(OUT/"moving_state_bias.png",bbox_inches="tight");fig.savefig(OUT/"moving_state_bias.pdf",bbox_inches="tight");plt.close(fig)

fig,axes=plt.subplots(3,1,figsize=(11,9),sharex=True,constrained_layout=True)
trace_choices={"original_midpoint_h2":("Old objective, β=(hω)²","#D55E00"),"quarter_b1":("Quarter predictor / potential","#0072B2"),"projection_verlet_E":("Verlet target + energy","#CC79A7"),"projection_verlet_EL":("Verlet target + energy + angular momentum","#E69F00"),"ls_dg_E":("DG mechanical residual + energy","#009E73")}
for variant,(label,color)in trace_choices.items():
    rr=read("traces/radial_orbit__"+variant+".csv");a={k:np.array([float(r[k])for r in rr])for k in rr[0]}
    if variant=="original_midpoint_h2":axes[0].plot(a["time"],a["ref_x0"],"k--",linewidth=1.7,label="Reference")
    axes[0].plot(a["time"],a["x0"],color=color,label=label,linewidth=1)
    axes[1].semilogy(a["time"],np.maximum(abs(a["energy"]-.82)/.82,1e-13),color=color)
    axes[2].semilogy(a["time"],np.maximum(abs(a["angular"]-.88),1e-13),color=color)
axes[0].set(ylabel="Horizontal position (m)",title="2D spring: identical 120 Hz, 10-second orbit case")
axes[1].set(ylabel="Relative energy error")
axes[2].set(ylabel="Angular momentum error",xlabel="Time (s)")
for ax in axes:ax.grid(alpha=.2)
handles,labels=axes[0].get_legend_handles_labels();fig.legend(handles,labels,loc="outside lower center",ncol=2,fontsize=9,frameon=False)
fig.suptitle("Energy preservation can hide incorrect central-force dynamics",fontsize=13)
fig.savefig(OUT/"radial_energy_angular.png",bbox_inches="tight");fig.savefig(OUT/"radial_energy_angular.pdf",bbox_inches="tight");plt.close(fig)
print(json.dumps(counts,indent=2))
