from pathlib import Path
import csv
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import NullLocator

ROOT = Path(__file__).resolve().parent
FIG = ROOT / "figures"
FIG.mkdir(exist_ok=True)
def read(name):
    return list(csv.DictReader((ROOT/name).open()))
one = read("one_step_results/one_step.csv")
key = lambda r: tuple(r[k] for k in ("case", "h", "gamma", "extension", "initial"))
full = {key(r): r for r in one if r["form"] == "full"}
summary = {}
for extension in ["balanced", "grounded_first", "grounded_last"]:
    rows = [r for r in one if r["form"] == "reduced" and r["budget"] == "128" and r["extension"] == extension]
    pairs = [r for r in rows if r["converged"] == "1" and full[key(r)]["converged"] == "1"]
    summary[extension] = dict(configurations=len(rows), converged=sum(r["converged"] == "1" for r in rows),
                             converged_pairs=len(pairs), max_position_gap=max(float(r["x_error_to_full"]) for r in pairs),
                             max_velocity_gap=max(float(r["v_error_to_full"]) for r in pairs),
                             max_full_stationarity=max(float(r["full_stationarity"]) for r in pairs))
tr = read("trajectory_results/trajectories.csv")
summary["trajectories"] = dict(runs=len(tr), unconverged_steps=sum(int(r["nonconverged"]) for r in tr),
                              max_reference_gap=max(float(r["reference_gap"]) for r in tr))
norm = read("normalization_results/trajectories.csv")
def trace(r):
    name = "_".join(r[k] for k in ("case", "target", "copies", "normalization")) + ".csv"
    return np.loadtxt(ROOT/"normalization_results"/name, delimiter=",", skiprows=1)[:,1:]
normalization = []
for r in norm:
    base = next(b for b in norm if b["case"] == r["case"] and b["target"] == r["target"] and b["copies"] == "1")
    normalization.append({**r, "motion_change_per_coordinate":float(np.sqrt(np.mean((trace(r)-trace(base))**2)))})
summary["normalization"] = dict(runs=len(norm),unconverged_steps=sum(int(r["nonconverged"]) for r in norm),
                                max_compensated_motion_change=max(r["motion_change_per_coordinate"] for r in normalization if r["normalization"] == "compensated"))
with (ROOT/"normalization_results"/"comparison.csv").open("w") as f:
    w=csv.DictWriter(f,fieldnames=normalization[0].keys());w.writeheader();w.writerows(normalization)
large = read("larger_results/one_step.csv")
summary["larger"] = {}
for form in ["full","reduced"]:
    rows=[r for r in large if r["form"]==form and r["budget"]=="128"]
    summary["larger"][form] = dict(configurations=len(rows),converged=sum(r["converged"]=="1" for r in rows))
if (ROOT/"nonlinear_results/results.csv").exists():
    controls=read("nonlinear_results/results.csv")
    summary["nonlinear_controls"]={m:dict(configurations=sum(r["method"]==m for r in controls),converged=sum(r["method"]==m and r["converged"]=="1" for r in controls)) for m in sorted(set(r["method"] for r in controls))}
if (ROOT/"polish_results/results.csv").exists():
    polish=read("polish_results/results.csv")
    summary["newton_polish"] = dict(configurations=len(polish),converged=sum(r["converged"]=="1" for r in polish))
if (ROOT/"disconnected_results/trajectories.csv").exists():
    dis=read("disconnected_results/trajectories.csv")
    summary["disconnected"] = dict(configurations=len(dis),fully_converged=sum(r["nonconverged"]=="0" for r in dis),max_reference_gap=max(float(r["reference_gap"]) for r in dis))
if (ROOT/"rollout_results/trajectories.csv").exists():
    rollouts=read("rollout_results/trajectories.csv")
    summary["larger_trajectories"]={m:dict(configurations=sum(r["method"]==m for r in rollouts),complete=sum(r["method"]==m and r["complete"]=="1" for r in rollouts)) for m in sorted(set(r["method"] for r in rollouts))}
(ROOT/"summary.json").write_text(json.dumps(summary,indent=2)+"\n")

plt.rcParams.update({"font.size":10,"axes.spines.top":False,"axes.spines.right":False})
fig,ax=plt.subplots(1,2,figsize=(12,4.2),layout="constrained")
names=["scalar_release","scalar_mixed","radial_orbit","elastic_mixed","chain_mixed"]
for label,form,ext,offset,color in [("Full work variables","full","balanced",-.18,"#3366aa"),("Eliminated work","reduced","balanced",.18,"#db8030")]:
    values=[float(next(r for r in tr if r["case"]==n and r["form"]==form and r["extension"]==ext and float(r["h"])>.01)["position_rms"]) for n in names]
    ax[0].bar(np.arange(5)+offset,values,.35,label=label,color=color)
ax[0].set(xticks=np.arange(5),xticklabels=["Scalar\nrelease","Scalar\nmoving","Radial\nspring","Elastic\npendulum","Two-spring\nchain"],ylabel="Position RMS to reference (m)",title="Same objective, same motion: 5 s at 60 Hz")
ax[0].legend(frameon=False)
for case,color in zip(["radial_orbit","elastic_mixed","chain_mixed"],["#3366aa","#db8030","#36835b"]):
    for target,style in [("frozen_power","-"),("implicit_midpoint_dg","--")]:
        rows=[r for r in normalization if r["case"]==case and r["target"]==target and r["normalization"]=="total"]
        ax[1].loglog([int(r["copies"]) for r in rows],[r["motion_change_per_coordinate"] for r in rows],style+"o",color=color,label=case.replace("_"," ") if target=="frozen_power" else None)
ax[1].set(xlabel="Independent identical copies (symmetry-reduced)",ylabel="Motion change from single copy (m / coordinate)",title="Total-energy normalization changes the objective")
ax[1].set_xticks([2,8,32],["2","8","32"])
ax[1].xaxis.set_minor_locator(NullLocator())
ax[1].set_ylim(1e-6,.2)
ax[1].legend(frameon=False,fontsize=8)
ax[1].text(.03,.04,"Solid: frozen power; dashed: predicted midpoint DG\nCompensated normalization: zero change in this test",transform=ax[1].transAxes,fontsize=8)
fig.savefig(FIG/"elimination_and_normalization.png",dpi=170)
fig.savefig(FIG/"elimination_and_normalization.pdf")
plt.close(fig)
if (ROOT/"disconnected_results/trajectories.csv").exists():
    fig,axes=plt.subplots(1,2,figsize=(11,3.8),layout="constrained")
    for ax,objective in zip(axes,["quarter","projection"]):
        data=np.loadtxt(ROOT/"disconnected_results"/("scalar_mixed_radial_orbit_"+objective+"_global_0.033333.csv"),delimiter=",",skiprows=1)
        for j,label,color in [(1,"Scalar spring","#3366aa"),(2,"Radial spring","#db8030"),(3,"Their sum","#444444")]:
            ax.plot(data[:,0],data[:,j],label=label,color=color,lw=1.7 if j!=3 else 1)
        ax.set(xlabel="Time (s)",ylabel="Change in component energy (J)",title=("Quarter objective" if objective=="quarter" else "Projection to Verlet target"))
        ax.legend(frameon=False)
    fig.suptitle("A single total-energy constraint couples two disconnected systems: 30 Hz")
    fig.savefig(FIG/"disconnected_energy_transfer.png",dpi=170)
    fig.savefig(FIG/"disconnected_energy_transfer.pdf")
    plt.close(fig)
print(json.dumps(summary,indent=2))
