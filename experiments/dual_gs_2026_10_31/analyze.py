"""Summarize the measured CSV files; no solver or trajectory repair here."""
import csv
import json
import math
import pathlib
import statistics
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

root = pathlib.Path(sys.argv[1] if len(sys.argv) > 1 else "results")
figdir = root / "figures"
figdir.mkdir(exist_ok=True)
def read(name):
    path = root / name
    if not path.exists(): return []
    rows = list(csv.DictReader(path.open()))
    for row in rows:
        for key, value in row.items():
            try: row[key] = float(value)
            except (ValueError, TypeError): pass
    return rows
def med(rows,key):
    vals=[r[key] for r in rows if math.isfinite(r[key])]
    return statistics.median(vals) if vals else None
methods=["refreshed_drop","refreshed_local_rx","frozen_full_rx"]
shortnames={"refreshed_drop":"Drop stationarity","refreshed_local_rx":"Retain local stationarity","frozen_full_rx":"Frozen GS + restoration"}
schedules=["material_storage_kinetic","kinetic_storage_material","storage_kinetic_material","edge","edge_storage_first","edge_reverse"]
labels=["material → storage → kinetic","kinetic → storage → material","storage → kinetic → material","edge interleaved","edge, storage first","edge, reverse traversal"]
one=read("one_step.csv")
summary={"one_step_rows":len(one),"one_step_by_method":{},"short":{},"long":{}}
for method in methods:
    selected=[r for r in one if r["method"]==method]
    final=[r for r in selected if r["budget"]/(3*r["N"]+(r["case"].startswith(("rod","elastic","chain")))*r["N"]+(r["case"].startswith("rod"))*r["N"])==128]
    summary["one_step_by_method"][method]={"rows":len(selected),"failures":sum(bool(r["failure"]) for r in selected),"final_rows":len(final),"full_kkt_below_1e-6":sum(r["full_residual"]<1e-6 for r in final),"median_final_full_residual":med(final,"full_residual"),"median_final_stationarity":med(final,"stationarity"),"median_final_energy":med(final,"energy")}

if one:
    fig, axes=plt.subplots(3,2,figsize=(12,10),sharex=True)
    for mi,method in enumerate(methods):
        for si,schedule in enumerate(schedules):
            rows=sorted([r for r in one if r["case"]=="chain3" and r["h"]==.08 and r["method"]==method and r["schedule"]==schedule],key=lambda r:r["budget"])
            for ci,key in enumerate(["full_residual","energy"]):
                axes[mi,ci].loglog([r["visits"] for r in rows],[max(1e-15,abs(r[key])) for r in rows],"o-",ms=3,label=labels[si])
        axes[mi,0].set_ylabel(shortnames[method]+"\nresidual")
        for ax in axes[mi]: ax.grid(True,alpha=.2)
    axes[0,0].set_title("Full KKT residual, including stationarity")
    axes[0,1].set_title("Relative one-step energy error")
    for ax in axes[-1]:ax.set_xlabel("Actual scalar multiplier visits")
    axes[0,0].legend(fontsize=8)
    fig.suptitle("Three moving masses, h = 0.08 s; identical objective and start")
    fig.tight_layout();fig.savefig(figdir/"chain3_convergence.png",dpi=180);fig.savefig(figdir/"chain3_convergence.pdf");plt.close(fig)

for prefix in ["short","long"]:
    rows=read(prefix+"_rollouts.csv")
    for method in methods+["dense_oracle"]:
        selected=[r for r in rows if r["method"]==method]
        complete=[r for r in selected if not r["failed"]]
        summary[prefix][method]={"rows":len(selected),"complete":len(complete),"median_max_energy_complete":med(complete,"max_energy"),"median_rms_position_complete":med(complete,"rms_position"),"worst_reference_gap":max([r["reference_gap"] for r in selected],default=0)}

traces=read("long_trajectories.csv") or read("short_trajectories.csv")
if traces:
    budget=32 if (root/"long_trajectories.csv").exists() else 16
    fig,axes=plt.subplots(3,2,figsize=(12,9))
    for rowi,case in enumerate(["scalar_mixed","rod_kick","chain3"]):
        for method in methods:
            data=[r for r in traces if r["case"]==case and r["h"]==.08 and r["method"]==method and r["schedule"]=="material_storage_kinetic" and r["budget_multiplier"]==budget]
            axes[rowi,0].plot([r["t"] for r in data],[r["x0"] for r in data],label=shortnames[method])
            axes[rowi,1].semilogy([r["t"] for r in data],[max(1e-15,r["energy"]) for r in data],label=shortnames[method])
        oracle=[r for r in traces if r["case"]==case and r["h"]==.08 and r["method"]=="dense_oracle"]
        axes[rowi,0].plot([r["t"] for r in oracle],[r["x0"] for r in oracle],"--",color="black",label="Same-objective dense oracle")
        axes[rowi,1].semilogy([r["t"] for r in oracle],[max(1e-15,r["energy"]) for r in oracle],"--",color="black")
        reference=read("long_references.csv" if (root/"long_trajectories.csv").exists() else "short_references.csv")
        ref=[r for r in reference if r["case"]==case and r["h"]==.08]
        axes[rowi,0].plot([r["t"] for r in ref],[r["x0"] for r in ref],":",color="gray",label="Fine RK4 physical reference")
        axes[rowi,0].set_ylabel(case+"\nx₀ (m)")
        axes[rowi,1].set_ylabel("Relative energy error")
        for ax in axes[rowi]:ax.grid(True,alpha=.2);ax.set_xlabel("Time (s)")
    axes[0,0].legend(fontsize=8);fig.suptitle(f"Material → storage → kinetic; h = 0.08 s; {budget} × distinct-row visits per step")
    fig.tight_layout();fig.savefig(figdir/"trajectories.png",dpi=180);plt.close(fig)

cost=read("fair_cost.csv")
if cost:
    fig,axes=plt.subplots(1,2,figsize=(11,4.5),sharey=True)
    for ax,case in zip(axes,["chain_mixed","chain3"]):
        for method in methods:
            data=sorted([r for r in cost if r["case"]==case and r["h"]==.08 and r["method"]==method and r["schedule"]=="material_storage_kinetic"],key=lambda r:r["evaluation_ceiling"])
            ax.loglog([r["row_evaluations"] for r in data],[max(1e-15,r["full_residual"]) for r in data],"o-",label=shortnames[method])
        ax.set_title("Two springs" if case=="chain_mixed" else "Three springs")
        ax.set_xlabel("Actual row-gradient evaluations")
        ax.grid(True,alpha=.2)
    axes[0].set_ylabel("Full KKT residual")
    axes[0].legend(fontsize=8)
    fig.suptitle("Matched evaluation budgets; h = 0.08 s; material → storage → kinetic")
    fig.tight_layout();fig.savefig(figdir/"matched_cost.png",dpi=180);fig.savefig(figdir/"matched_cost.pdf");plt.close(fig)

(root/"summary.json").write_text(json.dumps(summary,indent=2,allow_nan=False)+"\n")
print(json.dumps(summary,indent=2,allow_nan=False))
