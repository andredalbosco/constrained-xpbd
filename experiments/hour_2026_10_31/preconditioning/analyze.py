#!/usr/bin/env python3
"""Read measured CSVs and regenerate publication-style PNGs and RESULTS.md."""
import csv
import math
from pathlib import Path
from statistics import median
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parent
FIG=ROOT/"figures"
FIG.mkdir(exist_ok=True)
plt.rcParams.update({"font.size":10,"axes.spines.top":False,"axes.spines.right":False,
                     "axes.grid":True,"grid.alpha":.18,"figure.dpi":150,"savefig.dpi":180})

def read(name):
    with (ROOT/name).open() as f:return list(csv.DictReader(f))
def n(row,key):
    try:return float(row[key])
    except (ValueError,KeyError):return float("nan")
def select(rows,**keys):return [r for r in rows if all(str(r[k])==str(v) for k,v in keys.items())]
def fmt(v):
    if not math.isfinite(v):return "unavailable"
    if v==0:return "0"
    return f"{v:.3g}"
def save(fig,name):
    fig.tight_layout();fig.savefig(FIG/name,bbox_inches="tight");plt.close(fig)
def table(headers,rows):
    return "| "+" | ".join(headers)+" |\n|"+"|".join("---" for _ in headers)+"|\n"+"\n".join("| "+" | ".join(map(str,r))+" |" for r in rows)+"\n"

runs=read("runs.csv")
curvature=read("curvature_graph.csv")
for r in curvature:
    r["metric"]="curvature_graph";r["steps"]=str(round(.4/n(r,"h")))
allruns=runs+curvature
timing=read("timing.csv")
branches=read("branches.csv")
physical=read("physical.csv")
physical_refs=read("physical_references.csv")
ablation=read("krylov_ablation.csv")
forcing=read("forcing.csv")
cache=read("cache.csv")
refinement=read("phase_refinement_roundoff.csv")
strict=read("phase_refinement.csv")
labels={"none":"Identity","mass":"Mass","diagonal":"Diagonal","graph":"Isotropic\ngraph",
        "vertex_jacobi":"Vertex\nleft inverse","vertex_rows":"Vertex\nrow balance","vertex_right":"Vertex\nright inverse",
        "edge_rows":"Edge\nrow balance","edge_right":"Edge\nright balance","curvature_graph":"Curvature\ngraph","direct_block":"Direct block"}
order=list(labels)[:-1]
colors={"diagonal":"#777777","graph":"#e3932e","vertex_rows":"#2768ad","vertex_right":"#197b6b","edge_rows":"#a67dbc","edge_right":"#8b4aa0","curvature_graph":"#c44f48","direct_block":"#222222"}

fig,ax=plt.subplots(figsize=(13,4.8))
xx=np.arange(len(order));width=.37
counts=[]
for offset,budget in [(-width/2,"12"),(width/2,"48")]:
    values=[sum(r["converged"]==r["steps"] for r in select(allruns,metric=m,budget=budget)) for m in order]
    bars=ax.bar(xx+offset,values,width,label=f"{budget} Krylov products / step",color="#6ba3c8" if budget=="12" else "#284f7a")
    ax.bar_label(bars,padding=3,fontsize=9)
ax.axhline(42,color="black",ls=":",lw=1);ax.set(xticks=xx,xticklabels=[labels[m] for m in order],ylim=(0,47),ylabel="Fully converged trajectories / 42",title="Cheap factors help, but none makes every hard trajectory converge")
ax.legend(loc="upper left",ncol=2);save(fig,"convergence.png")

fig,axes=plt.subplots(1,2,figsize=(12,4.5))
selected=["diagonal","graph","vertex_right","edge_right","direct_block"]
for m in selected:
    ys=[]
    for N in [8,32,128,256]:
        rows=select(timing,N=N,regime="mass100",metric=m,budget=0 if m=="direct_block" else 48)
        ys.append(median(n(r,"cpu_ms") for r in rows))
    axes[0].plot([8,32,128,256],ys,"o-",color=colors[m],label=labels[m].replace("\n"," "))
axes[0].set(xscale="log",yscale="log",xlabel="Chain vertices",ylabel="Median CPU ms / ten steps",title="Mass contrast 100, h = 0.04, budget 48")
axes[0].set_xticks([8,32,128,256],labels=[8,32,128,256]);axes[0].legend(fontsize=8)
setup=[];rest=[]
for m in selected:
    rows=select(timing,N=128,regime="mass100",metric=m,budget=0 if m=="direct_block" else 48)
    st=median(n(r,"setup_us")/1000 for r in rows);total=median(n(r,"cpu_ms") for r in rows);setup.append(st);rest.append(total-st)
axes[1].bar(np.arange(len(selected)),setup,label="Metric setup / direct factor-and-solve",color="#e49a50")
axes[1].bar(np.arange(len(selected)),rest,bottom=setup,label="Remaining solver work",color="#467aa8")
axes[1].set(xticks=np.arange(len(selected)),xticklabels=[labels[m] for m in selected],ylabel="Median CPU ms / ten steps",title="N = 128: patch SVD setup dominates")
axes[1].legend(fontsize=8);save(fig,"isolated_cost.png")

fig,axes=plt.subplots(1,2,figsize=(11.5,4.4))
for policy,label,color in [("transformed","Transformed residual","#2768ad"),("physical_stop","Physical stopping only","#e3932e"),("physical_min","Physical QR minimum","#197b6b")]:
    rows=select(ablation,N=32,regime="mass100",metric="vertex_rows",budget=12,policy=policy)
    rows.sort(key=lambda r:n(r,"chunk"));axes[0].plot([n(r,"chunk") for r in rows],[n(r,"max_energy") for r in rows],"o-",label=label,color=color)
axes[0].set(yscale="log",xlabel="Maximum Krylov products before relinearization",ylabel="Max relative energy defect",title="N32 mass100: same 12-product budget")
axes[0].set_xticks([2,4,6,12]);axes[0].legend(fontsize=8)
policies=["fixed_transformed","fixed_physical","loose_physical","ew2_physical"]
policylabels=["Fixed 0.001\ntransformed","Fixed 0.001\nphysical","Loose 0.1\nphysical","Adaptive EW2\nphysical"]
rows=[select(forcing,N=32,regime="mass100",contrast=100,h=.04,metric="vertex_rows",budget=12,policy=p)[0] for p in policies]
axes[1].bar(np.arange(4),[n(r,"max_energy") for r in rows],color=["#777777","#e3932e","#6ca295","#197b6b"])
axes[1].set(yscale="log",xticks=np.arange(4),xticklabels=policylabels,ylabel="Max relative energy defect",title="Adaptive forcing spends more work on model refreshes")
for i,r in enumerate(rows):axes[1].annotate(f"{int(n(r,'outer'))} Newton\nrefreshes",(i,n(r,"max_energy")),xytext=(0,5),textcoords="offset points",ha="center",fontsize=8)
axes[1].set_ylim(1e-5,1);save(fig,"nonlinear_budget.png")

fig,axes=plt.subplots(1,3,figsize=(14,4.3))
preds=["ballistic","acceleration","secant","continuation"]
vals=[sum(r["completed"]==r["requested"] for r in select(branches,predictor=p)) for p in preds]
bars=axes[0].bar(preds,vals,color=["#777777","#aaa","#e3932e","#197b6b"]);axes[0].bar_label(bars,padding=3);axes[0].set(ylim=(0,40),ylabel="Completed, converged trajectories / 36",title="Following a root branch helps");axes[0].tick_params(axis="x",rotation=20)
traces=read("branch_traces.csv")
for pred,color in [("ballistic","#777777"),("continuation","#197b6b")]:
    rows=select(traces,N=32,regime="both100",contrast=100,motion="mixed",h=.04,predictor=pred)
    t=[n(r,"time") for r in rows]
    axes[1].plot(t,[max(1e-16,n(r,"energy_error")) for r in rows],"o-",label=pred,color=color)
    axes[2].plot(t,[n(r,"phase_error") for r in rows],"o-",label=pred,color=color)
axes[1].set(yscale="log",xlabel="Time (s)",ylabel="Relative energy defect",title="N32 both100, h = 0.04")
axes[2].set(xlabel="Time (s)",ylabel="Modal phase error (rad)",title="Conservation does not repair phase")
for ax in axes[1:]:ax.axvline(.16,color="#777777",ls=":",lw=1);ax.legend(fontsize=8)
save(fig,"continuation_and_phase.png")

fig,axes=plt.subplots(1,2,figsize=(11,4.4))
x=[n(r,"h_omega_max") for r in refinement];y=[n(r,"phase_error") for r in refinement]
axes[0].loglog(x,y,"o-",label="Recorded physical-reference phase error",color="#2768ad")
axes[0].loglog(x[-4:],[y[-1]*(xx/x[-1])**2 for xx in x[-4:]],"--",color="#444444",label="Second-order guide")
axes[0].set(xlabel="h × highest linearized frequency",ylabel="Endpoint modal phase error (rad)",title="Small-amplitude chain, physical time 0.02 s");axes[0].legend(fontsize=8)
axes[1].semilogx([n(r,"h") for r in strict],[n(r,"completed")/n(r,"steps") for r in strict],"o-",color="#c44f48",label="Strict residual tolerance 1e−10")
axes[1].semilogx([n(r,"h") for r in refinement],[n(r,"completed")/n(r,"steps") for r in refinement],"s--",color="#197b6b",label="Estimated cancellation floor")
axes[1].set(xlabel="Timestep (s)",ylabel="Fraction of requested steps converged",title="Absolute-coordinate subtraction reaches a residual floor",ylim=(-.05,1.15));axes[1].legend(fontsize=8)
save(fig,"phase_refinement.png")

followup=read("timing_followup.csv") if (ROOT/"timing_followup.csv").exists() else []
if followup:
    fig,axes=plt.subplots(1,2,figsize=(12,4.4))
    variants=["vertex_rows","vertex_chunk2","vertex_ew2","edge_right","edge_step_cache","curvature_graph","curvature_step_cache","curvature_ew2","direct_block"]
    for i,v in enumerate(variants):
        rows=select(followup,N=128,regime="mass100",variant=v,budget=0 if v=="direct_block" else 12)
        cpu=median(n(r,"cpu_ms") for r in rows);en=median(n(r,"max_energy") for r in rows);conv=median(n(r,"converged") for r in rows)
        axes[0].scatter(cpu,en,s=55,marker="o" if conv==10 else "x",color=plt.cm.tab10(i%10));axes[0].annotate(v.replace("_"," "),(cpu,en),xytext=(4,3),textcoords="offset points",fontsize=7)
    axes[0].set(xscale="log",yscale="log",xlabel="Isolated median CPU ms / ten steps",ylabel="Max relative energy defect",title="N128 mass100, budget 12; × = inexact horizon")
    variants=["edge_right","edge_step_cache","curvature_graph","curvature_step_cache"]
    total=[];setup=[]
    for v in variants:
        rows=select(followup,N=128,regime="mass100",variant=v,budget=48);total.append(median(n(r,"cpu_ms") for r in rows));setup.append(median(n(r,"setup_us")/1000 for r in rows))
    axes[1].bar(range(4),setup,label="Metric setup",color="#e49a50");axes[1].bar(range(4),np.array(total)-np.array(setup),bottom=setup,label="Other solver work",color="#467aa8")
    axes[1].set(xticks=range(4),xticklabels=["Edge fresh","Edge cached\nper step","Curvature fresh","Curvature cached\nper step"],ylabel="Isolated median CPU ms / ten steps",title="N128 mass100, budget 48: setup reuse")
    axes[1].legend(fontsize=8);save(fig,"followup_cost.png")

identities=read("identity_tests.csv")
other_audits=[read(n) for n in ["krylov_identity.csv","forcing_identity.csv","cache_identity.csv","curvature_graph_identity.csv"]]
report=["# Measured results: cheap preconditioning and continuation\n",
"Local 2×2 factors make the SVD-inspired direction balancing inexpensive. They improve convergence at the larger budget, but tight linear solves can waste a short budget before the nonlinear model is refreshed. Adaptive inexact-Newton forcing helps that tradeoff. For these serial chains, an exact O(N) block Newton solve remains an especially strong control.\n",
"Continuation solves more hard discrete steps; it does not correct their under-resolved physical phase. These are separate findings. Full equations, budgets, and verified primary-paper links are in [README.md](README.md).\n",
"## Validation and scope\n",
f"The main sweep contains {len(runs)} metric/direct trajectories over 42 scenarios; the curvature graph adds {len(curvature)}. All {len(runs)} main trajectories remained finite through 0.4 s, but only {sum(r['converged']==r['steps'] for r in runs)} converged at every step. Thirteen of the 42 same-timestep direct reference trajectories stop early; 247 main comparison rows therefore have only matched-prefix trajectory errors. Missing continuation/reference intervals are not counted as accurate motion.\n",
f"All {len(identities)} principal metric/root audits pass. Maximum directional Jacobian error is {fmt(max(n(r,'jacobian_error') for r in identities))}; block versus dense solve difference {fmt(max(n(r,'block_solve_error') for r in identities))}; off-root energy identity defect {fmt(max(n(r,'energy_identity') for r in identities))}; converged-root difference {fmt(max(n(r,'root_difference') for r in identities))}. Root difference here is mass-weighted position plus h times mass-weighted velocity difference. Additional audits cover {sum(map(len,other_audits))} Krylov-policy, forcing, cached-factor, and curvature-graph cases.\n",
"## Convergence at fixed Krylov budgets\n"]
rows=[]
for m in order:
    rows.append([m,*[f"{sum(r['converged']==r['steps'] for r in select(allruns,metric=m,budget=b))}/42" for b in [12,48]]])
report.append(table(["Metric","Budget 12","Budget 48"],rows))
report.append("The direct block control fully converges in 29/42 cases. Its allowance is larger and its actual cost is measured; it is not credited with an artificial Krylov budget.\n\n![Convergence](figures/convergence.png)\n")
report.append("## Isolated setup and solve cost\n")
rows=[]
for m in selected:
    rs=select(timing,N=128,regime="mass100",metric=m,budget=0 if m=="direct_block" else 48)
    rows.append([m,fmt(median(n(r,"cpu_ms") for r in rs)),fmt(median(n(r,"setup_us")/1000 for r in rs)),f"{int(n(rs[0],'converged'))}/10"])
report.append("Five serial repeats with rotated method order, all sibling numerical jobs paused/finished. Values below are N128, mass contrast 100, h = .04, ten steps; iterative budget 48. CPU includes setup.\n")
report.append(table(["Method","CPU ms","Setup ms","Converged steps"],rows))
report.append("For iterative methods, setup is preconditioner construction. For the direct control, that column includes the complete interleaved block factor-and-solve routine. Total CPU is comparable; the setup columns have these different scopes.\n")
report.append("At N256, vertex-right Jacobi takes "+fmt(median(n(r,"cpu_ms") for r in select(timing,N=256,regime="mass100",metric="vertex_right",budget=48)))+" ms and diagonal scaling "+fmt(median(n(r,"cpu_ms") for r in select(timing,N=256,regime="mass100",metric="diagonal",budget=48)))+" ms. Their convergence differs, so time alone would miss the useful result. The direct block control takes "+fmt(median(n(r,"cpu_ms") for r in select(timing,N=256,regime="mass100",metric="direct_block",budget=0)))+" ms. General meshes need a separate study.\n\n![CPU costs](figures/isolated_cost.png)\n")
report.append("## Why short budgets can get worse\n")
report.append("For N32 mass100 at h=.04, the first four-product vertex-row linear solve has physical relative residual 0.00872, compared with 0.0325 for diagonal scaling. Nevertheless, the propagated budget 12 trajectory is worse. Moving from four-product to two-product chunks reduces its maximum relative energy defect from 0.173 to 1.66e−4, while allowing more nonlinear model refreshes. It remains unconverged. Physical-residual stopping alone does not change this example. A physical QR minimum improves some errors but adds cost and is not uniformly better.\n")
rows=[]
for policy in ["fixed_transformed","fixed_physical","loose_physical","ew2_physical"]:
    rows.append([policy,*[f"{sum(r['converged']==r['requested'] for r in select(forcing,budget=b,policy=policy))}/252" for b in [12,48]]])
report.append(table(["Linear forcing / residual policy","Budget 12 fully converged","Budget 48 fully converged"],rows))
report.append("The EW2 example uses the same 120 total products over ten steps, but increases Newton refreshes from 30 to 50 and reduces maximum energy defect to 3.61e−5. Setup work grows. Across 42 cases, adaptive forcing is not uniformly better: the diagonal metric at budget 12 falls from 4 fully converged horizons to 3. This is evidence for distributing work across nonlinear models, not a universal cure or a claim that the resulting inexact states are valid solves.\n\n![Nonlinear budget](figures/nonlinear_budget.png)\n")
if followup:
    report.append("## Follow-up isolated cost and cached factors\n")
    rows=[]
    for v in ["vertex_rows","vertex_chunk2","vertex_ew2","edge_right","edge_step_cache","curvature_graph","curvature_step_cache","curvature_ew2","direct_block"]:
        rs=select(followup,N=128,regime="mass100",variant=v,budget=0 if v=="direct_block" else 12)
        rows.append([v,fmt(median(n(r,"cpu_ms") for r in rs)),fmt(median(n(r,"max_energy") for r in rs)),f"{int(median(n(r,'converged') for r in rs))}/10"])
    report.append("Seven rotated repeats in a second isolated window. N128 mass100, h=.04, ten steps, budget 12 unless direct.\n")
    report.append(table(["Variant","CPU ms","Max energy defect","Converged steps"],rows))
    report.append("At budget 12, curvature-aware graph scaling plus adaptive forcing converges all ten steps in 1.06 ms in this case, compared with 0.259 ms for the direct control. At budget 48, caching edge factors once per step reduces CPU from 4.99 ms to 1.52 ms while retaining ten converged steps. Caching holds only the metric; every nonlinear residual and Jacobian is still updated. The dedicated cache sweep includes regressions, for example per-step vertex-right caching changes 10/10 to 9/10 convergence on N128 mass100 at budget 48. No unconditional cache recommendation follows.\n\n![Follow-up costs](figures/followup_cost.png)\n")
report.append("## Branch continuation and physical phase\n")
report.append(table(["Guess / continuation","Fully converged horizons"],[[p,f"{sum(r['completed']==r['requested'] for r in select(branches,predictor=p))}/36"] for p in preds]))
report.append("Continuation agrees with the existing ballistic root on recovered prefixes to about 1e−12 in the nominal contrast 100 examples. Four hard cases find materially different roots on the shared prefix. Two continuation cases still fail: N32 both100, large motion, h=.1 (2/4 steps), and N32 both1000, mixed, h=.01 (23/40). Failed stages are included in cost and work.\n")
report.append(f"All nine physical RK4 references pass their paired-resolution checks: worst position gap {fmt(max(n(r,'position_gap') for r in physical_refs))} m and scaled velocity gap {fmt(max(n(r,'velocity_gap_scaled') for r in physical_refs))}. For N32 both100 at h=.04, continuation finishes 10/10 with about 7.94e−14 energy defect, yet physical RMS velocity error is 3.08 m/s and maximum modal phase error 2.80 rad. Comparing the completed continuation horizon with only the ballistic failure prefix would be misleading.\n\n![Continuation and phase](figures/continuation_and_phase.png)\n")
report.append("## Timestep and floating-point limits\n")
report.append("A separate 0.02 s probe scales the N8 both100 perturbation to 0.005 of its original amplitude. With the labeled cancellation-aware stopping floor, the last three endpoint phase errors are 0.07537, 0.01888, 0.004724 rad under timestep halving: approximately second-order physical phase convergence. The strict 1e−10 residual test stops early at h=1.25e−4 and produces zero accepted steps at smaller h, due to cancellation in absolute-position subtraction. Both versions are retained. The relaxed probe's largest actual residual is 8.69e−10 and maximum relative energy defect 3.86e−9; the estimated tolerance ceiling is looser than those achieved values.\n\n![Phase refinement and stopping](figures/phase_refinement.png)\n")
report.append("## What to carry forward\n\n- Use the 2×2 right block inverse as a cheap iterative control; preserve the direct block control on serial chains.\n- Couple preconditioning experiments to nonlinear refresh frequency and forcing tolerances. A fixed tiny linear tolerance can spend the budget in the wrong place.\n- Cache expensive patch factors only with explicit convergence checks; keep measuring setup.\n- Use fixed-old-state continuation to investigate a stalled discrete root, while checking branch changes and physical references separately.\n- Treat phase resolution and floating-point stopping as separate problems from energy preservation.\n")
(ROOT/"RESULTS.md").write_text("\n".join(report))
print(f"Wrote RESULTS.md and {len(list(FIG.glob('*.png')))} PNG figures.")
