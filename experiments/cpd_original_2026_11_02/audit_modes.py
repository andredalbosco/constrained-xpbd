"""Postprocess archived radius windows and angular-row ablation; no solves."""
from pathlib import Path
import csv
import numpy as np

ROOT=Path(__file__).resolve().parent


def read(path):
    with path.open(newline="") as f:return list(csv.DictReader(f))


def main():
    collected=[]
    for folder in ("results","angular_ablation","pd_control"):
        summary=ROOT/folder/"summary.csv"
        if not summary.exists():continue
        for row in read(summary):
            if row["hz"]!="60":continue
            if row["method"]=="CPD" and float(row["epsilon_ratio"])!=1e8:continue
            states=read(ROOT/folder/"trajectories"/(row["run_id"]+".csv"))
            accepted=[s for s in states if s["accepted"]=="1"]
            for start in range(5):
                selected=[s for s in accepted if start<=float(s["time"])<=start+1+1e-10]
                if len(selected)<2:continue
                radii=np.array([float(s["mode_length"]) for s in selected])
                collected.append(dict(source=folder,run_id=row["run_id"],case=row["case"],family=row["family"],
                                      method=row["method"],policy=row["policy"],angular_active=row["angular_active"],
                                      window_start=start,window_end=min(start+1,float(row["accepted_time"])),
                                      min_radius=min(radii),max_radius=max(radii),radius_range=np.ptp(radii),
                                      max_radial_kinetic=max(float(s["radial_kinetic"]) for s in selected),
                                      final_radial_kinetic=float(selected[-1]["radial_kinetic"]),
                                      final_tangential_kinetic=float(selected[-1]["tangential_kinetic"]),
                                      final_angular=float(selected[-1]["angular_global"]),
                                      max_angular_target_error=max(float(s["angular_target_error"]) for s in selected),
                                      max_actual_energy_error=max(abs(float(s["rel_energy"])) for s in selected)))
    with (ROOT/"results"/"mode_windows.csv").open("w",newline="") as f:
        writer=csv.DictWriter(f,fieldnames=list(collected[0]));writer.writeheader();writer.writerows(collected)
    lines=["# Radius and angular-row audit","",
           "All values below reuse measured trajectories at60Hz and epsilon ratio1e8. Relative radius of each free pair is its spring separation; no continuous-reference score is used.","",
           "| Full free-body case | Policy | 4–5s radius range | Final radius | Final radial KE | Angular drift over5s | Max actual E error |",
           "|---|---|---:|---:|---:|---:|---:|"]
    summary=read(ROOT/"results"/"summary.csv")
    for case in ("free_orbit","free_stiff"):
        for policy in ("paper","fixed","stationary"):
            row=next(r for r in summary if r["case"]==case and r["policy"]==policy and r["hz"]=="60" and float(r["epsilon_ratio"])==1e8)
            window=next(r for r in collected if r["run_id"]==row["run_id"] and r["source"]=="results" and r["window_start"]==4)
            lines.append(f"| {case} | {policy} | {window['radius_range']:.5g}m | {window['max_radius']:.9g}m | "
                         f"{window['final_radial_kinetic']:.5g}J | {float(row['max_angular_drift']):.5g} | "
                         f"{100*float(row['max_accepted_rel_energy']):.5g}% |")
    lines += ["","The orbiting masses keep rotating, but fixed16/full-KKT policies leave the stretched spring radius nearly constant. Correct energy and angular momentum therefore coexist with frozen radial motion.","",
              "## Angular-row removal","","Kstar remains unchanged in this ablation. This isolates removal of the angular constraint; it is not the full original CPD method.","",
              "| Case | L row | Initial angular momentum | Final angular momentum | Max actual E error | Final tangential KE |",
              "|---|---|---:|---:|---:|---:|"]
    if (ROOT/"angular_ablation"/"summary.csv").exists():
        ablated=read(ROOT/"angular_ablation"/"summary.csv")
        for off in ablated:
            on=next(r for r in summary if r["case"]==off["case"] and r["policy"]=="fixed" and r["hz"]=="60" and float(r["epsilon_ratio"])==1e8)
            for row in (on,off):
                lines.append(f"| {row['case']} | {row['angular_active']} | {float(row['initial_angular']):.6g} | "
                             f"{float(row['final_angular']):.6g} | {100*float(row['max_accepted_rel_energy']):.5g}% | "
                             f"{float(row['final_tangential_kinetic']):.6g}J |")
    lines += ["","Removing the angular row permits spin-up while maintaining energy. Keeping it suppresses that spin-up but does not restore spring extension oscillation.",""]
    (ROOT/"results"/"MODE_AUDIT.md").write_text("\n".join(lines))
    print(f"Wrote {len(collected)} measured mode windows")


if __name__=="__main__":main()
