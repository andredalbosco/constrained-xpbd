#!/usr/bin/env python3
"""Focused EL comparison from existing measured traces; no new simulations."""
import csv
import os
from collections import defaultdict
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/xpbd-implicit-matplotlib")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np

ROOT = Path(__file__).resolve().parent
CASES = ["spring_release", "pendulum_rest", "spring_orbit", "double_small"]
METHODS = ["xpbd", "EL_angular_10", "workL_s1_tau005"]
LABELS = ["Ordinary XPBD", "Compliant EL + angular balance", "Local work + angular balance + repayment"]
COLORS = ["#7d8993", "#c17c19", "#13734d"]
TITLES = ["Radial spring", "Rigid pendulum", "Spring orbit", "Double spring: small swing"]


def main():
    with (ROOT / "results/summary.csv").open() as f:
        stats = {(r["case"], r["method"]): r for r in csv.DictReader(f) if r["hz"] == "480"}
    rows = defaultdict(list)
    with (ROOT / "results/traces.csv").open() as f:
        reader = csv.reader(f)
        header = next(reader)
        col = {name: i for i, name in enumerate(header[3:])}
        for r in reader:
            if r[0] in CASES and r[1] in METHODS and r[2] == "480":
                rows[r[0], r[1]].append(list(map(float, r[3:])))
    traces = {k: np.asarray(v) for k, v in rows.items()}
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10,
                         "axes.spines.top": False, "axes.spines.right": False})
    fig, axes = plt.subplots(4, 3, figsize=(16, 14))
    fig.subplots_adjust(top=.88, bottom=.075, left=.075, right=.985, hspace=.47, wspace=.28)
    fig.suptitle("EL remains useful, but local work balance is more consistent", x=.075, y=.985,
                 ha="left", fontsize=20, fontweight="bold")
    fig.text(.075, .956, "Existing C++ experiments · 480 small steps/s · 10 seconds · backward-Euler objective", fontsize=11)
    handles = [Line2D([], [], color=c, lw=2, label=l) for c, l in zip(COLORS, LABELS)]
    handles.append(Line2D([], [], color="black", ls=":", label="Independent physical reference"))
    fig.legend(handles=handles, loc="upper left", bbox_to_anchor=(.07, .941), ncol=2, frameon=False)
    for i, case in enumerate(CASES):
        for method, color in zip(METHODS, COLORS):
            a = traces[case, method]
            s = stats[case, method]
            assert s["status"] == "completed"
            v = lambda name: a[:, col[name]]
            t = v("time")
            scale, initial = float(s["energy_scale"]), float(s["initial_energy"])
            axes[i, 0].plot(t, 100 * (1 + (v("energy") - initial) / scale), color=color, lw=1.4)
            if case == "pendulum_rest":
                motion = np.degrees(np.arctan2(v("x1"), -v("y1")))
                ref = np.degrees(np.arctan2(v("ref_x1"), -v("ref_y1")))
                ylabel = "Angle from vertical (degrees)"
            elif case == "spring_release":
                motion, ref = 100 * (v("x1") - 1), 100 * (v("ref_x1") - 1)
                ylabel = "Spring extension (cm)"
            else:
                body = "2" if case == "double_small" else "1"
                motion, ref = 100 * v("x" + body), 100 * v("ref_x" + body)
                ylabel = ("Second mass" if body == "2" else "Mass") + " x (cm)"
            axes[i, 1].plot(t, motion, color=color, lw=1.3)
            n = int(s["count"])
            error = 100 * np.sqrt(sum((v(f"x{k}") - v(f"ref_x{k}")) ** 2
                                     + (v(f"y{k}") - v(f"ref_y{k}")) ** 2
                                     for k in range(1, n + 1)) / n)
            measured_rms = np.sqrt(np.mean((error[1:] / 100) ** 2))
            assert np.isclose(measured_rms, float(s["rms_position"]), atol=1e-10)
            axes[i, 2].plot(t, error, color=color, lw=1.4)
        axes[i, 0].axhline(100, color="black", ls=":", lw=1.1)
        axes[i, 0].set_ylim(0, 105)
        axes[i, 0].set_ylabel("Physical excitation energy (%)")
        axes[i, 1].plot(t, ref, color="black", ls=":", lw=1.2)
        axes[i, 1].set_ylabel(ylabel)
        axes[i, 2].set_ylabel("Position error vs reference (cm)")
        for j, title in enumerate(["Energy", "Motion", "Error"]):
            axes[i, j].set_title(f"{TITLES[i]} · {title}", loc="left", fontweight="bold")
            axes[i, j].set_xlabel("Time (s)")
            axes[i, j].grid(alpha=.2)
    fig.text(.075, .025,
             "Energy excludes reservoirs; 100% is the initial energy above minimum potential. Error includes all coordinates and masses.\n"
             "All shown methods completed both 96-case stress screens. Corrections use caps and an energy ceiling; universal completion is unproved.",
             fontsize=10, color="#333")
    output = ROOT / "plots/el_status.png"
    fig.savefig(output, dpi=140, facecolor="white")
    plt.close(fig)
    print(output)


if __name__ == "__main__":
    main()
