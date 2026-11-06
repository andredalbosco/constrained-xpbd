#!/usr/bin/env python3
"""Plot saved XPBD behavior probes; requires matplotlib. No solver changes."""
import csv
from pathlib import Path
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection

ROOT = Path(sys.argv[1]) if len(sys.argv) > 1 else Path.cwd()
OUT = ROOT / "figures"
OUT.mkdir(exist_ok=True)
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 11,
                     "axes.spines.top": False, "axes.spines.right": False,
                     "axes.grid": True, "grid.alpha": .18, "axes.titleweight": "bold",
                     "figure.facecolor": "white", "savefig.facecolor": "white"})
CASES = [("xpbd_paired_b1000", "xpbd_b1000", "XPBD", "#536273"),
         ("active_sequential_b1000", "sequential_b1000", "Sequential: β = γ = 1000", "#087e8b"),
         ("active_sequential_b1", "sequential_b1", "Sequential: β = γ = 1", "#ca4a35")]


def read(path):
    with Path(path).open() as file:
        return [{key: float(value) for key, value in row.items()} for row in csv.DictReader(file)]


fig, axes = plt.subplots(2, 3, figsize=(15, 8.4), sharex="col")
titles = ["Spring: oscillation grows", "Pendulum: feedback changes the swing", "2D grid + kick: motion dies out"]
for column, scene in enumerate(["spring", "pendulum", "grid_kick"]):
    for simple, mesh, label, color in CASES:
        if scene == "grid_kick":
            rows = read(ROOT / "grid_cloth" / f"grid_kick_h240_{mesh}" / "trajectory.csv")
            time = [r["time"] for r in rows]
            coordinate = [r["middle_y"] - rows[0]["middle_y"] for r in rows]
            energy = [r["mechanical_energy"] - rows[0]["mechanical_energy"] for r in rows]
        else:
            rows = read(ROOT / "spring_pendulum" / "trajectories" / f"{scene}_h1-240_{simple}.csv")
            time = [r["t"] for r in rows]
            coordinate = [r["x"] - 1 if scene == "spring" else r["angle"] for r in rows]
            energy = [r["mechanical"] - rows[0]["mechanical"] for r in rows]
        axes[0, column].plot(time, coordinate, color=color, lw=1.25, label=label)
        axes[1, column].plot(time, energy, color=color, lw=1.8)
    axes[0, column].set_title(titles[column], loc="left", pad=12, fontsize=12)
    axes[0, column].set_ylabel(["Extension (m)", "Angle from downward vertical (rad)", "Middle-particle vertical displacement (m)"][column])
    axes[1, column].axhline(0, color="#333333", lw=.9, ls="--")
    axes[1, column].set_ylabel("Mechanical energy change (J)")
    axes[1, column].set_xlabel("Time (s)")
    axes[1, column].set_xlim(0, 10)
fig.suptitle("Sequential work feedback: motion and mechanical energy", x=.055, ha="left", fontsize=20, fontweight="bold", y=.99)
fig.text(.055, .94, "Previous-force work · h = 1/240 s · 20 sweeps · reservoirs initially zero · smaller β, γ means stronger feedback", fontsize=11, color="#46515c")
handles, labels = axes[0, 0].get_legend_handles_labels()
fig.legend(handles, labels, loc="upper center", bbox_to_anchor=(.5, .916), ncol=3, frameon=False)
fig.subplots_adjust(left=.075, right=.98, bottom=.1, top=.82, hspace=.30, wspace=.30)
fig.text(.055, .025, "Dashed zero line: mechanical-energy conservation. Reservoir bookkeeping closes separately, even when motion is wrong.", color="#46515c", fontsize=10)
fig.savefig(OUT / "motion_energy.png", dpi=160)
fig.savefig(OUT / "motion_energy.svg")
plt.close(fig)


def edges():
    for row in range(6):
        for col in range(6):
            i = row * 6 + col
            if col < 5:
                yield i, i + 1
            if row < 5:
                yield i, i + 6
            if row < 5 and col < 5:
                yield i, i + 7
                yield i + 1, i + 6


fig, axes = plt.subplots(1, 3, figsize=(13.6, 6.0), sharex=True, sharey=True)
for ax, (_, mesh, label, color) in zip(axes, CASES):
    folder = ROOT / "grid_cloth" / f"grid_kick_h240_{mesh}"
    rows = read(folder / "snapshots.csv")
    initial = {int(r["particle"]): r for r in rows if r["time"] == 0}
    final = {int(r["particle"]): r for r in rows if r["time"] == 10}
    for state, linestyle, alpha, width in [(initial, "--", .3, .7), (final, "-", 1, 1.3)]:
        segments = [[(state[a]["x"], state[a]["y"]), (state[b]["x"], state[b]["y"])] for a, b in edges()]
        ax.add_collection(LineCollection(segments, colors=color if state is final else "#536273", linewidths=width, linestyles=linestyle, alpha=alpha))
    pins = [r for r in final.values() if r["pinned"]]
    ax.scatter([r["x"] for r in pins], [r["y"] for r in pins], color="#20252a", s=30, zorder=5)
    with (folder / "summary.csv").open() as file:
        summary = next(csv.DictReader(file))
    ax.set_title(label.replace(": ", "\n"), color=color, fontsize=12)
    ax.set_aspect("equal")
    ax.set_xlim(-.83, .83)
    ax.set_ylim(.03, 1.68)
    ax.set_xlabel("x (m)")
    ax.text(.02, -.23, f"ΔE = {float(summary['final_energy_change']):+.3f} J\nFinal kinetic energy = {float(summary['final_kinetic_energy']):.2g} J", transform=ax.transAxes, fontsize=11)
axes[0].set_ylabel("y (m)")
fig.suptitle("Kicked 2D grid at 10 s: strong feedback leaves stored deformation", x=.055, ha="left", y=.98, fontsize=17, fontweight="bold")
fig.text(.055, .91, "6 × 6 particles · h = 1/240 s · 20 sweeps · initial kinetic energy 1.25 J · dotted mesh = initial shape", fontsize=11, color="#46515c")
fig.subplots_adjust(left=.055, right=.98, top=.80, bottom=.24, wspace=.18)
fig.savefig(OUT / "grid_shapes.png", dpi=160)
fig.savefig(OUT / "grid_shapes.svg")
plt.close(fig)


fig, axes = plt.subplots(1, 2, figsize=(12.8, 5.4), sharex=True)
for order, beta, label, color in [("xpbd", 1000, "XPBD", "#536273"),
                                  ("sequential", 1000, "Sequential: β = γ = 1000", "#087e8b"),
                                  ("sequential", 1, "Sequential: β = γ = 1", "#ca4a35")]:
    mode = "xpbd" if order == "xpbd" else "active-sequential"
    rows = read(ROOT / "grid_cloth" / f"cloth_h240_{order}_b{beta}" / f"cloth-{mode}.csv")
    time = [r["time"] for r in rows]
    axes[0].plot(time, [r["mechanical_energy"] - rows[0]["mechanical_energy"] for r in rows], label=label, color=color, lw=1.7)
    axes[1].plot(time, [100*r["relative_stretch"] for r in rows], color=color, lw=1.7)
axes[0].axhline(0, color="#333333", lw=.9, ls="--")
axes[0].set_ylabel("Mechanical energy change (J)")
axes[1].set_ylabel("Maximum relative spring-length error (%)")
for ax in axes:
    ax.set_xlim(0, 10)
    ax.set_xlabel("Time (s)")
fig.suptitle("Cloth: weaker damping can turn into excessive stored deformation", x=.065, ha="left", y=.99, fontsize=17, fontweight="bold")
fig.text(.065, .91, "9 × 9 particles · h = 1/240 s · 20 sweeps · two pinned corners · previous-force work", fontsize=11, color="#46515c")
handles, labels = axes[0].get_legend_handles_labels()
fig.legend(handles, labels, loc="upper center", bbox_to_anchor=(.5, .87), ncol=3, frameon=False)
fig.subplots_adjust(left=.085, right=.98, bottom=.14, top=.73, wspace=.26)
fig.savefig(OUT / "cloth_behavior.png", dpi=160)
fig.savefig(OUT / "cloth_behavior.svg")
plt.close(fig)
print(OUT)
