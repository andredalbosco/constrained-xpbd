"""Render archived trajectories using motion and physical-energy criteria.

No numerical simulation or continuous reference is evaluated here. All inputs
are read from the immutable block-orders campaign, and all outputs live here.
"""
from __future__ import annotations

import csv
from pathlib import Path
import argparse

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, FFMpegWriter


HERE = Path(__file__).resolve().parent
SOURCE = Path("/home/blacknoir/Projects/constrained-xpbd/experiments/block_orders_2026_11_01/results")
CASES = ("spring_release", "spring_orbit", "spring_stiff", "spring_gravity", "pendulum_release", "pendulum_kick")
COLORS = {"EL": "#d55e00", "LE": "#0072b2", "JOINT": "#7b2cbf"}
DISPLAY = {"EL": "Energy → material", "LE": "Material → energy", "JOINT": "Joint solve"}


def csv_rows(path):
    with path.open(newline="") as f:
        return list(csv.DictReader(f))


def numeric(row, field):
    return float(row[field])


def select(summary, case, order, hz, sweeps, ratio):
    found = [r for r in summary if r["case"] == case and r["order"] == order
             and int(r["hz"]) == hz and int(r["sweeps"]) == sweeps
             and float(r["epsilon_ratio"]) == ratio]
    if len(found) != 1:
        raise ValueError(f"Expected one run: {case} {order} {hz} {sweeps} {ratio}: {len(found)}")
    row = found[0]
    records = csv_rows(SOURCE / "trajectories" / (row["run_id"] + ".csv"))
    data = {key: np.array([float(r[key]) for r in records]) for key in records[0]}
    accepted = data["accepted"] > .5
    accepted_data = {k: v[accepted] for k, v in data.items()}
    length = np.hypot(accepted_data["x"], accepted_data["y"])
    speed = np.hypot(accepted_data["vx"], accepted_data["vy"])
    path = float(np.linalg.norm(np.diff(np.c_[accepted_data["x"], accepted_data["y"]], axis=0), axis=1).sum())
    return dict(row=row, data=data, accepted=accepted_data, length=length, speed=speed,
                path=path, initial_energy=data["energy"][0])


def index_at(run, time):
    return max(0, min(len(run["accepted"]["time"]) - 1,
                      np.searchsorted(run["accepted"]["time"], time + 1e-10, side="right") - 1))


def spring_points(p):
    length = np.linalg.norm(p)
    if length <= 1e-12:
        return [0, p[0]], [0, p[1]]
    direction = p / length
    perpendicular = np.array([-direction[1], direction[0]])
    along = np.r_[0., .08, np.linspace(.1, .9, 19), .92, 1.]
    side = np.zeros(len(along))
    side[2:-2] = .035 * np.where(np.arange(19) % 2, -1., 1.)
    points = along[:, None] * p + side[:, None] * perpendicular
    return points[:, 0], points[:, 1]


def extent(runs):
    radius = max(np.max(run["length"]) for run in runs)
    return max(1.1, min(radius * 1.12, 3.))


def motion_axes(ax, radius, *, labels=False):
    ax.set_xlim(-radius, radius)
    ax.set_ylim(-radius, radius)
    ax.set_aspect("equal")
    ax.plot(0, 0, "ks", ms=4, zorder=6)
    ax.axhline(0, color=".85", lw=.6, zorder=-2)
    ax.axvline(0, color=".85", lw=.6, zorder=-2)
    if labels:
        ax.set_xlabel("x (m)")
        ax.set_ylabel("y (m)")
    else:
        ax.set_xticks([])
        ax.set_yticks([])


def status(run, t):
    row = run["row"]
    if row["completed"] != "1" and t > float(row["accepted_time"]) + 1e-8:
        return f"FAILED at {float(row['accepted_time']):.3g} s · accepted prefix held", True
    if run["path"] < 1e-4:
        return "FROZEN motion · fails motion criterion", True
    i = index_at(run, t)
    j = index_at(run, max(0., t - .25))
    data = run["accepted"]
    travel = np.linalg.norm(np.diff(np.c_[data["x"][j:i + 1], data["y"][j:i + 1]], axis=0), axis=1).sum()
    if t >= .5 and travel < 1e-4:
        return "Nearly stationary over last 0.25 s", True
    return "", False


def contact_sheet(group, case, runs, out):
    times = (0., .25, .5, 1., 2., 3., 4., 5.)
    radius = extent(runs)
    fig, axes = plt.subplots(len(runs), len(times), figsize=(18, 3 * len(runs)),
                              squeeze=False, layout="constrained")
    for row_index, run in enumerate(runs):
        order = run["row"]["order"]
        data = run["accepted"]
        for col, t in enumerate(times):
            ax = axes[row_index, col]
            motion_axes(ax, radius)
            i = index_at(run, t)
            p = np.array([data["x"][i], data["y"][i]])
            xx, yy = spring_points(p) if case.startswith("spring") else ([0, p[0]], [0, p[1]])
            ax.plot(data["x"][:i + 1], data["y"][:i + 1], color=COLORS[order], alpha=.25, lw=.7)
            ax.plot(xx, yy, color=COLORS[order], lw=1.4)
            ax.plot(*p, "o", color=COLORS[order], ms=6)
            ax.set_title(f"t = {t:g} s", fontsize=10)
            if col == 0:
                ax.set_ylabel(DISPLAY[order], fontsize=10)
            text = (f"length {run['length'][i]:.4f} m\n"
                    f"speed {run['speed'][i]:.3g} m/s\n"
                    f"E {data['energy'][i]:.6g} J\n"
                    f"ΔE {100 * data['rel_energy'][i]:+.3g}%")
            warning, bad = status(run, t)
            if warning:
                text += "\n" + ("FAILED · prefix held" if warning.startswith("FAILED")
                                  else "FROZEN" if warning.startswith("FROZEN") else "Nearly stationary")
            ax.text(.02, .98, text, transform=ax.transAxes, va="top", fontsize=8,
                    color="#a40000" if bad else "#222222",
                    bbox=dict(facecolor="white", alpha=.8, edgecolor="none", pad=1.))
    row = runs[0]["row"]
    fig.suptitle(f"{case.replace('_', ' ').capitalize()} · {row['hz']} Hz · {row['sweeps']} sweeps · "
                 f"ε / E-scale = {float(row['epsilon_ratio']):g}\n"
                 "Snapshots of archived accepted states; black square is the fixed pivot", fontsize=13)
    fig.savefig(out / f"{group}_{case}_contacts.png", dpi=130)
    plt.close(fig)


def path_and_energy(group, selections, out):
    fig, axes = plt.subplots(6, 3, figsize=(14, 17), layout="constrained")
    for row_index, (case, runs) in enumerate(selections.items()):
        radius = extent(runs)
        motion_axes(axes[row_index, 0], radius, labels=True)
        for run in runs:
            row = run["row"]
            order = row["order"]
            data = run["accepted"]
            label = order + (" (failed prefix)" if row["completed"] != "1" else "")
            axes[row_index, 0].plot(data["x"], data["y"], color=COLORS[order], lw=1.1,
                                    alpha=.8, label=label)
            axes[row_index, 0].plot(data["x"][-1], data["y"][-1], "o", color=COLORS[order], ms=4)
            axes[row_index, 1].plot(data["time"], run["length"], color=COLORS[order], lw=1., label=label)
            axes[row_index, 2].plot(data["time"], 100 * data["rel_energy"], color=COLORS[order], lw=1., label=label)
            if row["completed"] != "1":
                for ax in axes[row_index, 1:]:
                    ax.axvline(float(row["accepted_time"]), ls=":", lw=.8, color=COLORS[order])
        axes[row_index, 2].axhline(0., color="black", ls="--", lw=.8,
                                  label="Initial physical energy")
        axes[row_index, 0].set_title(case.replace("_", " ").capitalize())
        axes[row_index, 1].set_ylabel("Length from fixed pivot (m)")
        axes[row_index, 1].ticklabel_format(axis="y", style="plain", useOffset=False)
        axes[row_index, 2].set_ylabel("Actual energy error (%)")
        axes[row_index, 2].ticklabel_format(axis="y", style="plain", useOffset=False)
        for ax in axes[row_index, 1:]:
            ax.set_xlabel("Time (s)")
            ax.grid(alpha=.2)
        axes[row_index, 2].legend(fontsize=7, loc="best")
    row = next(iter(selections.values()))[0]["row"]
    fig.suptitle(f"Archived accepted motion and actual energy · {row['hz']} Hz · {row['sweeps']} sweeps · "
                 f"ε / E-scale = {float(row['epsilon_ratio']):g}\n"
                 "Energy error = 100 × [E(t) − E(0)] / [E(0) − minimum U]; incomplete prefixes retain stop markers", fontsize=12)
    fig.savefig(out / f"{group}_paths_length_energy.png", dpi=130)
    plt.close(fig)


def movie(group, selections, out, fps):
    order_count = len(next(iter(selections.values())))
    columns = 2 if order_count == 3 else 3
    rows = 3 if order_count == 3 else 2
    fig = plt.figure(figsize=(18, 15 if order_count == 3 else 12), layout="constrained")
    outer = fig.add_gridspec(rows, columns)
    moving = []
    for case_index, (case, runs) in enumerate(selections.items()):
        inner = outer[case_index // columns, case_index % columns].subgridspec(2, len(runs), height_ratios=(3., 1.))
        radius = extent(runs)
        energy_ax = fig.add_subplot(inner[1, :])
        energy_ax.axhline(runs[0]["initial_energy"], color="black", ls="--", lw=.8)
        energy_ax.set_xlabel("Time (s)", fontsize=8)
        energy_ax.set_ylabel("Physical E (J)", fontsize=8)
        energy_ax.tick_params(labelsize=7)
        energy_ax.ticklabel_format(axis="y", style="plain", useOffset=False)
        energy_ax.grid(alpha=.2)
        energy_ax.set_xlim(0, 5.)
        energies = np.concatenate([run["accepted"]["energy"] for run in runs])
        e_range = max(np.ptp(energies), abs(runs[0]["initial_energy"]) * 1e-5, 1e-6)
        energy_ax.set_ylim(np.min(energies) - .1 * e_range, np.max(energies) + .1 * e_range)
        for order_index, run in enumerate(runs):
            row = run["row"]
            order = row["order"]
            color = COLORS[order]
            ax = fig.add_subplot(inner[0, order_index])
            motion_axes(ax, radius)
            ax.set_title(case.replace("_", " ") + "\n" + DISPLAY[order], fontsize=9)
            line, = ax.plot([], [], color=color, lw=1.4)
            trail, = ax.plot([], [], color=color, lw=.9, alpha=.3)
            dot, = ax.plot([], [], "o", color=color, ms=7)
            text = ax.text(.02, .98, "", transform=ax.transAxes, va="top", fontsize=8,
                           bbox=dict(facecolor="white", alpha=.8, edgecolor="none", pad=1.))
            energy_line, = energy_ax.plot([], [], color=color, lw=1.2, label=order)
            energy_dot, = energy_ax.plot([], [], "o", color=color, ms=3)
            moving.append(dict(run=run, case=case, line=line, trail=trail, dot=dot,
                               text=text, energy_line=energy_line, energy_dot=energy_dot))
        energy_ax.legend(loc="best", fontsize=7, ncol=len(runs))
    row = next(iter(selections.values()))[0]["row"]
    title = fig.suptitle("")
    def update(frame):
        t = frame / fps
        title.set_text(f"Actual motion and physical energy · {row['hz']} Hz · {row['sweeps']} sweeps · "
                       f"ε / E-scale = {float(row['epsilon_ratio']):g} · t = {t:.2f} s")
        for a in moving:
            run = a["run"]
            data = run["accepted"]
            i = index_at(run, t)
            p = np.array([data["x"][i], data["y"][i]])
            xx, yy = spring_points(p) if a["case"].startswith("spring") else ([0, p[0]], [0, p[1]])
            a["line"].set_data(xx, yy)
            a["dot"].set_data([p[0]], [p[1]])
            begin = index_at(run, max(0., t - .6))
            a["trail"].set_data(data["x"][begin:i + 1], data["y"][begin:i + 1])
            a["energy_line"].set_data(data["time"][:i + 1], data["energy"][:i + 1])
            a["energy_dot"].set_data([data["time"][i]], [data["energy"][i]])
            warning, bad = status(run, t)
            max_error = 100 * np.max(np.abs(data['rel_energy']))
            message = (f"length {run['length'][i]:.4f} m\nspeed {run['speed'][i]:.3g} m/s"
                       f"\nE {data['energy'][i]:.6g} J\nΔE {100 * data['rel_energy'][i]:+.3g}%"
                       f"\nmax |ΔE| {max_error:.3g}%")
            if warning:
                message += "\n" + ("FAILED · prefix held" if warning.startswith("FAILED")
                                    else "FROZEN motion" if warning.startswith("FROZEN") else "Nearly stationary")
            a["text"].set_text(message)
            a["text"].set_color("#a40000" if bad else "#222222")
        return []
    update(round(fps * 1.))
    fig.savefig(out / f"{group}_poster.png", dpi=100)
    # Fix the already verified geometry. Re-running constrained layout while
    # energy tick strings change is unnecessary for every animation frame.
    fig.set_layout_engine("none")
    animation = FuncAnimation(fig, update, frames=5 * fps + 1, interval=1000 / fps, blit=False)
    animation.save(out / f"{group}_motion.mp4", writer=FFMpegWriter(fps=fps, bitrate=3500,
                    extra_args=["-pix_fmt", "yuv420p"]), dpi=100)
    plt.close(fig)


def release_detail(summary, out):
    fig, axes = plt.subplots(2, 2, figsize=(11, 6), sharex="col", layout="constrained")
    for col, order in enumerate(("EL", "LE")):
        run = select(summary, "spring_release", order, 120, 16, 1e4)
        data = run["accepted"]
        changes = np.abs(np.diff(data["vx"]))
        if order == "LE":
            # The near-identical maximum at startup would conceal the
            # recurring turning-point restart. Show its first later event.
            candidates = np.flatnonzero((data["time"][1:] > .1) & (changes > .4))
            event = int(candidates[0]) + 1
        else:
            event = int(np.argmax(changes)) + 1
        center = data["time"][event]
        mask = (data["time"] >= max(0., center - .07)) & (data["time"] <= center + .07)
        axes[0, col].plot(data["time"][mask], data["x"][mask] - 1., "o-",
                          ms=4, lw=1.1, color=COLORS[order])
        axes[1, col].plot(data["time"][mask], data["vx"][mask], "o-",
                          ms=4, lw=1.1, color=COLORS[order])
        for ax in axes[:, col]:
            ax.axvline(center, ls=":", lw=.8, color="black")
            ax.grid(alpha=.2)
            ax.ticklabel_format(axis="y", style="plain", useOffset=False)
        change = data["vx"][event] - data["vx"][event - 1]
        axes[0, col].set_title(f"{DISPLAY[order]} · turn restart Δv = {change:+.3f} m/s")
        axes[0, col].set_ylabel("Spring extension (m)")
        axes[1, col].set_ylabel("Derived velocity (m/s)")
        axes[1, col].set_xlabel("Time (s)")
    fig.suptitle("Released spring · 120 Hz · 16 sweeps · ε / E-scale = 10000\n"
                 "Every timestep shown near EL's strongest turn restart and LE's first later turn restart", fontsize=12)
    fig.savefig(out / "release_velocity_detail.png", dpi=150)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--group", choices=("main", "high_penalty", "all"), default="all")
    parser.add_argument("--no-video", action="store_true")
    parser.add_argument("--fps", type=int, default=30)
    args = parser.parse_args()
    summary = csv_rows(SOURCE / "summary.csv")
    release_detail(summary, HERE)
    groups = dict(main=(120, 16, 1e4, ("EL", "LE")),
                  high_penalty=(60, 16, 1e8, ("EL", "LE", "JOINT")))
    manifest = []
    for group, (hz, sweeps, ratio, orders) in groups.items():
        if args.group != "all" and args.group != group:
            continue
        selections = {case: [select(summary, case, order, hz, sweeps, ratio) for order in orders]
                      for case in CASES}
        for case, runs in selections.items():
            contact_sheet(group, case, runs, HERE)
            for run in runs:
                row = run["row"]
                manifest.append(dict(group=group, case=case, order=row["order"], run_id=row["run_id"],
                                     completed=row["completed"], accepted_time=row["accepted_time"],
                                     max_actual_relative_energy=row["max_rel_energy"],
                                     accepted_path_length=run["path"],
                                     max_accepted_speed=float(np.max(run["speed"]))))
        path_and_energy(group, selections, HERE)
        if not args.no_video:
            movie(group, selections, HERE, args.fps)
        print(f"Finished {group}", flush=True)
    with (HERE / f"{args.group}_selected_runs.csv").open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(manifest[0]))
        writer.writeheader()
        writer.writerows(manifest)


if __name__ == "__main__":
    main()
