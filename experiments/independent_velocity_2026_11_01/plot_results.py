"""Inspect actual independent-velocity energy and displacement kinematics.

No numerical simulation, velocity replacement, or continuous-reference score
is performed. A chord-energy discrepancy is a diagnostic, not a requirement
that independent velocity equal backward position difference.
"""
from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, FFMpegWriter


ROOT = Path(__file__).resolve().parent
BASELINE = Path("/home/blacknoir/Projects/constrained-xpbd/experiments/block_orders_2026_11_01/results")
CASES = ("spring_release", "spring_orbit", "spring_stiff", "spring_gravity", "pendulum_release", "pendulum_kick")
COLORS = {"EL": "#d55e00", "LE": "#0072b2", "JOINT": "#7b2cbf"}
LABELS = {"EL": "Energy → material", "LE": "Material → energy", "JOINT": "Joint solve"}


def read(path):
    with path.open(newline="") as f:
        return list(csv.DictReader(f))


def write(path, rows):
    if not rows:
        return
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def select(rows, case, order, args, beta=None):
    found = [r for r in rows if r["case"] == case and r["order"] == order
             and float(r["hz"]) == args.hz and float(r["sweeps"]) == args.sweeps
             and float(r["epsilon_ratio"]) == args.epsilon_ratio
             and float(r["beta_v"]) == (args.beta_v if beta is None else beta)]
    return found[0] if found else None


def trajectory(results, row):
    if row is None:
        return None
    records = read(results / "trajectories" / (row["run_id"] + ".csv"))
    data = {key: np.array([float(r[key]) for r in records]) for key in records[0]}
    accepted = data["accepted"] > .5
    data = {key: value[accepted] for key, value in data.items()}
    data["factor"] = 1 + float(row["beta_v"]) * data["eta"]
    return data


def selected(rows, args):
    return {case: [(order, select(rows, case, order, args)) for order in args.orders.split(",")
                   if select(rows, case, order, args) is not None] for case in CASES}


def traces(rows, out, args):
    fig, axes = plt.subplots(6, 3, figsize=(15, 17), sharex="col", layout="constrained")
    for i, (case, entries) in enumerate(selected(rows, args).items()):
        for order, row in entries:
            data = trajectory(args.results, row)
            color = COLORS[order]
            t = data["time"]
            axes[i, 0].plot(t, 100 * data["rel_energy"], color=color, lw=1.2, label=order + " actual-v E")
            axes[i, 0].plot(t, 100 * data["chord_rel_energy"], color=color, lw=.9, ls="--", label=order + " chord E")
            axes[i, 1].plot(t, data["velocity_chord_mismatch"], color=color, lw=1.1, label=order)
            axes[i, 2].plot(t, data["factor"], color=color, lw=1.1, label=order)
            if row["completed"] != "1":
                for ax in axes[i]:
                    ax.axvline(float(row["accepted_time"]), color=color, ls=":", lw=.8)
        axes[i, 0].axhline(0, color="black", lw=.7, ls=":")
        axes[i, 2].axhline(1, color="black", lw=.6, ls=":")
        axes[i, 2].axhline(0, color="#a40000", lw=.7, ls="--")
        axes[i, 0].set_title(case.replace("_", " ").capitalize())
        axes[i, 0].set_ylabel("Actual / chord energy error (%)")
        axes[i, 1].set_ylabel("|v − (x − xₙ)/h| (m/s)")
        axes[i, 2].set_ylabel("1 + βᵥη")
        axes[i, 0].legend(fontsize=6.5, ncol=2)
        for ax in axes[i]:
            ax.grid(alpha=.22)
            ax.ticklabel_format(axis="y", style="plain", useOffset=False)
            ax.set_xlabel("Time (s)")
    fig.suptitle(f"Independent velocity · βᵥ={args.beta_v:g} · {args.hz:g} Hz · {args.sweeps} sweeps · "
                 f"ε / E-scale={args.epsilon_ratio:g}\n"
                 "Solid: actual E(x,v). Dashed: energy using chord velocity. Both normalized by initial E − minimum U.", fontsize=12)
    fig.savefig(out / "energy_coupling.png", dpi=130)
    plt.close(fig)


def paths(rows, out, args):
    fig, axes = plt.subplots(6, 3, figsize=(15, 17), layout="constrained")
    for i, (case, entries) in enumerate(selected(rows, args).items()):
        for order, row in entries:
            data = trajectory(args.results, row)
            color = COLORS[order]
            axes[i, 0].plot(data["x"], data["y"], color=color, lw=1.1, alpha=.8, label=order)
            axes[i, 0].plot(data["x"][-1], data["y"][-1], "o", color=color, ms=4)
            axes[i, 1].plot(data["time"], data["radius"], color=color, lw=1.1, label=order)
            axes[i, 2].plot(data["time"], data["speed"], color=color, lw=1.1, label=order + " actual |v|")
            axes[i, 2].plot(data["time"], data["chord_speed"], color=color, ls="--", lw=.9, label=order + " chord speed")
        axes[i, 0].plot(0, 0, "ks", ms=4)
        axes[i, 0].set_aspect("equal")
        axes[i, 0].set_title(case.replace("_", " ").capitalize())
        axes[i, 0].set_xlabel("x (m)")
        axes[i, 0].set_ylabel("y (m)")
        axes[i, 1].set_ylabel("Length from fixed pivot (m)")
        axes[i, 2].set_ylabel("Actual / chord speed (m/s)")
        axes[i, 2].legend(fontsize=6.5, ncol=2)
        for ax in axes[i, 1:]:
            ax.set_xlabel("Time (s)")
            ax.grid(alpha=.22)
            ax.ticklabel_format(axis="y", style="plain", useOffset=False)
    fig.suptitle("Accepted position paths and kinematics · independent v is retained between steps\n"
                 f"βᵥ={args.beta_v:g} · {args.hz:g} Hz · {args.sweeps} sweeps · ε / E-scale={args.epsilon_ratio:g}", fontsize=12)
    fig.savefig(out / "paths_lengths_speeds.png", dpi=130)
    plt.close(fig)


def contacts(rows, out, args):
    times = (0., .25, .5, 1., 2., 3., 4., 5.)
    for case, entries in selected(rows, args).items():
        if not entries:
            continue
        sources = [(order, row, trajectory(args.results, row)) for order, row in entries]
        radius = max(1.1, min(3., max(max(d["radius"]) for _, _, d in sources) * 1.1))
        fig, axes = plt.subplots(len(sources), len(times), figsize=(18, len(sources) * 3),
                                  squeeze=False, layout="constrained")
        for i, (order, row, data) in enumerate(sources):
            color = COLORS[order]
            for j, t in enumerate(times):
                index = min(max(0, np.searchsorted(data["time"], t + 1e-10, side="right") - 1), len(data["time"]) - 1)
                ax = axes[i, j]
                ax.set_xlim(-radius, radius)
                ax.set_ylim(-radius, radius)
                ax.set_aspect("equal")
                ax.set_xticks([])
                ax.set_yticks([])
                ax.plot(0, 0, "ks", ms=4)
                p = np.array([data["x"][index], data["y"][index]])
                xx, yy = spring_points(p) if case.startswith("spring") else ([0, p[0]], [0, p[1]])
                ax.plot(data["x"][:index + 1], data["y"][:index + 1], lw=.7, alpha=.25, color=color)
                ax.plot(xx, yy, color=color, lw=1.4)
                ax.plot(*p, "o", color=color, ms=6)
                ax.set_title(f"t={t:g} s", fontsize=9)
                if j == 0:
                    ax.set_ylabel(LABELS[order], fontsize=9)
                message = (f"length {data['radius'][index]:.4f} m\n"
                           f"E(v) {data['energy'][index]:.5g} J\n"
                           f"ΔE v/chord {100 * data['rel_energy'][index]:+.3g}/{100 * data['chord_rel_energy'][index]:+.3g}%\n"
                           f"|v−chord| {data['velocity_chord_mismatch'][index]:.3g}\n1+βη {data['factor'][index]:.4f}")
                bad = False
                if row["completed"] != "1" and t > float(row["accepted_time"]):
                    message += "\nFAILED · prefix held"
                    bad = True
                elif float(row["path_length"]) < 1e-4:
                    message += "\nFROZEN motion"
                    bad = True
                elif float(row["tail_path_length"]) < 1e-4 and t >= 4.:
                    message += "\nNearly stationary tail"
                    bad = True
                ax.text(.02, .98, message, transform=ax.transAxes, va="top", fontsize=7,
                        color="#a40000" if bad else "#222222",
                        bbox=dict(facecolor="white", alpha=.8, edgecolor="none", pad=1.))
        fig.suptitle(f"{case.replace('_', ' ').capitalize()} · independent velocity · βᵥ={args.beta_v:g} · "
                     f"{args.hz:g} Hz · {args.sweeps} sweeps · ε ratio={args.epsilon_ratio:g}\n"
                     "Actual-v energy and chord-energy diagnostic are reported separately; black square is fixed pivot", fontsize=12)
        fig.savefig(out / f"{case}_contacts.png", dpi=130)
        plt.close(fig)


def setting_table(rows, out):
    fields = ("order", "beta_v", "hz", "sweeps", "epsilon_ratio")
    groups = {}
    for row in rows:
        groups.setdefault(tuple(row[k] for k in fields), []).append(row)
    result = []
    for key, group in groups.items():
        result.append(dict(zip(fields, key)) | dict(
            cases=len(group), completed=sum(int(r["completed"]) for r in group),
            tail_stalled=sum(int(r["tail_stalled"]) for r in group),
            max_actual_energy_error=max(float(r["max_accepted_rel_energy"]) for r in group),
            max_chord_energy_error=max(float(r["max_accepted_chord_rel_energy"]) for r in group),
            max_velocity_chord_mismatch=max(float(r["max_velocity_chord_mismatch"]) for r in group),
            max_geometry_error=max(float(r["max_geometry_error"]) for r in group),
            max_tangency_error=max(float(r["max_tangency_error"]) for r in group),
            max_accepted_displacement=max(float(r["max_accepted_displacement"]) for r in group),
            min_tail_path=min(float(r["tail_path_length"]) for r in group),
            max_angular_error_central=max([float(r["max_angular_error"]) for r in group
                                         if r["angular_conservation_expected"] == "1"] or [0.])))
    result.sort(key=lambda r: (-r["completed"], r["tail_stalled"], r["max_actual_energy_error"]))
    write(out / "settings.csv", result)


def measured_table(rows, out, args):
    baseline = read(BASELINE / "summary.csv")
    lines = ["# Independent velocity: selected measured runs", "",
             f"βᵥ={args.beta_v:g}, {args.hz:g} Hz, {args.sweeps} sweeps, epsilon ratio={args.epsilon_ratio:g}.", "",
             "Energy percentages normalize by initial mechanical energy above minimum potential. Actual E uses returned independent velocity; chord E substitutes the backward position difference. Their difference is diagnostic, not an imposed equality. Incomplete runs are prefixes.", "",
             "| Case | Method | Complete? | Max actual E error | Max chord E error | Max v/chord mismatch | Tail travel | Min 1+βη |", "|---|---|---|---:|---:|---:|---:|---:|"]
    manifest = []
    for case, entries in selected(rows, args).items():
        for order, row in entries:
            data = trajectory(args.results, row)
            lines.append(f"| {case} | independent {order} | {row['completed']} | "
                         f"{100 * float(row['max_accepted_rel_energy']):.5g}% | "
                         f"{100 * float(row['max_accepted_chord_rel_energy']):.5g}% | "
                         f"{float(row['max_velocity_chord_mismatch']):.4g} m/s | "
                         f"{float(row['tail_path_length']):.4g} m | {min(data['factor']):.5g} |")
            manifest.append(dict(case=case, order=order, run_id=row["run_id"]))
        for order in ("EL", "LE"):
            old = [r for r in baseline if r["case"] == case and r["order"] == order
                   and float(r["hz"]) == args.hz and float(r["sweeps"]) == args.sweeps
                   and float(r["epsilon_ratio"]) == args.epsilon_ratio]
            if not old:
                continue
            row = old[0]
            old_records = read(BASELINE / "trajectories" / (row["run_id"] + ".csv"))
            accepted = [r for r in old_records if r["accepted"] == "1"]
            actual_error = max(abs(float(r["rel_energy"])) for r in accepted)
            tail = [r for r in accepted if float(r["time"]) > max(0., float(row["accepted_time"]) - 1.)]
            positions = np.array([[float(r["x"]), float(r["y"])] for r in tail])
            travel = np.linalg.norm(np.diff(positions, axis=0), axis=1).sum() if len(positions) > 1 else 0.
            lines.append(f"| {case} | derived-v {order} | {row['completed']} | {100 * actual_error:.5g}% | "
                         f"{100 * actual_error:.5g}% | 0 m/s | {travel:.4g} m | — |")
    lines += ["", "Derived-velocity baselines are reused archived trajectories; they were not rerun.", ""]
    (out / "SELECTED_RESULTS.md").write_text("\n".join(lines))
    write(out / "selected_runs.csv", manifest)


def spring_points(p):
    radius = np.linalg.norm(p)
    if radius <= 1e-12:
        return [0., p[0]], [0., p[1]]
    along = np.r_[0., .08, np.linspace(.1, .9, 19), .92, 1.]
    side = np.zeros(len(along))
    side[2:-2] = .035 * np.where(np.arange(19) % 2, -1., 1.)
    points = along[:, None] * p + side[:, None] * np.array([-p[1], p[0]]) / radius
    return points[:, 0], points[:, 1]


def movie(rows, out, args):
    selections = selected(rows, args)
    method_count = len(args.orders.split(","))
    columns = 2 if method_count >= 3 else 3
    nrows = 3 if columns == 2 else 2
    fig = plt.figure(figsize=(18, 15 if nrows == 3 else 12), layout="constrained")
    outer = fig.add_gridspec(nrows, columns)
    moving = []
    for case_index, (case, entries) in enumerate(selections.items()):
        inner = outer[case_index // columns, case_index % columns].subgridspec(2, max(1, len(entries)), height_ratios=(3., 1.))
        energy_ax = fig.add_subplot(inner[1, :])
        energy_ax.axhline(0., color="black", ls=":", lw=.7)
        energy_ax.set_xlim(0, 5.)
        energy_ax.set_ylabel("Energy error (%)", fontsize=8)
        energy_ax.set_xlabel("Time (s)", fontsize=8)
        energy_ax.grid(alpha=.2)
        energy_ax.tick_params(labelsize=7)
        energy_ax.ticklabel_format(axis="y", style="plain", useOffset=False)
        datasets = [(order, row, trajectory(args.results, row)) for order, row in entries]
        radius = max([max(data["radius"]) for _, _, data in datasets] or [1.]) * 1.1
        radius = max(1.1, min(radius, 3.))
        all_errors = np.concatenate([np.r_[100 * data["rel_energy"], 100 * data["chord_rel_energy"]]
                                     for _, _, data in datasets])
        spread = max(np.ptp(all_errors), 1e-4)
        energy_ax.set_ylim(min(all_errors) - .1 * spread, max(all_errors) + .1 * spread)
        for j, (order, row, data) in enumerate(datasets):
            color = COLORS[order]
            ax = fig.add_subplot(inner[0, j])
            ax.set_aspect("equal")
            ax.set_xlim(-radius, radius)
            ax.set_ylim(-radius, radius)
            ax.set_xticks([])
            ax.set_yticks([])
            ax.plot(0, 0, "ks", ms=4)
            ax.axhline(0, color=".85", lw=.6)
            ax.axvline(0, color=".85", lw=.6)
            ax.set_title(case.replace("_", " ") + "\n" + LABELS[order], fontsize=9)
            line, = ax.plot([], [], color=color, lw=1.4)
            dot, = ax.plot([], [], "o", color=color, ms=7)
            trail, = ax.plot([], [], color=color, lw=.8, alpha=.3)
            text = ax.text(.02, .98, "", transform=ax.transAxes, va="top", fontsize=7.5,
                           bbox=dict(facecolor="white", alpha=.8, edgecolor="none", pad=1.))
            actual_line, = energy_ax.plot([], [], color=color, lw=1.2, label=order + " actual")
            chord_line, = energy_ax.plot([], [], color=color, lw=.8, ls="--", label=order + " chord")
            moving.append(dict(case=case, row=row, data=data, line=line, dot=dot, trail=trail,
                               text=text, actual_line=actual_line, chord_line=chord_line))
        energy_ax.legend(fontsize=6, ncol=2)
    title = fig.suptitle("")
    fps = args.fps
    def update(frame):
        t = frame / fps
        title.set_text(f"Independent velocity retained · βᵥ={args.beta_v:g} · {args.hz:g} Hz · "
                       f"{args.sweeps} sweeps · ε ratio={args.epsilon_ratio:g} · t={t:.2f} s\n"
                       "Solid energy: actual independent v. Dashed energy: chord velocity diagnostic.")
        for a in moving:
            data, row = a["data"], a["row"]
            i = min(max(0, np.searchsorted(data["time"], t + 1e-9, side="right") - 1), len(data["time"]) - 1)
            p = np.array([data["x"][i], data["y"][i]])
            xx, yy = spring_points(p) if a["case"].startswith("spring") else ([0, p[0]], [0, p[1]])
            a["line"].set_data(xx, yy)
            a["dot"].set_data([p[0]], [p[1]])
            start = max(0, np.searchsorted(data["time"], t - .6))
            a["trail"].set_data(data["x"][start:i + 1], data["y"][start:i + 1])
            a["actual_line"].set_data(data["time"][:i + 1], 100 * data["rel_energy"][:i + 1])
            a["chord_line"].set_data(data["time"][:i + 1], 100 * data["chord_rel_energy"][:i + 1])
            message = (f"length {data['radius'][i]:.4f} m\n"
                       f"speed v/chord {data['speed'][i]:.3g}/{data['chord_speed'][i]:.3g}\n"
                       f"actual E {data['energy'][i]:.6g} J\n"
                       f"ΔE actual {100 * data['rel_energy'][i]:+.3g}%\n"
                       f"ΔE chord {100 * data['chord_rel_energy'][i]:+.3g}%\n"
                       f"|v−chord| {data['velocity_chord_mismatch'][i]:.3g} m/s\n"
                       f"1+βη {data['factor'][i]:.4f}")
            failed = row["completed"] != "1" and t > float(row["accepted_time"])
            frozen = float(row["tail_path_length"]) < 1e-4 and t >= 4.
            if failed:
                message += "\nFAILED · prefix held"
            elif frozen:
                message += "\nNearly stationary tail"
            a["text"].set_text(message)
            a["text"].set_color("#a40000" if failed or frozen else "#222222")
        return []
    update(fps)
    fig.savefig(out / "playback_poster.png", dpi=100)
    fig.set_layout_engine("none")
    animation = FuncAnimation(fig, update, frames=5 * fps + 1, interval=1000 / fps, blit=False)
    animation.save(out / "comparison.mp4", writer=FFMpegWriter(fps=fps, bitrate=3500,
                   extra_args=["-pix_fmt", "yuv420p"]), dpi=100)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", type=Path, default=ROOT / "results")
    parser.add_argument("--beta-v", type=float, default=1.)
    parser.add_argument("--hz", type=int, default=60)
    parser.add_argument("--sweeps", type=int, default=4)
    parser.add_argument("--epsilon-ratio", type=float, default=1e8)
    parser.add_argument("--orders", default="EL,LE,JOINT")
    parser.add_argument("--fps", type=int, default=30)
    parser.add_argument("--no-video", action="store_true")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    out = args.output or args.results / "figures"
    out.mkdir(parents=True, exist_ok=True)
    rows = read(args.results / "summary.csv")
    setting_table(rows, out)
    measured_table(rows, out, args)
    traces(rows, out, args)
    paths(rows, out, args)
    contacts(rows, out, args)
    if not args.no_video:
        movie(rows, out, args)
    print(f"Wrote {len(rows)} measured run summaries and selected visualizations to {out}")


if __name__ == "__main__":
    main()
