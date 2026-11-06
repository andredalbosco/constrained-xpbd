"""Inspect EL/LE block solves without hiding incomplete trajectories.

Run after run_experiments.py (or pass --results to another output directory).
This file reads measured CSV only; it does not rerun or alter a simulation.
"""
from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.animation import FFMpegWriter, FuncAnimation


ROOT = Path(__file__).resolve().parent
CASES = ("spring_release", "spring_orbit", "pendulum_release")
COLORS = {"EL": "#d55e00", "LE": "#0072b2", "JOINT": "#7b2cbf", "XPBD": "#777777"}
LABELS = {"EL": "Energy → material", "LE": "Material → energy", "JOINT": "Joint solve", "XPBD": "XPBD control"}


def read_csv(path):
    with Path(path).open(newline="") as f:
        return list(csv.DictReader(f))


def number(row, name, default=np.nan):
    try:
        return float(row.get(name, default))
    except (ValueError, TypeError):
        return float(default)


def truth(row, name):
    return str(row.get(name, "")).lower() in ("1", "true", "yes")


def select(summary, case, order, *, hz=60, sweeps=4, ratio=100):
    rows = [r for r in summary if r["case"] == case and r["order"] == order]
    if not rows:
        return None
    # Compare the specified experiment first. Controls may not carry the same
    # epsilon/sweep settings, so they use their nearest available metadata.
    def score(r):
        return (abs(number(r, "hz") - hz),
                abs(np.log10(max(number(r, "epsilon_ratio", ratio), 1e-300) / ratio)),
                abs(number(r, "sweeps") - sweeps), r["run_id"])
    return min(rows, key=score)


def trajectory(results, row):
    if row is None:
        return None
    path = results / "trajectories" / (row["run_id"] + ".csv")
    if not path.exists():
        return None
    records = read_csv(path)
    if not records:
        return None
    return {key: np.array([number(r, key) for r in records]) for key in records[0]}


def position_coordinate(case, x, y):
    if case == "spring_release":
        return x - 1.0
    if case == "spring_orbit":
        return np.unwrap(np.arctan2(y, x))
    return np.unwrap(np.arctan2(x, -y))


def failure_marker(ax, row, data):
    if truth(row, "completed"):
        return
    accepted_time = number(row, "accepted_time", data["time"][-1])
    ax.axvline(accepted_time, color="#a40000", lw=0.9, ls=":", alpha=.85)


def traces(results, summary, out, args):
    fig, axes = plt.subplots(3, 3, figsize=(14, 9), sharex="col", layout="constrained")
    metadata = []
    for col, case in enumerate(CASES):
        reference_drawn = False
        failures = []
        for order in ("XPBD", "JOINT", "EL", "LE"):
            row = select(summary, case, order, hz=args.hz, sweeps=args.sweeps, ratio=args.epsilon_ratio)
            data = trajectory(results, row)
            if data is None:
                continue
            t = data["time"]
            if not reference_drawn:
                axes[0, col].plot(t, position_coordinate(case, data["ref_x"], data["ref_y"]),
                                  "k--", lw=1.6, label="Continuous reference")
                reference_drawn = True
            label = LABELS[order]
            if not truth(row, "completed"):
                label += " (failed prefix)"
                failures.append(f"{order}: {number(row, 'accepted_time'):.3g} s")
            axes[0, col].plot(t, position_coordinate(case, data["x"], data["y"]),
                              color=COLORS[order], lw=1.3, label=label)
            axes[1, col].semilogy(t, np.maximum(np.abs(data["rel_energy"]), 1e-16),
                                  color=COLORS[order], lw=1.25)
            axes[2, col].semilogy(t, np.maximum(np.abs(data["residual"]), 1e-16),
                                  color=COLORS[order], lw=1.25)
            for ax in axes[:, col]:
                failure_marker(ax, row, data)
            metadata.append({"case": case, "order": order, "run_id": row["run_id"]})
        axes[0, col].set_title(case.replace("_", " ").capitalize())
        axes[0, col].set_ylabel("Extension (m)" if col == 0 else "Angle (rad)")
        axes[1, col].set_ylabel("Absolute normalized energy error")
        axes[2, col].set_ylabel("Scaled nonlinear residual")
        axes[2, col].set_xlabel("Time (s)")
        if failures:
            axes[0, col].text(.02, .02, "Failed: " + "; ".join(failures),
                               transform=axes[0, col].transAxes, fontsize=8, color="#a40000")
        for ax in axes[:, col]:
            ax.grid(alpha=.22)
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="outside lower center", ncol=3, fontsize=9)
    fig.suptitle(f"Block Gauss–Seidel: {args.hz:g} Hz, {args.sweeps:g} sweeps, "
                 f"ε / E-scale = {args.epsilon_ratio:g}\n"
                 "Dashed vertical marks retain the end of failed prefixes", fontsize=13)
    fig.savefig(out / "motion_energy_residual.png", dpi=160)
    plt.close(fig)
    return metadata


def heatmaps(summary, out, args):
    selected = [r for r in summary if r["order"] in ("EL", "LE")
                and abs(number(r, "hz") - args.hz) < 1e-8 and r["case"] in CASES]
    if not selected:
        return
    ratios = sorted({number(r, "epsilon_ratio") for r in selected if np.isfinite(number(r, "epsilon_ratio"))})
    sweeps = sorted({int(number(r, "sweeps")) for r in selected})
    if not ratios or not sweeps:
        return
    fig, axes = plt.subplots(3, 2, figsize=(12, 10), layout="constrained")
    images = []
    for i, case in enumerate(CASES):
        for j, order in enumerate(("EL", "LE")):
            ax = axes[i, j]
            values = np.full((len(sweeps), len(ratios)), np.nan)
            completion = np.full_like(values, np.nan)
            for r in selected:
                if r["case"] != case or r["order"] != order:
                    continue
                si = sweeps.index(int(number(r, "sweeps")))
                ei = ratios.index(number(r, "epsilon_ratio"))
                values[si, ei] = np.log10(max(number(r, "max_residual"), 1e-16))
                completion[si, ei] = truth(r, "completed")
            shown = np.ma.masked_invalid(values)
            im = ax.imshow(shown, origin="lower", aspect="auto", cmap="viridis", vmin=-10, vmax=4)
            images.append(im)
            for si in range(len(sweeps)):
                for ei in range(len(ratios)):
                    if not np.isfinite(values[si, ei]):
                        continue
                    label = f"{values[si, ei]:.1f}" + (" ×" if completion[si, ei] == 0 else "")
                    ax.text(ei, si, label, ha="center", va="center", fontsize=8,
                            color="black" if values[si, ei] > -3 else "white")
            ax.set_xticks(range(len(ratios)), [f"{v:g}" for v in ratios])
            ax.set_yticks(range(len(sweeps)), [str(v) for v in sweeps])
            ax.set_xlabel("ε / E-scale")
            ax.set_ylabel("Sweeps per timestep")
            ax.set_title(case.replace("_", " ").capitalize() + " · " + order)
    fig.colorbar(images[0], ax=axes.ravel().tolist(), label="log₁₀ maximum scaled nonlinear residual", shrink=.85)
    fig.suptitle(f"{args.hz:g} Hz · cell text gives log₁₀ residual; × marks an incomplete run", fontsize=13)
    fig.savefig(out / "residual_heatmap.png", dpi=160)
    plt.close(fig)


def high_penalty_release(results, summary, out):
    """Show branch failure and a converged frozen solution in one small panel."""
    case = "spring_release"
    chosen = {order: select(summary, case, order, hz=60, sweeps=16, ratio=1e8)
              for order in ("EL", "LE", "JOINT")}
    if not all(r is not None and abs(number(r, "epsilon_ratio") - 1e8) < 1 for r in chosen.values()):
        return []
    fig, axes = plt.subplots(1, 3, figsize=(14, 3.8), layout="constrained")
    source = trajectory(results, chosen["JOINT"])
    axes[0].plot(source["time"], source["ref_x"] - 1, "k--", lw=1.4, label="Continuous reference")
    failures = []
    manifest = []
    for order, row in chosen.items():
        data = trajectory(results, row)
        if data is None:
            continue
        accepted = data.get("accepted", np.ones_like(data["time"])) > .5
        label = LABELS[order]
        if not truth(row, "completed"):
            label += f" · fails at {number(row, 'accepted_time'):.3g} s"
            rejected = np.flatnonzero(~accepted)
            if len(rejected):
                index = rejected[-1]
                failures.append(f"{order}: rejected |x|={abs(data['x'][index]):.2g} m")
        axes[0].plot(data["time"][accepted], data["x"][accepted] - 1,
                     color=COLORS[order], lw=1.5, label=label)
        axes[1].semilogy(data["time"], np.maximum(np.abs(data["rel_energy"]), 1e-16),
                         color=COLORS[order], lw=1.3)
        axes[2].semilogy(data["time"], np.maximum(data["residual"], 1e-16),
                         color=COLORS[order], lw=1.3)
        if not truth(row, "completed"):
            for ax in axes:
                failure_marker(ax, row, data)
            for ax, values in zip(axes[1:], (np.abs(data["rel_energy"]), data["residual"])):
                ax.scatter(data["time"][~accepted], np.maximum(values[~accepted], 1e-16),
                           color=COLORS[order], marker="x", s=40, zorder=5)
        manifest.append({"case": case, "order": order, "run_id": row["run_id"]})
    if failures:
        axes[0].text(.02, .02, "\n".join(failures), transform=axes[0].transAxes,
                     fontsize=8, va="bottom", color="#a40000")
    axes[0].set_ylabel("Accepted extension (m)")
    axes[1].set_ylabel("Absolute normalized energy error")
    axes[2].set_ylabel("Scaled nonlinear residual")
    for ax in axes:
        ax.set_xlabel("Time (s)")
        ax.grid(alpha=.22)
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="outside lower center", ncol=2, fontsize=8)
    fig.suptitle("High energy penalty · released spring · 60 Hz · 16 sweeps · ε / E-scale = 10⁸\n"
                 "Solid motion is accepted prefix; × marks rejected trial in error panels", fontsize=12)
    fig.savefig(out / "high_penalty_release.png", dpi=160)
    plt.close(fig)
    return manifest


def spring_points(p, turns=9):
    length = np.linalg.norm(p)
    if not np.isfinite(length) or length == 0:
        return np.array([0., 0.]), np.array([0., 0.])
    direction = p / length
    perpendicular = np.array([-direction[1], direction[0]])
    along = np.concatenate(([0., .08], np.linspace(.1, .9, 2 * turns + 1), [.92, 1.]))
    side = np.zeros_like(along)
    side[2:-2] = .035 * np.where(np.arange(2 * turns + 1) % 2 == 0, 1., -1.)
    points = along[:, None] * p + side[:, None] * perpendicular
    return points[:, 0], points[:, 1]


def movie(results, summary, out, args):
    methods = ("REF", "EL", "LE", "JOINT", "XPBD")
    fig, axes = plt.subplots(3, 5, figsize=(16, 9), layout="constrained")
    artists = []
    max_time = 0.
    manifest = []
    for row_index, case in enumerate(CASES):
        candidates = [select(summary, case, order, hz=args.hz, sweeps=args.sweeps,
                             ratio=args.epsilon_ratio) for order in methods[1:]]
        available = [(r, trajectory(results, r)) for r in candidates]
        available = [(r, d) for r, d in available if d is not None]
        source = max(available, key=lambda item: item[1]["time"][-1])[0] if available else None
        reference = trajectory(results, source)
        if reference is None:
            continue
        # Use a common visible physical region in each row. Out-of-view states
        # are explicitly labeled during playback, not silently clipped.
        ref_pos = np.column_stack((reference["ref_x"], reference["ref_y"]))
        finite_ref = ref_pos[np.isfinite(ref_pos).all(axis=1)]
        radius = max(1.15, np.max(np.linalg.norm(finite_ref, axis=1)) * 1.15)
        max_time = max(max_time, reference["time"][-1])
        for col, method in enumerate(methods):
            ax = axes[row_index, col]
            chosen = source if method == "REF" else candidates[col - 1]
            data = reference if method == "REF" else trajectory(results, chosen)
            ax.set_aspect("equal")
            ax.set_xlim(-radius, radius)
            ax.set_ylim(-radius, radius)
            ax.set_xticks([])
            ax.set_yticks([])
            if row_index == 0:
                ax.set_title("Continuous reference" if method == "REF" else LABELS[method], fontsize=10)
            if col == 0:
                ax.set_ylabel(case.replace("_", " ").capitalize(), fontsize=10)
            ax.plot(0., 0., "ks", ms=4)
            if data is None:
                ax.text(.5, .5, "No trajectory", transform=ax.transAxes, ha="center")
                continue
            px = data["ref_x"] if method == "REF" else data["x"]
            py = data["ref_y"] if method == "REF" else data["y"]
            color = "#222222" if method == "REF" else COLORS[method]
            trail, = ax.plot([], [], color=color, lw=.8, alpha=.35)
            line, = ax.plot([], [], color=color, lw=1.2)
            dot, = ax.plot([], [], "o", color=color, ms=7)
            text = ax.text(.03, .97, "", transform=ax.transAxes, va="top", fontsize=8)
            accepted_time = number(chosen, "accepted_time", data["time"][-1])
            completed = method == "REF" or truth(chosen, "completed")
            if method != "REF":
                max_time = max(max_time, data["time"][-1])
                manifest.append({"case": case, "order": method, "run_id": chosen["run_id"]})
            artists.append(dict(line=line, dot=dot, trail=trail, text=text, data=data,
                                px=px, py=py, case=case, order=method, radius=radius,
                                completed=completed, accepted_time=accepted_time,
                                reason=chosen.get("reason", "")))
    title = fig.suptitle("")
    def update(frame):
        t = frame / args.fps
        title.set_text(f"Block Gauss–Seidel · {args.hz:g} Hz · {args.sweeps:g} sweeps · "
                       f"ε / E-scale = {args.epsilon_ratio:g} · t = {t:.2f} s")
        for a in artists:
            data = a["data"]
            index = int(np.searchsorted(data["time"], t, side="right") - 1)
            index = max(0, min(index, len(data["time"]) - 1))
            # On failed runs, retain the last accepted prefix endpoint.
            failed = not a["completed"] and t > a["accepted_time"] + 1e-8
            if failed:
                index = int(np.searchsorted(data["time"], a["accepted_time"] + 1e-8, side="right") - 1)
                index = max(0, min(index, len(data["time"]) - 1))
            p = np.array([a["px"][index], a["py"][index]])
            finite = np.isfinite(p).all()
            if finite:
                if "spring" in a["case"]:
                    xx, yy = spring_points(p)
                else:
                    xx, yy = [0, p[0]], [0, p[1]]
                a["line"].set_data(xx, yy)
                a["dot"].set_data([p[0]], [p[1]])
                begin = max(0, index - int(1.0 * args.hz))
                a["trail"].set_data(a["px"][begin:index + 1], a["py"][begin:index + 1])
            else:
                a["line"].set_data([], [])
                a["dot"].set_data([], [])
            if a["order"] == "REF":
                message = "Continuous ODE"
            else:
                message = f"E error: {data['rel_energy'][index]:+.2e}"
                message += f"\nρ: {data['rho'][index]:+.2e}"
            if failed:
                message += f"\nFAILED at {a['accepted_time']:.3g} s\nAccepted prefix held"
            if not finite or np.any(np.abs(p) > a["radius"]):
                message += f"\nOutside view: |x|={np.linalg.norm(p):.3g}"
            a["text"].set_text(message)
            a["text"].set_color("#a40000" if failed else "#222222")
        return [value for a in artists for key, value in a.items() if key in ("line", "dot", "trail", "text")]
    update(round(min(1., max_time) * args.fps))
    fig.savefig(out / "playback_poster.png", dpi=120)
    animation = FuncAnimation(fig, update, frames=round(max_time * args.fps) + 1,
                              interval=1000 / args.fps, blit=False)
    animation.save(out.parent / "comparison.mp4",
                   writer=FFMpegWriter(fps=args.fps, bitrate=2400,
                                      extra_args=["-pix_fmt", "yuv420p"]), dpi=110)
    plt.close(fig)
    return manifest


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", type=Path, default=ROOT / "results")
    parser.add_argument("--hz", type=float, default=60)
    parser.add_argument("--sweeps", type=int, default=4)
    parser.add_argument("--epsilon-ratio", type=float, default=100)
    parser.add_argument("--fps", type=int, default=30)
    parser.add_argument("--no-video", action="store_true")
    args = parser.parse_args()
    summary = read_csv(args.results / "summary.csv")
    out = args.results / "figures"
    out.mkdir(exist_ok=True, parents=True)
    manifest = traces(args.results, summary, out, args)
    heatmaps(summary, out, args)
    manifest += high_penalty_release(args.results, summary, out)
    if not args.no_video:
        manifest += movie(args.results, summary, out, args)
    with (out / "selected_runs.csv").open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=("case", "order", "run_id"))
        writer.writeheader()
        writer.writerows(manifest)
    print(f"Wrote measured figures to {out}")


if __name__ == "__main__":
    main()
