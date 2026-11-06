#!/usr/bin/env python3
"""Exhaustive scorecards for the original endpoint-energy gradient-policy study.

Read-only input: one summary CSV. Each of its 1,350 runs appears on exactly
one scorecard; incomplete trajectories have no displayed numerical score.
"""
import argparse
import collections
import csv
import hashlib
import json
import math
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/xpbd-scorecard-matplotlib")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.patches import Rectangle


CASES = ["spring_rest", "spring_radial", "spring_tangent",
         "pendulum_rest", "pendulum_tangent"]
CASE_TITLES = ["Spring\nreleased from rest", "Spring\nradial initial velocity",
               "Spring\ntangential initial velocity", "Pendulum\nreleased from rest",
               "Pendulum\ntangential initial velocity"]
POLICIES = ["fresh", "lag_sweep", "frozen_step", "previous_step"]
GAMMAS = {"fresh": [-1., .25, 1., 4., 16.],
          "lag_sweep": [-1., 1.], "frozen_step": [-1., 1.],
          "previous_step": [-1., 1.]}
SETTINGS = [(60, 1), (60, 5), (60, 20), (240, 1), (240, 5), (240, 20)]
KEY_FIELDS = ["case", "order", "gradient", "gamma", "retain_rx", "h", "sweeps"]
TEXT_FIELDS = {"case", "order", "gradient", "status"}
PALETTE = ["#f2f7f0", "#d7e7cd", "#f5edbf", "#f1d1a0", "#eaa998", "#d87a84"]
ENERGY_BOUNDS = [.1, 1., 5., 20., 100.]
POSITION_BOUNDS = [.01, .05, .1, .25, 1.]
STOP_COLOR = "#5e6670"


def configurations():
    result = [("baseline", "none", -1., 0)]
    for order in ["EL", "LE"]:
        for retain in [0, 1]:
            for policy in POLICIES:
                for gamma in GAMMAS[policy]:
                    result.append((order, policy, gamma, retain))
    assert len(result) == 45
    return result


def source_key(row):
    return tuple(row[f] for f in KEY_FIELDS)


def canonical_key(row):
    hz = round(1 / row["h"])
    if hz not in [60, 240] or not math.isclose(row["h"], 1/hz, rel_tol=1e-13):
        raise ValueError(f"Unexpected timestep: {row['h']}")
    return (row["case"], row["order"], row["gradient"], row["gamma"],
            int(row["retain_rx"]), hz, int(row["sweeps"]))


def read_summary(path):
    with path.open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    for row in rows:
        for field in row:
            if field not in TEXT_FIELDS:
                row[field] = float(row[field])
                if not math.isfinite(row[field]):
                    raise ValueError(f"Nonfinite {field}: {row}")
        if row["Escale"] <= 0 or row["max_abs_energy_error"] < 0 or row["rms_position_error"] < 0:
            raise ValueError("Invalid metric or normalization")
        if not math.isclose(row["requested_steps"] * row["h"], 10., abs_tol=1e-12):
            raise ValueError("Expected ten-second runs")
    index = {canonical_key(r): r for r in rows}
    if len(index) != len(rows):
        raise ValueError("Duplicate run key")
    expected = {(case, *cfg, hz, sweeps) for hz, sweeps in SETTINGS
                for cfg in configurations() for case in CASES}
    if set(index) != expected:
        raise ValueError(f"Matrix coverage mismatch: missing {expected-set(index)}, extra {set(index)-expected}")
    if len(rows) != 1350:
        raise ValueError("Expected 1,350 rows")
    return rows, index


def metric_color(value, bounds):
    return PALETTE[sum(value > b for b in bounds)]


def formatted(value):
    # Three significant figures; scientific notation keeps extreme runs readable.
    if value == 0:
        return "0"
    if value < .001 or value >= 10000:
        return f"{value:.2e}".replace("e+0", "e+").replace("e-0", "e−").replace("e-", "e−")
    return f"{value:.3g}"


def make_card(index, hz, sweeps):
    configs = configurations()
    fig = plt.figure(figsize=(20, 17.3), facecolor="white")
    fig.text(.026, .969, f"All gradient-policy results · {hz} Hz · {sweeps} outer sweep{'s' if sweeps != 1 else ''} per timestep",
             fontsize=22, weight="bold", color="#172635", va="top")
    fig.text(.026, .94, "Original endpoint-energy formulation  |  10 seconds  |  45 configurations × 5 initial conditions",
             fontsize=12, color="#465461", va="top")
    fig.text(.026, .919, "Smaller errors are better. Completing a run, or keeping energy constant, does not establish correct motion.",
             fontsize=12, weight="bold", color="#5c3c31", va="top")

    # Two synchronized axes: full settings at left and ten metric columns at right.
    bottom, height = .168, .680
    labels = fig.add_axes([.026, bottom, .307, height])
    ax = fig.add_axes([.337, bottom, .646, height])
    for axis in [labels, ax]:
        axis.set_ylim(44.5, -.5)
        axis.set_yticks([])
        for spine in axis.spines.values():
            spine.set_visible(False)
    labels.set_xlim(0, 1)
    labels.set_xticks([])
    ax.set_xlim(-.5, 9.5)
    ax.set_xticks([])
    colx = [.01, .145, .435, .755, .88]
    heads = ["#", "Order", "Gradient policy", "γ cutoff", "rₓ"]
    for x, text in zip(colx, heads):
        labels.text(x, -.95, text, ha="left", va="bottom", fontsize=10, weight="bold", clip_on=False)
    for c, title in enumerate(CASE_TITLES):
        ax.text(2*c+.5, -2.8, title, ha="center", va="bottom", fontsize=11.2,
                weight="bold", color="#172635", clip_on=False, linespacing=1.13)
        for j, heading in enumerate(["Max E error\n(%)", "Position RMS\n(m)"]):
            ax.text(2*c+j, -.97, heading, ha="center", va="bottom", fontsize=9.9,
                    weight="bold", clip_on=False)

    drawn_keys = []
    stopped = incompatible = 0
    for i, cfg in enumerate(configs):
        order, policy, gamma, retain = cfg
        if i % 2 == 0:
            labels.add_patch(Rectangle((0, i-.5), 1, 1, facecolor="#f3f5f7", edgecolor="none"))
        content = [str(i+1), "XPBD" if order == "baseline" else order, policy,
                   "off" if gamma < 0 else f"{gamma:g}", "drop" if not retain else "keep"]
        if order == "baseline":
            content[2:] = ["—", "—", "—"]
        for x, text in zip(colx, content):
            labels.text(x, i, text, ha="left", va="center", fontsize=10.3,
                        weight="bold" if order == "baseline" else "normal", color="#233140")
        for c, case in enumerate(CASES):
            key = (case, *cfg, hz, sweeps)
            row = index[key]
            drawn_keys.append(key)
            if row["completed_steps"] != row["requested_steps"]:
                stopped += 1
                ax.add_patch(Rectangle((2*c-.5, i-.5), 2, 1, facecolor=STOP_COLOR,
                                       edgecolor="white", linewidth=.4))
                ax.text(2*c+.5, i, "STOP", ha="center", va="center", fontsize=9.7,
                        color="white", weight="bold")
                continue
            bad_update = row["incompatible_updates"] > 0
            incompatible += int(bad_update)
            metrics = [100 * row["max_abs_energy_error"] / row["Escale"], row["rms_position_error"]]
            for j, (value, bounds) in enumerate(zip(metrics, [ENERGY_BOUNDS, POSITION_BOUNDS])):
                ax.add_patch(Rectangle((2*c+j-.5, i-.5), 1, 1,
                                       facecolor=metric_color(value, bounds), edgecolor="white", linewidth=.4))
                ax.text(2*c+j, i, formatted(value) + ("†" if bad_update else ""),
                        ha="center", va="center", fontsize=9.5, color="#22282d")

    # Group boundaries: baseline, then EL/drop, EL/keep, LE/drop, LE/keep.
    for y in [.5, 11.5, 22.5, 33.5, 44.5]:
        labels.axhline(y, color="#73808a", linewidth=.8)
        ax.axhline(y, color="#73808a", linewidth=.8)
    for x in [1.5, 3.5, 5.5, 7.5]:
        ax.axvline(x, color="#73808a", linewidth=.9)

    # Color bins stay fixed across all six pages and all five physical cases.
    for y, title, bins in [(.139, "Energy error (%)", ["≤ 0.1", "≤ 1", "≤ 5", "≤ 20", "≤ 100", "> 100"]),
                           (.119, "Position RMS (m)", ["≤ 0.01", "≤ 0.05", "≤ 0.1", "≤ 0.25", "≤ 1", "> 1"])]:
        fig.text(.026, y+.004, title, fontsize=10, weight="bold", va="center")
        for j, (color, label) in enumerate(zip(PALETTE, bins)):
            x = .148 + j*.065
            fig.add_artist(Rectangle((x, y-.003), .061, .014, transform=fig.transFigure,
                                     facecolor=color, edgecolor="white", linewidth=.5))
            fig.text(x+.0305, y+.004, label, ha="center", va="center", fontsize=9.5)
    fig.text(.563, .143, f"225 runs on this page: {225-stopped} complete, {stopped} STOP",
             fontsize=11, weight="bold", va="center")
    fig.text(.563, .122, f"† Completed with incompatible local updates: {incompatible} runs",
             fontsize=10, va="center")
    notes = [
        "Energy error = 100 maxₙ |Eₙ − E₀| / E_scale.  E_scale = initial energy above the minimum potential (J); E = endpoint kinetic + potential energy.",
        "Position RMS = root mean square distance to the RK4 reference (m), over the full 10 s. STOP = incomplete run; prefix metrics are deliberately hidden.",
        "EL: energy → material; LE: material → energy. rₓ keep/drop: retain/drop objective stationarity residual. γ off: no rejection; γ > 0: reject oversized energy corrections.",
        "fresh: current gradient; lag_sweep: previous sweep; frozen_step: timestep predictor; previous_step: preceding timestep. All plots use the same color bins.",
    ]
    for i, note in enumerate(notes):
        fig.text(.026, .092-i*.017, note, fontsize=9.5, color="#43505b", va="top")
    return fig, drawn_keys, {"hz": hz, "sweeps": sweeps, "runs": 225, "completed": 225-stopped,
                             "stopped": stopped, "completed_with_incompatible_updates": incompatible}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--summary", type=Path, default=Path("results/gradient_policies/summary.csv"))
    parser.add_argument("--output", type=Path, default=Path("results/gradient_policies/scorecards"))
    parser.add_argument("--dpi", type=int, default=180)
    args = parser.parse_args()
    rows, index = read_summary(args.summary)
    args.output.mkdir(parents=True, exist_ok=True)
    all_drawn, pages = [], []
    pdf_path = args.output / "all_gradient_policy_scorecards.pdf"
    with PdfPages(pdf_path, metadata={"Title": "All endpoint-energy gradient-policy results",
                                     "Subject": "1,350 runs; incomplete prefixes marked STOP"}) as pdf:
        for hz, sweeps in SETTINGS:
            fig, drawn, page = make_card(index, hz, sweeps)
            stem = f"scorecard_h{hz}_sweeps{sweeps}"
            fig.savefig(args.output / f"{stem}.png", dpi=args.dpi, facecolor="white")
            pdf.savefig(fig, facecolor="white")
            plt.close(fig)
            all_drawn.extend(drawn)
            pages.append(dict(page, png=f"{stem}.png"))
    if len(all_drawn) != 1350 or len(set(all_drawn)) != 1350 or set(all_drawn) != set(index):
        raise AssertionError("Rendered scorecard coverage is not exactly the input matrix")
    stopped = sum(r["completed_steps"] != r["requested_steps"] for r in rows)
    audit = {"source": str(args.summary.resolve()),
             "source_sha256": hashlib.sha256(args.summary.read_bytes()).hexdigest(),
             "source_runs": len(rows), "unique_source_keys": len(index),
             "rendered_run_entries": len(all_drawn), "unique_rendered_keys": len(set(all_drawn)),
             "metric_slots": len(all_drawn)*2, "numerical_metric_values": (len(rows)-stopped)*2,
             "suppressed_prefix_metric_values": stopped*2, "stopped_runs": stopped,
             "statuses": dict(collections.Counter(r["status"] for r in rows)),
             "key_fields": ["case", "order", "gradient", "gamma", "retain_rx", "hz", "sweeps"],
             "page_count": len(pages), "pages": pages,
             "energy_percent_formula": "100 * max_abs_energy_error / Escale",
             "Escale_definition": "initial energy above minimum potential, in joules",
             "position_units": "metres", "rounded_significant_digits": 3,
             "failed_prefix_metrics_displayed": False,
             "pdf": pdf_path.name, "rendered_keys": all_drawn}
    (args.output / "scorecard_audit.json").write_text(json.dumps(audit, indent=2)+"\n")
    print(json.dumps({k:v for k,v in audit.items() if k != "rendered_keys"}, indent=2))


if __name__ == "__main__":
    main()
