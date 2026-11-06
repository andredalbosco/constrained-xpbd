#!/usr/bin/env python3
"""Reproducible, completion-aware analysis of the composition experiment suite.

No trajectory with an incomplete prefix participates in motion comparisons.
Physical cohorts always retain case, timestep, duration and stress parameters.
Reported residual thresholds are diagnostics, not proofs of convergence or
visual quality. Compute timings collected during parallel runs are diagnostic.
"""
from __future__ import annotations
import argparse
import collections
import csv
import json
import math
import os
from pathlib import Path
import warnings

os.environ.setdefault("MPLCONFIGDIR", "/tmp/constrained-xpbd-composition-matplotlib")
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.lines import Line2D

COHORT = ["case", "hz", "duration", "stiffness_factor", "mass_ratio", "perturbation"]
CONFIG = ["family", "work", "treatment", "kinematics", "schedule", "curvature", "beta", "gamma", "iterations", "tangency", "keep_rx", "relinearize", "reverse"]
TIERS = ["core", "parameters", "solver", "long", "stress", "traces", "convergence"]
TEXT = {"id", "tier", "case", "family", "work", "treatment", "kinematics", "schedule", "curvature", "status", "reason", "source"}
WORKS = ["start", "end", "force_average", "separate_average", "dg_start", "dg_end", "dg_average", "frozen_power", "kinetic_start", "kinetic_end"]
KINEMATICS = ["hard_endpoint", "soft_endpoint", "soft_midpoint"]
CASES = ["scalar_release", "radial_orbit", "rod_release", "elastic_mixed", "chain_mixed"]
COLORS = {"independent": "#0072B2", "force_momentum": "#009E73", "xpbd": "#666666", "EL": "#D55E00", "fepr_projection": "#CC79A7", "shake": "#E69F00", "rattle": "#56B4E9"}
LABELS = {"independent": "Independent velocity", "force_momentum": "Momentum / work", "xpbd": "XPBD", "EL": "EL / LE", "fepr_projection": "Energy projection", "shake": "SHAKE", "rattle": "RATTLE"}
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10, "axes.titlesize": 11, "axes.labelsize": 10, "legend.fontsize": 9, "figure.dpi": 130, "savefig.dpi": 190, "axes.spines.top": False, "axes.spines.right": False})


def number(value, default=float("nan")):
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def read_csv(path):
    if not path.exists():
        return []
    result = []
    with path.open(newline="") as stream:
        for raw in csv.DictReader(stream):
            # Ignore a concurrently-written final partial row. It will appear
            # as pending until a subsequent analysis sees the complete record.
            if None in raw or any(v is None for v in raw.values()):
                continue
            row = {k: v if k in TEXT else number(v) for k, v in raw.items()}
            row["source"] = str(path)
            result.append(row)
    return result


def key(row, names=COHORT):
    return tuple(row.get(k) for k in names)


def canonical_key(row):
    values = dict(row)
    if row.get("family") == "independent" and row.get("treatment") in ("free", "global"):
        values.update(treatment="global", work="unused", schedule="joint", gamma=0., reverse=0.)
    if row.get("kinematics") == "hard_endpoint":
        values["beta"] = 0.
    if not str(row.get("case", "")).startswith("rod"):
        values["tangency"] = 0.
    return key(values, COHORT + CONFIG)


def annotate(row):
    r = dict(row)
    r["complete"] = r.get("status") == "completed" and r.get("completed", -1) == r.get("steps", -2) and r.get("steps", 0) > 0
    r["reference_ok"] = r.get("reference_valid") == 1
    r["motion_eligible"] = r["complete"] and r["reference_ok"] and all(math.isfinite(r.get(f, float("nan"))) for f in ["rms_position", "rms_velocity", "path_ratio", "max_energy_rel"])
    r["path_log_error"] = abs(math.log(max(r.get("path_ratio", 0), 1e-12)))
    r["path_nearly_frozen"] = r.get("path_ratio", 0) < .1
    r["balance_small"] = r.get("max_balance", float("inf")) <= 1e-6
    r["independent_stationarity_small"] = r.get("family") == "independent" and r.get("max_rx", float("inf")) <= 1e-5
    r["independent_residual_pass"] = r["complete"] and r["balance_small"] and r["independent_stationarity_small"] and (not str(r.get("case", "")).startswith("rod") or (r.get("max_rod_error", float("inf")) <= 1e-6 and (not r.get("tangency") or r.get("max_tangency", float("inf")) <= 1e-5)))
    # Explicit illustrative scales; score is only compared inside one cohort.
    # It does not certify visual quality or replace the individual metrics.
    r["balanced_score"] = max(r.get("rms_position", float("inf")) / .01, r.get("max_energy_rel", float("inf")) / .01, r["path_log_error"] / .1) if r["motion_eligible"] else float("inf")
    r["compute_us_per_step"] = 1000 * r.get("compute_ms", 0) / max(r.get("completed", 0), 1)
    r["visits_per_step"] = r.get("visits", 0) / max(r.get("completed", 0), 1)
    r["evaluations_per_step"] = r.get("evaluations", 0) / max(r.get("completed", 0), 1)
    return r


def deduplicate(rows):
    grouped = collections.defaultdict(list)
    for r in rows:
        grouped[canonical_key(r)].append(r)
    result = []
    for group in grouped.values():
        # Prefer a complete/reference-certified replicate, never a low-error
        # prefix. Equivalent aliases otherwise use a stable lexical choice.
        group.sort(key=lambda r: (not r["motion_eligible"], not r["complete"], r["id"]))
        representative = dict(group[0])
        representative["alias_replicates"] = len(group)
        representative["alias_ids"] = ";".join(r["id"] for r in group)
        result.append(representative)
    return result


def write_csv(path, rows, columns=None):
    if columns is None:
        columns = list(dict.fromkeys(k for row in rows for k in row))
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=columns, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def median(rows, field):
    values = [r[field] for r in rows if math.isfinite(r.get(field, float("nan")))]
    return float(np.median(values)) if values else ""


def aggregate(rows, fields):
    groups = collections.defaultdict(list)
    for row in rows:
        groups[key(row, fields)].append(row)
    output = []
    for k, group in sorted(groups.items(), key=lambda item: str(item[0])):
        eligible = [r for r in group if r["motion_eligible"]]
        done = [r for r in group if r["complete"]]
        out = dict(zip(fields, k))
        out.update(configurations=len(group), complete=len(done), incomplete=len(group)-len(done), reference_valid_complete=len(eligible), independent_residual_pass=sum(r["independent_residual_pass"] for r in group), nearly_frozen_complete=sum(r["path_nearly_frozen"] for r in eligible))
        for metric in ["max_energy_rel", "rms_position", "rms_velocity", "path_ratio", "max_rod_error", "max_tangency", "max_balance", "max_rx", "balanced_score", "compute_us_per_step", "visits_per_step"]:
            out["median_" + metric] = median(eligible, metric)
        output.append(out)
    return output


def pareto_mask(rows):
    if not rows:
        return []
    # Floors prevent machine-level energy differences from artificially making
    # every otherwise worse configuration nondominated.
    values = np.array([[max(r["rms_position"], 1e-8), max(r["rms_velocity"], 1e-7), max(r["max_energy_rel"], 1e-10), max(r["path_log_error"], 1e-6)] for r in rows])
    keep = np.ones(len(rows), dtype=bool)
    for i, v in enumerate(values):
        keep[i] = not np.any(np.all(values <= v, axis=1) & np.any(values < v, axis=1))
    return keep.tolist()


def selections(rows):
    groups = collections.defaultdict(list)
    for row in rows:
        if row["motion_eligible"]:
            groups[key(row)].append(row)
    output, evidence = [], []
    for cohort, group in sorted(groups.items(), key=lambda item: str(item[0])):
        group.sort(key=lambda r: (r["balanced_score"], r["rms_position"], r["id"]))
        frontier = pareto_mask(group)
        baselines = [r for r in group if r["family"] == "xpbd"]
        baseline = min(baselines, key=lambda r: (r["rms_position"], r["max_energy_rel"])) if baselines else None
        for i, (r, on_frontier) in enumerate(zip(group, frontier)):
            out = dict(r)
            out["balanced_rank_within_cohort"] = i + 1
            out["pareto_motion_energy_path"] = on_frontier
            if baseline:
                out["matched_xpbd_id"] = baseline["id"]
                for field in ["rms_position", "rms_velocity", "max_energy_rel", "path_log_error", "compute_us_per_step"]:
                    out["ratio_to_xpbd_" + field] = r[field] / max(baseline[field], 1e-12)
                motion_ok = r["rms_position"] <= 1.05 * baseline["rms_position"] and r["rms_velocity"] <= 1.05 * baseline["rms_velocity"] and r["path_log_error"] <= max(.01, 1.05 * baseline["path_log_error"])
                improves = r["rms_position"] < .8 * baseline["rms_position"] or r["max_energy_rel"] < .8 * baseline["max_energy_rel"]
                out["matched_candidate_gate"] = r["family"] != "xpbd" and motion_ok and improves and r["max_energy_rel"] <= max(1e-8, 1.05 * baseline["max_energy_rel"]) and not r["path_nearly_frozen"]
                if out["matched_candidate_gate"]:
                    evidence.append(out)
            if on_frontier or i < 5:
                output.append(out)
    return output, evidence


def select_cohort(rows, case, prefer_tier="core"):
    candidates = [r for r in rows if r["case"] == case]
    if not candidates:
        return []
    groups = collections.defaultdict(list)
    for r in candidates:
        groups[key(r)].append(r)
    # Prefer the planned nominal 120-Hz, 10-s slice even while it is partial.
    def preference(item):
        k, g = item
        return (k[1] == 120, k[2] == 10, k[3:] == (1., 1., 1.), any(r["tier"] == prefer_tier for r in g), len(g), str(k))
    return max(groups.items(), key=preference)[1]


def cohort_title(rows, fallback):
    if not rows:
        return fallback + " (pending)"
    r = rows[0]
    return f"{r['case'].replace('_', ' ')} · {r['hz']:g} Hz · {r['duration']:g} s\nk×{r['stiffness_factor']:g}, mass ratio {r['mass_ratio']:g}, perturbation {r['perturbation']:g}"


def save(fig, out, name):
    fig.savefig(out / (name + ".png"), bbox_inches="tight", facecolor="white")
    fig.savefig(out / (name + ".pdf"), bbox_inches="tight", facecolor="white")
    plt.close(fig)


def plots(rows, out):
    fig, axes = plt.subplots(2, 3, figsize=(14, 8), constrained_layout=True)
    for ax, case in zip(axes.flat, CASES):
        cohort = select_cohort(rows, case)
        eligible = [r for r in cohort if r["motion_eligible"]]
        for family, color in COLORS.items():
            group = [r for r in eligible if r["family"] == family]
            ax.scatter([max(r["rms_position"], 1e-9) for r in group], [max(r["max_energy_rel"], 1e-12) for r in group], c=color, s=18, alpha=.55, linewidths=0, label=LABELS[family])
        frozen = [r for r in eligible if r["path_nearly_frozen"]]
        ax.scatter([max(r["rms_position"], 1e-9) for r in frozen], [max(r["max_energy_rel"], 1e-12) for r in frozen], facecolors="none", edgecolors="black", s=45, linewidths=.6)
        ax.set(xscale="log", yscale="log", xlabel="Position RMS error (m)", ylabel="Maximum relative energy error", title=cohort_title(cohort, case))
        ax.grid(alpha=.18, which="both")
        ax.text(.98, .03, f"complete / certified: {len(eligible)}\nincomplete excluded: {sum(not r['complete'] for r in cohort)}", ha="right", va="bottom", transform=ax.transAxes, fontsize=8)
    axes.flat[-1].axis("off")
    handles = [Line2D([], [], marker="o", ls="", color=c, label=LABELS[f]) for f, c in COLORS.items()]
    handles.append(Line2D([], [], marker="o", ls="", color="black", markerfacecolor="none", label="Path length < 10% of reference"))
    axes.flat[-1].legend(handles=handles, loc="upper left", frameon=False)
    axes.flat[-1].text(.02, .30, "Each panel is one matched physical cohort.\nPrefixes from failed runs never enter this plot.\nLow energy error alone does not imply good motion.\nMachine-level values are floored for log axes.", transform=axes.flat[-1].transAxes, va="top", fontsize=10)
    fig.suptitle("Energy preservation and trajectory accuracy", fontsize=15)
    save(fig, out, "energy_vs_motion")

    fig, axes = plt.subplots(2, 5, figsize=(18, 12), constrained_layout=True)
    image = None
    for col, case in enumerate(CASES):
        cohort = select_cohort(rows, case)
        for rownum, treatment in enumerate(["prescribed", "penalized"]):
            ax = axes[rownum, col]
            mat = np.full((len(WORKS), 3), np.nan)
            states = {}
            for y, work in enumerate(WORKS):
                for x, kin in enumerate(KINEMATICS):
                    matches = [r for r in cohort if r["family"] == "independent" and r["work"] == work and r["treatment"] == treatment and r["kinematics"] == kin and r["schedule"] == "joint" and r["curvature"] == "drop" and r["beta"] == 1 and r["gamma"] == 1 and r["iterations"] == 4 and r["keep_rx"] == 1 and r["relinearize"] == 1 and r["reverse"] == 0 and r["tangency"] == int(case.startswith("rod"))]
                    if matches:
                        r = matches[0]
                        if r["motion_eligible"]:
                            mat[y, x] = math.log10(max(r["balanced_score"], 1e-3))
                            states[y, x] = "F" if r["path_nearly_frozen"] else ("!" if not r["independent_residual_pass"] else "")
                        else:
                            states[y, x] = "×" if not r["complete"] else "R"
            image = ax.imshow(mat, vmin=-2, vmax=3, cmap="viridis_r", aspect="auto")
            ax.set_xticks(range(3), ["Hard\nendpoint", "Soft\nendpoint", "Soft\nmidpoint"])
            ax.set_yticks(range(len(WORKS)), [w.replace("_", " ") for w in WORKS] if col == 0 else [""] * len(WORKS))
            detail = f"{cohort[0]['hz']:g} Hz · {cohort[0]['duration']:g} s" if cohort else "pending"
            ax.set_title(case.replace("_", " ") + "\n" + detail + " · " + treatment)
            for (y, x), status in states.items():
                if status:
                    ax.text(x, y, status, ha="center", va="center", color="black" if not math.isfinite(mat[y, x]) or mat[y, x] < 0 else "white", fontsize=10)
    fig.colorbar(image, ax=axes, label="log₁₀ max(position RMS / 0.01 m, energy error / 0.01, |log path ratio| / 0.1)", shrink=.7)
    fig.suptitle("Work × velocity target: fixed solver slice\nJoint solve; 4 iterations; β=γ=1; dropped multiplier curvature. F: nearly frozen, !: residual threshold missed, ×: incomplete, R: reference uncertified; white cell: pending.", fontsize=12)
    save(fig, out, "work_kinematics_quality")

    fig, axes = plt.subplots(2, 2, figsize=(12, 8), constrained_layout=True)
    schedules = ["joint", "constraint_block", "separate", "symmetric"]
    schedule_colors = ["#0072B2", "#D55E00", "#009E73", "#CC79A7"]
    for rownum, case in enumerate(["elastic_mixed", "chain_mixed"]):
        cohort = select_cohort([r for r in rows if r["tier"] == "solver"], case, "solver")
        for col, metric in enumerate(["max_balance", "max_rx"]):
            ax = axes[rownum, col]
            for schedule, color in zip(schedules, schedule_colors):
                for curvature, style in [("drop", "--"), ("exact", "-")]:
                    group = [r for r in cohort if r["complete"] and r["family"] == "independent" and r["work"] == "dg_average" and r["treatment"] == "penalized" and r["kinematics"] == "soft_endpoint" and r["schedule"] == schedule and r["curvature"] == curvature and r["beta"] == 1 and r["gamma"] == 1 and r["keep_rx"] == 1 and r["relinearize"] == 1 and r["reverse"] == 0]
                    group.sort(key=lambda r: r["iterations"])
                    if group:
                        ax.plot([r["iterations"] for r in group], [max(r[metric], 1e-14) for r in group], color=color, linestyle=style, marker="o", markersize=4, label=schedule.replace("_", " ") + ", " + curvature)
            ax.set(xscale="log", yscale="log", xlabel="Iteration budget per timestep", ylabel="Maximum normalized balance residual" if col == 0 else "Maximum scaled stationarity residual", title=cohort_title(cohort, case))
            ax.grid(alpha=.2, which="both")
    handles, labels = axes[0, 0].get_legend_handles_labels()
    if handles:
        fig.legend(handles, labels, loc="outside lower center", ncol=4, frameon=False)
    fig.suptitle("Solver schedules and iteration budgets\nIndependent velocity; DG-average work; penalized work; soft endpoint target; β=γ=1; retained stationarity; fresh derivatives.", fontsize=12)
    save(fig, out, "solver_residuals")

    fig, axes = plt.subplots(1, 2, figsize=(12, 5), constrained_layout=True)
    for ax, case in zip(axes, ["rod_release", "rod_kick"]):
        cohort = select_cohort(rows, case)
        eligible = [r for r in cohort if r["motion_eligible"]]
        for family, color in COLORS.items():
            group = [r for r in eligible if r["family"] == family]
            ax.scatter([max(r["max_tangency"], 1e-12) for r in group], [max(r["max_energy_rel"], 1e-12) for r in group], s=24, alpha=.65, c=color, linewidths=0, label=LABELS[family])
        frozen = [r for r in eligible if r["path_nearly_frozen"]]
        ax.scatter([max(r["max_tangency"], 1e-12) for r in frozen], [max(r["max_energy_rel"], 1e-12) for r in frozen], s=55, facecolors="none", edgecolors="black", linewidths=.7)
        ax.set(xscale="log", yscale="log", xlabel="Maximum radial endpoint velocity (m/s)", ylabel="Maximum relative energy error", title=cohort_title(cohort, case))
        ax.grid(alpha=.2, which="both")
    handles = [Line2D([], [], marker="o", ls="", color=c, label=LABELS[f]) for f, c in COLORS.items()]
    fig.legend(handles=handles, loc="outside lower center", ncol=4, frameon=False)
    fig.suptitle("Rigid pendulum: energy and endpoint tangency\nBlack rings flag path length below 10% of reference; completed, certified runs only.", fontsize=13)
    save(fig, out, "rigid_energy_tangency")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    direct = (args.results / "summary.csv").exists()
    tiers = TIERS + sorted(p.name for p in args.results.iterdir() if p.is_dir() and p.name not in TIERS and not p.name.startswith("smoke") and (p / "summary.csv").exists())
    paths = [args.results / "summary.csv"] if direct else [args.results / tier / "summary.csv" for tier in tiers]
    rows = [annotate(r) for path in paths for r in read_csv(path)]
    if not rows:
        print("No complete CSV records yet; analysis awaits results.")
        return
    unique = deduplicate(rows)
    manifest_dir = (args.results.parent.parent if direct else args.results.parent) / "manifests"
    planned = []
    manifest_tiers = sorted(set(r["tier"] for r in rows)) if direct else tiers
    for tier in manifest_tiers:
        planned.extend(read_csv(manifest_dir / (tier + ".csv")))
    observed = {r["id"]: r for r in rows}
    coverage_groups = collections.defaultdict(list)
    for r in planned or rows:
        coverage_groups[key(r, ["tier", "family"] + COHORT)].append(r)
    coverage = []
    for k, group in sorted(coverage_groups.items(), key=lambda item: str(item[0])):
        seen = [observed[r["id"]] for r in group if r["id"] in observed]
        out = dict(zip(["tier", "family"] + COHORT, k))
        out.update(planned=len(group), observed=len(seen), pending=len(group)-len(seen), complete=sum(r["complete"] for r in seen), incomplete=sum(not r["complete"] for r in seen), reference_valid_complete=sum(r["motion_eligible"] for r in seen), independent_residual_pass=sum(r["independent_residual_pass"] for r in seen))
        coverage.append(out)
    write_csv(args.out / "coverage.csv", coverage)
    family = aggregate(unique, COHORT + ["family"])
    write_csv(args.out / "family_summary.csv", family)
    composition = aggregate([r for r in unique if r["family"] == "independent"], COHORT + CONFIG[1:])
    write_csv(args.out / "work_composition.csv", composition)
    best, evidence = selections(unique)
    best_columns = ["id", "tier"] + COHORT + CONFIG + ["complete", "reference_ok", "independent_residual_pass", "pareto_motion_energy_path", "balanced_rank_within_cohort", "balanced_score", "max_energy_rel", "rms_position", "rms_velocity", "path_ratio", "max_rod_error", "max_tangency", "max_balance", "max_rx", "compute_us_per_step", "visits_per_step", "matched_xpbd_id", "matched_candidate_gate", "ratio_to_xpbd_rms_position", "ratio_to_xpbd_rms_velocity", "ratio_to_xpbd_max_energy_rel", "ratio_to_xpbd_path_log_error", "ratio_to_xpbd_compute_us_per_step", "alias_replicates", "alias_ids"]
    write_csv(args.out / "best_combinations.csv", best, best_columns)
    write_csv(args.out / "candidate_evidence.csv", evidence, best_columns)
    failures = []
    for r in rows:
        categories = []
        if not r["complete"]:
            categories.append("incomplete_execution")
        if not r["reference_ok"]:
            categories.append("uncertified_reference")
        if r["complete"] and r["family"] == "independent" and not r["independent_residual_pass"]:
            categories.append("completed_residual_threshold_missed")
        if r["complete"] and r["path_nearly_frozen"]:
            categories.append("nearly_frozen_path")
        if categories:
            failures.append(dict(r, failure_categories=";".join(categories)))
    write_csv(args.out / "failures.csv", failures, ["id", "tier"] + COHORT + CONFIG + ["failure_categories", "status", "reason", "completed", "steps", "reference_valid", "max_energy_rel", "rms_position", "rms_velocity", "path_ratio", "max_balance", "max_rx", "max_rod_error", "max_tangency"])
    plots(unique, args.out)
    data = {
        "results_root": str(args.results.resolve()), "loaded_files": [str(p) for p in paths if p.exists()],
        "planned": len(planned), "observed": len(rows), "pending": max(0, len(planned)-len(observed)), "unique_compositions": len(unique),
        "complete": sum(r["complete"] for r in rows), "incomplete": sum(not r["complete"] for r in rows), "certified_complete": sum(r["motion_eligible"] for r in rows),
        "statuses": dict(collections.Counter(r["status"] for r in rows)), "failure_reasons": dict(collections.Counter(r["reason"] for r in rows if not r["complete"])),
        "independent_residual_pass": sum(r["independent_residual_pass"] for r in rows), "nearly_frozen_complete": sum(r["complete"] and r["path_nearly_frozen"] for r in rows),
        "matched_candidate_evidence_rows": len(evidence), "matched_candidate_families": dict(collections.Counter(r["family"] for r in evidence)),
        "cohort_fields": COHORT,
        "rules": {
            "completion": "status=completed AND completed=steps AND steps>0; failed prefixes excluded from all rankings and figures",
            "reference": "RMS comparisons require reference_valid=1; certified by paired reference resolutions in study.cpp",
            "residuals": "independent solver pass: max_balance<=1e-6 and max_rx<=1e-5; rods additionally length<=1e-6 and requested tangency<=1e-5; completed does not imply this pass",
            "pareto": "position RMS, velocity RMS, maximum relative energy error, absolute log(path ratio); floors 1e-8,1e-7,1e-10,1e-6 respectively",
            "balanced_score": "max(position RMS/0.01 m, energy error/0.01, absolute log(path ratio)/0.1); illustrative disclosed scales, used only within identical physical cohorts",
            "nearly_frozen": "path_ratio<0.1; diagnostic, not a visual-quality assessment",
            "xpbd_comparator": "complete/certified XPBD with smallest position RMS in the same case,hz,duration,stiffness,mass-ratio,perturbation cohort",
            "candidate_gate": "complete/certified; position and velocity RMS <=1.05x matched XPBD; log-path error <=max(0.01,1.05x baseline); energy <=max(1e-8,1.05x baseline); either position or energy improves at least20%; path ratio>=0.1; residual pass is reported separately",
            "timing": "compute-only, parallel-contention affected; relative timings are diagnostic until a serial benchmark",
            "aliases": "explicit free/global-energy aliases deduplicated; ignored hard-endpoint beta and non-rigid tangency deduplicated",
            "pending": "manifest jobs without a CSV record are pending, not failed executions",
        },
        "coverage": coverage,
        "top_balanced_by_cohort": [{k: r[k] for k in best_columns if k in r} for r in best if r["balanced_rank_within_cohort"] <= 3],
    }
    def json_safe(value):
        if isinstance(value, dict):
            return {k: json_safe(v) for k, v in value.items()}
        if isinstance(value, list):
            return [json_safe(v) for v in value]
        if isinstance(value, float) and not math.isfinite(value):
            return None
        return value
    landscape_path = args.results / "landscape.csv"
    if landscape_path.exists():
        # Preserve original column types here: landscape experiment fields are
        # not part of the trajectory CSV schema.
        with landscape_path.open(newline="") as stream:
            data["scalar_landscape"] = list(csv.DictReader(stream))
    (args.out / "report_data.json").write_text(json.dumps(json_safe(data), indent=2, allow_nan=False) + "\n")
    notes = ["# Numerical findings", "", f"Observed **{len(rows):,} / {len(planned):,}** planned runs. **{data['complete']:,}** completed their full requested duration; **{data['incomplete']:,}** did not. Completion is not convergence: **{data['independent_residual_pass']:,}** independent-velocity runs pass the stated balance, stationarity and requested geometric thresholds. **{data['nearly_frozen_complete']:,}** completed trajectories have less than 10% of the reference path length.", "", "All selections below require full completion and a certified reference. Each comparison holds case, timestep, duration and stress factors fixed. The balanced score is the maximum of position RMS / 0.01 m, relative energy error / 0.01, and absolute log(path ratio) / 0.1. These are disclosed illustrative scales, not a visual-quality assessment.", "", "## Nominal 120 Hz, 10 second comparisons", "", "| Case | Family | Configuration | Work / velocity target | Position RMS (m) | Max relative energy error | Path ratio | Balance residual | Stationarity residual |", "|---|---|---|---|---:|---:|---:|---:|---:|"]
    for case in CASES:
        cohort = [r for r in unique if r["motion_eligible"] and r["case"] == case and r["hz"] == 120 and r["duration"] == 10 and r["stiffness_factor"] == r["mass_ratio"] == r["perturbation"] == 1]
        for family_name in ["xpbd", "independent", "force_momentum", "rattle"]:
            group = [r for r in cohort if r["family"] == family_name]
            if not group:
                continue
            r = min(group, key=lambda r: (r["balanced_score"], r["rms_position"], r["id"]))
            notes.append(f"| {case} | {family_name} | {r['id']} | {r['work']} / {r['kinematics']} | {r['rms_position']:.6g} | {r['max_energy_rel']:.6g} | {r['path_ratio']:.6g} | {r['max_balance']:.6g} | {r['max_rx']:.6g} |")
    notes += ["", "Residual columns have method-specific meanings and normalizations; compare stationarity within a family. The independent solver's threshold-pass flag is provided separately in the CSVs. Its small energy error may coexist with a poor trajectory or nonzero optimization residual. An exact work identity does not imply the numerical solve reached its constraints.", "", "## Solver evidence", "", "The dedicated solver tier uses its actual duration from the manifest (0.1 seconds here); it diagnoses short-horizon iteration behavior and is not a long-horizon stability test. The convergence tier uses 64 iterations and its own one-second horizon. The long tier separately tests 60 seconds. Curvature, scheduling and iteration comparisons must retain those distinctions.", "", "True edge blocks include both endpoint kinetic rows and keep shared vertex multipliers. In the fixed DG-average, penalized, soft-endpoint slice shown in solver_residuals.png, joint and overlapping edge blocks reduce balance residuals much more than storage/kinetic split passes on the two-spring chain. More iterations do not consistently repair the split schedule at the budgets tested.", "", "## Scalar objective bias, independent of solve convergence", "", "For release from rest with positive spring extension C_n, write omega²=k/m and v=−h omega² C_n s. Physical acceleration has s→1. Exact energy conservation gives x−x_n=−h² omega² C_n s²/2+O(h⁴). Substituting into the displayed independent-velocity objective gives the following leading nonzero stationary equations for the moving branch:", "", "- Soft endpoint: (beta+1)s²−3s+2(1−beta)=0; its moving minimum has s=[3+sqrt(1+8 beta²)]/[2(1+beta)]. At beta=1, s=1.5.", "- Soft midpoint, using the implemented unit-v residual v+v_n−2(x−x_n)/h: (beta+4)s²−6s+2(1−beta)=0; the moving minimum has s=[3+sqrt(1+6 beta+2 beta²)]/(beta+4). At beta=1, s=1.2; at beta=0.1, s≈1.042144.", "", "The landscape experiment confirms this distinction: at 120→240 Hz, scalar-release velocity/reference ratios approach 1.49397→1.49848 for endpoint beta=1 and 1.19740→1.19935 for midpoint beta=1, despite approximately solved objectives and energy balances. Midpoint beta=0.1 gives 1.04118→1.04190. Therefore this particular finite-compliance objective carries a release-acceleration bias, even when solved accurately. These statements apply to the chosen global-energy objective and moving branch; they do not establish a universal limitation of every independent-velocity formulation. The frozen state can also be a local minimum for beta<1, so initialization and branch selection remain relevant.", "", "## Interpretation limits", "", "- Momentum-closed midpoint spring methods and RATTLE are reference candidates supported by matched trajectory evidence, not automatic replacements for the current solver.", "- The work_kinematics_quality figure fixes four iterations, joint scheduling, beta=gamma=1, dropped multiplier curvature, and rigid tangency when relevant. It illustrates composition effects at that finite budget; it is not a ranking of converged optimization objectives.", "- Candidate evidence is screened against the best certified XPBD position RMS in the same physical cohort. The candidate gate requires both position and velocity accuracy to stay within 5%, path agreement to remain comparable, and at least 20% improvement in either position or energy. Residual passage is reported independently.", "- Family medians summarize heterogeneous tested configurations within one physical cohort. They do not establish a universal property of every member of a family.", "- Failed prefixes remain in failures.csv and never enter rankings. Pending jobs are identified separately in coverage.csv.", "- Timings are compute-only and affected by parallel contention; use a serial benchmark before making performance claims.", "- No visual-quality claim is inferred solely from RMS, energy, path length, or residual statistics.", ""]
    (args.out / "findings.md").write_text("\n".join(notes))
    print(json.dumps({k: data[k] for k in ["planned", "observed", "pending", "complete", "incomplete", "unique_compositions", "matched_candidate_evidence_rows"]}, indent=2))


if __name__ == "__main__":
    main()
