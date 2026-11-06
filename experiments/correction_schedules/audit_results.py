#!/usr/bin/env python3
"""Read-only independent counter, filter, trajectory, and reference audit.

Usage: python3 audit_results.py [experiment_directory] [report_json]
The optional report path writes only the audit result, never simulation data.
"""
import collections
import csv
import json
import math
from pathlib import Path
import sys


def audit_summaries(root):
    reports = {}
    issues = []
    for directory in ("sweep", "selected"):
        with (root / directory / "summary.csv").open() as source:
            rows = list(csv.DictReader(source))
        maximum_cap_excess = 0.0
        for row in rows:
            key = (directory, row["case"], row["method"], row["substeps"])
            n = int(row["evaluations"])
            singular = int(row["singular_updates"])
            nonzero = n - singular
            schedule = row["schedule"]
            fallback = int(row["guarded_fallbacks"])
            completed = int(row["completed_steps"])
            aborted = row["status"].startswith("aborted")

            def check(condition, message):
                if not condition:
                    issues.append([*key, message])

            check(n == completed + int(aborted), "one evaluation per attempted substep")
            check(int(row["linearizations"]) == 1, "one linearization")
            check(int(row["requested_steps"]) == round(float(row["duration"]) / float(row["h"])), "duration/substeps")
            if schedule == "baseline":
                scalar, block, spectral, attempts = n, 0, 0, 0
            elif schedule in ("jacobi", "EL", "LE", "alternating"):
                scalar, block, spectral, attempts = 2 * nonzero + singular, 0, 0, n
            elif schedule == "LEL":
                scalar, block, spectral, attempts = 3 * nonzero + singular, 0, 0, n
            elif schedule == "ELE":
                scalar, block, spectral, attempts = 3 * nonzero + singular, 0, 0, 2 * nonzero + singular
            elif schedule == "joint":
                scalar, block, spectral, attempts = singular, nonzero, 0, n
            elif schedule == "guarded_joint":
                scalar, block, spectral, attempts = 2 * fallback + singular, nonzero - fallback, 0, n
            elif schedule == "spectral":
                scalar, block, spectral, attempts = singular, 0, nonzero, n
            else:
                raise AssertionError(schedule)
            check(int(row["scalar_solves"]) == scalar and int(row["block_solves"]) == block and int(row["spectral_solves"]) == spectral, "schedule operation accounting")
            check(int(row["energy_attempts"]) == attempts, "energy-row visits")
            policy = row["policy"]
            if "cap" in policy:
                excess = float(row["max_accepted_ratio"]) - float(row["cap"])
                maximum_cap_excess = max(maximum_cap_excess, excess)
                check(excess <= 1e-12 * max(1, float(row["cap"])), "reported displacement cap")
            else:
                check(int(row["capped"]) == 0, "disabled cap counter")
            if "direction" not in policy:
                check(int(row["direction_skips"]) == 0, "disabled direction counter")
            if "repair" not in policy:
                check(int(row["repair_skips"]) == 0 and int(row["repair_scaled"]) == 0, "disabled repair counters")
            if schedule != "baseline":
                check(int(row["accepted"]) + int(row["direction_skips"]) + int(row["repair_skips"]) == n, "acceptance partition")
        reports[directory] = {
            "rows": len(rows),
            "statuses": dict(collections.Counter(row["status"] for row in rows)),
            "evaluated_substeps": sum(int(row["evaluations"]) for row in rows),
            "largest_reported_cap_overshoot": maximum_cap_excess,
        }
    return reports, issues


def audit_traces(root):
    previous = None
    issues = []
    counts = collections.Counter()
    maxima = collections.defaultdict(float)
    worsening = collections.Counter()
    with (root / "selected" / "traces.csv").open() as source:
        for row in csv.DictReader(source):
            step = int(row["step"])
            key = (row["case"], row["method"], row["substeps"])
            current = tuple(float(row[field]) for field in ("x", "y", "vx", "vy"))
            if step == 0:
                previous = (key, current)
                continue
            assert previous[0] == key
            x, y, vx, vy = previous[1]
            mass, length, alpha, h, gravity = (float(row[field]) for field in ("mass", "length", "alpha", "h", "gravity"))
            rod = int(row["pendulum"])
            policy = row["policy"]
            px, py = x + h * vx, y + h * vy - h * h * gravity
            radius = math.hypot(px, py)
            nx, ny = px / radius, py / radius
            constraint = radius - length
            at = 0 if rod else alpha / (h * h)
            baseline_lambda = -constraint / (1 / mass + at)
            bx, by = px + nx * baseline_lambda / mass, py + ny * baseline_lambda / mass
            ex, ey = current[0] - bx, current[1] - by
            actual_extra = math.hypot(ex, ey)
            H0, change, blend, predicted = (float(row[field]) for field in ("H0", "candidate_change", "blend", "predicted_H"))
            maxima["extra_diagnostic_discrepancy"] = max(maxima["extra_diagnostic_discrepancy"], abs(actual_extra - float(row["extra_displacement"])))
            maxima["prediction_identity_discrepancy"] = max(maxima["prediction_identity_discrepancy"], abs(predicted - H0 - blend * change))

            def energy(xx, yy, ux, uy):
                extension = math.hypot(xx, yy) - length
                return .5 * mass * (ux * ux + uy * uy) + (0 if rod else extension * extension / (2 * alpha)) + mass * gravity * yy

            actual_H0 = energy(bx, by, (bx - x) / h, (by - y) / h) - float(row["target_energy"])
            maxima["actual_baseline_H_discrepancy"] = max(maxima["actual_baseline_H_discrepancy"], abs(actual_H0 - float(row["actual_H0"])))
            if "cap" in policy:
                old_radius = math.hypot(x, y)
                spring = 0 if rod else -(old_radius - length) / (alpha * mass * old_radius)
                scale = h * math.hypot(vx, vy) + h * h * math.hypot(spring * x, -gravity + spring * y)
                bound = float(row["cap"]) * scale
                excess = actual_extra - bound
                maxima["actual_cap_excess"] = max(maxima["actual_cap_excess"], excess)
                counts["capped_policy_updates"] += 1
                if excess > 2e-12 * (1 + bound):
                    issues.append([*key, step, "actual displacement cap", excess])
            if "repair" in policy:
                counts["repair_policy_updates"] += 1
                delta = blend * change
                fraction = float(row["fraction"])
                if abs(delta) > fraction * abs(H0) + 2e-12 * (1 + abs(H0)):
                    issues.append([*key, step, "repair fraction"])
                if abs(H0 + delta) > abs(H0) + 2e-12 * (1 + abs(H0)):
                    issues.append([*key, step, "predicted residual worsening"])
                actual_H = energy(*current) - float(row["target_energy"])
                if blend > 0:
                    counts["accepted_repair_updates"] += 1
                    if abs(actual_H) > abs(actual_H0) + 1e-10:
                        counts["accepted_repairs_worsening_actual_residual"] += 1
                        worsening[row["case"] + "/" + row["method"]] += 1
            if "direction" in policy:
                q_factor = 0 if rod else constraint / alpha
                work = q_factor * (nx * ex + ny * ey)
                if work > 2e-10 * (1 + abs(q_factor) * actual_extra):
                    issues.append([*key, step, "accepted uphill objective direction", work])
            if row["case"] in ("vertical_oscillation", "spring_release"):
                time = float(row["time"])
                omega = math.sqrt(1 / (mass * alpha))
                expected = (0, -(length + mass * gravity * alpha) - .2 * math.cos(omega * time)) if row["case"] == "vertical_oscillation" else (length + .25 * math.cos(omega * time), 0)
                error = math.hypot(float(row["ref_x"]) - expected[0], float(row["ref_y"]) - expected[1])
                maxima["reference_error_vs_analytic_oscillator"] = max(maxima["reference_error_vs_analytic_oscillator"], error)
            previous = (key, current)
            counts["trace_updates"] += 1
    assert maxima["reference_error_vs_analytic_oscillator"] < 1e-9
    return {"counts": dict(counts), "maxima": dict(maxima), "actual_repair_worsening_by_case_method": dict(worsening)}, issues


def main():
    root = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(__file__).parent
    summaries, summary_issues = audit_summaries(root)
    traces, trace_issues = audit_traces(root)
    report = {"root": str(root), "summaries": summaries, "traces": traces, "issues": summary_issues + trace_issues,
              "interpretation": "Repair constrains the frozen linear prediction. Actual nonlinear energy residual may worsen relative to the material-only endpoint even when every filter invariant passes."}
    encoded = json.dumps(report, indent=2) + "\n"
    if len(sys.argv) > 2:
        Path(sys.argv[2]).write_text(encoded)
    print(encoded)
    if report["issues"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
