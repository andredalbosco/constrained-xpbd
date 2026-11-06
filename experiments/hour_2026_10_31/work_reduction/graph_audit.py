"""Independent weighted-work elimination checks, including cycles/components.

These are graph-QP checks, not mechanical trajectory benchmarks.
"""
from pathlib import Path
import csv
import json
import numpy as np

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "graph_results"
OUT.mkdir(exist_ok=True)
rng = np.random.default_rng(20261031)


def incidence(n, edges):
    b = np.zeros((n, len(edges)))
    for j, (u, v) in enumerate(edges):
        b[u, j], b[v, j] = -1, 1
    return b


graphs = {
    "tree8": (8, [(i, i + 1) for i in range(7)]),
    "cycle8": (8, [(i, (i + 1) % 8) for i in range(8)]),
    "complete4": (4, [(i, j) for i in range(4) for j in range(i + 1, 4)]),
    "two_trees": (8, [(i, i + 1) for i in [0, 1, 2, 4, 5, 6]]),
    "two_cycles": (8, [(4*k+i, 4*k+(i+1) % 4) for k in range(2) for i in range(4)]),
}
rows = []
for name, (n, edges) in graphs.items():
    b = incidence(n, edges)
    m = b.shape[1]
    rank = np.linalg.matrix_rank(b)
    # Remove one redundant row per component using an independent row-rank scan.
    kept = []
    for i in range(n):
        candidate = b[kept + [i]]
        if np.linalg.matrix_rank(candidate) > len(kept):
            kept.append(i)
    br = b[kept]
    cycles = np.eye(m) - np.linalg.pinv(b) @ b
    for trial in range(48):
        q = np.diag(np.exp(rng.uniform(-3, 3, m)))
        if trial % 2:
            factor = rng.normal(size=(m, m)) / np.sqrt(m)
            q += factor.T @ factor
        qi = np.linalg.inv(q)
        target = rng.normal(size=m)
        delta = b @ rng.normal(size=m)
        d = delta - b @ target
        lap = b @ qi @ b.T
        dual = np.linalg.pinv(lap, rcond=1e-13) @ d
        w = target + qi @ b.T @ dual
        kkt = np.block([[q, -br.T], [br, np.zeros((rank, rank))]])
        reference = np.linalg.solve(kkt, np.r_[q @ target, delta[kept]])[:m]
        rp, cp = rng.permutation(n), rng.permutation(m)
        bp, qp = b[rp][:, cp], qi[cp][:, cp]
        wp = target[cp] + qp @ bp.T @ np.linalg.pinv(bp @ qp @ bp.T, rcond=1e-13) @ d[rp]
        errors = {
            "balance_error": np.linalg.norm(b @ w - delta),
            "full_kkt_error": np.linalg.norm(w - reference),
            "cycle_stationarity": np.linalg.norm(cycles @ q @ (w-target)),
            "permutation_error": np.linalg.norm(wp - w[cp]),
        }
        passed = all(e < 2e-9 for e in errors.values())
        rows.append(dict(graph=name, trial=trial, components=n-rank, cycle_rank=m-rank,
                         **errors, passed=int(passed)))
        assert passed, rows[-1]

# Total energy conservation alone is insufficient for disconnected components.
b = incidence(*graphs["two_trees"])
delta = np.r_[np.ones(4)/4, -np.ones(4)/4]
remaining = delta - b @ np.linalg.pinv(b) @ delta
disconnected = {
    "total_energy_change": float(delta.sum()),
    "component_energy_changes": [float(delta[:4].sum()), float(delta[4:].sum())],
    "unrepresentable_balance_norm": float(np.linalg.norm(remaining)),
}
assert abs(delta.sum()) < 1e-15 and np.linalg.norm(remaining) > .5
with (OUT / "validation.csv").open("w") as f:
    writer = csv.DictWriter(f, fieldnames=rows[0].keys())
    writer.writeheader()
    writer.writerows(rows)
report = {
    "checks": len(rows), "passed": sum(r["passed"] for r in rows),
    "maximum_errors": {key: max(r[key] for r in rows) for key in errors},
    "disconnected_counterexample": disconnected,
}
(OUT / "summary.json").write_text(json.dumps(report, indent=2) + "\n")
print(json.dumps(report, indent=2))
