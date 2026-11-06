"""Linear-time work-flow and transpose operators for a forest energy graph."""
from pathlib import Path
import csv
import numpy as np


class ForestFlow:
    def __init__(self, n, edges):
        self.n, self.edges = n, edges
        adj = [[] for _ in range(n)]
        for j, (u, v) in enumerate(edges):
            adj[u].append((v, j, 1.0))  # child v has +1 in incidence column
            adj[v].append((u, j, -1.0))
        self.parent = [-1] * n
        self.edge = [-1] * n
        self.sign = [0.0] * n
        self.components, self.order = [], []
        seen = set()
        for root in range(n):
            if root in seen:
                continue
            seen.add(root)
            component = [root]
            for u in component:
                self.order.append(u)
                for v, j, sign in adj[u]:
                    if j == self.edge[u]:
                        continue
                    if v in seen:
                        raise ValueError("energy graph contains a cycle")
                    seen.add(v)
                    self.parent[v], self.edge[v], self.sign[v] = u, j, sign
                    component.append(v)
            self.components.append(component)

    def project(self, x):
        x = np.array(x, dtype=float, copy=True)
        for c in self.components:
            x[c] -= x[c].mean()
        return x

    def apply(self, delta):
        subtotal = self.project(delta)
        work = np.zeros(len(self.edges))
        for u in reversed(self.order):
            if self.parent[u] >= 0:
                work[self.edge[u]] = subtotal[u] / self.sign[u]
                subtotal[self.parent[u]] += subtotal[u]
        return work

    def transpose(self, edge_gradient):
        potential = np.zeros(self.n)
        for u in self.order:
            if self.parent[u] >= 0:
                potential[u] = potential[self.parent[u]] + edge_gradient[self.edge[u]] / self.sign[u]
        return self.project(potential)


if __name__ == "__main__":
    rng = np.random.default_rng(98765)
    out = Path(__file__).resolve().parent / "graph_results"
    out.mkdir(exist_ok=True)
    rows = []
    for n in [4, 8, 32, 128, 1024, 16384]:
        for components in [1, 2]:
            edges = []
            for begin, end in [(i*n//components, (i+1)*n//components) for i in range(components)]:
                for v in range(begin+1, end):
                    u = int(rng.integers(begin, v))
                    edges.append((u, v) if rng.integers(2) else (v, u))
            flow = ForestFlow(n, edges)
            delta, y = rng.normal(size=n), rng.normal(size=len(edges))
            w = flow.apply(delta)
            divergence = np.zeros(n)
            for j, (u, v) in enumerate(edges):
                divergence[u] -= w[j]
                divergence[v] += w[j]
            balance = np.linalg.norm(divergence-flow.project(delta))
            adjoint = abs(y @ w-delta @ flow.transpose(y)) / max(1, abs(y @ w))
            oracle = float("nan")
            if n <= 128:
                b = np.zeros((n, len(edges)))
                for j, (u, v) in enumerate(edges):
                    b[u, j], b[v, j] = -1, 1
                p = np.linalg.pinv(b)
                oracle = np.linalg.norm(w-p@delta) + np.linalg.norm(flow.transpose(y)-p.T@y)
            passed = balance < 1e-9 and adjoint < 1e-10 and (n > 128 or oracle < 1e-9)
            rows.append(dict(nodes=n, components=components, balance_error=balance,
                             adjoint_error=adjoint, dense_oracle_error=oracle, passed=int(passed)))
            assert passed, rows[-1]
    with (out / "tree_operator.csv").open("w") as f:
        writer = csv.DictWriter(f, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)
    print(f"{len(rows)}/{len(rows)} tree flow and adjoint checks passed; largest graph {max(r['nodes'] for r in rows)} nodes")
