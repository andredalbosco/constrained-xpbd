"""Independent numerical audits of the two block-GS spring solvers.

The independent reference below assembles the full linearized equations and
solves them with numpy.linalg.solve. It does not call the production increment
formula when constructing the reference solution.
"""
from __future__ import annotations

import argparse
from dataclasses import replace
import importlib
import json
from pathlib import Path

import numpy as np


def potential(x, mass, stiffness, length, gravity):
    return 0.5 * stiffness * (np.linalg.norm(x) - length) ** 2 - mass * gravity @ x


def independent(x, lam, rho, eta, xn, vn, h, mass, stiffness, length,
                gravity, epsilon, energy_target, kinetic_floor, rigid=False):
    """Return residual, approximation, exact derivative, and physical quantities."""
    dim = len(x)
    radius = np.linalg.norm(x)
    normal = x / radius
    constraint = radius - length
    normal_derivative = (np.eye(dim) - np.outer(normal, normal)) / radius
    compliance = 0.0 if rigid else 1.0 / stiffness
    potential_gradient = stiffness * constraint * normal - mass * gravity
    potential_hessian = stiffness * (
        np.outer(normal, normal) + constraint * normal_derivative)
    energy_gradient = mass * (x - xn) / h ** 2 + potential_gradient
    position_metric = mass * np.eye(dim) / h ** 2
    gap = energy_target - kinetic_floor
    energy = (0.5 * mass * np.dot(x - xn, x - xn) / h ** 2
              + potential(x, mass, stiffness, length, gravity))
    residual = np.r_[
        mass * (x - xn - h * vn) / h ** 2 - lam * normal - mass * gravity
        + eta * energy_gradient,
        constraint + compliance * lam,
        epsilon * rho + eta * gap,
        energy - (1 - rho) * energy_target - rho * kinetic_floor,
    ]
    matrix = np.zeros((dim + 3, dim + 3))
    matrix[:dim, :dim] = position_metric
    matrix[:dim, dim] = -normal
    matrix[dim, :dim] = normal
    matrix[dim, dim] = compliance
    matrix[:dim, dim + 2] = energy_gradient
    matrix[dim + 1, dim + 1] = epsilon
    matrix[dim + 1, dim + 2] = gap
    matrix[dim + 2, :dim] = energy_gradient
    matrix[dim + 2, dim + 1] = gap
    exact = matrix.copy()
    exact[:dim, :dim] += (
        -lam * normal_derivative + eta * (position_metric + potential_hessian))
    return dict(residual=residual, matrix=matrix, exact=exact,
                normal=normal, constraint=constraint,
                gradient=potential_gradient, hessian=potential_hessian,
                energy_gradient=energy_gradient, energy=energy,
                compliance=compliance, gap=gap)


def finite_difference(function, point):
    point = np.asarray(point, dtype=float)
    value = np.asarray(function(point))
    derivative = np.empty((value.size, point.size))
    for column in range(point.size):
        step = 2e-6 * max(1.0, abs(point[column]))
        plus, minus = point.copy(), point.copy()
        plus[column] += step
        minus[column] -= step
        derivative[:, column] = (function(plus) - function(minus)) / (2 * step)
    return derivative


def relative_error(actual, expected):
    return float(np.linalg.norm(actual - expected) / max(1.0, np.linalg.norm(expected)))


class Checks:
    def __init__(self):
        self.count = 0
        self.maxima = {}

    def close(self, name, actual, expected, tolerance=2e-10):
        self.count += 1
        error = relative_error(np.asarray(actual), np.asarray(expected))
        self.maxima[name] = max(error, self.maxima.get(name, 0.0))
        if not np.isfinite(error) or error > tolerance:
            raise AssertionError(f"{name}: relative error {error:.3e} > {tolerance:.3e}")

    def true(self, name, condition):
        self.count += 1
        if not condition:
            raise AssertionError(name)


def independent_math_checks(checks):
    rng = np.random.default_rng(20261101)
    max_contraction = 0.0
    for sample in range(60):
        mass = 10 ** rng.uniform(-0.5, 0.5)
        stiffness = 10 ** rng.uniform(1.0, 3.0)
        length = rng.uniform(0.5, 1.5)
        h = rng.uniform(0.02, 0.2)
        epsilon = 10 ** rng.uniform(-1, 2)
        angle = rng.uniform(-np.pi, np.pi)
        x = rng.uniform(0.7, 1.5) * length * np.array([np.cos(angle), np.sin(angle)])
        xn = x + rng.normal(size=2) * 0.05
        vn = rng.normal(size=2)
        gravity = np.array([0.0, -1.7]) if sample % 2 else np.zeros(2)
        energy_target = rng.uniform(1, 20)
        kinetic_floor = rng.uniform(-1, 0)
        lam, rho, eta = rng.normal(size=3) * np.array([10.0, 0.2, 0.3])
        args = (xn, vn, h, mass, stiffness, length, gravity, epsilon,
                energy_target, kinetic_floor)
        state = independent(x, lam, rho, eta, *args)
        residual, matrix = state["residual"], state["matrix"]
        point = np.r_[x, lam, rho, eta]
        derivative = finite_difference(
            lambda q: independent(q[:2], q[2], q[3], q[4], *args)["residual"], point)
        checks.close("full residual finite differences", derivative, state["exact"], 3e-8)
        energy_derivative = finite_difference(
            lambda q: np.array([independent(q, lam, rho, eta, *args)["energy"]]), x)
        checks.close("physical energy finite differences", energy_derivative[0],
                     state["energy_gradient"], 3e-8)
        potential_derivative = finite_difference(
            lambda q: independent(q, lam, rho, eta, *args)["gradient"], x)
        checks.close("physical potential Hessian finite differences", potential_derivative,
                     state["hessian"], 3e-8)

        # Elimination must recover the earlier, undualized approximate system,
        # including its right side when lambda is not constitutively feasible.
        undual = matrix[np.ix_([0, 1, 3, 4], [0, 1, 3, 4])].copy()
        undual[:2, :2] += stiffness * np.outer(state["normal"], state["normal"])
        undual_residual = residual[[0, 1, 3, 4]].copy()
        undual_residual[:2] += stiffness * state["normal"] * residual[2]
        full_delta = np.linalg.solve(matrix, -residual)
        undual_delta = np.linalg.solve(undual, -undual_residual)
        checks.close("dualization elimination including r_lambda", full_delta[[0, 1, 3, 4]],
                     undual_delta, 3e-9)

        normal, energy_gradient = state["normal"], state["energy_gradient"]
        inverse_mass_metric = h * h / mass
        material_diagonal = 1 / stiffness + inverse_mass_metric * (normal @ normal)
        energy_diagonal = (inverse_mass_metric * (energy_gradient @ energy_gradient)
                           + state["gap"] ** 2 / epsilon)
        cross = inverse_mass_metric * (normal @ energy_gradient)
        multiplier_matrix = np.array([[material_diagonal, -cross],
                                      [-cross, energy_diagonal]])
        multiplier_rhs = np.array([
            -residual[2] + inverse_mass_metric * normal @ residual[:2],
            residual[4] - inverse_mass_metric * energy_gradient @ residual[:2]
            - state["gap"] * residual[3] / epsilon,
        ])
        multiplier_delta = np.linalg.solve(multiplier_matrix, multiplier_rhs)
        checks.close("independent multiplier Schur solve", multiplier_delta,
                     full_delta[[2, 4]], 3e-9)
        contraction = cross ** 2 / (material_diagonal * energy_diagonal)
        max_contraction = max(max_contraction, contraction)
        checks.true("frozen GS contraction below one", 0 <= contraction < 1)
        for row_order in ((0, 1), (1, 0)):
            propagator = np.eye(2)
            for row in row_order:
                coordinate_update = np.eye(2)
                coordinate_update[row, :] -= (multiplier_matrix[row, :]
                                              / multiplier_matrix[row, row])
                propagator = coordinate_update @ propagator
            checks.close("frozen GS exact spectral radius",
                         np.sort(np.abs(np.linalg.eigvals(propagator))),
                         [0.0, contraction])

        # Check both frozen schedules against a dense joint solve. These tests
        # deliberately keep the linearization frozen while updating its residual.
        for order in ("EL", "LE"):
            total = np.zeros(5)
            for sweep in range(5000):
                for block in order:
                    rows = [0, 1, 2] if block == "L" else [0, 1, 3, 4]
                    live_residual = residual + matrix @ total
                    increment = np.linalg.solve(matrix[np.ix_(rows, rows)],
                                                -live_residual[rows])
                    total[rows] += increment
                    checks.close("frozen block equations after immediate push",
                                 (residual + matrix @ total)[rows], np.zeros(len(rows)),
                                 2e-8)
                if relative_error(total, full_delta) < 3e-9:
                    break
            checks.close(f"frozen {order} versus joint solution", total, full_delta, 3e-8)
    return {"random_states": 60, "max_frozen_GS_contraction": max_contraction}


def production_checks(checks):
    """Production adapter is filled once solver.py's public API is available."""
    module = importlib.import_module("solver")
    return validate_production_api(checks, module)


def validate_production_api(checks, module):
    rng = np.random.default_rng(20261102)
    for sample in range(80):
        mass = 10 ** rng.uniform(-0.5, 0.5)
        rigid = sample % 3 == 0
        stiffness = 0.0 if rigid else 10 ** rng.uniform(1.0, 3.0)
        length = rng.uniform(0.5, 1.5)
        h = rng.uniform(0.02, 0.15)
        epsilon = 10 ** rng.uniform(-2, 8)
        angle = rng.uniform(-np.pi, np.pi)
        x = rng.uniform(0.7, 1.5) * length * np.array([np.cos(angle), np.sin(angle)])
        xn = x + rng.normal(size=2) * 0.05
        vn = rng.normal(size=2)
        gravity = np.array([0.0, -3.0]) if sample % 2 else np.zeros(2)
        energy_target = rng.uniform(1.0, 20.0)
        kinetic_floor = rng.uniform(-2.0, 0.0)
        lam, rho, eta = rng.normal(size=3) * np.array([10.0, 0.2, 0.3])
        scene = module.Scene("validation", mass, stiffness, length, gravity, xn, vn, rigid)
        args = (xn, vn, h, mass, stiffness, length, gravity, epsilon,
                energy_target, kinetic_floor)
        reference = independent(x, lam, rho, eta, *args, rigid=rigid)
        evaluation = module.evaluate_state(scene, xn, vn, h, x, lam, rho, eta,
                                           epsilon, target=energy_target,
                                           kstar=kinetic_floor)
        implementation_residual = np.r_[evaluation.r_x, evaluation.r_lambda,
                                        evaluation.r_rho, evaluation.r_eta]
        checks.close("production analytical residuals", implementation_residual,
                     reference["residual"])
        checks.close("production physical energy", evaluation.physical_energy,
                     reference["energy"])
        checks.close("production energy gradient", evaluation.energy_gradient,
                     reference["energy_gradient"])
        checks.close("production potential gradient", evaluation.potential_gradient,
                     reference["gradient"])
        checks.close("production normal", evaluation.normal, reference["normal"])
        joint = module.joint_increment(scene, h, epsilon, evaluation)
        checks.close("production joint increment versus independent full matrix",
                     np.r_[joint.dx, joint.dlam, joint.drho, joint.deta],
                     np.linalg.solve(reference["matrix"], -reference["residual"]), 3e-9)
        for block in "LE":
            rows = [0, 1, 2] if block == "L" else [0, 1, 3, 4]
            expected = np.zeros(5)
            expected[rows] = np.linalg.solve(reference["matrix"][np.ix_(rows, rows)],
                                            -reference["residual"][rows])
            increment = (module.material_increment(scene, h, evaluation) if block == "L"
                         else module.energy_increment(scene, h, epsilon, evaluation))
            actual = np.r_[increment.dx, increment.dlam, increment.drho, increment.deta]
            checks.close(f"production {block} increment versus independent dense block",
                         actual, expected, 3e-9)
            checks.close(f"production {block} immediate included-row residual",
                         (reference["residual"] + reference["matrix"] @ actual)[rows],
                         np.zeros(len(rows)), 3e-8)
            if block == "L":
                checks.close("material holds rho and eta", actual[[3, 4]], [0.0, 0.0])
            else:
                checks.close("energy holds lambda", actual[2], 0.0)
                checks.close("energy restores rho-eta equation", reference["residual"][3]
                             + epsilon * increment.drho
                             + reference["gap"] * increment.deta, 0.0, 3e-8)
    replayed_visits = 0
    for scene in module.scenes():
        force = scene.m * float(np.linalg.norm(scene.g))
        expected_floor = (-force * scene.l if scene.rigid
                          else -force * scene.l - force ** 2 / (2 * scene.k))
        checks.close("gravity-adapted reference floor formula", module.energy_floor(scene),
                     expected_floor)
        if force > 0:
            minimum_radius = scene.l if scene.rigid else scene.l + force / scene.k
            minimum_position = minimum_radius * scene.g / np.linalg.norm(scene.g)
        else:
            minimum_position = np.array([scene.l, 0.0])
        checks.close("reference floor attained at physical potential minimum",
                     module.energy(scene, minimum_position, np.zeros(2)), expected_floor)
        for epsilon in (1.0, 100.0, 1e4):
            for order in ("EL", "LE"):
                for sweeps in (1, 2, 4):
                    h = 1.0 / 120
                    actual = module.block_step(scene, scene.x0, scene.v0, h, epsilon,
                                               order, sweeps=sweeps, trace=True)
                    checks.true("replayed fixed-budget step accepted", actual.accepted)
                    checks.true("every visit retained in trace", len(actual.history) == 2*sweeps)
                    checks.true("visit count matches budget", actual.visits == 2*sweeps)
                    # Cold start only once per physical timestep. Thereafter
                    # all three auxiliary quantities accumulate independently.
                    point = np.r_[scene.x0+h*scene.v0, 0.0, 0.0, 0.0]
                    target = module.energy(scene, scene.x0, scene.v0)
                    floor = module.energy_floor(scene)
                    args = (scene.x0, scene.v0, h, scene.m, scene.k, scene.l,
                            scene.g, epsilon, target, floor)
                    visit = 0
                    for sweep in range(sweeps):
                        for block in order:
                            reference = independent(point[:2], point[2], point[3], point[4],
                                                    *args, rigid=scene.rigid)
                            rows = [0, 1, 2] if block == "L" else [0, 1, 3, 4]
                            change = np.linalg.solve(reference["matrix"][np.ix_(rows, rows)],
                                                     -reference["residual"][rows])
                            point[rows] += change
                            snapshot = actual.history[visit]
                            checks.true("correct sequential block order", snapshot["block"] == block)
                            checks.true("correct sequential sweep number", snapshot["sweep"] == sweep)
                            checks.close("refreshed dense-block replay versus trace",
                                         np.r_[snapshot["x"], snapshot["lam"],
                                               snapshot["rho"], snapshot["eta"]], point, 3e-8)
                            refreshed = independent(point[:2], point[2], point[3], point[4],
                                                    *args, rigid=scene.rigid)
                            checks.close("post-visit refreshed residual trace",
                                         np.r_[snapshot["r_x"], snapshot["r_lambda"],
                                               snapshot["r_rho"], snapshot["r_eta"]],
                                         refreshed["residual"], 3e-8)
                            checks.close("post-visit physical energy trace",
                                         snapshot["physical_energy"], refreshed["energy"], 3e-8)
                            visit += 1
                            replayed_visits += 1
                    checks.close("final state versus dense refreshed block replay",
                                 np.r_[actual.x, actual.lam, actual.rho, actual.eta], point, 3e-8)
                    checks.close("velocity derived from displacement",
                                 actual.v, (actual.x-scene.x0)/h)
                    if not actual.converged:
                        checks.true("fixed budget residual explicitly reported",
                                    actual.reason == "fixed_budget_residual")

    equilibrium = module.Scene("zero-energy equilibrium", 1.0, 100.0, 1.0,
                               np.zeros(2), np.array([1.0, 0.0]), np.zeros(2))
    evaluation = module.evaluate_state(equilibrium, equilibrium.x0, equilibrium.v0,
                                       1/60, equilibrium.x0, 0.0, 0.0, 0.0, 100.0)
    increment = module.energy_increment(equilibrium, 1/60, 100.0, evaluation)
    checks.true("compatible zero energy row explicitly flagged", increment.compatible_rank_defect)
    checks.close("compatible zero energy row returns zero increment",
                 np.r_[increment.dx, increment.dlam, increment.drho, increment.deta], np.zeros(5))
    try:
        module.energy_increment(equilibrium, 1/60, 100.0, replace(evaluation, r_eta=1.0))
    except FloatingPointError as error:
        checks.true("incompatible zero energy row explicit failure",
                    str(error) == "incompatible_zero_energy_row")
    else:
        raise AssertionError("incompatible zero energy row was silently accepted")
    for order in ("EL", "LE"):
        result = module.block_step(equilibrium, equilibrium.x0, equilibrium.v0,
                                   1/60, 100.0, order, sweeps=4, trace=True)
        checks.true("equilibrium remains accepted and converged", result.accepted and result.converged)
        checks.close("equilibrium stays at rest", np.r_[result.x, result.v],
                     np.r_[equilibrium.x0, equilibrium.v0])
        checks.true("no invented equilibrium energy multiplier", result.eta == 0.0)
        checks.true("fixed budget does not stop early", result.visits == 8)
    try:
        module.block_step(equilibrium, equilibrium.x0, equilibrium.v0,
                          1/60, 100.0, sweeps=0)
    except ValueError:
        checks.true("zero sweep budget rejected explicitly", True)
    else:
        raise AssertionError("zero sweep budget silently accepted")
    return {"production_arbitrary_states": 80, "replayed_nonlinear_visits": replayed_visits}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--independent-only", action="store_true")
    parser.add_argument("--write-report", action="store_true")
    args = parser.parse_args()
    checks = Checks()
    details = independent_math_checks(checks)
    if not args.independent_only:
        details.update(production_checks(checks))
    result = dict(status="passed", checks=checks.count, details=details,
                  largest_relative_errors=checks.maxima)
    print(json.dumps(result, indent=2, sort_keys=True))
    if args.write_report:
        report = Path(__file__).with_name("VALIDATION.md")
        report.write_text("""# Independent validation

The audit uses independently assembled dense linear systems rather than a
second copy of the production closed-form increment formulas.

## Coverage

- Sixty deterministic random spring states: physical energy and potential
  derivatives, the exact residual Jacobian, elimination of the carried material
  multiplier including nonzero `r_lambda`, and frozen EL/LE versus a full joint
  solve.
- Eighty arbitrary spring/rod states, with nonzero accumulated multipliers and
  relaxation variables: production residuals, both full blocks, and the joint
  increment versus independent dense solves. Gravity contributes to the
  objective outside the energy multiplier as well as to physical energy.
- Full refreshed replay of 504 visits across all six scenes, both orders,
  epsilon 1/100/10000, and one/two/four sweeps. Every post-visit state and residual
  is compared, checking multiplier accumulation, immediate position response,
  refresh timing, block ordering, and derived velocity.
- Gravity-adapted energy floors are checked against independently derived
  potential minima. Compatible zero energy rows are flagged; incompatible rows
  fail explicitly. Equilibrium stays at rest in both orders. A zero-sweep budget
  is explicitly rejected by the public API.

## Interpretation

For a frozen spring linearization the reduced two-multiplier system has positive
diagonals and the GS factor is the squared off-diagonal divided by their product.
The audited random states all have factor below one; both orders converge to the
same dense solution. This does not imply convergence of the refreshed nonlinear
iteration or stability of the physical rollout. Finite sweep budgets can retain
nonzero residuals, and the implementation reports that separately from a finite
accepted step.

The checks validate the stated adapted equations and their implementation. They
do not establish equivalence to the original CPD authors' software, optimal
motion selection, or reliable long-time stability. The hard rod is an extension
with alpha zero, not a claim about the original CPD spring model.

Run from this directory:

```bash
OPENBLAS_NUM_THREADS=1 /tmp/xpbd-nonlinear-venv/bin/python validate.py --write-report
```

## Numerical receipt

```json
"""
                          + json.dumps(result, indent=2, sort_keys=True)
                          + "\n```\n")


if __name__ == "__main__":
    main()
