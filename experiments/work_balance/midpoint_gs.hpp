#pragma once
// Midpoint discrete-gradient equations solved by curvature-dropped Newton + dual GS.
// No Hessians or dense linear solver are used. These equations are not generally
// the KKT equations of the original implicit-Euler minimization.
#include "scenes.hpp"
#include <Eigen/Geometry>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

using namespace cxpbd;

struct Row {
    Vec3 discrete_gradient, endpoint_gradient;
    double residual = 0.0, compliance = 0.0, scale = 1.0;
};
struct Residual {
    std::vector<Row> rows;
    std::vector<Vec3> x;
    double position = 0.0, material = 0.0;
};

inline Residual residual(const Model& old, const std::vector<Vec3>& q,
                  const std::vector<Vec3>& predictor, const std::vector<double>& lambda, double h) {
    Residual r;
    r.x.resize(q.size());
    for (std::size_t i = 0; i < q.size(); ++i)
        r.x[i] = old.particles[i].mass * (q[i] - predictor[i]);
    for (std::size_t j = 0; j < old.springs.size(); ++j) {
        const auto& s = old.springs[j];
        const Vec3 before = old.particles[s.a].x - old.particles[s.b].x;
        const Vec3 middle = q[s.a] - q[s.b];
        const Vec3 after = 2.0 * middle - before;
        Row row;
        row.compliance = 4.0 * s.compliance / (h * h);
        if (s.compliance > 0.0) {
            const double sum_length = before.norm() + after.norm();
            if (after.norm() < 1e-12 || sum_length < 1e-12)
                throw std::runtime_error("collapsed distance gradient");
            row.discrete_gradient = 2.0 * middle / sum_length;
            row.endpoint_gradient = after.normalized();
            row.residual = 0.5 * sum_length - s.rest_length + row.compliance * lambda[j];
        } else {
            // C=.5*(||r||^2-L^2), discrete gradient=(r_new+r_old)/2.
            row.discrete_gradient = middle;
            row.endpoint_gradient = after;
            row.residual = 0.25 * (before.squaredNorm() + after.squaredNorm()
                                  - 2.0 * s.rest_length * s.rest_length);
            row.scale = s.rest_length;
        }
        r.x[s.a] -= row.discrete_gradient * lambda[j];
        r.x[s.b] += row.discrete_gradient * lambda[j];
        if (!row.discrete_gradient.allFinite() || !row.endpoint_gradient.allFinite()
            || !std::isfinite(row.residual))
            throw std::runtime_error("nonfinite work-balance row");
        r.material += std::pow(row.residual / row.scale, 2);
        r.rows.push_back(row);
    }
    for (std::size_t i = 0; i < q.size(); ++i)
        if (!old.particles[i].pinned)
            r.position += (r.x[i] / old.particles[i].mass).squaredNorm();
    if (!std::isfinite(r.position) || !std::isfinite(r.material))
        throw std::runtime_error("nonfinite work-balance residual");
    return r;
}

struct StepReport {
    double rx = 0.0, constitutive = 0.0, kinetic_work_defect = 0.0, identity_error = 0.0;
    int iterations = 0, backtracks = 0;
};

inline StepReport advance(Model& model, double h, int outer_sweeps, int linear_sweeps,
                   bool keep_rx, bool symmetric_jacobian) {
    const Model old = model;
    std::vector<Vec3> q, predictor;
    for (const auto& p : old.particles) {
        Vec3 value = p.x;
        if (!p.pinned) value += 0.5 * h * p.v + 0.25 * h * h * old.gravity;
        q.push_back(value); predictor.push_back(value);
    }
    std::vector<double> lambda(old.springs.size(), 0.0);
    StepReport report;
    for (int outer = 0; outer < outer_sweeps; ++outer) {
        const Residual r = residual(old, q, predictor, lambda, h);
        const double merit = r.material + (keep_rx ? r.position : 0.0);
        if (merit < 1e-28) break;
        std::vector<Vec3> dq(q.size(), Vec3::Zero());
        std::vector<double> dlambda(lambda.size(), 0.0);
        if (keep_rx)
            for (std::size_t i = 0; i < q.size(); ++i)
                dq[i] = -inverse_mass(old.particles[i]) * r.x[i];
        // Solve the frozen multiplier system by Gauss-Seidel. Its rows are
        // J_end M^-1 J_discrete^T + 4 alpha/h^2, generally nonsymmetric.
        for (int inner = 0; inner < linear_sweeps; ++inner) {
            for (std::size_t j = 0; j < lambda.size(); ++j) {
                const auto& s = old.springs[j];
                const Row& row = r.rows[j];
                const Vec3 left = symmetric_jacobian ? row.discrete_gradient : row.endpoint_gradient;
                const double wa = inverse_mass(old.particles[s.a]);
                const double wb = inverse_mass(old.particles[s.b]);
                const double denominator = (wa + wb) * left.dot(row.discrete_gradient) + row.compliance;
                if (std::abs(denominator) < 1e-14) throw std::runtime_error("singular GS diagonal");
                const double increment = -(row.residual + left.dot(dq[s.a] - dq[s.b])
                                           + row.compliance * dlambda[j]) / denominator;
                dlambda[j] += increment;
                dq[s.a] += wa * row.discrete_gradient * increment;
                dq[s.b] -= wb * row.discrete_gradient * increment;
            }
        }
        bool accepted = false;
        double scale = 1.0;
        for (int attempt = 0; attempt < 20; ++attempt, scale *= 0.5) {
            std::vector<Vec3> next_q = q;
            std::vector<double> next_lambda = lambda;
            for (std::size_t i = 0; i < q.size(); ++i) next_q[i] += scale * dq[i];
            for (std::size_t j = 0; j < lambda.size(); ++j) next_lambda[j] += scale * dlambda[j];
            const Residual next = residual(old, next_q, predictor, next_lambda, h);
            const double next_merit = next.material + (keep_rx ? next.position : 0.0);
            if (std::isfinite(next_merit) &&
                (next_merit < 1e-28 || next_merit < merit * (1.0 - 1e-4 * scale))) {
                q = next_q; lambda = next_lambda; accepted = true; break;
            }
            ++report.backtracks;
        }
        ++report.iterations;
        if (!accepted) break; // Leave a visible residual, rather than claiming convergence.
    }
    const Residual final = residual(old, q, predictor, lambda, h);
    report.rx = std::sqrt(final.position);
    report.constitutive = std::sqrt(final.material);
    std::vector<Vec3> impulses(q.size(), Vec3::Zero());
    for (std::size_t i = 0; i < q.size(); ++i) {
        auto& p = model.particles[i];
        if (!p.pinned) {
            p.x = 2.0 * q[i] - old.particles[i].x;
            p.v = 2.0 * (p.x - old.particles[i].x) / h - old.particles[i].v;
            impulses[i] = h * p.mass * old.gravity;
        }
        if (!p.x.allFinite() || !p.v.allFinite())
            throw std::runtime_error("nonfinite midpoint state");
    }
    for (std::size_t j = 0; j < lambda.size(); ++j) {
        const auto& s = model.springs[j];
        const Vec3 impulse = (4.0 / h) * final.rows[j].discrete_gradient * lambda[j];
        impulses[s.a] += impulse; impulses[s.b] -= impulse;
    }
    double delta_t = 0.0, work = 0.0, expected_defect = 0.0;
    for (std::size_t i = 0; i < q.size(); ++i) {
        const auto& p = model.particles[i];
        const auto& before = old.particles[i];
        const Vec3 average = 0.5 * (p.v + before.v);
        delta_t += 0.5 * p.mass * (p.v.squaredNorm() - before.v.squaredNorm());
        work += average.dot(impulses[i]);
        expected_defect += (4.0 / h) * average.dot(final.x[i]);
    }
    report.kinetic_work_defect = delta_t - work;
    report.identity_error = std::abs(report.kinetic_work_defect - expected_defect);
    if (!std::isfinite(report.kinetic_work_defect) || !std::isfinite(report.identity_error))
        throw std::runtime_error("nonfinite midpoint work diagnostic");
    return report;
}
