#pragma once

#include <Eigen/Core>
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>

namespace cxpbd {

using Vec3 = Eigen::Vector3d;
enum class Mode { XPBD, Passive, Active };
enum class ActiveOrder { Paired, PairedLocal, Sequential };

struct Particle {
    Vec3 x = Vec3::Zero();
    Vec3 v = Vec3::Zero();
    double mass = 1.0;
    bool pinned = false;
    double reservoir = 0.0; // Signed R_i, carried between timesteps.
};

struct Spring {
    int a = 0, b = 1;
    double rest_length = 1.0;
    double compliance = 0.001; // Zero denotes a hard distance constraint.
    double lambda = 0.0;      // Current timestep-scaled material multiplier.
    double previous_lambda = 0.0;
    double reservoir = 0.0;   // Signed S_j, carried between timesteps.
};

struct Model {
    std::vector<Particle> particles;
    std::vector<Spring> springs;
    Vec3 gravity = Vec3::Zero();
    double previous_dt = 0.0;
};

struct Settings {
    double dt = 1.0 / 240.0;
    int iterations = 20;
    Mode mode = Mode::XPBD;
    ActiveOrder active_order = ActiveOrder::Paired;
    double kinetic_reservoir_compliance = 1000.0; // beta, in energy units; larger = weaker feedback.
    double elastic_reservoir_compliance = 1000.0; // gamma; conservative demo default, not a stability bound.
};

struct Diagnostics {
    double mechanical_energy = 0.0;
    double reservoir_energy = 0.0;
    double reservoir_norm = 0.0;
    double energy_balance_error = 0.0;
    double max_balance_residual = 0.0;
    double max_feedback_residual = 0.0;
    double max_material_residual = 0.0;
    double position_residual_norm = 0.0;
    int backtracks = 0;
    int rejected_updates = 0;
};

inline double inverse_mass(const Particle& p) {
    return p.pinned ? 0.0 : 1.0 / p.mass;
}

inline double spring_energy(const Model& model, const Spring& s) {
    if (s.compliance == 0.0) return 0.0;
    const double c = (model.particles[s.a].x - model.particles[s.b].x).norm()
                   - s.rest_length;
    return 0.5 * c * c / s.compliance;
}

inline double mechanical_energy(const Model& model) {
    double energy = 0.0;
    for (const auto& p : model.particles)
        energy += 0.5 * p.mass * p.v.squaredNorm() - p.mass * model.gravity.dot(p.x);
    for (const auto& s : model.springs) energy += spring_energy(model, s);
    return energy;
}

inline double total_reservoir(const Model& model) {
    double value = 0.0;
    for (const auto& p : model.particles) value += p.reservoir;
    for (const auto& s : model.springs) value += s.reservoir;
    return value;
}

inline double maximum_relative_stretch(const Model& model) {
    double value = 0.0;
    for (const auto& s : model.springs) {
        const double length = (model.particles[s.a].x - model.particles[s.b].x).norm();
        value = std::max(value, std::abs(length - s.rest_length) / s.rest_length);
    }
    return value;
}

// Validate inputs before advancing any state. Stationary pins have zero velocity.
inline void validate(const Model& model, const Settings& settings) {
    if (!(settings.dt > 0.0) || !std::isfinite(settings.dt) || settings.iterations < 1)
        throw std::invalid_argument("dt must be finite and positive; iterations must be positive");
    if (!model.gravity.allFinite() || !std::isfinite(model.previous_dt) || model.previous_dt < 0.0)
        throw std::invalid_argument("invalid gravity or force-history timestep");
    if (settings.mode == Mode::Active &&
        (!(settings.kinetic_reservoir_compliance > 0.0) ||
         !(settings.elastic_reservoir_compliance > 0.0) ||
         !std::isfinite(settings.kinetic_reservoir_compliance) ||
         !std::isfinite(settings.elastic_reservoir_compliance)))
        throw std::invalid_argument("active reservoir compliances must be finite and positive");
    for (const auto& p : model.particles) {
        if (!(p.mass > 0.0) || !std::isfinite(p.mass) || !p.x.allFinite() ||
            !p.v.allFinite() || !std::isfinite(p.reservoir) ||
            (p.pinned && p.v.squaredNorm() != 0.0))
            throw std::invalid_argument("invalid particle or moving pin");
    }
    for (const auto& s : model.springs) {
        if (s.a < 0 || s.b < 0 || s.a == s.b ||
            s.a >= static_cast<int>(model.particles.size()) ||
            s.b >= static_cast<int>(model.particles.size()) ||
            !(s.rest_length > 0.0) || !std::isfinite(s.rest_length) ||
            s.compliance < 0.0 || !std::isfinite(s.compliance) ||
            !std::isfinite(s.previous_lambda) || !std::isfinite(s.reservoir))
            throw std::invalid_argument("invalid spring");
    }
}

inline Diagnostics step(Model& model, const Settings& settings) {
    validate(model, settings);
    Diagnostics out;
    const double h = settings.dt, h2 = h * h;
    const double beta = settings.kinetic_reservoir_compliance;
    const double gamma = settings.elastic_reservoir_compliance;
    const double before = mechanical_energy(model) + total_reservoir(model);
    const std::size_t n = model.particles.size(), m = model.springs.size();
    std::vector<Vec3> old_x(n), old_v(n), predicted(n), particle_force(n, Vec3::Zero());
    std::vector<Vec3> old_force(m, Vec3::Zero()); // Force on spring endpoint a.
    std::vector<double> old_u(m), old_r(n), old_s(m);
    std::vector<double> q_mu(n, 0.0), q_nu(m, 0.0); // q = h^2 times physical multiplier.
    std::vector<bool> incident(n, false);

    for (std::size_t i = 0; i < n; ++i) {
        const auto& p = model.particles[i];
        old_x[i] = p.x; old_v[i] = p.v; old_r[i] = p.reservoir;
    }
    // Freeze the work estimate for the entire timestep, before predicting positions.
    for (std::size_t j = 0; j < m; ++j) {
        auto& s = model.springs[j];
        old_s[j] = s.reservoir;
        old_u[j] = spring_energy(model, s);
        const Vec3 edge = old_x[s.a] - old_x[s.b];
        const double length = edge.norm();
        if (length > 1e-14 * s.rest_length) {
            double force = 0.0;
            if (model.previous_dt > 0.0)
                force = s.previous_lambda / (model.previous_dt * model.previous_dt);
            else if (s.compliance > 0.0)
                force = -(length - s.rest_length) / s.compliance;
            old_force[j] = force * edge / length;
        }
        particle_force[s.a] += old_force[j];
        particle_force[s.b] -= old_force[j];
        incident[s.a] = incident[s.b] = true;
        s.lambda = 0.0;
    }
    for (std::size_t i = 0; i < n; ++i) {
        auto& p = model.particles[i];
        if (!p.pinned) p.x += h * p.v + h2 * model.gravity;
        predicted[i] = p.x;
    }

    auto kinetic_balance = [&](std::size_t i) {
        const auto& p = model.particles[i];
        const Vec3 dx = p.x - old_x[i];
        const double delta_t = 0.5 * p.mass * (dx.squaredNorm() / h2 - old_v[i].squaredNorm());
        return delta_t - (particle_force[i] + p.mass * model.gravity).dot(dx);
    };
    auto elastic_balance = [&](std::size_t j) {
        const auto& s = model.springs[j];
        const Vec3 dx = (model.particles[s.a].x - old_x[s.a])
                      - (model.particles[s.b].x - old_x[s.b]);
        return spring_energy(model, s) - old_u[j] + old_force[j].dot(dx);
    };

    // One particle-owned multiplier, using work from ALL its incident springs.
    auto solve_kinetic = [&](std::size_t i) {
        auto& p = model.particles[i];
        const double w = inverse_mass(p);
        const Vec3 gradient = p.mass * (p.x - old_x[i]) / h2
                            - particle_force[i] - p.mass * model.gravity;
        const double denominator = w * gradient.squaredNorm() + beta / h2;
        const double rhs = old_r[i] - kinetic_balance(i) - (beta / h2) * q_mu[i];
        const double dmu = rhs / denominator;
        const Vec3 x = p.x;
        const double mu = q_mu[i];
        bool accepted = false;
        double scale = 1.0;
        for (int attempt = 0; attempt < 16; ++attempt, scale *= 0.5) {
            p.x = x + scale * w * gradient * dmu;
            q_mu[i] = mu + scale * dmu;
            const double residual = kinetic_balance(i) + (beta / h2) * q_mu[i] - old_r[i];
            if (std::isfinite(residual) && std::abs(residual) <= std::abs(rhs) * (1.0 + 1e-12)) {
                accepted = true; break;
            }
            ++out.backtracks;
        }
        if (!accepted) {
            p.x = x; q_mu[i] = mu; ++out.rejected_updates;
        }
    };
    auto solve_endpoints = [&](const Spring& s) {
        solve_kinetic(s.a);
        solve_kinetic(s.b);
    };
    // The scalar potential-work row is relinearized AFTER the material update.
    auto solve_potential = [&](std::size_t j) {
        const auto& s = model.springs[j];
        auto& pa = model.particles[s.a];
        auto& pb = model.particles[s.b];
        const double wa = inverse_mass(pa), wb = inverse_mass(pb);
        const Vec3 edge = pa.x - pb.x;
        const double length = edge.norm();
        if (length <= 1e-14 * s.rest_length) { ++out.rejected_updates; return; }
        Vec3 gradient = old_force[j];
        if (s.compliance > 0.0)
            gradient += ((length - s.rest_length) / s.compliance) * edge / length;
        const double denominator = (wa + wb) * gradient.squaredNorm() + gamma / h2;
        const double rhs = old_s[j] - elastic_balance(j) - (gamma / h2) * q_nu[j];
        const double dnu = rhs / denominator;
        const Vec3 xa = pa.x, xb = pb.x;
        const double nu = q_nu[j];
        bool accepted = false;
        double scale = 1.0;
        for (int attempt = 0; attempt < 16; ++attempt, scale *= 0.5) {
            pa.x = xa + scale * wa * gradient * dnu;
            pb.x = xb - scale * wb * gradient * dnu;
            q_nu[j] = nu + scale * dnu;
            const double residual = elastic_balance(j) + (gamma / h2) * q_nu[j] - old_s[j];
            if (std::isfinite(residual) && std::abs(residual) <= std::abs(rhs) * (1.0 + 1e-12)) {
                accepted = true; break;
            }
            ++out.backtracks;
        }
        if (!accepted) {
            pa.x = xa; pb.x = xb; q_nu[j] = nu; ++out.rejected_updates;
        }
    };

    for (int iteration = 0; iteration < settings.iterations; ++iteration) {
        for (std::size_t j = 0; j < m; ++j) {
            auto& s = model.springs[j];
            auto& pa = model.particles[s.a];
            auto& pb = model.particles[s.b];
            const double wa = inverse_mass(pa), wb = inverse_mass(pb), w = wa + wb;
            const Vec3 edge = pa.x - pb.x;
            const double length = edge.norm();
            if (w == 0.0) {
                // Fixed geometry still has solvable compliant material/work multipliers.
                if (s.compliance > 0.0)
                    s.lambda = -h2 * (length - s.rest_length) / s.compliance;
                if (settings.mode == Mode::Active)
                    q_nu[j] = h2 * (old_s[j] - elastic_balance(j)) / gamma;
                if (settings.mode == Mode::Active && settings.active_order != ActiveOrder::Paired)
                    solve_endpoints(s);
                continue;
            }
            if (length <= 1e-14 * s.rest_length) {
                ++out.rejected_updates;
                if (settings.mode == Mode::Active && settings.active_order != ActiveOrder::Paired)
                    solve_endpoints(s);
                continue;
            }
            const Vec3 g = edge / length;
            const double c = length - s.rest_length;
            const double alpha_tilde = s.compliance / h2;
            const double a = w + alpha_tilde;
            const double rhs_c = -c - alpha_tilde * s.lambda;

            if (settings.mode != Mode::Active || settings.active_order == ActiveOrder::Sequential) {
                // Unmodified ordinary XPBD, also used by passive bookkeeping.
                const double dlambda = rhs_c / a;
                pa.x += wa * g * dlambda;
                pb.x -= wb * g * dlambda;
                s.lambda += dlambda;
                if (settings.mode == Mode::Active) {
                    solve_potential(j);
                    solve_endpoints(s);
                }
                continue;
            }

            // Pair the material row with d_j(x) + gamma/h^2 q_nu - S_old = 0.
            // Keep the analytical U(x) gradient: no off-manifold substitution of lambda.
            Vec3 z = old_force[j];
            if (s.compliance > 0.0) z += (c / s.compliance) * g;
            const double b = w * g.dot(z);
            const double d = w * z.squaredNorm() + gamma / h2;
            const double rhs_d = old_s[j] - elastic_balance(j) - (gamma / h2) * q_nu[j];
            // Positive-term Schur complement avoids cancellation in d - b*b/a.
            const double parallel = g.dot(z);
            const Vec3 perpendicular = z - parallel * g;
            const double schur = gamma / h2 + w * perpendicular.squaredNorm()
                               + (w * alpha_tilde / a) * parallel * parallel;
            const double dnu = (rhs_d - (b / a) * rhs_c) / schur;
            const double dlambda = (rhs_c - b * dnu) / a;
            const Vec3 correction = g * dlambda + z * dnu;
            const Vec3 xa = pa.x, xb = pb.x;
            const double lambda = s.lambda, nu = q_nu[j];
            const double merit = rhs_c * rhs_c / a + rhs_d * rhs_d / d;
            bool accepted = false;
            double scale = 1.0;
            for (int attempt = 0; attempt < 16; ++attempt, scale *= 0.5) {
                pa.x = xa + scale * wa * correction;
                pb.x = xb - scale * wb * correction;
                s.lambda = lambda + scale * dlambda;
                q_nu[j] = nu + scale * dnu;
                const double next_c = (pa.x - pb.x).norm() - s.rest_length + alpha_tilde * s.lambda;
                const double next_d = elastic_balance(j) + (gamma / h2) * q_nu[j] - old_s[j];
                const double next_merit = next_c * next_c / a + next_d * next_d / d;
                if (std::isfinite(next_merit) && next_merit <= merit * (1.0 + 1e-12)) {
                    accepted = true; break;
                }
                ++out.backtracks;
            }
            if (!accepted) {
                pa.x = xa; pb.x = xb; s.lambda = lambda; q_nu[j] = nu;
                ++out.rejected_updates;
            }
            if (settings.active_order == ActiveOrder::PairedLocal) solve_endpoints(s);
        }

        if (settings.mode == Mode::Active) {
            for (std::size_t i = 0; i < n; ++i) {
                if (settings.active_order == ActiveOrder::Paired || !incident[i]) solve_kinetic(i);
            }
        }
    }

    // Close the local work accounts from the ACCEPTED positions, even at finite iterations.
    // R=beta*mu and S=gamma*nu need not yet hold: report their residual separately.
    std::vector<Vec3> position_residual(n, Vec3::Zero());
    for (std::size_t i = 0; i < n; ++i) {
        auto& p = model.particles[i];
        p.v = p.pinned ? Vec3::Zero().eval() : ((p.x - old_x[i]) / h).eval();
        const double b = kinetic_balance(i);
        if (settings.mode != Mode::XPBD) p.reservoir = old_r[i] - b;
        out.max_balance_residual = std::max(out.max_balance_residual,
            std::abs(b + p.reservoir - old_r[i]));
        if (settings.mode == Mode::Active)
            out.max_feedback_residual = std::max(out.max_feedback_residual,
                std::abs(p.reservoir - (beta / h2) * q_mu[i]));
        const Vec3 gradient = p.mass * (p.x - old_x[i]) / h2
                            - particle_force[i] - p.mass * model.gravity;
        position_residual[i] = p.mass * (p.x - predicted[i]) - gradient * q_mu[i];
        out.reservoir_norm += p.reservoir * p.reservoir;
    }
    for (std::size_t j = 0; j < m; ++j) {
        auto& s = model.springs[j];
        const double d = elastic_balance(j);
        if (settings.mode != Mode::XPBD) s.reservoir = old_s[j] - d;
        out.max_balance_residual = std::max(out.max_balance_residual,
            std::abs(d + s.reservoir - old_s[j]));
        if (settings.mode == Mode::Active)
            out.max_feedback_residual = std::max(out.max_feedback_residual,
                std::abs(s.reservoir - (gamma / h2) * q_nu[j]));
        const Vec3 edge = model.particles[s.a].x - model.particles[s.b].x;
        const double length = edge.norm(), c = length - s.rest_length;
        out.max_material_residual = std::max(out.max_material_residual,
            std::abs(c + (s.compliance / h2) * s.lambda));
        if (length > 1e-14 * s.rest_length) {
            const Vec3 g = edge / length;
            Vec3 z = old_force[j];
            if (s.compliance > 0.0) z += (c / s.compliance) * g;
            const Vec3 contribution = g * s.lambda + z * q_nu[j];
            position_residual[s.a] -= contribution;
            position_residual[s.b] += contribution;
        }
        s.previous_lambda = s.lambda;
        out.reservoir_norm += s.reservoir * s.reservoir;
    }
    for (std::size_t i = 0; i < n; ++i) {
        const auto& p = model.particles[i];
        if (!p.x.allFinite() || !p.v.allFinite() || !std::isfinite(p.reservoir))
            throw std::runtime_error("non-finite simulation state");
        if (!p.pinned) out.position_residual_norm += position_residual[i].squaredNorm();
    }
    model.previous_dt = h;
    out.position_residual_norm = std::sqrt(out.position_residual_norm) / h2;
    out.reservoir_norm = std::sqrt(out.reservoir_norm);
    out.mechanical_energy = mechanical_energy(model);
    out.reservoir_energy = total_reservoir(model);
    out.energy_balance_error = out.mechanical_energy + out.reservoir_energy - before;
    if (!std::isfinite(out.mechanical_energy) || !std::isfinite(out.reservoir_energy))
        throw std::runtime_error("non-finite simulation energy");
    return out;
}

} // namespace cxpbd
