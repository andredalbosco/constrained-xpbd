#include "scenes.hpp"

#include <Eigen/Geometry>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

using namespace cxpbd;

namespace {

void require(bool condition, const std::string& message) {
    if (!condition) throw std::runtime_error(message);
}

double potential(const Model& model, const Spring& spring) {
    if (spring.compliance == 0.0) return 0.0;
    const double extension = (model.particles[spring.a].x - model.particles[spring.b].x).norm()
                           - spring.rest_length;
    return extension * extension / (2.0 * spring.compliance);
}

double energy(const Model& model) {
    double value = 0.0;
    for (const auto& particle : model.particles)
        value += particle.mass * (0.5 * particle.v.squaredNorm()
                                - model.gravity.dot(particle.x));
    for (const auto& spring : model.springs) value += potential(model, spring);
    return value;
}

double reservoirs(const Model& model) {
    double value = 0.0;
    for (const auto& particle : model.particles) value += particle.reservoir;
    for (const auto& spring : model.springs) value += spring.reservoir;
    return value;
}

struct Audit {
    double local_error = 0.0;
    double extended_error = 0.0;
};

// Reconstruct the rows independently from the accepted states and the old force history.
Audit audit(const Model& before, const Model& after, double h) {
    std::vector<Vec3> forces(before.particles.size(), Vec3::Zero());
    Audit result;
    for (std::size_t j = 0; j < before.springs.size(); ++j) {
        const auto& spring = before.springs[j];
        const Vec3 edge = before.particles[spring.a].x - before.particles[spring.b].x;
        const double length = edge.norm();
        Vec3 force = Vec3::Zero();
        if (length > 1e-14 * spring.rest_length) {
            const double magnitude = before.previous_dt > 0.0
                ? spring.previous_lambda / (before.previous_dt * before.previous_dt)
                : (spring.compliance > 0.0
                    ? -(length - spring.rest_length) / spring.compliance : 0.0);
            force = magnitude * edge / length;
        }
        forces[spring.a] += force;
        forces[spring.b] -= force;
        const Vec3 relative_step = (after.particles[spring.a].x - before.particles[spring.a].x)
                                - (after.particles[spring.b].x - before.particles[spring.b].x);
        const double d = potential(after, after.springs[j]) - potential(before, spring)
                       + force.dot(relative_step);
        result.local_error = std::max(result.local_error,
            std::abs(d + after.springs[j].reservoir - spring.reservoir));
    }
    for (std::size_t i = 0; i < before.particles.size(); ++i) {
        const auto& old = before.particles[i];
        const auto& next = after.particles[i];
        require(next.x.allFinite() && next.v.allFinite() && std::isfinite(next.reservoir),
                "non-finite particle state");
        const Vec3 dx = next.x - old.x;
        const Vec3 expected_velocity = old.pinned ? Vec3::Zero().eval() : (dx / h).eval();
        require((next.v - expected_velocity).norm() < 1e-11, "velocity reconstruction mismatch");
        if (old.pinned) require(dx.squaredNorm() == 0.0, "stationary pin moved");
        const double kinetic_change = 0.5 * old.mass * (next.v.squaredNorm() - old.v.squaredNorm());
        const double b = kinetic_change - forces[i].dot(dx) - old.mass * before.gravity.dot(dx);
        result.local_error = std::max(result.local_error,
            std::abs(b + next.reservoir - old.reservoir));
    }
    result.extended_error = std::abs(energy(after) + reservoirs(after)
                                  - energy(before) - reservoirs(before));
    const double scale = 1.0 + std::abs(energy(before)) + std::abs(reservoirs(before));
    require(result.local_error < 2e-11 * scale, "local accepted-state work account failed");
    require(result.extended_error < 2e-11 * scale, "extended energy account failed");
    return result;
}

void finite_diagnostics(const Diagnostics& d) {
    const double fields[] = {d.mechanical_energy, d.reservoir_energy, d.reservoir_norm,
        d.energy_balance_error, d.max_balance_residual, d.max_feedback_residual,
        d.max_material_residual, d.position_residual_norm};
    for (double field : fields) require(std::isfinite(field), "non-finite solver diagnostic");
}

Model free_pair(double seed) {
    Model model;
    Particle left, right;
    left.x = Vec3(-1.0, 0.0, 0.0); left.v = Vec3(-1.0, 0.0, 0.0);
    right.x = -left.x; right.v = -left.v;
    left.reservoir = right.reservoir = seed;
    model.particles = {left, right};
    return model;
}

void test_signed_transfer_and_carry() {
    Settings settings;
    settings.mode = Mode::Active;
    settings.kinetic_reservoir_compliance = 1.0;
    settings.elastic_reservoir_compliance = 1.0;
    for (double seed : {0.2, -0.2}) {
        Model model = free_pair(seed);
        const Model old = model;
        finite_diagnostics(step(model, settings));
        audit(old, model, settings.dt);
        const double energy_change = energy(model) - energy(old);
        require(seed * energy_change > 1e-5, "seeded reservoir did not transfer energy with expected sign");
        require(std::abs(reservoirs(model)) < std::abs(reservoirs(old)), "signed reservoir debt did not shrink");
        Vec3 momentum = Vec3::Zero(), angular = Vec3::Zero();
        for (const auto& p : model.particles) {
            momentum += p.mass * p.v;
            angular += p.x.cross(p.mass * p.v);
        }
        require(momentum.norm() < 1e-11 && angular.norm() < 1e-11,
                "symmetric free-particle transfer broke momentum symmetry");
        const Model carried = model;
        Model reset = model;
        for (auto& p : reset.particles) p.reservoir = 0.0;
        finite_diagnostics(step(model, settings));
        finite_diagnostics(step(reset, settings));
        audit(carried, model, settings.dt);
        require((model.particles[0].x - reset.particles[0].x).norm() > 1e-7,
                "carried reservoir had no effect on the next timestep");
        std::cout << "signed seed=" << seed << " first mechanical transfer=" << energy_change
                  << " remaining reservoir=" << reservoirs(model) << '\n';
    }
}

void test_variable_timestep() {
    Model model = make_spring();
    Settings settings;
    settings.mode = Mode::Active;
    settings.dt = 1.0 / 80.0;
    finite_diagnostics(step(model, settings));
    const Model old = model;
    settings.dt = 1.0 / 320.0;
    finite_diagnostics(step(model, settings));
    audit(old, model, settings.dt);
    const auto& spring = old.springs[0];
    const Vec3 normal = (old.particles[spring.a].x - old.particles[spring.b].x).normalized();
    const Vec3 relative_step = (model.particles[spring.a].x - old.particles[spring.a].x)
                            - (model.particles[spring.b].x - old.particles[spring.b].x);
    const double work_error_if_current_h_used = spring.previous_lambda * normal.dot(relative_step)
        * (1.0 / (settings.dt * settings.dt) - 1.0 / (old.previous_dt * old.previous_dt));
    require(std::abs(work_error_if_current_h_used) > 1e-4,
            "variable-timestep test did not distinguish old and current force scaling");
    require(model.previous_dt == settings.dt, "accepted timestep history was not updated");
    std::cout << "variable timestep: wrong force denominator would change work by "
              << work_error_if_current_h_used << '\n';
}

void test_finite_iteration_disclosure() {
    Model model = make_cloth();
    const Model old = model;
    Settings settings;
    settings.mode = Mode::Active;
    settings.dt = 1.0 / 60.0;
    settings.iterations = 1;
    const Diagnostics diagnostic = step(model, settings);
    finite_diagnostics(diagnostic);
    audit(old, model, settings.dt);
    require(diagnostic.max_feedback_residual > 1e-10,
            "one-sweep test did not expose the separate feedback residual");
    require(diagnostic.position_residual_norm > 1e-10,
            "one-sweep test did not expose the separate position residual");
    std::cout << "one sweep: ledger=" << diagnostic.max_balance_residual
              << " feedback=" << diagnostic.max_feedback_residual
              << " material=" << diagnostic.max_material_residual
              << " position=" << diagnostic.position_residual_norm << '\n';
}

void test_fixed_geometry_multipliers() {
    Model model = make_spring();
    model.particles[1].pinned = true;
    model.particles[0].reservoir = 0.1;
    model.springs[0].reservoir = -0.2;
    const Model old = model;
    Settings settings;
    settings.mode = Mode::Active;
    const Diagnostics diagnostic = step(model, settings);
    audit(old, model, settings.dt);
    require(diagnostic.max_material_residual < 1e-12,
            "fixed compliant spring did not update its material multiplier");
    require(diagnostic.max_feedback_residual < 1e-12,
            "fixed reservoirs did not update their feedback multipliers");
    require(model.springs[0].previous_lambda != 0.0,
            "fixed compliant spring lost its force history");
}

void test_scene(const std::string& name) {
    Model model = make_scene(name);
    const Model initial = model;
    const double initial_extended = energy(model) + reservoirs(model);
    Settings settings;
    settings.mode = Mode::Active;
    double max_motion = 0.0, max_stretch = 0.0, max_feedback = 0.0;
    double max_position = 0.0, max_ledger = 0.0, max_material = 0.0;
    int rejected = 0;
    for (int frame = 0; frame < 2400; ++frame) {
        const Model old = model;
        const Diagnostics diagnostic = step(model, settings);
        finite_diagnostics(diagnostic);
        const Audit check = audit(old, model, settings.dt);
        max_ledger = std::max(max_ledger, check.extended_error);
        max_feedback = std::max(max_feedback, diagnostic.max_feedback_residual);
        max_position = std::max(max_position, diagnostic.position_residual_norm);
        max_material = std::max(max_material, diagnostic.max_material_residual);
        rejected += diagnostic.rejected_updates;
        max_stretch = std::max(max_stretch, maximum_relative_stretch(model));
        for (std::size_t i = 0; i < model.particles.size(); ++i) {
            const auto& particle = model.particles[i];
            max_motion = std::max(max_motion, (particle.x - initial.particles[i].x).norm());
            require(particle.x.norm() < 20.0 && particle.v.norm() < 100.0,
                    name + ": state escaped the scene's metre/second scale");
        }
        if (name == "pendulum") {
            require(potential(model, model.springs[0]) == 0.0 && spring_energy(model, model.springs[0]) == 0.0,
                    "hard pendulum acquired artificial elastic potential");
        }
    }
    std::cout << name << ": motion=" << max_motion << " stretch=" << max_stretch
              << " mechanical delta=" << energy(model) - energy(initial)
              << " reservoir=" << reservoirs(model) << " ledger=" << max_ledger
              << " feedback=" << max_feedback << " material=" << max_material
              << " position=" << max_position << " rejected=" << rejected << '\n';
    require(max_motion > 1e-3, name + ": active scene did not move");
    const double stretch_bound = name == "pendulum" ? 1e-3 : 0.5;
    require(max_stretch < stretch_bound, name + ": spring deformation exceeded the scene bound");
    require(std::abs(energy(model) + reservoirs(model) - initial_extended)
                < 2e-9 * (1.0 + std::abs(initial_extended)),
            name + ": cumulative extended-energy drift");
}

} // namespace

int main() {
    try {
        std::cout << std::setprecision(12);
        test_signed_transfer_and_carry();
        test_variable_timestep();
        test_finite_iteration_disclosure();
        test_fixed_geometry_multipliers();
        for (const std::string name : {"spring", "pendulum", "grid", "cloth"}) test_scene(name);
        std::cout << "Active reservoir checks passed. Extended balance is bookkeeping; mechanical conservation is not asserted.\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "FAIL: " << error.what() << '\n';
        return 1;
    }
}
