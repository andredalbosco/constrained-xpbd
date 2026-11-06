#include "scenes.hpp"

#include <algorithm>
#include <cmath>
#include <exception>
#include <iostream>
#include <stdexcept>
#include <string>

using namespace cxpbd;

namespace {

void require(bool condition, const std::string& message) {
    if (!condition) throw std::runtime_error(message);
}

void near(double actual, double expected, double tolerance, const std::string& message) {
    if (!std::isfinite(actual) || std::abs(actual - expected) > tolerance) {
        throw std::runtime_error(message + ": actual=" + std::to_string(actual)
                                 + ", expected=" + std::to_string(expected));
    }
}

void finite_state(const Model& model) {
    for (const Particle& particle : model.particles) {
        require(particle.x.allFinite() && particle.v.allFinite(), "non-finite particle state");
    }
    require(std::isfinite(mechanical_energy(model)), "non-finite mechanical energy");
}

void unchanged_pins(const Model& model, const Model& initial) {
    for (std::size_t i = 0; i < model.particles.size(); ++i) {
        if (!initial.particles[i].pinned) continue;
        near((model.particles[i].x - initial.particles[i].x).norm(), 0.0, 1e-14,
             "pinned position changed");
        near(model.particles[i].v.norm(), 0.0, 1e-14, "pinned velocity changed");
    }
}

Settings baseline_settings() {
    Settings settings;
    settings.mode = Mode::XPBD;
    settings.dt = 1.0 / 240.0;
    settings.iterations = 20;
    return settings;
}

void spring_single_step() {
    Model model = make_spring();
    const Model initial = model;
    Settings settings = baseline_settings();
    settings.dt = 0.03;
    settings.iterations = 1;
    model.particles[1].v.x() = -0.2;
    const double omega_squared = 1.0 / (model.particles[1].mass * model.springs[0].compliance);
    const double displacement = model.particles[1].x.x() - model.springs[0].rest_length;
    const double velocity = model.particles[1].v.x();
    const double denominator = 1.0 + settings.dt * settings.dt * omega_squared;
    const double expected_displacement = (displacement + settings.dt * velocity) / denominator;
    const double expected_velocity = (velocity - settings.dt * omega_squared * displacement)
                                     / denominator;
    step(model, settings);
    near(model.particles[1].x.x(), 1.0 + expected_displacement, 1e-13,
         "spring single-step position differs from implicit Euler");
    near(model.particles[1].v.x(), expected_velocity, 1e-12,
         "spring single-step velocity differs from implicit Euler");
    unchanged_pins(model, initial);
}

void spring_oscillation() {
    Model model = make_spring();
    const Settings settings = baseline_settings();
    double displacement = 0.25;
    double velocity = 0.0;
    double previous_energy = mechanical_energy(model);
    const double omega_squared = 1.0 / (model.particles[1].mass * model.springs[0].compliance);
    const double denominator = 1.0 + settings.dt * settings.dt * omega_squared;
    int crossings = 0;
    for (int frame = 0; frame < 240; ++frame) {
        const double next_displacement = (displacement + settings.dt * velocity) / denominator;
        const double next_velocity = (velocity - settings.dt * omega_squared * displacement)
                                     / denominator;
        if (displacement * next_displacement < 0.0) ++crossings;
        displacement = next_displacement;
        velocity = next_velocity;
        step(model, settings);
        near(model.particles[1].x.x(), 1.0 + displacement, 2e-10,
             "spring trajectory differs from implicit Euler recurrence");
        near(model.particles[1].v.x(), velocity, 2e-9,
             "spring velocity differs from implicit Euler recurrence");
        const double energy = mechanical_energy(model);
        near(energy, previous_energy / denominator, 2e-10,
             "spring energy does not match implicit Euler dissipation");
        previous_energy = energy;
    }
    require(crossings >= 3, "spring did not complete the expected oscillations");
    near(total_reservoir(model), 0.0, 1e-14, "base XPBD unexpectedly changed reservoirs");
}

void pendulum_length_and_motion() {
    Model model = make_pendulum();
    const Model initial = model;
    const Settings settings = baseline_settings();
    double greatest_motion = 0.0;
    for (int frame = 0; frame < 240; ++frame) {
        step(model, settings);
        finite_state(model);
        near((model.particles[1].x - model.particles[0].x).norm(), 1.0, 1e-12,
             "hard pendulum changed length");
        near(model.particles[1].x.z(), 0.0, 1e-14, "pendulum left its plane");
        greatest_motion = std::max(greatest_motion,
            (model.particles[1].x - initial.particles[1].x).norm());
    }
    unchanged_pins(model, initial);
    require(greatest_motion > 0.5, "pendulum did not swing");
}

void grid_motion_and_stretch() {
    Model model = make_grid();
    const Model initial = model;
    const Settings settings = baseline_settings();
    double greatest_motion = 0.0;
    double greatest_stretch = 0.0;
    for (int frame = 0; frame < 240; ++frame) {
        step(model, settings);
        finite_state(model);
        greatest_stretch = std::max(greatest_stretch, maximum_relative_stretch(model));
        for (std::size_t i = 0; i < model.particles.size(); ++i) {
            near(model.particles[i].x.z(), 0.0, 1e-14, "2D grid left its plane");
            greatest_motion = std::max(greatest_motion,
                (model.particles[i].x - initial.particles[i].x).norm());
        }
    }
    unchanged_pins(model, initial);
    require(greatest_motion > 1e-3, "grid remained stationary under gravity");
    require(greatest_stretch < 0.2, "grid spring stretch exceeded 20 percent");
}

void cloth_sag_and_stretch() {
    Model model = make_cloth();
    const Model initial = model;
    const Settings settings = baseline_settings();
    double greatest_sag = 0.0;
    double greatest_stretch = 0.0;
    for (int frame = 0; frame < 240; ++frame) {
        step(model, settings);
        finite_state(model);
        greatest_stretch = std::max(greatest_stretch, maximum_relative_stretch(model));
        for (std::size_t i = 0; i < model.particles.size(); ++i) {
            greatest_sag = std::max(greatest_sag,
                initial.particles[i].x.y() - model.particles[i].x.y());
        }
    }
    unchanged_pins(model, initial);
    require(greatest_sag > 0.2, "cloth did not sag under gravity");
    require(greatest_stretch < 0.25, "cloth spring stretch exceeded 25 percent");
}

Vec3 linear_momentum(const Model& model) {
    Vec3 momentum = Vec3::Zero();
    for (const Particle& particle : model.particles) momentum += particle.mass * particle.v;
    return momentum;
}

void free_spring_linear_momentum() {
    Model model = make_spring();
    model.particles[0].pinned = false;
    model.particles[0].mass = 2.0;
    model.particles[1].mass = 3.0;
    model.particles[0].v = Vec3(0.4, -0.1, 0.2);
    model.particles[1].v = Vec3(-0.2, 0.3, -0.1);
    const Vec3 initial_momentum = linear_momentum(model);
    const Settings settings = baseline_settings();
    for (int frame = 0; frame < 240; ++frame) {
        step(model, settings);
        finite_state(model);
        near((linear_momentum(model) - initial_momentum).norm(), 0.0, 1e-10,
             "free spring changed total linear momentum");
    }
}

void passive_trajectory_and_energy() {
    Settings xpbd_settings = baseline_settings();
    Settings passive_settings = xpbd_settings;
    passive_settings.mode = Mode::Passive;
    for (const char* name : {"spring", "pendulum", "grid", "cloth"}) {
        Model xpbd = make_scene(name);
        Model passive = xpbd;
        const double initial_energy = mechanical_energy(passive) + total_reservoir(passive);
        const double tolerance = 1e-9 * std::max(1.0, std::abs(initial_energy));
        for (int frame = 0; frame < 120; ++frame) {
            step(xpbd, xpbd_settings);
            step(passive, passive_settings);
            finite_state(passive);
            for (std::size_t i = 0; i < xpbd.particles.size(); ++i) {
                require((xpbd.particles[i].x.array() == passive.particles[i].x.array()).all(),
                        std::string(name) + ": passive mode changed the XPBD positions");
                require((xpbd.particles[i].v.array() == passive.particles[i].v.array()).all(),
                        std::string(name) + ": passive mode changed the XPBD velocities");
            }
            near(mechanical_energy(passive) + total_reservoir(passive), initial_energy,
                 tolerance, std::string(name) + ": passive augmented energy changed");
        }
    }
}

} // namespace

int main() {
    int failures = 0;
    const auto run = [&failures](const char* name, void (*test)()) {
        try {
            test();
            std::cout << "[pass] " << name << '\n';
        } catch (const std::exception& error) {
            ++failures;
            std::cerr << "[fail] " << name << ": " << error.what() << '\n';
        }
    };
    run("spring analytic single step", spring_single_step);
    run("spring implicit Euler oscillation and dissipation", spring_oscillation);
    run("hard pendulum length, pins, and motion", pendulum_length_and_motion);
    run("2D grid motion, plane, and stretch", grid_motion_and_stretch);
    run("3D cloth sag, pins, and stretch", cloth_sag_and_stretch);
    run("free spring linear momentum", free_spring_linear_momentum);
    run("passive exact trajectory parity and augmented energy", passive_trajectory_and_energy);
    return failures == 0 ? 0 : 1;
}
