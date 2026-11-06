#pragma once

#include "xpbd.hpp"

#include <cmath>
#include <stdexcept>
#include <string>

namespace cxpbd {

inline void add_spring(Model& model, int a, int b, double compliance) {
    Spring spring;
    spring.a = a;
    spring.b = b;
    spring.rest_length = (model.particles[a].x - model.particles[b].x).norm();
    spring.compliance = compliance;
    model.springs.push_back(spring);
}

inline Model make_spring() {
    Model model;
    Particle anchor;
    anchor.x = Vec3::Zero();
    anchor.v = Vec3::Zero();
    anchor.pinned = true;
    Particle bob;
    bob.x = Vec3(1.25, 0.0, 0.0);
    bob.v = Vec3::Zero();
    model.particles = {anchor, bob};
    add_spring(model, 0, 1, 0.01);
    model.springs[0].rest_length = 1.0;
    return model;
}

inline Model make_pendulum() {
    Model model;
    model.gravity = Vec3(0.0, -9.81, 0.0);
    Particle anchor;
    anchor.x = Vec3::Zero();
    anchor.v = Vec3::Zero();
    anchor.pinned = true;
    Particle bob;
    bob.x = Vec3(std::sin(0.65), -std::cos(0.65), 0.0);
    bob.v = Vec3::Zero();
    model.particles = {anchor, bob};
    add_spring(model, 0, 1, 0.0);
    return model;
}

inline Model make_grid() {
    Model model;
    model.gravity = Vec3(0.0, -9.81, 0.0);
    constexpr int columns = 6;
    constexpr int rows = 6;
    constexpr double spacing = 0.25;
    constexpr double compliance = 1e-4;
    for (int row = 0; row < rows; ++row) {
        for (int column = 0; column < columns; ++column) {
            Particle particle;
            particle.x = Vec3((column - 0.5 * (columns - 1)) * spacing,
                              1.5 - row * spacing, 0.0);
            particle.v = Vec3::Zero();
            particle.mass = 0.1;
            particle.pinned = row == 0 && (column == 0 || column == columns - 1);
            model.particles.push_back(particle);
        }
    }
    for (int row = 0; row < rows; ++row) {
        for (int column = 0; column < columns; ++column) {
            const int i = row * columns + column;
            if (column + 1 < columns) add_spring(model, i, i + 1, compliance);
            if (row + 1 < rows) add_spring(model, i, i + columns, compliance);
            if (column + 1 < columns && row + 1 < rows) {
                add_spring(model, i, i + columns + 1, compliance);
                add_spring(model, i + 1, i + columns, compliance);
            }
        }
    }
    return model;
}

inline Model make_cloth() {
    Model model;
    model.gravity = Vec3(0.0, -9.81, 0.0);
    constexpr int columns = 9;
    constexpr int rows = 9;
    constexpr double spacing = 0.2;
    constexpr double compliance = 2e-4;
    for (int row = 0; row < rows; ++row) {
        for (int column = 0; column < columns; ++column) {
            Particle particle;
            particle.x = Vec3((column - 0.5 * (columns - 1)) * spacing,
                              1.5, row * spacing);
            particle.v = Vec3::Zero();
            particle.mass = 0.05;
            particle.pinned = row == 0 && (column == 0 || column == columns - 1);
            model.particles.push_back(particle);
        }
    }
    for (int row = 0; row < rows; ++row) {
        for (int column = 0; column < columns; ++column) {
            const int i = row * columns + column;
            if (column + 1 < columns) add_spring(model, i, i + 1, compliance);
            if (row + 1 < rows) add_spring(model, i, i + columns, compliance);
            if (column + 1 < columns && row + 1 < rows) {
                add_spring(model, i, i + columns + 1, compliance);
                add_spring(model, i + 1, i + columns, compliance);
            }
        }
    }
    return model;
}

inline Model make_scene(const std::string& name) {
    if (name == "spring") return make_spring();
    if (name == "pendulum") return make_pendulum();
    if (name == "grid") return make_grid();
    if (name == "cloth") return make_cloth();
    throw std::invalid_argument("unknown scene: " + name);
}

} // namespace cxpbd
