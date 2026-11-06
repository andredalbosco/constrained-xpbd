#include "scenes.hpp"
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

using namespace cxpbd;

struct WorkAudit {
    double kinetic_change = 0.0, potential_change = 0.0, gravity_change = 0.0;
    double frozen_internal = 0.0, frozen_external = 0.0;
    double midpoint_internal = 0.0, midpoint_external = 0.0;
    double sequential_internal = 0.0, sequential_external = 0.0;
    double material_internal = 0.0;
    double impulse_defect = 0.0;
};

double kinetic(const Model& m) {
    double value = 0.0;
    for (const auto& p : m.particles) value += 0.5 * p.mass * p.v.squaredNorm();
    return value;
}

// Instrument the baseline update, and compare its positions with the production step.
WorkAudit audit_step(Model& model, double h, int sweeps) {
    const Model old = model;
    WorkAudit result;
    std::vector<Vec3> impulse_a(model.springs.size(), Vec3::Zero());
    std::vector<Vec3> material_impulses(model.particles.size(), Vec3::Zero());
    std::vector<Vec3> actual_impulses(model.particles.size(), Vec3::Zero());
    for (std::size_t i = 0; i < model.particles.size(); ++i) {
        auto& p = model.particles[i];
        if (p.pinned) continue;
        p.x += h * p.v + h * h * model.gravity;
        const Vec3 velocity = (p.x - old.particles[i].x) / h;
        const Vec3 impulse = h * p.mass * model.gravity;
        result.sequential_external += 0.5 * (p.v + velocity).dot(impulse);
        actual_impulses[i] += impulse;
        material_impulses[i] += impulse;
    }
    for (auto& spring : model.springs) spring.lambda = 0.0;
    for (int sweep = 0; sweep < sweeps; ++sweep) {
        for (std::size_t j = 0; j < model.springs.size(); ++j) {
            auto& spring = model.springs[j];
            auto& a = model.particles[spring.a];
            auto& b = model.particles[spring.b];
            const double wa = inverse_mass(a), wb = inverse_mass(b);
            if (wa + wb == 0.0) continue;
            const Vec3 edge = a.x - b.x;
            const double length = edge.norm();
            const Vec3 gradient = edge / length;
            const double alpha = spring.compliance / (h * h);
            const double dlambda = -(length - spring.rest_length + alpha * spring.lambda)
                                 / (wa + wb + alpha);
            const Vec3 va = (a.x - old.particles[spring.a].x) / h;
            const Vec3 vb = (b.x - old.particles[spring.b].x) / h;
            a.x += wa * gradient * dlambda;
            b.x -= wb * gradient * dlambda;
            spring.lambda += dlambda;
            const Vec3 next_va = (a.x - old.particles[spring.a].x) / h;
            const Vec3 next_vb = (b.x - old.particles[spring.b].x) / h;
            const Vec3 impulse = gradient * (dlambda / h);
            result.sequential_internal += 0.5 * (va + next_va - vb - next_vb).dot(impulse);
            impulse_a[j] += impulse;
        }
    }
    for (std::size_t i = 0; i < model.particles.size(); ++i) {
        auto& p = model.particles[i];
        const auto& before = old.particles[i];
        p.v = (p.x - before.x) / h;
        const Vec3 delta = p.x - before.x;
        const Vec3 mean_velocity = 0.5 * (p.v + before.v);
        result.midpoint_external += h * p.mass * model.gravity.dot(mean_velocity);
        result.frozen_external += p.mass * model.gravity.dot(delta);
        result.gravity_change -= p.mass * model.gravity.dot(delta);
    }
    result.kinetic_change = kinetic(model) - kinetic(old);
    for (std::size_t j = 0; j < model.springs.size(); ++j) {
        auto& s = model.springs[j];
        const auto& a = model.particles[s.a];
        const auto& b = model.particles[s.b];
        const auto& old_a = old.particles[s.a];
        const auto& old_b = old.particles[s.b];
        const Vec3 old_edge = old_a.x - old_b.x;
        const Vec3 delta = (a.x - old_a.x) - (b.x - old_b.x);
        const Vec3 mean_velocity = 0.5 * (a.v + old_a.v - b.v - old_b.v);
        const double old_force = old.previous_dt > 0.0
            ? old.springs[j].previous_lambda / (old.previous_dt * old.previous_dt)
            : (s.compliance > 0.0 ? -(old_edge.norm() - s.rest_length) / s.compliance : 0.0);
        result.frozen_internal += (old_force * old_edge.normalized()).dot(delta);
        result.midpoint_internal += impulse_a[j].dot(mean_velocity);
        const Vec3 final_impulse = (s.lambda / h) * (a.x - b.x).normalized();
        result.material_internal += final_impulse.dot(mean_velocity);
        actual_impulses[s.a] += impulse_a[j];
        actual_impulses[s.b] -= impulse_a[j];
        material_impulses[s.a] += final_impulse;
        material_impulses[s.b] -= final_impulse;
        result.potential_change += spring_energy(model, s) - spring_energy(old, old.springs[j]);
        s.previous_lambda = s.lambda;
    }
    model.previous_dt = h;
    for (std::size_t i = 0; i < model.particles.size(); ++i) {
        const auto& p = model.particles[i];
        if (p.pinned) continue;
        const Vec3 mean = 0.5 * (p.v + old.particles[i].v);
        result.impulse_defect += mean.dot(p.mass * (p.v - old.particles[i].v) - material_impulses[i]);
        const double exact_defect = (p.mass * (p.v - old.particles[i].v) - actual_impulses[i]).norm();
        if (exact_defect > 1e-8) throw std::runtime_error("incremental impulse audit failed");
    }
    Model reference = old;
    Settings settings;
    settings.dt = h; settings.iterations = sweeps;
    step(reference, settings);
    for (std::size_t i = 0; i < model.particles.size(); ++i)
        if ((model.particles[i].x - reference.particles[i].x).norm() > 1e-13)
            throw std::runtime_error("instrumented update differs from baseline");
    return result;
}

Model two_springs() {
    Model m;
    Particle a, b, c;
    a.x = Vec3(-1.0, 0.0, 0.0); a.pinned = true;
    b.x = Vec3(0.8, 0.5, 0.0); b.pinned = true;
    c.x = Vec3(0.0, 0.8, 0.0); c.v = Vec3(0.7, -0.2, 0.0);
    m.particles = {a, b, c};
    add_spring(m, 0, 2, 0.001);
    add_spring(m, 1, 2, 0.001);
    m.springs[0].rest_length *= 0.9;
    m.springs[1].rest_length *= 0.9;
    return m;
}

int main() {
    try {
        std::cout << std::setprecision(15)
            << "scene,h,sweeps,mechanical_drift,max_frozen_kinetic_defect,max_exact_kinetic_defect,"
               "max_sequential_kinetic_defect,max_material_kinetic_defect,max_impulse_identity_error,"
               "max_elastic_defect_with_kinetic_work,max_gravity_defect_with_kinetic_work,"
               "max_internal_allocation_difference\n";
        for (const std::string scene : {"spring", "pendulum", "two_springs"}) {
            for (double h : {1.0 / 240.0, 1.0 / 60.0, 1.0 / 15.0}) {
                for (int sweeps : {1, 20}) {
                    Model m = scene == "two_springs" ? two_springs() : make_scene(scene);
                    const double initial = mechanical_energy(m);
                    double frozen = 0, exact = 0, sequential = 0, material = 0;
                    double residual_identity = 0, elastic = 0, gravity = 0, allocation = 0;
                    for (int t = 0; t < static_cast<int>(std::lround(10.0 / h)); ++t) {
                        const WorkAudit a = audit_step(m, h, sweeps);
                        frozen = std::max(frozen, std::abs(a.kinetic_change-a.frozen_internal-a.frozen_external));
                        exact = std::max(exact, std::abs(a.kinetic_change-a.midpoint_internal-a.midpoint_external));
                        sequential = std::max(sequential, std::abs(a.kinetic_change-a.sequential_internal-a.sequential_external));
                        material = std::max(material, std::abs(a.kinetic_change-a.material_internal-a.midpoint_external));
                        residual_identity = std::max(residual_identity,
                            std::abs(a.kinetic_change-a.material_internal-a.midpoint_external-a.impulse_defect));
                        elastic = std::max(elastic, std::abs(a.potential_change+a.midpoint_internal));
                        gravity = std::max(gravity, std::abs(a.gravity_change+a.midpoint_external));
                        allocation = std::max(allocation, std::abs(a.midpoint_internal-a.sequential_internal));
                    }
                    if (exact > 1e-9 || sequential > 1e-9 || residual_identity > 1e-9)
                        throw std::runtime_error("kinetic work identity did not close");
                    std::cout << scene << ',' << h << ',' << sweeps << ',' << mechanical_energy(m)-initial
                        << ',' << frozen << ',' << exact << ',' << sequential << ',' << material
                        << ',' << residual_identity << ',' << elastic << ',' << gravity << ',' << allocation << '\n';
                }
            }
        }
        // Same net impulse and kinetic change, different work attribution by visit order.
        const Vec3 a(1.0,0.0,0.0), b(1.0,1.0,0.0);
        const double a_first = 0.5*a.squaredNorm();
        const double a_second = b.dot(a)+0.5*a.squaredNorm();
        std::cerr << "order example m=1,v0=0: W_A(A then B)=" << a_first
                  << ", W_A(B then A)=" << a_second << ", allocation shift=" << a.dot(b) << '\n';
        return 0;
    } catch (const std::exception& e) {
        std::cerr << e.what() << '\n'; return 1;
    }
}
