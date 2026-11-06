// Reproducible spring work / kinematics study. No project headers are used.
// Build: c++ -std=c++17 -O2 -I/path/to/eigen spring_study.cpp -o spring_study
// Run from the desired output directory: ./spring_study [duration_seconds]
#include <Eigen/Dense>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
constexpr double mass = 1.0, stiffness = 100.0, amplitude = 0.25;
constexpr double omega = 10.0, pi = 3.14159265358979323846;

void check(bool condition, const std::string& message) {
    if (!condition) throw std::runtime_error(message);
}
struct State { double q, v; };
double kinetic(State s) { return 0.5 * mass * s.v * s.v; }
double elastic(State s) { return 0.5 * stiffness * s.q * s.q; }
double energy(State s) { return kinetic(s) + elastic(s); }

State advance(State s, double h, const std::string& method) {
    if (method == "backward_euler") {
        const double q = (mass * (s.q + h * s.v)) / (mass + stiffness * h * h);
        return {q, (q - s.q) / h};
    }
    if (method == "discrete_gradient_backward_kinematics") {
        const double a = stiffness * h * h / (2.0 * mass);
        const double q = ((1.0 - a) * s.q + h * s.v) / (1.0 + a);
        return {q, (q - s.q) / h};
    }
    if (method == "implicit_midpoint") {
        const double a = stiffness * h * h / (4.0 * mass);
        const double q = ((1.0 - a) * s.q + h * s.v) / (1.0 + a);
        return {q, 2.0 * (q - s.q) / h - s.v};
    }
    if (method == "midpoint_xpbd") {
        // The unknown is midpoint position y. C(y)=y and alpha=1/k.
        // Predict over half a step; the compliance scaling is 4 alpha / h^2.
        double y = s.q + 0.5 * h * s.v;
        const double alpha_tilde = 4.0 / (stiffness * h * h);
        const double lambda = -y / (1.0 / mass + alpha_tilde);
        y += lambda / mass;
        return {2.0 * y - s.q, s.v + 4.0 * lambda / (mass * h)};
    }
    throw std::runtime_error("unknown linear method");
}

bool midpoint_method(const std::string& method) {
    return method == "implicit_midpoint" || method == "midpoint_xpbd";
}

void linear_study(double duration) {
    std::ofstream csv("spring_linear.csv");
    csv.exceptions(std::ios::failbit | std::ios::badbit);
    csv << std::setprecision(17);
    csv << "method,h,steps,time,initial_energy,final_energy,relative_final_energy_error,"
           "max_relative_energy_error,phase_error_rad,normalized_endpoint_error,"
           "max_kinetic_minus_frozen_work,max_elastic_plus_frozen_work,"
           "max_kinetic_minus_impulse_work,max_elastic_plus_impulse_work,"
           "max_kinetic_minus_discrete_gradient_work,max_elastic_plus_discrete_gradient_work,"
           "max_kinetic_minus_actual_force_spatial_work,max_dynamics_residual,"
           "max_kinematics_residual,max_expected_dissipation_identity_error\n";
    const std::vector<std::string> methods = {"backward_euler",
        "discrete_gradient_backward_kinematics", "implicit_midpoint", "midpoint_xpbd"};
    for (const auto& method : methods) {
        for (int hz : {30, 60, 120, 240, 480, 960}) {
            const double h = 1.0 / hz;
            const int steps = static_cast<int>(std::lround(duration / h));
            const double time = steps * h;
            State state{amplitude, 0.0};
            const double initial_energy = energy(state);
            double max_energy = 0.0, phase = 0.0, last_angle = 0.0;
            double old_t = 0.0, old_u = 0.0, impulse_t = 0.0, impulse_u = 0.0;
            double dg_t = 0.0, dg_u = 0.0, actual_t = 0.0, dynamics = 0.0;
            double kinematics = 0.0, dissipation_identity = 0.0;
            for (int i = 0; i < steps; ++i) {
                const State old = state;
                state = advance(old, h, method);
                const double dq = state.q - old.q, dv = state.v - old.v;
                const double qbar = 0.5 * (old.q + state.q), vbar = 0.5 * (old.v + state.v);
                const double dt = kinetic(state) - kinetic(old), du = elastic(state) - elastic(old);
                const double frozen_work = -stiffness * old.q * dq;
                const double impulse_work = mass * dv * vbar;
                const double dg_work = -stiffness * qbar * dq;
                const double force = -stiffness * (method == "backward_euler" ? state.q : qbar);
                old_t = std::max(old_t, std::abs(dt - frozen_work));
                old_u = std::max(old_u, std::abs(du + frozen_work));
                impulse_t = std::max(impulse_t, std::abs(dt - impulse_work));
                impulse_u = std::max(impulse_u, std::abs(du + impulse_work));
                dg_t = std::max(dg_t, std::abs(dt - dg_work));
                dg_u = std::max(dg_u, std::abs(du + dg_work));
                actual_t = std::max(actual_t, std::abs(dt - force * dq));
                dynamics = std::max(dynamics, std::abs(mass * dv - h * force));
                kinematics = std::max(kinematics, std::abs(dq - h * (midpoint_method(method) ? vbar : state.v)));
                double expected_change = 0.0;
                if (method == "backward_euler")
                    expected_change = -0.5 * mass * dv * dv - 0.5 * stiffness * dq * dq;
                if (method == "discrete_gradient_backward_kinematics")
                    expected_change = -0.5 * mass * dv * dv;
                dissipation_identity = std::max(dissipation_identity, std::abs(dt + du - expected_change));
                max_energy = std::max(max_energy, std::abs(energy(state) / initial_energy - 1.0));
                const double angle = std::atan2(-state.v / omega, state.q);
                phase += std::remainder(angle - last_angle, 2.0 * pi);
                last_angle = angle;
            }
            const State exact{amplitude * std::cos(omega * time), -amplitude * omega * std::sin(omega * time)};
            const double endpoint_error = std::hypot(state.q - exact.q, (state.v - exact.v) / omega) / amplitude;
            const double phase_error = phase - omega * time;
            if (method == "backward_euler" || midpoint_method(method)) {
                const double increment = method == "backward_euler"
                    ? std::atan(omega * h) : 2.0 * std::atan(0.5 * omega * h);
                check(std::abs(phase - steps * increment) < 1e-8, "analytic oscillator phase check failed");
            }
            check(impulse_t < 1e-11 && dg_u < 1e-11, "work identity check failed");
            check(dissipation_identity < 1e-10, "derived dissipation identity failed");
            if (midpoint_method(method)) check(max_energy < 1e-9, "quadratic midpoint energy preservation failed");
            csv << method << ',' << h << ',' << steps << ',' << time << ',' << initial_energy << ','
                << energy(state) << ',' << energy(state) / initial_energy - 1.0 << ',' << max_energy << ','
                << phase_error << ',' << endpoint_error << ',' << old_t << ',' << old_u << ','
                << impulse_t << ',' << impulse_u << ',' << dg_t << ',' << dg_u << ','
                << actual_t << ',' << dynamics << ',' << kinematics << ',' << dissipation_identity << '\n';
            std::cout << method << " h=1/" << hz << " E/E0-1=" << energy(state) / initial_energy - 1.0
                      << " phase_error=" << phase_error << " max_dT-Wimp=" << impulse_t
                      << " max_dU+Wimp=" << impulse_u << '\n';
        }
    }
    csv.close();
}

void degeneracy_study() {
    std::ofstream csv("spring_degeneracy.csv");
    csv.exceptions(std::ios::failbit | std::ios::badbit);
    csv << std::setprecision(17);
    csv << "case,h,q_new,v_new,energy_error,backward_euler_objective,kinetic_work,elastic_work,backward_euler_momentum_residual\n";
    const double h = 1.0 / 60.0;
    const State old{amplitude, 0.0};
    auto emit = [&](const std::string& name, State next) {
        const double objective = 0.5 * mass * std::pow((next.q - old.q - h * old.v) / h, 2) + elastic(next);
        csv << name << ',' << h << ',' << next.q << ',' << next.v << ',' << energy(next) - energy(old) << ','
            << objective << ',' << kinetic(next) - kinetic(old) << ',' << elastic(old) - elastic(next) << ','
            << mass * (next.v - old.v) + h * stiffness * next.q << '\n';
    };
    emit("energy_constrained_backward_no_motion", old);
    const double dq = -2.0 * stiffness * old.q * h * h / (mass + stiffness * h * h);
    emit("energy_constrained_backward_moving", {old.q + dq, dq / h});
    emit("midpoint_turning_point", advance(old, h, "implicit_midpoint"));
    emit("exact_turning_point", {amplitude * std::cos(omega * h), -amplitude * omega * std::sin(omega * h)});
    const State at_equilibrium = advance({0.0, 0.0}, h, "midpoint_xpbd");
    check(at_equilibrium.q == 0.0 && at_equilibrium.v == 0.0, "zero-state degeneration failed");
    // A step can have dq=0 while reversing velocity. The polynomial discrete gradient remains finite.
    const State chord_start{amplitude, 0.5 * h * stiffness * amplitude / mass};
    const State chord_end = advance(chord_start, h, "midpoint_xpbd");
    check(std::abs(chord_end.q - chord_start.q) < 1e-15, "zero-displacement turning chord failed");
    check(std::abs(chord_end.v + chord_start.v) < 1e-14, "turning chord velocity reversal failed");
    std::cout << "zero-displacement midpoint chord: q=" << chord_end.q << " v_old=" << chord_start.v
              << " v_new=" << chord_end.v << " (no energy/work division used)\n";
    csv.close();
}

using Vec3 = Eigen::Vector3d;
using Mat3 = Eigen::Matrix3d;
constexpr double m1 = 2.0, m2 = 1.0, reduced_mass = m1 * m2 / (m1 + m2), rest = 1.0;
struct RelativeState { Vec3 r, u; };
struct Gradient { Vec3 value; Mat3 derivative; };

Gradient radial_gradient(const Vec3& old_r, const Vec3& r, const std::string& method) {
    if (method == "midpoint_discrete_gradient") {
        const double sum_length = old_r.norm() + r.norm();
        check(sum_length > 1e-10 && r.norm() > 1e-10, "radial DG collapse");
        const Vec3 midpoint = 0.5 * (old_r + r);
        const double a = 1.0 - 2.0 * rest / sum_length;
        return {stiffness * a * midpoint,
            stiffness * (0.5 * a * Mat3::Identity()
                + (2.0 * rest / (sum_length * sum_length)) * midpoint * (r / r.norm()).transpose())};
    }
    const Vec3 point = method == "backward_euler" ? r : (0.5 * (old_r + r)).eval();
    const double length = point.norm();
    check(length > 1e-10, "radial force collapse");
    const double scale = method == "backward_euler" ? 1.0 : 0.5;
    return {stiffness * (1.0 - rest / length) * point,
        scale * stiffness * ((1.0 - rest / length) * Mat3::Identity()
            + rest / (length * length * length) * point * point.transpose())};
}

RelativeState radial_step(RelativeState old, double h, const std::string& method,
                          double& max_newton_residual, int& max_iterations) {
    const bool backward = method == "backward_euler";
    const double factor = h * h / ((backward ? 1.0 : 2.0) * reduced_mass);
    Vec3 r = old.r + h * old.u;
    int iteration = 0;
    for (; iteration < 30; ++iteration) {
        const Gradient gradient = radial_gradient(old.r, r, method);
        const Vec3 residual = r - old.r - h * old.u + factor * gradient.value;
        if (residual.norm() < 2e-14 * (1.0 + r.norm())) break;
        const Vec3 correction = (Mat3::Identity() + factor * gradient.derivative).fullPivLu().solve(-residual);
        bool accepted = false;
        double scale = 1.0;
        for (int attempt = 0; attempt < 20; ++attempt, scale *= 0.5) {
            const Vec3 candidate = r + scale * correction;
            const Vec3 next_residual = candidate - old.r - h * old.u
                + factor * radial_gradient(old.r, candidate, method).value;
            if (next_residual.norm() < residual.norm()) {
                r = candidate; accepted = true; break;
            }
        }
        check(accepted, "radial Newton line search failed");
    }
    const Vec3 residual = r - old.r - h * old.u + factor * radial_gradient(old.r, r, method).value;
    max_newton_residual = std::max(max_newton_residual, residual.norm());
    max_iterations = std::max(max_iterations, iteration);
    check(iteration < 30 && residual.norm() < 1e-11, "radial nonlinear solve failed");
    const Vec3 u = backward ? ((r - old.r) / h).eval() : (2.0 * (r - old.r) / h - old.u).eval();
    return {r, u};
}

double radial_elastic(const Vec3& r) { return 0.5 * stiffness * std::pow(r.norm() - rest, 2); }
double radial_kinetic(const Vec3& u) { return 0.5 * reduced_mass * u.squaredNorm(); }

void radial_study(double duration) {
    std::ofstream csv("spring_radial.csv");
    csv.exceptions(std::ios::failbit | std::ios::badbit);
    csv << std::setprecision(17);
    csv << "method,h,steps,time,initial_internal_energy,relative_final_energy_error,max_relative_energy_error,"
           "max_absolute_linear_momentum_error,max_absolute_angular_momentum_error,"
           "max_kinetic_minus_impulse_work,max_elastic_plus_impulse_work,"
           "max_elastic_plus_actual_force_spatial_work,max_newton_residual,max_newton_iterations\n";
    for (const std::string method : {"backward_euler", "implicit_midpoint", "midpoint_discrete_gradient"}) {
        for (int hz : {30, 60, 120, 240}) {
            const double h = 1.0 / hz;
            const int steps = static_cast<int>(std::lround(duration / h));
            RelativeState state{Vec3(1.25, 0.0, 0.0), Vec3(0.0, 1.0, 0.0)};
            const Vec3 initial_c(0.4, 0.7, -0.1), translation(0.3, -0.2, 0.1);
            const Vec3 initial_p = (m1 + m2) * translation;
            const Vec3 initial_l = initial_c.cross(initial_p) + reduced_mass * state.r.cross(state.u);
            const double initial_energy = radial_elastic(state.r) + radial_kinetic(state.u);
            double max_energy = 0.0, max_p = 0.0, max_l = 0.0, impulse_t = 0.0, impulse_u = 0.0;
            double spatial_u = 0.0, nonlinear_residual = 0.0;
            int nonlinear_iterations = 0;
            for (int step = 0; step < steps; ++step) {
                const RelativeState old = state;
                state = radial_step(old, h, method, nonlinear_residual, nonlinear_iterations);
                const Vec3 center = initial_c + (step + 1) * h * translation;
                const Vec3 x1 = center + (m2 / (m1 + m2)) * state.r;
                const Vec3 x2 = center - (m1 / (m1 + m2)) * state.r;
                const Vec3 v1 = translation + (m2 / (m1 + m2)) * state.u;
                const Vec3 v2 = translation - (m1 / (m1 + m2)) * state.u;
                const Vec3 p = m1 * v1 + m2 * v2;
                const Vec3 l = x1.cross(m1 * v1) + x2.cross(m2 * v2);
                max_p = std::max(max_p, (p - initial_p).norm());
                max_l = std::max(max_l, (l - initial_l).norm());
                const double dt = radial_kinetic(state.u) - radial_kinetic(old.u);
                const double du = radial_elastic(state.r) - radial_elastic(old.r);
                const double wi = reduced_mass * (state.u - old.u).dot(0.5 * (state.u + old.u));
                const double ws = -radial_gradient(old.r, state.r, method).value.dot(state.r - old.r);
                impulse_t = std::max(impulse_t, std::abs(dt - wi));
                impulse_u = std::max(impulse_u, std::abs(du + wi));
                spatial_u = std::max(spatial_u, std::abs(du + ws));
                max_energy = std::max(max_energy,
                    std::abs((radial_kinetic(state.u) + radial_elastic(state.r)) / initial_energy - 1.0));
            }
            const double final_error = (radial_kinetic(state.u) + radial_elastic(state.r)) / initial_energy - 1.0;
            if (method != "backward_euler") check(max_l < 1e-8, "midpoint angular momentum preservation failed");
            if (method == "midpoint_discrete_gradient") check(max_energy < 1e-8, "radial DG energy preservation failed");
            check(max_p < 1e-12, "center-of-mass reconstruction failed");
            csv << method << ',' << h << ',' << steps << ',' << steps * h << ',' << initial_energy << ','
                << final_error << ',' << max_energy << ',' << max_p << ',' << max_l << ','
                << impulse_t << ',' << impulse_u << ',' << spatial_u << ',' << nonlinear_residual << ','
                << nonlinear_iterations << '\n';
            std::cout << "radial " << method << " h=1/" << hz << " E/E0-1=" << final_error
                      << " max_E_error=" << max_energy << " max_L_error=" << max_l << '\n';
        }
    }
    csv.close();
}

struct LinearizedStep {
    RelativeState state;
    Vec3 position_residual;
    Vec3 material_impulse;
    double material_residual;
};

LinearizedStep radial_linearized_step(RelativeState old, double h, int sweeps,
                                     bool keep_position_residual, bool symmetric) {
    const Vec3 predicted = old.r + 0.5 * h * old.u;
    Vec3 midpoint = predicted;
    double lambda = 0.0;
    const double alpha_tilde = 4.0 / (stiffness * h * h);
    const double old_length = old.r.norm();
    for (int sweep = 0; sweep < sweeps; ++sweep) {
        const Vec3 endpoint = 2.0 * midpoint - old.r;
        const double length = endpoint.norm();
        check(length > 1e-10 && old_length > 1e-10, "linearized spring collapse");
        const Vec3 j_bar = (endpoint + old.r) / (length + old_length);
        const Vec3 j_end = endpoint / length;
        const Vec3 derivative = symmetric ? j_bar : j_end;
        const double constraint = 0.5 * (length + old_length) - rest;
        const Vec3 actual_rx = reduced_mass * (midpoint - predicted) - j_bar * lambda;
        const Vec3 used_rx = keep_position_residual ? actual_rx : Vec3::Zero().eval();
        const double r_lambda = constraint + alpha_tilde * lambda;
        const double denominator = derivative.dot(j_bar) / reduced_mass + alpha_tilde;
        check(std::isfinite(denominator) && denominator > 0.0, "linearized denominator invalid");
        const double d_lambda = (-r_lambda + derivative.dot(used_rx) / reduced_mass) / denominator;
        midpoint += (j_bar * d_lambda - used_rx) / reduced_mass;
        lambda += d_lambda;
    }
    const Vec3 endpoint = 2.0 * midpoint - old.r;
    const double length = endpoint.norm();
    const Vec3 j_bar = (endpoint + old.r) / (length + old_length);
    const Vec3 rx = reduced_mass * (midpoint - predicted) - j_bar * lambda;
    const double r_lambda = 0.5 * (length + old_length) - rest + alpha_tilde * lambda;
    const Vec3 velocity = 4.0 * (midpoint - old.r) / h - old.u;
    check(endpoint.allFinite() && velocity.allFinite(), "non-finite linearized trajectory");
    return {{endpoint, velocity}, rx, (4.0 * lambda / h) * j_bar, r_lambda};
}

void linearized_study(double duration) {
    std::ofstream csv("spring_linearized.csv");
    csv.exceptions(std::ios::failbit | std::ios::badbit);
    csv << std::setprecision(17);
    csv << "variant,sweeps,h,steps,time,relative_final_energy_error,max_relative_energy_error,"
           "max_absolute_angular_momentum_error,max_position_residual,"
           "max_force_residual,max_material_residual,max_kinetic_minus_material_impulse_work,"
           "max_elastic_plus_material_impulse_work,max_kinetic_minus_actual_impulse_work,"
           "max_impulse_residual_identity_error,final_position_distance_to_converged_DG,"
           "final_velocity_distance_to_converged_DG\n";
    const double h = 1.0 / 60.0;
    const int steps = static_cast<int>(std::lround(std::min(duration, 10.0) / h));
    for (bool keep_rx : {true, false}) {
        for (bool symmetric : {false, true}) {
            const std::string variant = std::string(keep_rx ? "retain_rx" : "drop_rx")
                                      + (symmetric ? "_Jbar" : "_Jend");
            for (int sweeps : {1, 5, 20, 100}) {
                RelativeState state{Vec3(1.25, 0.0, 0.0), Vec3(0.0, 1.0, 0.0)};
                RelativeState reference = state;
                const Vec3 initial_l = reduced_mass * state.r.cross(state.u);
                const double initial_energy = radial_elastic(state.r) + radial_kinetic(state.u);
                double max_energy = 0.0, max_l = 0.0, max_rx = 0.0, max_force = 0.0, max_r_lambda = 0.0;
                double kinetic_defect = 0.0, elastic_defect = 0.0, actual_impulse_defect = 0.0;
                double impulse_identity = 0.0, reference_residual = 0.0;
                int reference_iterations = 0;
                for (int step = 0; step < steps; ++step) {
                    const RelativeState old = state;
                    const LinearizedStep result = radial_linearized_step(old, h, sweeps, keep_rx, symmetric);
                    state = result.state;
                    reference = radial_step(reference, h, "midpoint_discrete_gradient",
                                            reference_residual, reference_iterations);
                    const Vec3 actual_impulse = reduced_mass * (state.u - old.u);
                    const Vec3 average_velocity = 0.5 * (state.u + old.u);
                    const double dt = radial_kinetic(state.u) - radial_kinetic(old.u);
                    const double du = radial_elastic(state.r) - radial_elastic(old.r);
                    const double material_work = result.material_impulse.dot(average_velocity);
                    kinetic_defect = std::max(kinetic_defect, std::abs(dt - material_work));
                    elastic_defect = std::max(elastic_defect, std::abs(du + material_work));
                    actual_impulse_defect = std::max(actual_impulse_defect,
                        std::abs(dt - actual_impulse.dot(average_velocity)));
                    impulse_identity = std::max(impulse_identity,
                        (actual_impulse - result.material_impulse - (4.0 / h) * result.position_residual).norm());
                    max_rx = std::max(max_rx, result.position_residual.norm());
                    max_force = std::max(max_force, (4.0 / (h * h)) * result.position_residual.norm());
                    max_r_lambda = std::max(max_r_lambda, std::abs(result.material_residual));
                    max_energy = std::max(max_energy,
                        std::abs((radial_kinetic(state.u) + radial_elastic(state.r)) / initial_energy - 1.0));
                    max_l = std::max(max_l, (reduced_mass * state.r.cross(state.u) - initial_l).norm());
                }
                const double final_error = (radial_kinetic(state.u) + radial_elastic(state.r)) / initial_energy - 1.0;
                check(actual_impulse_defect < 1e-11 && impulse_identity < 1e-10,
                      "linearized work/residual identity failed");
                if (keep_rx && sweeps == 100) {
                    check(max_energy < 1e-8 && max_l < 1e-8, "retained-residual solver did not converge conservatively");
                    check((state.r - reference.r).norm() < 1e-7, "retained-residual solver differs from DG reference");
                }
                csv << variant << ',' << sweeps << ',' << h << ',' << steps << ',' << steps * h << ','
                    << final_error << ',' << max_energy << ',' << max_l << ',' << max_rx << ',' << max_force << ','
                    << max_r_lambda << ',' << kinetic_defect << ',' << elastic_defect << ','
                    << actual_impulse_defect << ',' << impulse_identity << ','
                    << (state.r - reference.r).norm() << ',' << (state.u - reference.u).norm() << '\n';
                std::cout << "linearized " << variant << " sweeps=" << sweeps
                          << " max_E_error=" << max_energy << " max_L_error=" << max_l
                          << " force_residual=" << max_force << " material_residual=" << max_r_lambda << '\n';
            }
        }
    }
    csv.close();
}
} // namespace

int main(int argc, char** argv) {
    try {
        const double duration = argc > 1 ? std::stod(argv[1]) : 20.0;
        check(duration > 0.0 && duration < 10000.0, "duration must be in (0, 10000)");
        std::cout << std::setprecision(12);
        linear_study(duration);
        degeneracy_study();
        radial_study(duration);
        linearized_study(duration);
        std::cout << "All independent work / analytic checks passed. CSV files written.\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "FAIL: " << error.what() << '\n';
        return 1;
    }
}
