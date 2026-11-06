#include "scenes.hpp"

#include <chrono>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

using namespace cxpbd;

namespace {
Mode parse_mode(const std::string& name) {
    if (name == "xpbd") return Mode::XPBD;
    if (name == "passive") return Mode::Passive;
    if (name == "active") return Mode::Active;
    throw std::invalid_argument("unknown mode: " + name);
}

ActiveOrder parse_order(const std::string& name) {
    if (name == "paired") return ActiveOrder::Paired;
    if (name == "paired-local") return ActiveOrder::PairedLocal;
    if (name == "sequential") return ActiveOrder::Sequential;
    throw std::invalid_argument("unknown order: " + name);
}

void run(const std::string& scene, const std::string& mode, Settings settings,
         int steps, const std::filesystem::path& output) {
    settings.mode = parse_mode(mode);
    const std::string method = mode != "active" || settings.active_order == ActiveOrder::Paired
        ? mode : mode + (settings.active_order == ActiveOrder::Sequential ? "-sequential" : "-paired-local");
    Model model = make_scene(scene);
    const Model initial = model;
    const double initial_energy = mechanical_energy(model);
    std::ofstream csv;
    if (!output.empty()) {
        std::filesystem::create_directories(output);
        csv.open(output / (scene + "-" + method + ".csv"));
        if (!csv) throw std::runtime_error("cannot open diagnostics CSV");
        csv << std::setprecision(17)
            << "step,time,mechanical_energy,particle_reservoir,spring_reservoir,reservoir_norm,"
               "extended_energy_drift,relative_stretch,max_displacement,material_residual,"
               "feedback_residual,position_residual,backtracks,rejected_updates\n";
        csv << "0,0," << initial_energy << ",0,0,0,0,"
            << maximum_relative_stretch(model) << ",0,0,0,0,0,0\n";
    }
    double max_drift = 0.0, max_extended_drift = 0.0;
    double max_stretch = maximum_relative_stretch(model), max_displacement = 0.0, max_feedback = 0.0;
    long long backtracks = 0, rejected = 0;
    Diagnostics diagnostics;
    const auto start = std::chrono::steady_clock::now();
    for (int index = 1; index <= steps; ++index) {
        diagnostics = step(model, settings);
        double particle_reservoir = 0.0, spring_reservoir = 0.0, displacement = 0.0;
        for (std::size_t i = 0; i < model.particles.size(); ++i) {
            particle_reservoir += model.particles[i].reservoir;
            displacement = std::max(displacement,
                (model.particles[i].x - initial.particles[i].x).norm());
        }
        for (const auto& s : model.springs) spring_reservoir += s.reservoir;
        const double drift = diagnostics.mechanical_energy - initial_energy;
        const double extended_drift = drift + diagnostics.reservoir_energy;
        const double stretch = maximum_relative_stretch(model);
        max_drift = std::max(max_drift, std::abs(drift));
        max_extended_drift = std::max(max_extended_drift, std::abs(extended_drift));
        max_stretch = std::max(max_stretch, stretch);
        max_displacement = std::max(max_displacement, displacement);
        max_feedback = std::max(max_feedback, diagnostics.max_feedback_residual);
        backtracks += diagnostics.backtracks;
        rejected += diagnostics.rejected_updates;
        if (csv) csv << index << ',' << index * settings.dt << ','
            << diagnostics.mechanical_energy << ',' << particle_reservoir << ','
            << spring_reservoir << ',' << diagnostics.reservoir_norm << ',' << extended_drift << ','
            << stretch << ',' << displacement << ',' << diagnostics.max_material_residual << ','
            << diagnostics.max_feedback_residual << ',' << diagnostics.position_residual_norm << ','
            << diagnostics.backtracks << ',' << diagnostics.rejected_updates << '\n';
    }
    const double milliseconds = std::chrono::duration<double, std::milli>(
        std::chrono::steady_clock::now() - start).count();
    if (!output.empty()) {
        std::ofstream state(output / (scene + "-" + method + "-particles.csv"));
        if (!state) throw std::runtime_error("cannot open particle state CSV");
        state << std::setprecision(17) << "particle,x,y,z,vx,vy,vz,mass,pinned,reservoir\n";
        for (std::size_t i = 0; i < model.particles.size(); ++i) {
            const auto& p = model.particles[i];
            state << i << ',' << p.x.x() << ',' << p.x.y() << ',' << p.x.z() << ','
                  << p.v.x() << ',' << p.v.y() << ',' << p.v.z() << ','
                  << p.mass << ',' << p.pinned << ',' << p.reservoir << '\n';
        }
        if (!csv || !state) throw std::runtime_error("failed writing output CSV");
    }
    std::cout << scene << ',' << method << ',' << model.particles.size() << ',' << model.springs.size()
              << ',' << initial_energy << ',' << diagnostics.mechanical_energy
              << ',' << diagnostics.reservoir_energy << ',' << max_drift
              << ',' << max_extended_drift << ',' << max_stretch << ',' << max_displacement
              << ',' << max_feedback << ',' << diagnostics.reservoir_norm
              << ',' << milliseconds / steps << ',' << backtracks << ',' << rejected << '\n';
}
} // namespace

int main(int argc, char** argv) {
    try {
        Settings settings;
        std::string scene = "all", mode = "all";
        std::filesystem::path output;
        int steps = 1200;
        for (int i = 1; i < argc; ++i) {
            const std::string option = argv[i];
            if (option == "--help") {
                std::cout << "Usage: constrained-xpbd [options]\n"
                    "  --scene spring|pendulum|grid|cloth|all (default all)\n"
                    "  --mode xpbd|passive|active|all         (default all)\n"
                    "  --order paired|paired-local|sequential (active mode; default paired)\n"
                    "  --steps N                            (default 1200)\n"
                    "  --dt H                               (default 1/240)\n"
                    "  --iterations N                       (default 20)\n"
                    "  --beta B --gamma G                    (kinetic/potential compliance; default 1000)\n"
                    "                                        Smaller is stronger; both must be positive.\n"
                    "  --output DIRECTORY                   (optional per-step CSV + final positions)\n";
                return 0;
            }
            if (i + 1 >= argc) throw std::invalid_argument("missing value for " + option);
            const std::string value = argv[++i];
            std::size_t parsed = 0;
            if (option == "--scene") scene = value;
            else if (option == "--mode") mode = value;
            else if (option == "--order") settings.active_order = parse_order(value);
            else if (option == "--output") output = value;
            else if (option == "--steps" || option == "--iterations") {
                const int number = std::stoi(value, &parsed);
                if (parsed != value.size() || number < 1) throw std::invalid_argument("invalid " + option);
                if (option == "--steps") steps = number;
                else settings.iterations = number;
            } else {
                const double number = std::stod(value, &parsed);
                if (parsed != value.size() || !std::isfinite(number) || number <= 0.0)
                    throw std::invalid_argument("invalid " + option);
                if (option == "--dt") settings.dt = number;
                else if (option == "--beta") settings.kinetic_reservoir_compliance = number;
                else if (option == "--gamma") settings.elastic_reservoir_compliance = number;
                else throw std::invalid_argument("unknown option: " + option);
            }
        }
        const std::vector<std::string> scenes = scene == "all"
            ? std::vector<std::string>{"spring", "pendulum", "grid", "cloth"}
            : std::vector<std::string>{scene};
        const std::vector<std::string> modes = mode == "all"
            ? std::vector<std::string>{"xpbd", "passive", "active"}
            : std::vector<std::string>{mode};
        std::cout << std::setprecision(10)
            << "scene,mode,particles,springs,initial_energy,final_energy,final_reservoir,"
               "max_mechanical_drift,max_extended_drift,max_relative_stretch,max_displacement,"
               "max_feedback_residual,reservoir_norm,ms_per_step_including_output,backtracks,rejected_updates\n";
        for (const auto& name : scenes)
            for (const auto& method : modes) run(name, method, settings, steps, output);
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << '\n';
        return 1;
    }
}
