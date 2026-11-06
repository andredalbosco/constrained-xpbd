#include "/home/blacknoir/Projects/constrained-xpbd/experiments/traditional_cases/solver.hpp"
#include <cmath>
#include <iomanip>
#include <iostream>

// The original radial solver remains on its positive-x affine branch during
// the first two release steps. Check that exact historical event independently.
int main() {
    using namespace traditional;
    std::cout << "hz,step,original_x,original_v,original_energy,affine_x,affine_v,affine_energy,position_difference,velocity_difference\n";
    for (int hz : {60, 120, 240, 480}) {
        Model model;
        model.mass[0] = 1;
        model.alpha[0] = .01;
        model.length[0] = 1;
        State old;
        old.x[0] = {1.25, 0};
        double q = .25, v = 0;
        const double h = 1.0 / hz, target = 3.125;
        Options options;
        options.schedule = Schedule::EL;
        options.policy = Policy::Raw;
        for (int step = 1; step <= 2; ++step) {
            auto ctx = context(model, old, h, target);
            Iterate iterate{ctx.predictor};
            Diagnostics stats;
            if (!update(ctx, iterate, options, stats)) return 2;
            State next;
            next.x = iterate.x;
            next.v[0] = (next.x[0] - old.x[0]) / h;
            const double prediction = q + h * v;
            const double residual = .5 * v * v + 50 * prediction * prediction - target;
            const double gradient = v / h + 100 * prediction;
            const double q_next = (prediction - residual / gradient) / (1 + 100 * h * h);
            const double v_next = (q_next - q) / h;
            const double dx = next.x[0].x - (1 + q_next);
            const double dv = next.v[0].x - v_next;
            std::cout << std::setprecision(17) << hz << ',' << step << ','
                      << next.x[0].x << ',' << next.v[0].x << ',' << energy(model, next) << ','
                      << 1 + q_next << ',' << v_next << ',' << .5 * v_next * v_next + 50 * q_next * q_next << ','
                      << dx << ',' << dv << '\n';
            // Near-zero gradients amplify representation differences; test
            // their resulting displacement, rather than pretending bit parity.
            if (std::abs(dx) > 1e-7 || std::abs(dv) > 1e-4) return 3;
            old = next;
            q = q_next;
            v = v_next;
        }
    }
}
