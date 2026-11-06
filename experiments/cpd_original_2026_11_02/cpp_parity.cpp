// Source cross-check against the unchanged local 3D CPD baseline.
#include "baselines.hpp"
#include <iomanip>
#include <iostream>

int main() {
    std::cout << std::setprecision(17);
    std::cout << "case,cap,epsilon,iterations,rho,rank,status,x0,y0,vx0,vy0,x1,y1,vx1,vy1\n";
    for (int kind = 0; kind < 3; ++kind) {
        const char* name = kind == 0 ? "free_release" : kind == 1 ? "free_orbit" : "free_stiff";
        const double radius = kind == 0 ? 1.25 : 1.1;
        const double speed = kind == 0 ? 0.0 : 0.8;
        const double stiffness = kind == 2 ? 1000.0 : 100.0;
        for (double epsilon : {1e4, 1e8}) {
            for (int cap : {1, 4, 16}) {
                System system;
                system.particles = {{{radius/2,0,0},{0,speed/2,0},.5},
                                    {{-radius/2,0,0},{0,-speed/2,0},.5}};
                system.constraints = {{0,1,1.0,1.0/stiffness}};
                work_balance::Options options;
                options.maxIterations = cap;
                options.tolerance = 1e-30;
                options.cpRegularization = epsilon;
                const auto result = work_balance::step(system, 1.0/60.0,
                                                       work_balance::Baseline::Cpd, options);
                const auto& a = system.particles[0];
                const auto& b = system.particles[1];
                std::cout << name << ',' << cap << ',' << epsilon << ',' << result.iterations
                          << ',' << result.energyRelaxation << ',' << result.constraintRank
                          << ',' << result.status << ',' << a.position.x() << ',' << a.position.y()
                          << ',' << a.velocity.x() << ',' << a.velocity.y() << ',' << b.position.x()
                          << ',' << b.position.y() << ',' << b.velocity.x() << ',' << b.velocity.y() << '\n';
            }
        }
    }
}
