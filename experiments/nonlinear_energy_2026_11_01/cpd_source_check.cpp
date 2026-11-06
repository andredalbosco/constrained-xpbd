#include "baselines.hpp"
#include <iomanip>
#include <iostream>

// Compile against the existing, unchanged dynamics-research baseline.
int main() {
    using namespace work_balance;
    std::cout << "case,h,epsilon,relax,cap,x,y,vx,vy,rho,converged,iterations,stationarity,constraint,status\n";
    for (int c=0;c<4;++c) for (double h : {.02,.08})
    for (double eps : {1e4,1e8}) for (bool relax : {false,true})
    for (int cap : {1,5,200}) {
        System s;
        double mass = c<2 ? 2. : 1.4;
        double k = c<2 ? 100. : 80.;
        Vec3 x = c<2 ? Vec3(1.3,0,0) : Vec3(1.2,-.4,0);
        Vec3 v = c==0 ? Vec3::Zero() : (c==1 ? Vec3(.3,0,0) : Vec3(.3,.7,0));
        s.particles = {{{0,0,0},Vec3::Zero(),0},{x,v,1/mass}};
        s.constraints = {{1,0,1,1/k}};
        s.gravity = c==3 ? Vec3(0,-9.81,0) : Vec3::Zero();
        Options o; o.maxIterations=cap; o.cpRegularization=eps; o.cpRelaxEnergy=relax;
        auto r=step(s,h,Baseline::Cpd,o);
        auto p=s.particles[1];
        std::cout << std::setprecision(17) << c << ',' << h << ',' << eps << ',' << relax << ',' << cap
                  << ',' << p.position.x() << ',' << p.position.y()
                  << ',' << p.velocity.x() << ',' << p.velocity.y()
                  << ',' << r.energyRelaxation << ',' << r.converged << ',' << r.iterations
                  << ',' << r.stationarityResidual << ',' << r.constraintResidual << ',' << r.status << '\n';
    }
}
