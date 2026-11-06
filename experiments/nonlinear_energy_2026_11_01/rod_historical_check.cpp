#include "/home/blacknoir/Projects/constrained-xpbd/experiments/traditional_cases/relinearized.hpp"
#include <iomanip>
#include <iostream>
int main() {
    using namespace traditional;
    std::cout << "hz,theta,omega,mass,length,step,x,y,vx,vy,energy,singular\n";
    for (int hz : {60,480}) for (double theta : {.25,.65,1.2}) for (double omega : {0.,.7}) {
        Model model;
        const double mass=theta<.5?1.:1.7, length=theta<1?.7:1.4, h=1./hz;
        model.mass[0]=mass;model.length[0]=length;model.alpha[0]=0;model.gravity={0,-9.81};
        State old;
        old.x[0]={length*std::sin(theta),-length*std::cos(theta)};
        old.v[0]={length*omega*std::cos(theta),length*omega*std::sin(theta)};
        const double target=energy(model,old);
        Options options;options.schedule=Schedule::EL;options.policy=Policy::Raw;
        for (int step=1;step<=2;++step) {
            const auto ctx=context(model,old,h,target);
            Iterate iterate{ctx.predictor};Diagnostics stats;
            // planar.raw_EL refreshes the rod normal after moving E, so its
            // matching original control is the refreshed EL implementation.
            if (!update_relinearized(ctx,iterate,options,stats)) return 2;
            State next;next.x=iterate.x;next.v[0]=(next.x[0]-old.x[0])/h;
            std::cout<<std::setprecision(17)<<hz<<','<<theta<<','<<omega<<','<<mass<<','<<length<<','<<step<<','
                     <<next.x[0].x<<','<<next.x[0].y<<','<<next.v[0].x<<','<<next.v[0].y<<','<<energy(model,next)<<','<<stats.singular<<'\n';
            old=next;
        }
    }
}
