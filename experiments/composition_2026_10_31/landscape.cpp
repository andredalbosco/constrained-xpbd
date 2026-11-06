#include "independent.hpp"
#include <fstream>
#include <iomanip>
#include <iostream>
using namespace compose;
// Independent dense search of the scalar constant-energy ellipse. This checks
// the objective itself without relying on the KKT optimizer's local branch.
int main(int argc,char**argv){std::ofstream out(argc>1?argv[1]:"landscape.csv");out<<std::setprecision(17)<<"case,hz,kinematics,beta,minimum_x,minimum_v,minimum_objective,numerical_x,numerical_v,numerical_objective,objective_gap,energy_residual,stationarity,reference_x,reference_v\n";
    for(const auto&name:{"scalar_release","scalar_mixed"}){for(int hz:{60,120,240}){for(const auto&kin:{"soft_endpoint","soft_midpoint"}){for(double beta:{.1,1.,10.}){
        auto c=make_case(name);auto&m=c.model;auto&s=c.initial;double h=1./hz,E=energy(m,s),A=std::sqrt(2*E/m.stiffness[0]),B=std::sqrt(2*E/m.mass[0]);
        auto state=[&](double t){State z{V(1),V(1)};z.x[0]=m.length[0]+A*std::cos(t);z.v[0]=B*std::sin(t);return z;};
        auto obj=[&](double t){auto z=state(t);double u=(z.x[0]-s.x[0])/h,r=std::string(kin)=="soft_midpoint"?z.v[0]+s.v[0]-2*u:z.v[0]-u;return (.5*m.mass[0]*(u-s.v[0])*(u-s.v[0])+potential(m,z.x)+.5*m.mass[0]*r*r/beta)/E;};
        const int n=65536;const double twopi=2*std::acos(-1.0),dt=twopi/n;int best=0;double val=obj(0);for(int k=1;k<n;++k){double v=obj(k*dt);if(v<val){val=v;best=k;}}
        double lo=(best-1)*dt,hi=(best+1)*dt;for(int k=0;k<100;++k){double a=lo+(hi-lo)/3,b=hi-(hi-lo)/3;if(obj(a)<obj(b))hi=b;else lo=a;}double t=.5*(lo+hi);auto z=state(t);Config cfg;cfg.family="independent";cfg.treatment="global";cfg.kinematics=kin;cfg.beta=beta;cfg.iterations=128;cfg.curvature="exact";auto r=independent_step(m,s,h,cfg,E);auto ref=reference_step(m,s,h,256);double u=(r.state.x[0]-s.x[0])/h,rr=std::string(kin)=="soft_midpoint"?r.state.v[0]+s.v[0]-2*u:r.state.v[0]-u;double numeric=(.5*m.mass[0]*(u-s.v[0])*(u-s.v[0])+potential(m,r.state.x)+.5*m.mass[0]*rr*rr/beta)/E;
        out<<name<<','<<hz<<','<<kin<<','<<beta<<','<<z.x[0]<<','<<z.v[0]<<','<<obj(t)<<','<<r.state.x[0]<<','<<r.state.v[0]<<','<<numeric<<','<<numeric-obj(t)<<','<<r.balance<<','<<r.rx<<','<<ref.x[0]<<','<<ref.v[0]<<'\n';
    }}}}
    std::cout<<"36 scalar objective landscapes searched with 65536 samples plus local refinement\n";
}
