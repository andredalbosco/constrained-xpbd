#include "relinearized.hpp"
#include <iostream>
using namespace traditional;
void require(bool ok,const char* message){if(!ok)throw std::runtime_error(message);}
void near(double a,double b,double t,const char* message){require(std::isfinite(a)&&std::isfinite(b)&&std::abs(a-b)<=t,message);}
void near(V2 a,V2 b,double t,const char* message){near(norm(a-b),0,t,message);}
void spring_affine_and_zero_gradient(){
    Model m;State s;s.x[0]={1.25,0};const double h=1./480,E0=energy(m,s);
    for(int step=0;step<2;++step){
        const auto c=context(m,s,h,E0);Iterate frozen{c.predictor},fresh=frozen;
        Options o;o.schedule=Schedule::EL;Diagnostics a,b;
        require(update(c,frozen,o,a)&&update_relinearized(c,fresh,o,b),"EL update failed");
        near(fresh.x[0],frozen.x[0],1e-11,"affine spring EL is unchanged by rebuilding material row");
        near(fresh.lambda[0],frozen.lambda[0],1e-12,"affine spring material multiplier agrees");
        near(fresh.eta,frozen.eta,1e-10,"affine spring energy multiplier agrees");
        s.v[0]=(fresh.x[0]-s.x[0])/h;s.x=fresh.x;
    }
    State old;old.x[0]={1.25,0};auto c=context(m,old,h,E0);Iterate fresh{c.predictor},baseline=fresh;
    Options o;o.schedule=Schedule::LE;Diagnostics d,b;
    require(update_relinearized(c,fresh,o,d)&&fresh_material(c,baseline,b),"fresh LE release failed");
    require(d.singular==1,"zero energy gradient after first material solve must be skipped within resolution");
    near(fresh.x[0],baseline.x[0],5e-16,"first refreshed LE step takes baseline motion instead of freezing");
    require(norm(fresh.x[0]-old.x[0])>1e-5,"refreshed LE moves a released spring");
}
void independent_pendulum_visit_oracles(){
    Model m;m.alpha[0]=0;m.gravity={0,-9.81};State old;
    old.x[0]={std::sin(.4),-std::cos(.4)};old.v[0]=tangent(old.x[0])*1.2;
    const double h=1./120,target=energy(m,old)+.1;const auto c=context(m,old,h,target);
    auto E=[&](V2 x){const V2 v=(x-old.x[0])/h;return .5*dot(v,v)-dot(m.gravity,x);};
    auto gradient=[&](V2 x){return (x-old.x[0])/(h*h)-m.gravity;};
    for(auto schedule:{Schedule::EL,Schedule::LE}){
        V2 x=c.predictor[0];double lambda=0,eta=0;
        auto material=[&](){const double radius=norm(x);lambda=-(radius-1);x=x/radius;};
        auto energy_visit=[&](){const V2 g=gradient(x);eta=-(E(x)-target)/dot(g,g);x+=g*eta;};
        if(schedule==Schedule::EL){energy_visit();material();}else{material();energy_visit();}
        Iterate z{c.predictor};Options o;o.schedule=schedule;Diagnostics d;
        require(update_relinearized(c,z,o,d),"pendulum fresh update failed");
        near(z.x[0],x,3e-14,"fresh positions match independently reconstructed visit sequence");
        near(z.lambda[0],lambda,3e-15,"fresh lambda uses post-energy geometry");
        near(z.eta,eta,3e-15,"fresh nu uses post-material energy when in LE order");
        if(schedule==Schedule::EL)near(norm(z.x[0]),1,5e-16,"last rebuilt material row preserves rod length exactly");
    }
}
void repair_cap_properties(){
    Model spring,rod;rod.alpha[0]=0;rod.gravity={0,-9.81};
    State ss;ss.x[0]={1.2,0};ss.v[0]={-.6,.2};
    State ps;ps.x[0]={std::sin(.4),-std::cos(.4)};ps.v[0]=tangent(ps.x[0])*1.2;
    for(int k=0;k<2;++k)for(auto schedule:{Schedule::EL,Schedule::LE}){
        const Model m=k?rod:spring;const State old=k?ps:ss;
        const auto c=context(m,old,1./240,energy(m,old)+.01);const Iterate initial{c.predictor};
        const auto r=evaluate(c,initial);const auto base=material_pass(m,r,{});const auto dx0=displacement(m,r,base);
        Options raw;raw.schedule=schedule;Options bounded=raw;bounded.policy=Policy::RepairCap;
        Iterate proposal=initial,accepted=initial;Diagnostics rd,bd;
        require(update_relinearized(c,proposal,raw,rd)&&update_relinearized(c,accepted,bounded,bd),"guarded fresh update failed");
        Vec extra{},actual{};for(int i=0;i<m.count;++i){extra[i]=proposal.x[i]-initial.x[i]-dx0[i];actual[i]=accepted.x[i]-initial.x[i]-dx0[i];}
        require(mass_norm(m,actual)<=.1*motion_scale(c)+2e-15,"actual refreshed extra displacement obeys cap");
        near(actual[0],extra[0]*bd.last_blend,2e-15,"accepted fresh positions blend actual proposal, not frozen reconstruction");
        const double H0=r.H+inner(r.g,dx0,m.count),change=inner(r.g,actual,m.count);
        require(H0*change<=1e-12&&std::abs(change)<=.5*std::abs(H0)+1e-12,"fresh guard bounds repair and improves predicted baseline energy");
        near(accepted.lambda[0],base[0]+bd.last_blend*(proposal.lambda[0]-base[0]),2e-15,"material multiplier follows same accepted blend");
        near(accepted.eta,bd.last_blend*proposal.eta,2e-15,"energy multiplier follows same accepted blend");
    }
}
int main(){try{spring_affine_and_zero_gradient();independent_pendulum_visit_oracles();repair_cap_properties();
    std::cout<<"PASS: refreshed affine equivalence; zero-gradient release; independent pendulum visit oracles; actual position and multiplier guard bounds\n";
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
