#pragma once
#include "solver.hpp"
namespace agenda {
struct Case {std::string name;Model model;State initial;};
V2 direction(double angle){return {std::sin(angle),-std::cos(angle)};}
std::vector<Case> cases(){
    Model spring,rod,chain;rod.alpha[0]=0;rod.gravity={0,-9.81};chain.count=2;chain.gravity=rod.gravity;
    State release,orbit,rest,kick;release.x[0]={1.25,0};orbit.x[0]={1.1,0};orbit.v[0]={0,.8};
    rest.x[0]=direction(.65);kick.x[0]=direction(.4);kick.v[0]=tangent(kick.x[0])*1.2;
    auto double_state=[](double a,double b,double dr1,double dr2,double vr1,double vt1,double vr2,double vt2){
        State s;const auto n1=direction(a),n2=direction(b);
        s.x[0]=n1*(1.1962+dr1);s.x[1]=s.x[0]+n2*(1.0981+dr2);
        s.v[0]=n1*vr1+tangent(n1)*vt1;s.v[1]=s.v[0]+n2*vr2+tangent(n2)*vt2;return s;};
    return {{"spring_release",spring,release},{"spring_orbit",spring,orbit},{"pendulum_rest",rod,rest},{"pendulum_kick",rod,kick},
        {"double_small",chain,double_state(.15,-.1,0,0,0,0,0,0)},
        {"double_large",chain,double_state(.65,-.45,.1,.2,0,0,0,0)},
        {"double_mixed",chain,double_state(.45,-.3,.05,.1,-.2,.5,.3,.8)}};
}
State reference_step(const Model& m,State s,double h){
    if(m.count==1&&m.alpha[0]==0){
        const double l=m.length[0],a=std::atan2(s.x[0].x,-s.x[0].y),w=dot(tangent(s.x[0]),s.v[0])/(l*l);
        auto f=[&](double t){return -norm(m.gravity)*std::sin(t)/l;};
        const double a1=w,w1=f(a),a2=w+h*w1/2,w2=f(a+h*a1/2),a3=w+h*w2/2,w3=f(a+h*a2/2),a4=w+h*w3,w4=f(a+h*a3);
        const double an=a+h*(a1+2*a2+2*a3+a4)/6,wn=w+h*(w1+2*w2+2*w3+w4)/6;
        s.x[0]=direction(an)*l;s.v[0]=tangent(s.x[0])*wn;return s;
    }
    const Vec a1=acceleration(m,s.x);State s2=s,s3=s,s4=s;
    for(int i=0;i<m.count;++i){s2.x[i]+=s.v[i]*(h/2);s2.v[i]+=a1[i]*(h/2);}const Vec a2=acceleration(m,s2.x);
    for(int i=0;i<m.count;++i){s3.x[i]+=s2.v[i]*(h/2);s3.v[i]+=a2[i]*(h/2);}const Vec a3=acceleration(m,s3.x);
    for(int i=0;i<m.count;++i){s4.x[i]+=s3.v[i]*h;s4.v[i]+=a3[i]*h;}const Vec a4=acceleration(m,s4.x);
    for(int i=0;i<m.count;++i){s.x[i]+=(s.v[i]+s2.v[i]*2+s3.v[i]*2+s4.v[i])*(h/6);s.v[i]+=(a1[i]+a2[i]*2+a3[i]*2+a4[i])*(h/6);}return s;
}
std::vector<State> reference(const Case& c,double h,int steps,int subdivisions){
    std::vector<State> out;out.reserve(steps+1);out.push_back(c.initial);State s=c.initial;
    for(int j=0;j<steps;++j){for(int k=0;k<subdivisions;++k)s=reference_step(c.model,s,h/subdivisions);out.push_back(s);}return out;
}
}
