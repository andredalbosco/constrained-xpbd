#include "gradient_policies.hpp"
#include "frozen_schedules.hpp"
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>

using namespace energy_ordering;
namespace gp=gradient_policies;
namespace {
void require(bool ok,const std::string& message){if(!ok)throw std::runtime_error(message);}
void near(double actual,double expected,double tolerance,const std::string& message){
    if(!std::isfinite(actual)||!std::isfinite(expected)||std::abs(actual-expected)>tolerance)
        throw std::runtime_error(message+"; actual="+std::to_string(actual)+", expected="+std::to_string(expected));
}
void near(V2 actual,V2 expected,double tolerance,const std::string& message){near(norm(actual-expected),0,tolerance,message);}
void same(Iterate a,Iterate b,double tolerance,const std::string& message){
    near(a.x,b.x,tolerance,message+" x");near(a.lambda,b.lambda,tolerance,message+" lambda");near(a.eta,b.eta,tolerance,message+" eta");
}
Model elastic_pendulum(){Model m;m.pendulum=false;m.mass=1;m.length=1;m.alpha=.01;m.acceleration={0,-9.81};return m;}
double physical_energy(const Model& m,const State& s){
    const double r=std::hypot(s.x.x,s.x.y),extension=r-m.length;
    return .5*m.mass*(s.v.x*s.v.x+s.v.y*s.v.y)+extension*extension/(2*m.alpha)
        -m.mass*(m.acceleration.x*s.x.x+m.acceleration.y*s.x.y);
}
V2 physical_force(const Model& m,V2 x){
    const double r=std::hypot(x.x,x.y),coefficient=-(r-m.length)/(m.alpha*r);
    return {coefficient*x.x+m.mass*m.acceleration.x,coefficient*x.y+m.mass*m.acceleration.y};
}
double physical_motion_scale(const Context& c){
    const auto f=physical_force(c.model,c.old.x);
    return c.h*std::hypot(c.old.v.x,c.old.v.y)+c.h*c.h*std::hypot(f.x,f.y)/c.model.mass;
}
struct IndependentRows {V2 n,g,R;double rC,rE,A,B,D;};
IndependentRows independent_rows(const Context& c,Iterate z){
    const auto& m=c.model;IndependentRows r;
    const double radius=std::hypot(z.x.x,z.x.y),extension=radius-m.length,h2=c.h*c.h;
    r.n=z.x/radius;
    r.g=(z.x-c.old.x)*(m.mass/h2)-physical_force(m,z.x);
    r.rC=extension+m.alpha*z.lambda/h2;
    r.rE=physical_energy(m,{z.x,(z.x-c.old.x)/c.h})-c.target+c.beta*z.eta/h2;
    const V2 y=c.old.x+c.old.v*c.h+m.acceleration*h2;
    r.R=(z.x-y)*m.mass-r.n*z.lambda-r.g*z.eta;
    r.A=1/m.mass+m.alpha/h2;r.B=dot(r.n,r.g)/m.mass;r.D=dot(r.g,r.g)/m.mass+c.beta/h2;
    return r;
}

void energy_gradient_and_gravity(){
    Model m=elastic_pendulum();m.mass=1.7;m.alpha=.025;
    const State old{{.72,-1.13},{.31,-.24}};
    const Context c=context(m,old,.037,0,physical_energy(m,old));
    const Iterate z{{.69,-1.16},-.013,.002};
    const Rows actual=evaluate(c,z);const auto expected=independent_rows(c,z);
    require(actual.valid,"gravity derivative fixture invalid");
    near(actual.H,expected.rE,2e-13,"gravity and spring physical energy");
    near(actual.g,expected.g,4e-13,"gravity and spring energy gradient");
    near(actual.R,expected.R,2e-14,"gravity predictor and current stationarity");
    for(V2 direction:{V2{1,0},V2{0,1},V2{.6,.8}}){
        constexpr double eps=1e-6;
        auto H=[&](V2 x){return physical_energy(m,{x,(x-c.old.x)/c.h})-c.target;};
        const double finite_difference=(H(z.x+direction*eps)-H(z.x-direction*eps))/(2*eps);
        near(dot(actual.g,direction),finite_difference,3e-8,"full energy gradient finite difference");
    }
    // Gravitational potential decreases for a downward displacement.
    Model no_gravity=m;no_gravity.acceleration={};
    const double gravity_contribution=physical_energy(m,{z.x,{}})-physical_energy(no_gravity,{z.x,{}});
    near(gravity_contribution,m.mass*9.81*z.x.y,2e-14,"gravity potential sign");
}

void total_force_motion_scale(){
    Model m=elastic_pendulum();const State equilibrium{{0,-1.0981},{}};
    const Context eq=context(m,equilibrium,1./60,0,physical_energy(m,equilibrium));
    near(physical_force(m,equilibrium.x),{0,0},2e-14,"spring and gravity forces cancel at equilibrium");
    require(gp::motion_scale(eq)<1e-16,"equilibrium cutoff must use net force, not separate force magnitudes");
    const State unstretched{{0,-1},{}};const Context g=context(m,unstretched,1./60,0,physical_energy(m,unstretched));
    near(gp::motion_scale(g),(1./60)*(1./60)*9.81,2e-18,"gravity-only motion scale");
    m.mass=2.5;m.alpha=.02;
    const State oblique{{.72,-.96},{.3,-.4}};
    const Context c=context(m,oblique,.04,0,physical_energy(m,oblique));
    near(gp::motion_scale(c),physical_motion_scale(c),2e-16,"cutoff uses old total physical force and mass");
    const double separate=.04*.5+.04*.04*(norm(m.acceleration)+.2/(m.alpha*m.mass));
    require(std::abs(gp::motion_scale(c)-separate)>1e-3,"force fixture must detect summing magnitudes before cancellation");
}

void gravity_el_oracle(){
    Model m=elastic_pendulum();m.mass=1.3;m.alpha=.02;
    const State old{{.6,-1.05},{.4,-.2}};
    const Context c=context(m,old,.03,0,physical_energy(m,old));
    const Iterate start{c.predictor+V2{.001,-.002},-.006,.0004};
    const IndependentRows r=independent_rows(c,start);
    const double proposed=-r.rE/r.D;
    const double displacement=norm(r.g*(proposed/m.mass));
    require(displacement>1e-6,"EL oracle needs a nonzero energy correction");
    const double ratio=displacement/physical_motion_scale(c);
    for(double gamma:{-1.,ratio/2,ratio*2}){
        const bool reject=gamma>=0&&displacement>gamma*physical_motion_scale(c);
        const double de=reject?0:proposed;
        // Solve the independent upper-triangular 2x2 system: energy, then material.
        const double dl=(-r.rC-r.B*de)/r.A;
        const Iterate expected{start.x+(r.n*dl+r.g*de)/m.mass,start.lambda+dl,start.eta+de};
        Iterate actual=start;gp::History history;history.begin(c,actual);gp::Diagnostics d;
        gp::Options o;o.order=gp::Order::EL;o.gradient=gp::Gradient::Fresh;o.retain_rx=false;o.gamma=gamma;
        require(gp::update(c,actual,o,history,d),"gravity EL update failed");
        same(actual,expected,3e-13,"independent gravity EL scalar oracle");
        require(d.oversized==(reject?1:0)&&d.accepted==(reject?0:1),"gravity cutoff oracle decision");
        const V2 dx=actual.x-start.x;
        near(r.rC+dot(r.n,dx)+(m.alpha/(c.h*c.h))*dl,0,3e-13,"EL closes the frozen elastic material row");
        near(dx*m.mass-r.n*dl-r.g*de,{0,0},3e-13,"drop-rx increment equation");
    }
}

void static_equilibrium(){
    const Model m=elastic_pendulum();const State equilibrium{{0,-1.0981},{}};
    const double target=physical_energy(m,equilibrium);
    for(double h:{1./60,1./240}){
        const Context c=context(m,equilibrium,h,0,target);
        const Iterate root{equilibrium.x,-h*h*9.81,0};const Rows r=evaluate(c,root);
        // The gravity-centered implementation subtracts nearby positions before
        // dividing by h^2, so its roundoff grows at the smaller timestep.
        near(r.g,{0,0},1e-11,"static energy gradient");near(r.rC,0,2e-16,"static compliant material equation");
        near(r.H,0,0,"static physical energy equation");near(r.R,{0,0},3e-16,"static force stationarity");
        for(bool cutoff:{false,true}){
            State s=equilibrium;long rejects=0,incompatible=0;
            for(int step=0;step<static_cast<int>(std::lround(10/h));++step){
                const Context ctx=context(m,s,h,0,target);Iterate z{ctx.predictor};
                gp::History hist;hist.begin(ctx,z);gp::Diagnostics diag;Counters count;
                for(int sweep=0;sweep<5;++sweep){
                    if(cutoff){
                        gp::Options o;o.order=gp::Order::EL;o.gradient=gp::Gradient::Fresh;o.retain_rx=false;o.gamma=1;
                        require(gp::update(ctx,z,o,hist,diag),"equilibrium cutoff solve failed");
                    }else require(material_update(ctx,z,count),"equilibrium baseline solve failed");
                }
                hist.end(ctx,z);s={z.x,(z.x-s.x)/h};rejects+=diag.oversized;incompatible+=diag.core.incompatible;
                near(s.x,equilibrium.x,2e-12,"equilibrium position drift");near(s.v,{0,0},2e-10,"equilibrium velocity drift");
                near(physical_energy(m,s),target,2e-12,"equilibrium energy drift");
            }
            if(cutoff)require(rejects>0&&incompatible==0,"equilibrium uses zero-force cutoff without incompatible energy rows");
        }
    }
}

void vertical_baseline_formula(){
    const Model m=elastic_pendulum();const double radius_eq=1.0981,q0=.05,w0=-.2,frequency=10;
    const State eq{{0,-radius_eq},{}};const double Eeq=physical_energy(m,eq),excitation=.5*w0*w0+.5*100*q0*q0;
    for(double h:{1./60,1./240}){
        State s{{0,-(radius_eq+q0)},{0,-w0}};
        const int steps=static_cast<int>(std::lround(10/h));const double d=1+frequency*frequency*h*h,theta=std::atan(frequency*h);
        for(int n=1;n<=steps;++n){
            const Context c=context(m,s,h,0,physical_energy(m,s));Iterate z{c.predictor};Counters count;
            require(material_update(c,z,count),"vertical baseline update failed");
            s={z.x,(z.x-s.x)/h};
            const double attenuation=std::pow(d,-.5*n),angle=n*theta;
            const double q=attenuation*(q0*std::cos(angle)+(w0/frequency)*std::sin(angle));
            const double w=attenuation*(w0*std::cos(angle)-frequency*q0*std::sin(angle));
            near(s.x,{0,-(radius_eq+q)},3e-12,"gravity-shifted backward-Euler position formula");
            near(s.v,{0,-w},3e-11,"gravity-shifted backward-Euler velocity formula");
            near(physical_energy(m,s)-Eeq,excitation*std::pow(d,-n),3e-12,"gravity-shifted backward-Euler energy damping");
        }
    }
}

void no_gravity_parity(){
    Model m;m.mass=1.2;m.alpha=.015;const State old{{1.25,.2},{-.7,.3}};
    const Context c=context(m,old,1./60,0,physical_energy(m,old));const Iterate start{c.predictor};
    Iterate a=start,b=start;gp::History hist;hist.begin(c,a);gp::Diagnostics diag;Counters count;
    gp::Options o;o.order=gp::Order::EL;o.gradient=gp::Gradient::Fresh;o.retain_rx=false;o.gamma=-1;
    require(gp::update(c,a,o,hist,diag),"zero-gravity cutoff-policy reference failed");
    require(schedule_update(c,b,Schedule::EL,false,count),"original zero-gravity EL failed");
    same(a,b,2e-13,"no-gravity original EL parity");
}
} // namespace
int main(){
    try{
        energy_gradient_and_gravity();total_force_motion_scale();gravity_el_oracle();static_equilibrium();vertical_baseline_formula();no_gravity_parity();
        std::cout<<"PASS: elastic-gravity energy derivatives; static equilibrium; old net-force cutoff scale; independent EL oracle; "
                   "vertical backward-Euler damping; no-gravity parity\n";
    }catch(const std::exception& e){std::cerr<<"FAIL: "<<e.what()<<'\n';return 1;}
}
