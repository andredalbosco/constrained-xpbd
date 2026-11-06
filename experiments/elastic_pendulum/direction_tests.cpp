#include "gradient_policies.hpp"
#include "frozen_schedules.hpp"
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>

using namespace energy_ordering;
namespace gp=gradient_policies;
namespace {
void require(bool ok,const std::string& message){if(!ok)throw std::runtime_error(message);}
void near(double a,double b,double tolerance,const std::string& message){
    if(!std::isfinite(a)||!std::isfinite(b)||std::abs(a-b)>tolerance)
        throw std::runtime_error(message+"; actual="+std::to_string(a)+", expected="+std::to_string(b));
}
void near(V2 a,V2 b,double tolerance,const std::string& message){near(norm(a-b),0,tolerance,message);}
void same(Iterate a,Iterate b,double tolerance,const std::string& message){
    near(a.x,b.x,tolerance,message+" x");near(a.lambda,b.lambda,tolerance,message+" lambda");near(a.eta,b.eta,tolerance,message+" eta");
}
double physical_energy(const Model& m,const State& s){
    const double c=std::hypot(s.x.x,s.x.y)-m.length;
    return .5*m.mass*(s.v.x*s.v.x+s.v.y*s.v.y)+(m.pendulum?0:c*c/(2*m.alpha))
        -m.mass*(m.acceleration.x*s.x.x+m.acceleration.y*s.x.y);
}
double objective(const Context& c,V2 x){
    const V2 y=c.old.x+c.old.v*c.h+c.model.acceleration*(c.h*c.h),d=x-y;
    const double extension=std::hypot(x.x,x.y)-c.model.length;
    return .5*c.model.mass*dot(d,d)/(c.h*c.h)+(c.model.pendulum?0:extension*extension/(2*c.model.alpha));
}
struct Independent {V2 n,g,q,R;double rc,H,A,B,D;};
Independent rows(const Context& c,Iterate z){
    const auto& m=c.model;const double radius=std::hypot(z.x.x,z.x.y),h2=c.h*c.h,extension=radius-m.length;
    const V2 y=c.old.x+c.old.v*c.h+m.acceleration*h2;
    Independent r;r.n=z.x/radius;
    const V2 elastic=m.pendulum?V2{}:r.n*(extension/m.alpha);
    r.q=(z.x-y)*(m.mass/h2)+elastic;
    r.g=(z.x-c.old.x)*(m.mass/h2)+elastic-m.acceleration*m.mass;
    r.R=(z.x-y)*m.mass-r.n*z.lambda-r.g*z.eta;
    const double at=m.pendulum?0:m.alpha/h2;
    r.rc=extension+at*z.lambda;r.H=physical_energy(m,{z.x,(z.x-c.old.x)/c.h})-c.target+c.beta*z.eta/h2;
    r.A=1/m.mass+at;r.B=dot(r.n,r.g)/m.mass;r.D=dot(r.g,r.g)/m.mass+c.beta/h2;
    return r;
}
gp::Options options(bool gate=true){gp::Options o;o.order=gp::Order::EL;o.gradient=gp::Gradient::Fresh;o.retain_rx=false;o.gamma=-1;o.objective_gate=gate;return o;}
Iterate apply(const Context& c,Iterate start,const gp::Options& o,gp::Diagnostics& d){
    gp::History h;h.begin(c,start);require(gp::update(c,start,o,h,d),"direction policy update failed");return start;
}
Context rest_context(){Model m;m.alpha=1./64;return context(m,{{1.25,0},{}},1./8,0,2);}

void objective_derivatives(){
    for(bool hard:{false,true}){
        Model m;m.pendulum=hard;m.mass=1.7;m.length=1.2;m.alpha=.025;m.acceleration={.3,-9.81};
        const State old{{.72,-1.13},{.31,-.24}};const Context c=context(m,old,.037,0,physical_energy(m,old));
        const Iterate z{{.69,-1.16},-.013,.002};const auto independent=rows(c,z);
        const V2 actual=gp::objective_gradient(c,z,evaluate(c,z));
        near(actual,independent.q,2e-13,"full implicit objective gradient");
        near(actual,independent.g-old.v*(m.mass/c.h),4e-13,"objective/energy gradient identity");
        for(V2 direction:{V2{1,0},V2{0,1},V2{.6,.8}}){
            constexpr double eps=1e-6;
            const double fd=(objective(c,z.x+direction*eps)-objective(c,z.x-direction*eps))/(2*eps);
            near(dot(actual,direction),fd,3e-8,"objective finite-difference derivative");
        }
        Iterate changed=z;changed.lambda+=10;changed.eta-=3;
        near(gp::objective_gradient(c,changed,evaluate(c,changed)),actual,0,
             "objective gradient must not contain constraint or energy multipliers");
    }
}

void exact_uphill_and_downhill(){
    const Context c=rest_context();
    for(double x:{1.0625,1.1875,1.375}){
        const Iterate start{{x,0}};const auto r=rows(c,start);const double de=-r.H/r.D;
        const V2 dxE=r.g*de;const double work=dot(r.q,dxE);const bool reject=work>0;
        require(std::abs(work)>.1,"binary direction fixture must have a nonzero slope");
        const double accepted_de=reject?0:de,dl=(-r.rc-r.B*accepted_de)/r.A;
        const Iterate expected{start.x+r.n*dl+r.g*accepted_de,dl,accepted_de};
        gp::Diagnostics d;const Iterate actual=apply(c,start,options(),d);
        same(actual,expected,0,"binary uphill/downhill energy gate");
        require(d.directional_skips==(reject?1:0)&&d.accepted==(reject?0:1)&&d.oversized==0,
                "direction rejection must have its own counter");
        require(d.core.scalar_solves==(reject?1:2),"material solve always runs after directional decision");
        near(r.rc+dot(r.n,actual.x-start.x)+(c.model.alpha/(c.h*c.h))*dl,0,0,
             "EL material row uses accepted energy multiplier");
        if(x==1.0625){
            const V2 elastic=r.n*((x-1)/c.model.alpha);
            require(dot(elastic,dxE)<0&&work>0,
                    "fixture must distinguish full objective from elastic force alone");
        }
        if(!reject){
            near(actual.eta,de,0,"downhill energy increment must be accepted whole, without amplitude clipping");
            require(dot(r.q,r.g*actual.eta)<0,"accepted downhill direction");
        }
        near(d.max_accepted_objective_work,0,0,"accepted direction cannot have positive objective work");
    }
}

void tangent_strict_boundary(){
    Model m;m.pendulum=true;m.alpha=0;
    for(double vy:{std::nextafter(1.,0.),1.,std::nextafter(1.,2.)}){
        const State old{{1,0},{-1,vy}};const Context c=context(m,old,1,0,1);const Iterate start{{1,1}};
        const auto r=rows(c,start);const double de=-r.H/r.D;const V2 dxE=r.g*de;
        near(dxE,{0,.5},0,"nonzero tangent energy direction");
        const double slope=dot(r.q,dxE);const bool reject=vy<1;
        require((slope>0)==reject,"strict directional boundary fixture");
        gp::Diagnostics d;const Iterate actual=apply(c,start,options(),d);
        require(d.directional_skips==(reject?1:0)&&d.accepted==(reject?0:1),
                "positive tiny slope rejected; equality and negative slope accepted");
        near(actual.eta,reject?0:.5,0,"strict tangent decision multiplier");
        if(vy==1){
            near(r.q,{1,0},0,"tangent fixture has nonzero objective gradient");
            near(slope,0,0,"tangent fixture exact zero slope");
            require(objective(c,start.x+dxE)>objective(c,start.x),
                    "linearized gate deliberately permits tangent moves with second-order objective increase");
        }
    }
}

void zero_objective_gradient(){
    Model m;m.alpha=1./64;const Context c=context(m,{{1.125,0},{-1,0}},1./8,0,2);const Iterate start{c.predictor};
    const auto r=rows(c,start);near(r.q,{0,0},0,"second-step objective minimum");near(r.g,{-8,0},0,"second-step nonzero energy gradient");
    gp::Diagnostics d;const Iterate actual=apply(c,start,options(),d);
    same(actual,{{.90625,0},.09375,.0234375},0,"zero-gradient gate permits full energy increment and subsequent material response");
    near(d.max_candidate,.1875,0,"zero-gradient case retains the large energy candidate");
    require(d.accepted==1&&d.directional_skips==0&&d.oversized==0,"zero objective gradient is allowed under strict positive-slope gate");
    require(objective(c,actual.x)>objective(c,start.x),"first-order neutral acceptance is not a nonlinear objective line search");
    // Zero ENERGY correction is also accepted, regardless of the objective slope.
    const Context rest=rest_context();d={};const Iterate no_energy=apply(rest,{rest.predictor},options(),d);
    same(no_energy,{{1.125,0},-.125,0},0,"zero energy correction retains the material solve");
    require(d.accepted==1&&d.directional_skips==0,"zero energy increment is neutral");
}

void gravity_vector_gate(){
    Model m;m.mass=1.3;m.alpha=.02;m.acceleration={0,-9.81};
    const State old{{.6,-1.05},{.4,-.2}};Context c=context(m,old,.03,0,physical_energy(m,old));
    const Iterate start{c.predictor+V2{.001,-.002},-.006,.0004};const auto base_rows=rows(c,start);
    const double qg=dot(base_rows.q,base_rows.g);require(std::abs(qg)>1,"gravity fixture needs a well-defined direction");
    for(bool retain:{false,true})for(bool uphill:{false,true}){
        const V2 base=retain?base_rows.R*(-1/m.mass):V2{};
        const double numerator=(qg>0?1.:-1.)*(uphill?1.:-1.);
        const double current_energy=physical_energy(m,{start.x,(start.x-c.old.x)/c.h});
        c.target=current_energy+dot(base_rows.g,base)+numerator;
        const auto r=rows(c,start);const double raw_de=(-r.H-dot(r.g,base))/r.D;
        const V2 raw=r.g*(raw_de/m.mass);require((dot(r.q,raw)>0)==uphill,"gravity gate fixture direction");
        auto o=options();o.retain_rx=retain;gp::Diagnostics d;
        const Iterate actual=apply(c,start,o,d);
        const double de=uphill?0:raw_de,dl=(-r.rc-dot(r.n,base)-r.B*de)/r.A;
        const Iterate expected{start.x+base+(r.n*dl+r.g*de)/m.mass,start.lambda+dl,start.eta+de};
        same(actual,expected,3e-13,"independent gravity vector gate and material solve");
        require(d.directional_skips==(uphill?1:0)&&d.accepted==(uphill?0:1)&&d.oversized==0,"gravity directional diagnostics");
        near((actual.x-start.x)*m.mass-r.n*dl-r.g*de,retain?r.R*(-1):V2{},4e-13,"retained base applied once after directional decision");
        if(!uphill)near((actual.eta-start.eta),raw_de,2e-16,"downhill gravity energy correction is unscaled");
    }
}

void gate_off_legacy_parity(){
    const Context c=rest_context();
    require(!gp::Options{}.objective_gate,"legacy default must leave directional gate off");
    for(Iterate start:{Iterate{{1.0625,0}},Iterate{{1.375,0},-.01,.002}})
    for(bool retain:{false,true})for(gp::Order order:{gp::Order::LE,gp::Order::EL}){
        auto o=options(false);o.retain_rx=retain;o.order=order;gp::Diagnostics d;
        const Iterate actual=apply(c,start,o,d);Iterate expected=start;Counters count;
        require(schedule_update(c,expected,order==gp::Order::EL?Schedule::EL:Schedule::LE,retain,count),"legacy frozen comparison failed");
        same(actual,expected,3e-13,"gate-off fresh unbounded parity");
        require(d.directional_skips==0&&d.oversized==0,"gate-off path must not reject by direction or size");
    }
    gp::Diagnostics d;apply(c,{{1.0625,0}},options(false),d);
    require(d.max_accepted_objective_work>0,"ungated diagnostic must expose accepted uphill directions");
}
} // namespace
int main(){
    try{
        objective_derivatives();exact_uphill_and_downhill();tangent_strict_boundary();zero_objective_gradient();gravity_vector_gate();gate_off_legacy_parity();
        std::cout<<"PASS: full objective finite differences; uphill rejection/downhill acceptance; strict tangent/zero boundaries; "
                   "material response; gravity vector oracle; legacy gate-off parity\n";
    }catch(const std::exception& e){std::cerr<<"FAIL: "<<e.what()<<'\n';return 1;}
}
