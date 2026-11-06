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
void near(double actual,double expected,double tolerance,const std::string& message){
    if(!std::isfinite(actual)||!std::isfinite(expected)||std::abs(actual-expected)>tolerance)
        throw std::runtime_error(message+"; actual="+std::to_string(actual)+", expected="+std::to_string(expected));
}
void near(V2 actual,V2 expected,double tolerance,const std::string& message){near(norm(actual-expected),0,tolerance,message);}
void same(Iterate a,Iterate b,double tolerance,const std::string& message){
    near(a.x,b.x,tolerance,message+" position");near(a.lambda,b.lambda,tolerance,message+" lambda");near(a.eta,b.eta,tolerance,message+" eta");
}
double physical_energy(const Model& m,const State& s){
    const double extension=std::hypot(s.x.x,s.x.y)-m.length;
    return .5*m.mass*(s.v.x*s.v.x+s.v.y*s.v.y)
        +(m.pendulum?0:extension*extension/(2*m.alpha))
        -m.mass*(m.acceleration.x*s.x.x+m.acceleration.y*s.x.y);
}
struct Physical {V2 n,g,R;double rc,re,at,bt;};
Physical physical(const Context& c,const Iterate& z){
    const auto& m=c.model;const double radius=std::hypot(z.x.x,z.x.y),extension=radius-m.length,h2=c.h*c.h;
    Physical r;
    r.n=z.x/radius;r.g=(z.x-c.old.x)*(m.mass/h2)-m.acceleration*m.mass;
    if(!m.pendulum)r.g+=r.n*(extension/m.alpha);
    r.at=m.pendulum?0:m.alpha/h2;r.bt=c.beta/h2;
    r.rc=extension+r.at*z.lambda;r.re=physical_energy(m,{z.x,(z.x-c.old.x)/c.h})-c.target+r.bt*z.eta;
    const V2 y=c.old.x+c.old.v*c.h+m.acceleration*h2;
    r.R=(z.x-y)*m.mass-r.n*z.lambda-r.g*z.eta;
    return r;
}
Iterate oracle(const Context& c,Iterate initial,V2 used,gp::Order order,bool retain){
    const Physical r=physical(c,initial);const double m=c.model.mass;
    const V2 base=retain?r.R*(-1/m):V2{};
    const double matrix[2][2]={{dot(r.n,r.n)/m+r.at,dot(r.n,used)/m},
                               {dot(r.n,used)/m,dot(used,used)/m+r.bt}};
    const double rhs[2]={-r.rc-dot(r.n,base),-r.re-dot(used,base)};
    double d[2]={0,0};
    for(int visit=0;visit<2;++visit){
        const int i=order==gp::Order::LE?visit:1-visit;
        const double residual=rhs[i]-matrix[i][0]*d[0]-matrix[i][1]*d[1];
        d[i]+=residual/matrix[i][i];
    }
    initial.x+=base+(r.n*d[0]+used*d[1])/m;initial.lambda+=d[0];initial.eta+=d[1];
    return initial;
}
Context fixture(){
    Model m;m.mass=1.7;m.length=1.2;m.alpha=.035;m.acceleration={.2,-.3};
    const State old{{1.1,.7},{.15,-.12}};
    return context(m,old,.08,0,physical_energy(m,old)-.03);
}
const Iterate initial{{1.06,.76},-.014,.003};
Context exact_rest(){Model m;m.alpha=1./64;return context(m,{{1.25,0},{}},1./8,0,2);}

void fresh_parity(){
    const Context c=fixture();
    for(bool retain:{false,true})for(gp::Order order:{gp::Order::LE,gp::Order::EL}){
        Iterate a=initial,b=initial;gp::History h;h.begin(c,a);gp::Diagnostics d;Counters count;
        gp::Options o;o.order=order;o.retain_rx=retain;o.gradient=gp::Gradient::Fresh;o.gamma=-1;
        require(gp::update(c,a,o,h,d),"fresh policy update failed");
        require(schedule_update(c,b,order==gp::Order::LE?Schedule::LE:Schedule::EL,retain,count),"existing frozen update failed");
        same(a,b,2e-13,"unbounded fresh parity");
        same(a,oracle(c,initial,physical(c,initial).g,order,retain),2e-13,"independent fresh oracle");
        require(d.core.evaluations==1&&d.attempts==1&&d.accepted==1&&d.oversized==0,"fresh diagnostics");
    }
}

void stale_rows_and_true_residual(){
    const Context c=fixture();const Physical p=physical(c,initial);
    for(gp::Gradient policy:{gp::Gradient::LagSweep,gp::Gradient::FrozenStep,gp::Gradient::PreviousStep})
    for(bool retain:{false,true})for(gp::Order order:{gp::Order::LE,gp::Order::EL}){
        gp::History history;history.begin(c,initial);
        history.previous_sweep={3,-4};history.frozen={-2,5};history.previous_step={7,1};history.has_previous_step=true;
        const V2 chosen=policy==gp::Gradient::LagSweep?history.previous_sweep:
            (policy==gp::Gradient::FrozenStep?history.frozen:history.previous_step);
        const auto before_history=history;
        gp::Options o;o.order=order;o.gradient=policy;o.retain_rx=retain;o.gamma=-1;
        Iterate actual=initial;gp::Diagnostics d;
        require(gp::update(c,actual,o,history,d),"stale update failed");
        same(actual,oracle(c,initial,chosen,order,retain),3e-13,"stale approximate-column oracle");
        const V2 top=(actual.x-initial.x)*c.model.mass-p.n*(actual.lambda-initial.lambda)-chosen*(actual.eta-initial.eta);
        near(top,retain?p.R*(-1):V2{},5e-13,"top residual must use CURRENT physical energy gradient");
        near(history.previous_sweep,p.g,2e-14,"lagged snapshot must contain the physical PRE-update gradient");
        near(history.frozen,before_history.frozen,0,"frozen cache changed within timestep");
        near(history.previous_step,before_history.previous_step,0,"previous-step cache changed within timestep");
        require(d.stale_visits==1&&d.core.evaluations==1,"stale update count");
        if(retain)require(norm((p.g-chosen)*initial.eta)>1e-3,"fixture must detect substituting stale g in physical residual");
    }
}

void cache_lifecycle(){
    const Context c=fixture();gp::History h;
    h.begin(c,initial);const V2 start=physical(c,initial).g;
    near(h.frozen,start,2e-14,"begin frozen gradient");near(h.previous_sweep,start,2e-14,"begin lagged gradient");
    require(!h.has_previous_step,"begin must not invent previous-timestep data");
    Iterate changed=initial;changed.x+=V2{.01,-.02};const Rows changed_rows=evaluate(c,changed);
    near(gp::used_gradient(changed_rows,h,gp::Gradient::PreviousStep),changed_rows.g,0,"first-step previous-gradient fallback must be fresh");
    gp::Options o;o.order=gp::Order::EL;o.gradient=gp::Gradient::LagSweep;o.retain_rx=true;
    for(int visit=0;visit<3;++visit){
        const Iterate before=changed;const V2 used=h.previous_sweep;
        gp::Diagnostics d;
        require(gp::update(c,changed,o,h,d),"cache lifecycle update failed");
        same(changed,oracle(c,before,used,o.order,o.retain_rx),5e-12,"cache lifecycle oracle");
        near(h.previous_sweep,physical(c,before).g,2e-12,"lag sweep must retain pre-update gradient");
        near(h.frozen,start,2e-14,"frozen predictor gradient persisted");
    }
    h.end(c,changed);const V2 previous=physical(c,changed).g;
    require(h.has_previous_step,"end must mark previous timestep available");near(h.previous_step,previous,2e-12,"end must snapshot FINAL gradient");
    const State next{changed.x,(changed.x-c.old.x)/c.h};
    const Context next_c=context(c.model,next,c.h,0,c.target);const Iterate next_z{next_c.predictor};
    h.begin(next_c,next_z);const V2 next_prediction=physical(next_c,next_z).g;
    near(h.previous_step,previous,2e-12,"begin next step must preserve prior final gradient");
    near(h.frozen,next_prediction,2e-12,"next timestep resets frozen gradient");
    near(h.previous_sweep,next_prediction,2e-12,"next timestep resets lagged gradient");
    near(gp::used_gradient(evaluate(next_c,next_z),h,gp::Gradient::PreviousStep),previous,2e-12,"previous-step policy reads prior final cache");
}

void rejected_energy_material_base(){
    const Context c=fixture();const Physical r=physical(c,initial);
    for(bool retain:{false,true})for(gp::Order order:{gp::Order::LE,gp::Order::EL})
    for(gp::Gradient policy:{gp::Gradient::Fresh,gp::Gradient::FrozenStep}){
        gp::History h;h.begin(c,initial);h.frozen={3,-4};gp::Diagnostics d;
        gp::Options o;o.order=order;o.gradient=policy;o.retain_rx=retain;o.gamma=0;
        Iterate actual=initial;require(gp::update(c,actual,o,h,d),"rejection update failed");
        const V2 base=retain?r.R*(-1/c.model.mass):V2{};
        const double dl=(-r.rc-dot(r.n,base))/(1/c.model.mass+r.at);
        const Iterate expected{initial.x+base+r.n*(dl/c.model.mass),initial.lambda+dl,initial.eta};
        same(actual,expected,3e-13,"reject only energy; apply material and position base once");
        require(d.oversized==1&&d.accepted==0&&d.core.scalar_solves==1,"rejected energy must not count as accepted/solved");
        near(d.max_accepted_ratio,0,0,"rejected ratio must not enter accepted diagnostics");
    }
}

void exact_boundaries(){
    const Context c=exact_rest();const Iterate start{c.predictor};
    near(gp::motion_scale(c),.25,0,"binary old-force motion scale");
    for(double gamma:{std::nextafter(.5,0.),.5,std::nextafter(.5,1.)}){
        gp::Options o;o.order=gp::Order::LE;o.gradient=gp::Gradient::Fresh;o.retain_rx=true;o.gamma=gamma;
        gp::History h;h.begin(c,start);gp::Diagnostics d;Iterate actual=start;
        require(gp::update(c,actual,o,h,d),"bound boundary update failed");
        const bool rejected=gamma<.5;
        near(d.max_candidate,.125,0,"exact candidate energy displacement");near(d.max_ratio,.5,0,"exact candidate ratio");
        require(d.oversized==(rejected?1:0)&&d.accepted==(rejected?0:1),"strict greater-than boundary policy");
        near(actual.x,{rejected?1.125:1.25,0},0,"boundary accepted/rejected position");
        near(actual.eta,rejected?0:1./128,0,"boundary accepted/rejected multiplier");
    }
    // A zero displacement is accepted even with gamma=0.
    gp::Options o;o.order=gp::Order::EL;o.gamma=0;gp::History h;h.begin(c,start);gp::Diagnostics d;Iterate z=start;
    require(gp::update(c,z,o,h,d),"zero candidate update failed");
    require(d.accepted==1&&d.oversized==0,"zero candidate must be accepted at zero bound");
    near(d.max_candidate,0,0,"zero candidate displacement");
    // With no old velocity or force, the bound is zero for every finite gamma.
    Model m;m.alpha=1./64;const Context still=context(m,{{1,0},{}},1./8,0,0);const Iterate displaced{{1.125,0}};
    near(gp::motion_scale(still),0,0,"zero physical motion scale");
    o.gamma=100;o.order=gp::Order::EL;h={};h.begin(still,displaced);d={};z=displaced;
    require(gp::update(still,z,o,h,d),"zero scale update failed");
    require(d.oversized==1&&d.accepted==0,"nonzero candidate at zero scale must reject");
    near(d.max_ratio,1e300,0,"diagnostic cap for zero scale");
}

void dangerous_second_step(){
    Model m;m.alpha=1./64;const State old{{1.125,0},{-1,0}};const Context c=context(m,old,1./8,0,2);
    const Iterate start{c.predictor};const Physical p=physical(c,start);
    near(start.x,{1,0},0,"second-step predictor");near(p.g,{-8,0},0,"second-step energy gradient");near(p.re,-1.5,0,"second-step energy residual");
    for(gp::Order order:{gp::Order::LE,gp::Order::EL})for(bool retain:{false,true}){
        gp::Options o;o.order=order;o.retain_rx=retain;o.gamma=-1;gp::History h;h.begin(c,start);gp::Diagnostics d;Iterate free=start;
        require(gp::update(c,free,o,h,d),"unbounded second step failed");
        near(d.max_candidate,1.5*.125,0,"second-step unbounded -1.5q1 energy contribution");
        near(free.x.x-start.x.x,order==gp::Order::LE?-.1875:-.09375,0,
             "EL material-last correction partly cancels the energy contribution");
        o.gamma=.5;h={};h.begin(c,start);d={};Iterate bounded=start;
        require(gp::update(c,bounded,o,h,d),"bounded second step failed");
        same(bounded,start,0,"reject second-step energy while zero material row remains unchanged");
        require(d.oversized==1&&d.accepted==0,"second-step dangerous increment must reject");
        near(d.max_ratio,.75,0,"second-step candidate is 75 percent of physical scale");
    }
}

void singular_gradient(){
    const Context c=exact_rest();Iterate minimum{c.predictor};Counters temp;
    require(material_update(c,minimum,temp),"minimum fixture material update failed");
    near(physical(c,minimum).g,{0,0},0,"minimum fixture has zero fresh energy gradient");
    gp::Options o;o.order=gp::Order::EL;o.gradient=gp::Gradient::Fresh;o.gamma=.5;
    gp::History h;h.begin(c,minimum);gp::Diagnostics d;Iterate result=minimum;
    require(gp::update(c,result,o,h,d),"singular gradient must explicitly skip");
    same(result,minimum,0,"zero-gradient minimum should remain unchanged");
    require(d.attempts==1&&d.core.singular==1&&d.core.incompatible==1&&d.accepted==0&&d.oversized==0,
            "singular rows must be distinguished from displacement rejection");
    // A cached zero gradient is skipped even when the current physical gradient is nonzero.
    const Context general=fixture();h={};h.begin(general,initial);h.frozen={};d={};result=initial;
    o.gradient=gp::Gradient::FrozenStep;o.retain_rx=true;o.gamma=-1;
    require(gp::update(general,result,o,h,d),"cached-zero gradient update failed");
    require(d.core.singular==1&&d.accepted==0&&d.oversized==0,"cached-zero gradient policy");
    near(result.eta,initial.eta,0,"cached-zero energy row must preserve multiplier");
}
} // namespace
int main(){
    try{
        fresh_parity();stale_rows_and_true_residual();cache_lifecycle();rejected_energy_material_base();
        exact_boundaries();dangerous_second_step();singular_gradient();
        std::cout<<"PASS: original endpoint-gradient policy parity; stale linearization/current residual; cache lifecycle; rejection-only energy; "
                   "exact displacement boundaries; second-step overshoot guard; singular gradient classification\n";
    }catch(const std::exception& e){std::cerr<<"FAIL: "<<e.what()<<'\n';return 1;}
}
