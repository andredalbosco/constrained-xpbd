#include "small_steps.hpp"
#include <complex>
#include <cstring>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>

using namespace energy_ordering;
namespace ss=small_steps_el;
namespace gp=gradient_policies;
namespace {
void require(bool ok,const std::string& message){if(!ok)throw std::runtime_error(message);}
void near(double a,double b,double tolerance,const std::string& message){
    if(!std::isfinite(a)||!std::isfinite(b)||std::abs(a-b)>tolerance)
        throw std::runtime_error(message+"; actual="+std::to_string(a)+", expected="+std::to_string(b));
}
void near(V2 a,V2 b,double tolerance,const std::string& message){near(norm(a-b),0,tolerance,message);}
void same(Iterate a,Iterate b,double tolerance,const std::string& message){
    near(a.x,b.x,tolerance,message+" position");near(a.lambda,b.lambda,tolerance,message+" lambda");near(a.eta,b.eta,tolerance,message+" eta");
}
double physical_energy(const Model& m,const State& s){
    const double extension=std::hypot(s.x.x,s.x.y)-m.length;
    return .5*m.mass*dot(s.v,s.v)+(m.pendulum?0:extension*extension/(2*m.alpha))-m.mass*dot(m.acceleration,s.x);
}
double physical_scale(const Context& c){
    const auto& m=c.model;V2 acceleration=m.acceleration;
    if(!m.pendulum){const double r=std::hypot(c.old.x.x,c.old.x.y);acceleration+=c.old.x*(-(r-m.length)/(m.mass*m.alpha*r));}
    return c.h*norm(c.old.v)+c.h*c.h*norm(acceleration);
}
struct Independent {V2 n,g,q;double rc,re,at,bt,A,B,D;};
Independent rows(const Context& c,const Iterate& z){
    const auto& m=c.model;const double h2=c.h*c.h,radius=std::hypot(z.x.x,z.x.y),extension=radius-m.length;
    Independent r;r.n=z.x/radius;
    const V2 elastic=m.pendulum?V2{}:r.n*(extension/m.alpha);
    r.q=(z.x-(c.old.x+c.old.v*c.h+m.acceleration*h2))*(m.mass/h2)+elastic;
    r.g=(z.x-c.old.x)*(m.mass/h2)+elastic-m.acceleration*m.mass;
    r.at=m.pendulum?0:m.alpha/h2;r.bt=c.beta/h2;
    r.rc=extension+r.at*z.lambda;r.re=physical_energy(m,{z.x,(z.x-c.old.x)/c.h})-c.target+r.bt*z.eta;
    r.A=1/m.mass+r.at;r.B=dot(r.n,r.g)/m.mass;r.D=dot(r.g,r.g)/m.mass+r.bt;
    return r;
}
Iterate apply(const Context& c,Iterate z,const ss::Options& o,ss::Diagnostics& d){
    require(ss::update(c,z,o,d),"small-step update failed");return z;
}
void check_increment(const Context& c,Iterate start,Iterate actual,double de,double tolerance,const std::string& message){
    const auto r=rows(c,start);const double dl=(-r.rc-r.B*de)/r.A;
    const Iterate expected{start.x+(r.n*dl+r.g*de)/c.model.mass,start.lambda+dl,start.eta+de};
    same(actual,expected,tolerance,message);
    const V2 dx=actual.x-start.x;
    near(dx*c.model.mass-r.n*(actual.lambda-start.lambda)-r.g*(actual.eta-start.eta),{},tolerance*4,message+" dropped top row");
    near(r.rc+dot(r.n,dx)+r.at*(actual.lambda-start.lambda),0,tolerance*4,message+" linearized material row closes");
}
Context rest_context(){Model m;m.alpha=1./64;return context(m,{{1.25,0},{}},1./8,0,2);}

void raw_legacy_parity(){
    require(ss::Options{}.gate==ss::Gate::Raw&&ss::Options{}.relaxation==1&&ss::Options{}.cap<0,"default controls preserve raw gate");
    for(bool gravity:{false,true})for(double beta:{0.,.005}){
        Model m;m.alpha=.02;m.mass=1.3;if(gravity)m.acceleration={0,-9.81};
        const State old{{.6,-1.05},{.4,-.2}};const Context c=context(m,old,.03,beta,physical_energy(m,old));
        for(Iterate start:{Iterate{c.predictor},Iterate{c.predictor+V2{.003,-.002},-.006,.0004}}){
            for(ss::Gate gate:{ss::Gate::None,ss::Gate::Raw}){
                ss::Options o;o.gate=gate;ss::Diagnostics d;const Iterate actual=apply(c,start,o,d);
                gp::Options legacy;legacy.order=gp::Order::EL;legacy.gradient=gp::Gradient::Fresh;legacy.retain_rx=false;legacy.gamma=-1;legacy.objective_gate=gate==ss::Gate::Raw;
                gp::History history;history.begin(c,start);gp::Diagnostics old_d;Iterate expected=start;
                require(gp::update(c,expected,legacy,history,old_d),"legacy comparison update failed");
                same(actual,expected,0,"unrelaxed uncapped legacy EL parity");
                require(d.core.evaluations==1&&d.attempts==1,"one frozen evaluation and energy attempt per substep");
                require(d.accepted==old_d.accepted&&d.directional_skips==old_d.directional_skips,"legacy gate decisions preserved");
                require(d.core.scalar_solves==2,"count candidate division even when directional gate rejects");
            }
        }
    }
}

void relaxation_and_compliance(){
    // A deliberately nonzero eta tests beta/h^2 in BOTH residual and denominator.
    Model m;m.alpha=1./64;
    for(double h:{1./8,1./16}){
        const Context c=context(m,{{1.25,0},{-.5,0}},h,1./32,2);
        const Iterate start{{1.1875,0},-.015625,.03125};const auto r=rows(c,start);
        const double raw=-r.re/r.D;
        require(std::abs(r.bt*start.eta)>.01,"compliance fixture must exercise accumulated energy multiplier");
        near(r.bt,c.beta/(h*h),0,"physical beta is fixed across substep sizes");
        for(double omega:{0.,.125,.5,1.}){
            ss::Options o;o.gate=ss::Gate::None;o.relaxation=omega;ss::Diagnostics d;
            const Iterate actual=apply(c,start,o,d);const double de=omega*raw;
            check_increment(c,start,actual,de,2e-14,"relaxed compliant EL multiplier consistency");
            near(r.re+r.D*(actual.eta-start.eta),(1-omega)*r.re,3e-14,"relaxation leaves expected first energy-row residual");
            near(norm(r.g*((actual.eta-start.eta)/m.mass)),omega*norm(r.g*(raw/m.mass)),2e-14,"eta and raw energy displacement use same relaxation");
            require(d.accepted==1&&d.capped==0&&d.directional_skips==0,"relaxation acceptance diagnostics");
        }
    }
}

void cap_is_scaling_not_rejection(){
    const Context c=rest_context();const Iterate start{{1.375,0}};const auto r=rows(c,start);
    const double raw=-r.re/r.D,raw_size=norm(r.g*(raw/c.model.mass)),scale=physical_scale(c);
    near(scale,.25,0,"cap scale uses OLD physical elastic force");
    require(dot(r.q,r.g*raw)<0,"cap fixture must pass raw gate");
    for(double fraction:{0.,.25,.5,1.,2.}){
        ss::Options o;o.cap=fraction*raw_size/scale;ss::Diagnostics d;
        const Iterate actual=apply(c,start,o,d);const double de=raw*std::min(1.,fraction);
        check_increment(c,start,actual,de,2e-14,"capped energy correction and material response");
        near(norm(r.g*(actual.eta-start.eta)),std::min(raw_size,o.cap*scale),2e-15,"cap displacement magnitude");
        require(d.accepted==1&&d.directional_skips==0&&d.capped==(fraction<1?1:0),"clipping is accepted scaling, equality is uncapped");
        if(fraction>0)require((actual.eta-start.eta)*raw>0,"cap preserves energy direction and nonzero multiplier");
    }
    // Cap acts AFTER relaxation: .25 raw fits below .5 raw cap without clipping.
    ss::Options o;o.relaxation=.25;o.cap=.5*raw_size/scale;ss::Diagnostics d;
    const Iterate relaxed=apply(c,start,o,d);check_increment(c,start,relaxed,.25*raw,2e-14,"relaxation precedes cap");
    require(d.capped==0,"cap compares the relaxed energy displacement");
    // Gate acts BEFORE both controls: even omega=0/cap=0 retains its decision.
    const Iterate uphill{{1.0625,0}};o.relaxation=0;o.cap=0;d={};
    const Iterate rejected=apply(c,uphill,o,d);check_increment(c,uphill,rejected,0,2e-14,"raw uphill gate before damping and cap");
    require(d.directional_skips==1&&d.accepted==0&&d.capped==0,"raw direction tested before zero controls");
}

void effective_gate_and_response(){
    // Hard rod removes the energy correction's normal component from the final
    // displacement. Raw and effective slopes can therefore have opposite signs.
    Model m;m.pendulum=true;m.alpha=0;
    for(bool raw_downhill:{false,true}){
        const State old{{0,-1},raw_downhill?V2{3,0}:V2{-1,2}};
        const Context c=context(m,old,1,0,3);const Iterate start{{1,0}};const auto r=rows(c,start);
        const double raw=-r.re/r.D;near(raw,1,0,"opposite-sign gate fixture multiplier");
        near(r.g,{1,1},0,"opposite-sign gate fixture energy gradient");
        const V2 raw_dx=r.g*raw,effective=(r.g-r.n*(r.B/r.A))*raw;
        require((dot(r.q,raw_dx)<0)==raw_downhill,"raw objective slope fixture");
        require((dot(r.q,effective)>0)==raw_downhill,"effective objective slope reverses raw fixture");
        ss::Options zero;zero.gate=ss::Gate::None;zero.relaxation=0;ss::Diagnostics d0;
        const Iterate material_only=apply(c,start,zero,d0);
        for(ss::Gate gate:{ss::Gate::None,ss::Gate::Raw,ss::Gate::Effective}){
            const bool rejected=(gate==ss::Gate::Raw&&!raw_downhill)||(gate==ss::Gate::Effective&&raw_downhill);
            ss::Options o;o.gate=gate;ss::Diagnostics d;const Iterate actual=apply(c,start,o,d);
            check_increment(c,start,actual,rejected?0:raw,0,"raw versus effective direction decision");
            near(actual.x-material_only.x,rejected?V2{}:effective,0,"effective direction equals response relative to material-only update");
            require(d.directional_skips==(rejected?1:0)&&d.accepted==(rejected?0:1),"effective gate counters");
            if(gate!=ss::Gate::None)near(d.max_accepted_objective_work,0,0,"accepted tested direction nonpositive");
        }
    }
    // Exact zero objective gradient is accepted by either gate.
    Model spring;spring.alpha=1./64;const Context c=context(spring,{{1.125,0},{-1,0}},1./8,0,2);
    for(ss::Gate gate:{ss::Gate::Raw,ss::Gate::Effective}){
        ss::Options o;o.gate=gate;ss::Diagnostics d;const Iterate actual=apply(c,{c.predictor},o,d);
        same(actual,{{.90625,0},.09375,.0234375},0,"neutral objective gradient allows energy correction");
        require(d.accepted==1&&d.directional_skips==0,"neutral objective slope accepted");
    }
}

void predictor_reset_gate_equivalence(){
    // At the actual reset predictor, q = n*C/alpha. Thus removing the
    // material row's response scales the objective slope by at/A > 0.
    const State states[]={{{1.25,0},{}},{{0,-1.2},{.7,0}},{{.8,-1.05},{-.4,.3}},{{.72,-.54},{.2,-.5}}};
    long compared=0,uphill=0,downhill=0;
    auto bit_equal=[](double a,double b){return std::memcmp(&a,&b,sizeof(double))==0;};
    for(double mass:{.25,1.,3.})for(double alpha:{.001,.01,.08})
    for(double h:{1./8,1./240,1./480,1./960})for(bool gravity:{false,true})
    for(double beta:{0.,.002})for(const State& old:states){
        Model m;m.mass=mass;m.alpha=alpha;if(gravity)m.acceleration={0,-9.81};
        const Context c=context(m,old,h,beta,physical_energy(m,old));const Iterate start{c.predictor};
        require(start.lambda==0&&start.eta==0,"actual substep multipliers reset");
        const auto r=rows(c,start);const double de=-r.re/r.D;
        const V2 raw_dx=r.g*(de/m.mass),effective=(r.g-r.n*(r.B/r.A))*(de/m.mass);
        const double raw_slope=dot(r.q,raw_dx),effective_slope=dot(r.q,effective),factor=r.at/r.A;
        const double magnitude=std::max(1.,norm(r.q)*norm(raw_dx));
        near(r.q,r.n*((norm(start.x)-m.length)/m.alpha),0,"reset predictor objective gradient is radial");
        require(factor>0&&factor<1,"compliant material response gives strictly positive slope factor");
        near(effective_slope,factor*raw_slope,128*std::numeric_limits<double>::epsilon()*magnitude,
             "predictor/reset effective slope equals alpha_tilde/A times raw slope");
        if(std::abs(raw_slope)<1e-9*magnitude)continue;
        require((raw_slope>0)==(effective_slope>0),"raw and effective predictor gate decisions agree away from zero");
        for(double cap:{-1.,.01}){
            ss::Options raw;raw.gate=ss::Gate::Raw;raw.cap=cap;raw.relaxation=cap<0?1:.25;
            ss::Options effective_option=raw;effective_option.gate=ss::Gate::Effective;
            ss::Diagnostics dr,deff;const Iterate a=apply(c,start,raw,dr),b=apply(c,start,effective_option,deff);
            require(bit_equal(a.x.x,b.x.x)&&bit_equal(a.x.y,b.x.y)&&bit_equal(a.lambda,b.lambda)&&bit_equal(a.eta,b.eta),
                    "actual predictor/reset raw/effective outputs are bit-exact");
            require(dr.accepted==deff.accepted&&dr.directional_skips==deff.directional_skips&&dr.capped==deff.capped,
                    "predictor/reset gate and cap counters identical");
            require(dr.core.evaluations==1&&deff.core.evaluations==1,"one frozen update in predictor equivalence fixture");
            ++compared;if(raw_slope>0)++uphill;else ++downhill;
        }
    }
    require(compared>500&&uphill>0&&downhill>0,"predictor equivalence covers substantial accepted and rejected cases");
}

void static_equilibrium(){
    // Binary fixture isolates real equilibrium drift from decimal roundoff.
    Model m;m.alpha=1./64;m.acceleration={0,-8};const State old{{0,-1.125},{}};
    const Context c=context(m,old,1./8,0,physical_energy(m,old));
    near(physical_scale(c),0,0,"old equilibrium has zero velocity and zero NET force");
    near(gp::motion_scale(c),0,0,"production scale uses gravity plus elastic force");
    for(ss::Gate gate:{ss::Gate::None,ss::Gate::Raw,ss::Gate::Effective}){
        for(double cap:{0.,.1,1.}){
            ss::Options o;o.gate=gate;o.cap=cap;ss::Diagnostics d;const Iterate actual=apply(c,{c.predictor},o,d);
            near(actual.x,old.x,0,"finite cap preserves static equilibrium");near(actual.eta,0,0,"zero motion cap zeroes eta correction");
            require(d.capped==1&&d.accepted==1,"zero-motion cap is a fully scaled accepted correction");
        }
        ss::Options o;o.gate=gate;ss::Diagnostics d;const Iterate actual=apply(c,{c.predictor},o,d);
        near(actual.x,{0,-1.09375},0,"uncapped EL equilibrium drift is explicit");
        require(d.accepted==1&&d.directional_skips==0,"both direction gates allow equilibrium's candidate");
    }
}

void zero_gradient_policy(){
    Model m;m.pendulum=true;m.alpha=0;const State old{{1,0},{}};
    for(double beta:{0.,.125}){
        const Context c=context(m,old,.25,beta,2);const Iterate start{old.x,0,.5};ss::Options o;o.gate=ss::Gate::None;ss::Diagnostics d;
        const Iterate actual=apply(c,start,o,d);
        near(actual.x,old.x,0,"zero gradient cannot move endpoint");
        if(beta==0){near(actual.eta,start.eta,0,"incompatible hard zero-gradient row skips energy");require(d.core.singular==1&&d.core.incompatible==1&&d.core.scalar_solves==1,"hard zero-gradient policy diagnostics");}
        else{near(actual.eta,1,0,"compliance balances zero-gradient scalar row through accumulated eta");require(d.core.singular==0&&d.accepted==1&&d.core.scalar_solves==2,"compliant zero-gradient row remains solvable");}
    }
}

void small_step_baseline_formula(){
    // With eta correction disabled, one radial material row is the exact BE
    // oscillator update. Test 4/8/16 substeps over 60 frames, not extra sweeps.
    Model m;m.alpha=1./64;m.acceleration={0,-8};constexpr double equilibrium=-1.125,omega=8,frame_h=1./60;
    const State initial{{0,equilibrium+.0625},{0,.25}};
    const double excitation0=.5*initial.v.y*initial.v.y+.5*omega*omega*.0625*.0625;
    for(int substeps:{4,8,16}){
        const double h=frame_h/substeps;State state=initial;ss::Diagnostics d;ss::Options o;o.gate=ss::Gate::None;o.relaxation=0;
        const int total=60*substeps;
        for(int step=0;step<total;++step){
            const Context c=context(m,state,h,0,physical_energy(m,initial));const Iterate next=apply(c,{c.predictor},o,d);
            const double q=state.x.y-equilibrium,denominator=1+h*h*omega*omega;
            const double expected_q=(q+h*state.v.y)/denominator,expected_v=(state.v.y-h*omega*omega*q)/denominator;
            const State updated{next.x,(next.x-state.x)/h};
            near(updated.x.y-equilibrium,expected_q,6e-16,"one material sweep per substep BE position");
            near(updated.v.y,expected_v,4e-13,"one material sweep per substep BE velocity");
            state=updated;
        }
        const std::complex<double> initial_phase(.0625,.25/omega);
        const auto factor=std::complex<double>(1,-h*omega)/(1+h*h*omega*omega);
        const auto expected=initial_phase*std::pow(factor,total);
        near(state.x.y-equilibrium,expected.real(),5e-12,"small-step closed-form oscillator phase");
        near(state.v.y,expected.imag()*omega,5e-11,"small-step closed-form oscillator velocity");
        const double excitation=.5*state.v.y*state.v.y+.5*omega*omega*(state.x.y-equilibrium)*(state.x.y-equilibrium);
        near(excitation/excitation0,std::pow(1+h*h*omega*omega,-total),2e-11,"baseline energy damping per substep");
        require(d.core.evaluations==total&&d.attempts==total,"exactly one EL attempt per substep");
    }
}
} // namespace
int main(int argc,char** argv){
    try{
        if(argc==2&&std::string(argv[1])=="predictor-equivalence"){
            predictor_reset_gate_equivalence();
            std::cout<<"PASS: predictor/reset raw-effective slope identity and bit-exact update equivalence\n";
            return 0;
        }
        raw_legacy_parity();relaxation_and_compliance();cap_is_scaling_not_rejection();effective_gate_and_response();static_equilibrium();zero_gradient_policy();small_step_baseline_formula();
        predictor_reset_gate_equivalence();
        std::cout<<"PASS: legacy parity; relaxation and compliance; cap scaling/order; raw/effective gates; static equilibrium audit; zero-gradient rows; 4/8/16-substep baseline formula; predictor/reset gate equivalence\n";
    }catch(const std::exception& e){std::cerr<<"FAIL: "<<e.what()<<'\n';return 1;}
}
