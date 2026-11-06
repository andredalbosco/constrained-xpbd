#include "correction.hpp"
#include "frozen_schedules.hpp"
#include "refreshed.hpp"
#include <array>
#include <cmath>
#include <complex>
#include <cstdint>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>

using namespace energy_ordering;
namespace cs=correction_schedules;
namespace {
void require(bool ok,const std::string& message){if(!ok)throw std::runtime_error(message);}
void near(double actual,double expected,double tolerance,const std::string& message){
    if(!std::isfinite(actual)||!std::isfinite(expected)||std::abs(actual-expected)>tolerance)
        throw std::runtime_error(message+"; actual="+std::to_string(actual)+", expected="+std::to_string(expected));
}
void near(V2 actual,V2 expected,double tolerance,const std::string& message){near(norm(actual-expected),0,tolerance,message);}
void same(Iterate actual,Iterate expected,double tolerance,const std::string& message){
    near(actual.x,expected.x,tolerance,message+" position");
    near(actual.lambda,expected.lambda,tolerance,message+" lambda");
    near(actual.eta,expected.eta,tolerance,message+" eta");
}
double physical_energy(const Model& m,const State& s){
    const double extension=std::hypot(s.x.x,s.x.y)-m.length;
    return .5*m.mass*dot(s.v,s.v)+(m.pendulum?0:extension*extension/(2*m.alpha))-m.mass*dot(m.acceleration,s.x);
}
double physical_scale(const Context& c){
    V2 acceleration=c.model.acceleration;
    if(!c.model.pendulum){const double radius=norm(c.old.x);acceleration+=c.old.x*(-(radius-c.model.length)/(c.model.mass*c.model.alpha*radius));}
    return c.h*norm(c.old.v)+c.h*c.h*norm(acceleration);
}
struct IndependentRows {V2 n,g,q;double C,H,at,bt,A,B,D;};
IndependentRows independent_rows(const Context& c,const Iterate& z){
    const auto& m=c.model;const double radius=norm(z.x),h2=c.h*c.h,extension=radius-m.length;
    IndependentRows r;r.n=z.x/radius;
    const V2 elastic=m.pendulum?V2{}:r.n*(extension/m.alpha);
    r.q=(z.x-(c.old.x+c.old.v*c.h+m.acceleration*h2))*(m.mass/h2)+elastic;
    r.g=(z.x-c.old.x)*(m.mass/h2)+elastic-m.acceleration*m.mass;
    r.at=m.pendulum?0:m.alpha/h2;r.bt=c.beta/h2;
    r.C=extension+r.at*z.lambda;
    r.H=physical_energy(m,{z.x,(z.x-c.old.x)/c.h})-c.target+r.bt*z.eta;
    r.A=1/m.mass+r.at;r.B=dot(r.n,r.g)/m.mass;r.D=dot(r.g,r.g)/m.mass+r.bt;
    return r;
}
struct Pair {double l=0,e=0;};
Pair baseline(const IndependentRows& r){return {-r.C/r.A,0};}
V2 displacement(const Context& c,const IndependentRows& r,Pair p){return (r.n*p.l+r.g*p.e)/c.model.mass;}
Pair scalar_oracle(const IndependentRows& r,const std::string& visits){
    // Compute TOTAL increments. A repeated row replaces its previous increment.
    long double l=0,e=0;
    for(char row:visits){
        if(row=='L')l=(-static_cast<long double>(r.C)-r.B*e)/r.A;
        else if(row=='E')e=(-static_cast<long double>(r.H)-r.B*l)/r.D;
        else throw std::runtime_error("unknown oracle row");
    }
    return {static_cast<double>(l),static_cast<double>(e)};
}
Pair joint_oracle(const IndependentRows& r){
    // Independent long-double Cramer's rule, not the production Schur solve.
    const long double A=r.A,B=r.B,D=r.D,C=r.C,H=r.H,det=A*D-B*B;
    require(det>0,"joint oracle needs a positive definite block");
    return {static_cast<double>((-C*D+B*H)/det),static_cast<double>((B*C-A*H)/det)};
}
Context coupled_fixture(double H=.5){
    Model m;m.alpha=1./64;const State old{{1.46875,-.046875},{-1.75,.375}};
    Context c=context(m,old,1./8,0,0);
    c.target=physical_energy(m,{c.predictor,(c.predictor-old.x)/c.h})-H;
    return c;
}
void fixture_checks(){
    const Context c=coupled_fixture();const Iterate start{c.predictor};const auto r=independent_rows(c,start);
    near(start.x,{1.25,0},0,"coupled fixture predictor");near(r.g,{2,3},0,"coupled fixture energy gradient");
    near(r.q,{16,0},0,"coupled fixture objective gradient");near(r.A,2,0,"coupled fixture A");
    near(r.B,2,0,"coupled fixture B");near(r.D,13,0,"coupled fixture D");near(r.C,.25,0,"coupled fixture material residual");
    near(r.H,.5,0,"coupled fixture energy residual");
    near(r.H+dot(r.g,displacement(c,r,baseline(r))),.25,0,"coupled fixture baseline energy prediction");
}

Pair spectral_oracle(const IndependentRows& r,double power,double cap){
    const long double sa=std::sqrt(static_cast<long double>(r.A)),sd=std::sqrt(static_cast<long double>(r.D));
    const long double rho=r.B/(sa*sd),u=-r.C/sa,v=-r.H/sd;
    auto mode=[&](long double eigenvalue,long double rhs){
        require(eigenvalue>1e-12L,"full-rank spectral oracle fixture");
        long double gain=std::pow(eigenvalue,-static_cast<long double>(power));
        if(cap>0)gain=std::min(gain,static_cast<long double>(cap));
        return gain*rhs/2;
    };
    const long double plus=mode(1+rho,u+v),minus=mode(1-rho,u-v);
    return {static_cast<double>((plus+minus)/sa),static_cast<double>((plus-minus)/sd)};
}
Pair schedule_oracle(const IndependentRows& r,cs::Schedule schedule,long step=0,double power=.5,double gain_cap=4){
    switch(schedule){
    case cs::Schedule::Baseline:return baseline(r);
    case cs::Schedule::Jacobi:return {-r.C/r.A,-r.H/r.D};
    case cs::Schedule::EL:return scalar_oracle(r,"EL");
    case cs::Schedule::LE:return scalar_oracle(r,"LE");
    case cs::Schedule::LEL:return scalar_oracle(r,"LEL");
    case cs::Schedule::ELE:return scalar_oracle(r,"ELE");
    case cs::Schedule::Alternating:return scalar_oracle(r,step%2==0?"EL":"LE");
    case cs::Schedule::Joint:case cs::Schedule::GuardedJoint:return joint_oracle(r);
    case cs::Schedule::Spectral:return spectral_oracle(r,power,gain_cap);
    }
    throw std::runtime_error("unknown oracle schedule");
}
const std::array<cs::Schedule,10> schedules={cs::Schedule::Baseline,cs::Schedule::Jacobi,cs::Schedule::EL,cs::Schedule::LE,
    cs::Schedule::Joint,cs::Schedule::GuardedJoint,cs::Schedule::Alternating,cs::Schedule::LEL,cs::Schedule::ELE,cs::Schedule::Spectral};
Iterate apply(const Context& c,Iterate z,const cs::Options& options,long step,cs::Diagnostics& d){
    require(cs::update(c,z,options,step,d),"correction update failed");return z;
}
Iterate expected_iterate(const Context& c,const IndependentRows& r,Iterate start,Pair increment){
    return {start.x+displacement(c,r,increment),start.lambda+increment.l,start.eta+increment.e};
}

void schedule_oracles(){
    for(double H:{-.5,.125,.5,2.})for(double beta:{0.,.03125}){
        Context c=coupled_fixture(H);c.beta=beta;
        for(Iterate start:{Iterate{c.predictor},Iterate{c.predictor+V2{.015625,-.03125},-.0078125,.001953125}}){
            const auto r=independent_rows(c,start);
            for(auto schedule:schedules)for(long step:{0,1}){
                cs::Options o;o.schedule=schedule;cs::Diagnostics d;
                const Iterate actual=apply(c,start,o,step,d);
                const Pair expected=schedule_oracle(r,schedule,step,o.power,o.gain_cap);
                same(actual,expected_iterate(c,r,start,expected),2e-14,"independent two-row schedule oracle");
                require(d.core.evaluations==1,"each schedule freezes exactly one row evaluation");
                const V2 dx=actual.x-start.x;
                near(dx*c.model.mass-r.n*(actual.lambda-start.lambda)-r.g*(actual.eta-start.eta),{},5e-14,
                     "position is the sum of final multiplier increments");
                if(schedule==cs::Schedule::EL||schedule==cs::Schedule::LEL||schedule==cs::Schedule::Joint||schedule==cs::Schedule::GuardedJoint)
                    near(r.C+dot(r.n,dx)+r.at*(actual.lambda-start.lambda),0,5e-14,"material-ending schedule closes material row");
                if(schedule==cs::Schedule::LE||schedule==cs::Schedule::ELE||schedule==cs::Schedule::Joint||schedule==cs::Schedule::GuardedJoint)
                    near(r.H+dot(r.g,dx)+r.bt*(actual.eta-start.eta),0,8e-14,"energy-ending schedule closes energy row");
                const long expected_energy=schedule==cs::Schedule::Baseline?0:(schedule==cs::Schedule::ELE?2:1);
                require(d.energy_attempts==expected_energy,"schedule counts actual energy-row visits");
            }
        }
    }
}

void immediate_deferred_equivalence(){
    const Context c=coupled_fixture();const Iterate start{c.predictor};const auto r=independent_rows(c,start);
    for(const std::string visits:{"EL","LE","LEL","ELE"}){
        V2 immediate=start.x;Pair p;V2 incorrectly_accumulated=start.x;
        for(char row:visits){
            if(row=='L'){
                const double next=(-r.C-r.B*p.e)/r.A;
                immediate+=r.n*((next-p.l)/c.model.mass);incorrectly_accumulated+=r.n*(next/c.model.mass);p.l=next;
            }else{
                const double next=(-r.H-r.B*p.l)/r.D;
                immediate+=r.g*((next-p.e)/c.model.mass);incorrectly_accumulated+=r.g*(next/c.model.mass);p.e=next;
            }
        }
        const auto schedule=visits=="EL"?cs::Schedule::EL:(visits=="LE"?cs::Schedule::LE:(visits=="LEL"?cs::Schedule::LEL:cs::Schedule::ELE));
        cs::Options o;o.schedule=schedule;cs::Diagnostics d;const Iterate actual=apply(c,start,o,0,d);
        near(actual.x,immediate,3e-16,"immediate differences and deferred final sum agree on frozen rows");
        near(actual.lambda,p.l,0,"final material increment is replaced, not accumulated");near(actual.eta,p.e,0,"final energy increment is replaced, not accumulated");
        if(visits.size()==3)require(norm(actual.x-incorrectly_accumulated)>.01,"symmetric fixture detects double counting revisited multipliers");
        require(d.core.evaluations==1,"immediate oracle never requires derivative recomputation");
    }
}

void old_el_parity(){
    for(double beta:{0.,.02})for(double H:{-.5,.5}){
        Context c=coupled_fixture(H);c.beta=beta;
        for(Iterate start:{Iterate{c.predictor},Iterate{c.predictor+V2{.01,.02},-.01,.002}}){
            cs::Options o;o.schedule=cs::Schedule::EL;cs::Diagnostics d;const Iterate actual=apply(c,start,o,0,d);
            Iterate expected=start;Counters count;
            require(schedule_update(c,expected,energy_ordering::Schedule::EL,false,count),"legacy frozen EL failed");
            same(actual,expected,0,"raw candidate preserves old EL exactly");
        }
    }
}

Rows matrix_rows(double A,double B,double D,double C=.25,double H=.5){
    Rows r;r.A=A;r.B=B;r.D=D;r.rC=C;r.rE=H;r.C=C;r.H=H;
    r.g={std::sqrt(D),0};r.n={1,0};r.gradient_scale=1;
    r.schur=static_cast<double>(static_cast<long double>(D)-static_cast<long double>(B)*B/A);
    return r;
}
IndependentRows independent_matrix(const Rows& r){return {r.n,r.g,{},r.rC,r.rE,r.at,r.bt,r.A,r.B,r.D};}
void rank_and_guard_handling(){
    const Context c=coupled_fixture();
    for(bool compatible:{false,true}){
        Rows r=matrix_rows(1,2,4,.25,compatible?.5:1.);
        cs::Options o;o.schedule=cs::Schedule::Joint;cs::Diagnostics d;const auto z=cs::propose(c,r,o,0,d);
        near(z.lambda,-.25,0,"rank-deficient joint falls back to baseline");near(z.eta,0,0,"rank fallback leaves energy multiplier alone");
        require(d.rank_fallbacks==1&&d.core.rank_deficient_blocks==1&&d.core.incompatible==(compatible?0:1),"joint rank compatibility diagnostics");
        o.schedule=cs::Schedule::Spectral;o.power=1;o.gain_cap=-1;d={};const auto pseudo=cs::propose(c,r,o,0,d);
        require(d.core.rank_deficient_blocks==1&&d.core.incompatible==(compatible?0:1),"spectral nullspace compatibility diagnostics");
        near(pseudo.lambda,compatible?-.125:-.1875,0,"normalized spectral pseudoinverse material multiplier");
        near(pseudo.eta,compatible?-.0625:-.09375,0,"normalized spectral pseudoinverse energy multiplier");
    }
    for(auto schedule:schedules){
        Rows r=matrix_rows(2,0,0,.25,1);r.g={};r.schur=0;
        cs::Options o;o.schedule=schedule;cs::Diagnostics d;const auto z=cs::propose(c,r,o,0,d);
        near(z.lambda,-.125,0,"zero gradient keeps material baseline");near(z.eta,0,0,"zero gradient skips energy");
        require(d.core.singular==(schedule==cs::Schedule::Baseline?0:1),"zero-gradient rank policy is explicit");
    }
    // A positive but small normalized eigenvalue makes GuardedJoint select EL.
    Rows r=matrix_rows(2,(1-5e-4)*std::sqrt(26.),13);cs::Options o;o.schedule=cs::Schedule::GuardedJoint;cs::Diagnostics d;
    const auto actual=cs::propose(c,r,o,0,d);const auto expected=scalar_oracle(independent_matrix(r),"EL");
    near(actual.lambda,expected.l,2e-16,"guarded joint EL fallback lambda");near(actual.eta,expected.e,2e-16,"guarded joint EL fallback eta");
    require(d.guarded_fallbacks==1&&d.core.block_solves==0&&d.core.scalar_solves==2,"guard fallback costs two scalar rows");
    // Exact guard boundary is strict, while the rank boundary is inclusive.
    r=matrix_rows(1,.5,1);o.guard_tol=.5;d={};cs::propose(c,r,o,0,d);require(d.guarded_fallbacks==0&&d.core.block_solves==1,"guard equality permits joint");
    o.guard_tol=std::nextafter(.5,1.);d={};cs::propose(c,r,o,0,d);require(d.guarded_fallbacks==1,"just above guard boundary falls back");
    o.schedule=cs::Schedule::Joint;o.rank_tol=.5;d={};cs::propose(c,r,o,0,d);require(d.rank_fallbacks==1,"rank equality falls back");
}

void normalized_spectral_limits(){
    const Context c=coupled_fixture();
    for(double A:{.125,2.,16.})for(double D:{.25,13.,64.})for(double rho:{-.9,-.25,0.,.5,.95}){
        Rows r=matrix_rows(A,rho*std::sqrt(A*D),D,.375,-.625);const auto oracle=independent_matrix(r);
        cs::Options o;o.schedule=cs::Schedule::Spectral;o.gain_cap=-1;
        for(double power:{0.,.25,.5,1.}){
            o.power=power;cs::Diagnostics d;const auto actual=cs::propose(c,r,o,0,d);
            const auto expected=spectral_oracle(oracle,power,-1);
            near(actual.lambda,expected.l,4e-12,"normalized spectral power lambda");near(actual.eta,expected.e,4e-12,"normalized spectral power eta");
            if(power==0){near(actual.lambda,-r.rC/r.A,4e-15,"p=0 Jacobi material limit");near(actual.eta,-r.rE/r.D,4e-15,"p=0 Jacobi energy limit");}
            if(power==1){const auto joint=joint_oracle(oracle);near(actual.lambda,joint.l,5e-12,"p=1 exact joint material limit");near(actual.eta,joint.e,5e-12,"p=1 exact joint energy limit");}
            require(d.spectral_solves==1&&d.energy_attempts==1,"one normalized block operation per spectral proposal");
        }
        // Diagonal scaling defines the norm in which the eigenmode gain is bounded.
        for(double cap:{.5,1.,4.}){
            o.power=1;o.gain_cap=cap;cs::Diagnostics d;const auto actual=cs::propose(c,r,o,0,d);
            const auto expected=spectral_oracle(oracle,1,cap);
            near(actual.lambda,expected.l,4e-13,"bounded spectral gain material multiplier");near(actual.eta,expected.e,4e-13,"bounded spectral gain energy multiplier");
            const double input=std::hypot(r.rC/std::sqrt(A),r.rE/std::sqrt(D));
            const double output=std::hypot(actual.lambda*std::sqrt(A),actual.eta*std::sqrt(D));
            require(output<=cap*input*(1+1e-14),"spectral cap bounds normalized multiplier norm amplification");
        }
    }
}

void filter_properties(){
    const std::array<cs::Policy,6> policies={cs::Policy::Raw,cs::Policy::Direction,cs::Policy::Cap,cs::Policy::DirectionCap,cs::Policy::Repair,cs::Policy::RepairCap};
    for(double H:{-.5,.125,.25,.5,2.})for(double beta:{0.,.03125}){
        Context c=coupled_fixture(H);c.beta=beta;
        const Iterate start{c.predictor,0,.001953125};const auto r=independent_rows(c,start);const Pair base=baseline(r);
        for(auto schedule:schedules)for(auto policy:policies){
            cs::Options o;o.schedule=schedule;o.policy=policy;o.fraction=.25;o.cap=.04;
            const Pair raw=schedule_oracle(r,schedule,0,o.power,o.gain_cap);
            const V2 dx0=displacement(c,r,base),extra=displacement(c,r,{raw.l-base.l,raw.e});
            const double H0=r.H+dot(r.g,dx0),change=dot(r.g,extra)+r.bt*raw.e,scale=physical_scale(c);
            double blend=1;const bool direction=policy==cs::Policy::Direction||policy==cs::Policy::DirectionCap;
            const bool repair=policy==cs::Policy::Repair||policy==cs::Policy::RepairCap;
            const bool cap=policy==cs::Policy::Cap||policy==cs::Policy::DirectionCap||policy==cs::Policy::RepairCap;
            if(schedule!=cs::Schedule::Baseline){
                if(direction&&dot(r.q,extra)>0)blend=0;
                if(repair){
                    const double tolerance=128*std::numeric_limits<double>::epsilon()*(1+std::abs(c.target)+std::abs(r.H-r.bt*start.eta)+std::abs(r.bt));
                    if(std::abs(H0)<=tolerance||H0*change>=0)blend=0;
                    else blend=std::min(blend,o.fraction*std::abs(H0)/std::abs(change));
                }
                if(cap&&norm(extra)>0)blend=std::min(blend,o.cap*scale/norm(extra));
            }
            const Pair accepted{base.l+blend*(raw.l-base.l),blend*raw.e};
            cs::Diagnostics d;const Iterate actual=apply(c,start,o,0,d);
            same(actual,expected_iterate(c,r,start,accepted),4e-14,"filter blends BOTH multipliers toward material-only baseline");
            near(d.last_blend,blend,5e-14,"independent trust-filter blend");
            require(d.last_blend>=0&&d.last_blend<=1,"filter blend is convex");
            const V2 actual_extra=actual.x-start.x-dx0;
            near(actual_extra,extra*blend,5e-14,"position extra matches blended paired multiplier response");
            near(d.last_H0,H0,2e-14,"energy prediction starts after material baseline");
            near(d.last_candidate_change,change,4e-14,"candidate energy change includes compliance eta term");
            near(d.last_predicted_H,H0+blend*change,4e-14,"accepted linear energy prediction");
            if(schedule!=cs::Schedule::Baseline){
                if(direction)require(dot(r.q,actual_extra)<=1e-12,"direction filter rejects positive objective slope of extra");
                if(cap)require(norm(actual_extra)<=o.cap*scale+1e-14,"physical cap bounds actual EXTRA displacement");
                if(repair){
                    require(std::abs(blend*change)<=o.fraction*std::abs(H0)+1e-14,"repair corrects at most requested fraction");
                    require(std::abs(H0+blend*change)<=std::abs(H0)+1e-14,"repair never worsens predicted residual");
                    require(H0*(H0+blend*change)>=-1e-14,"fraction under one cannot overshoot zero");
                }
            }
            const double candidate_material=r.C+r.A*raw.l+r.B*raw.e;
            near(r.C+dot(r.n,actual.x-start.x)+r.at*(actual.lambda-start.lambda),blend*candidate_material,7e-14,
                 "material residual interpolates from the closed baseline row");
            require(d.core.evaluations==1,"filters preserve one frozen linearization");
        }
    }
}

void targeted_filter_boundaries(){
    // Repair must not use the predictor's H in place of the baseline's H0.
    for(double H:{.125,.25}){
        const Context c=coupled_fixture(H);const Iterate start{c.predictor};const auto r=independent_rows(c,start);
        cs::Options o;o.schedule=cs::Schedule::EL;o.policy=cs::Policy::Repair;cs::Diagnostics d;
        const Iterate actual=apply(c,start,o,0,d);
        same(actual,expected_iterate(c,r,start,baseline(r)),0,"nonimproving or zero-baseline-residual correction skips to baseline");
        require(d.repair_skips==1&&d.last_blend==0,"repair rejection explicitly recorded");
    }
    // Capping the combined correction must scale material's compensating lambda.
    const Context c=coupled_fixture();const Iterate start{c.predictor};const auto r=independent_rows(c,start);const Pair base=baseline(r),raw=joint_oracle(r);
    const V2 extra=displacement(c,r,{raw.l-base.l,raw.e});
    cs::Options o;o.schedule=cs::Schedule::Joint;o.policy=cs::Policy::Cap;o.cap=.25*norm(extra)/physical_scale(c);cs::Diagnostics d;
    const Iterate actual=apply(c,start,o,0,d);
    near(actual.eta,.25*raw.e,2e-17,"cap scales eta");near(actual.lambda,base.l+.25*(raw.l-base.l),3e-17,"cap blends lambda toward baseline");
    near(r.C+dot(r.n,actual.x-start.x)+r.at*actual.lambda,0,3e-16,"joint cap preserves closed material row");
    require(std::abs(r.C+r.A*raw.l+r.B*.25*raw.e)>.01,"fixture detects scaling eta while leaving candidate lambda unchanged");
    require(d.capped==1&&d.accepted==1,"partial cap counts accepted scaling");
    // Exact neutral objective gradient is permitted, without amplitude filtering.
    Model m;m.alpha=1./64;const Context neutral=context(m,{{1.125,0},{-1,0}},1./8,0,2);
    o={};o.schedule=cs::Schedule::EL;o.policy=cs::Policy::Direction;d={};
    const Iterate unbounded=apply(neutral,{neutral.predictor},o,0,d);
    same(unbounded,{{.90625,0},.09375,.0234375},0,"neutral objective gradient permits full correction");
    require(d.direction_skips==0&&d.last_blend==1,"zero objective slope accepted");
}

void physical_unit_scaling(){
    // Rescale mass, length, time units while retaining the same dimensionless
    // problem. Position scales by L, lambda by M*L, and eta by T^2.
    const Context original=coupled_fixture();const Iterate start{original.predictor};
    for(auto schedule:schedules)for(auto policy:{cs::Policy::Raw,cs::Policy::DirectionCap,cs::Policy::RepairCap}){
        cs::Options o;o.schedule=schedule;o.policy=policy;o.cap=.04;o.fraction=.25;
        cs::Diagnostics original_d;const Iterate original_result=apply(original,start,o,0,original_d);
        for(double mass_scale:{.125,1.,8.})for(double length_scale:{.125,1.,8.})for(double time_scale:{.125,1.,8.}){
            Model model=original.model;
            model.mass*=mass_scale;model.length*=length_scale;
            model.alpha*=time_scale*time_scale/mass_scale;
            model.acceleration=model.acceleration*(length_scale/(time_scale*time_scale));
            const State old{original.old.x*length_scale,original.old.v*(length_scale/time_scale)};
            const double energy_scale=mass_scale*length_scale*length_scale/(time_scale*time_scale);
            const Context c=context(model,old,original.h*time_scale,original.beta*energy_scale,original.target*energy_scale);
            cs::Diagnostics d;const Iterate actual=apply(c,{c.predictor},o,0,d);
            near(actual.x/length_scale,original_result.x,3e-14,"unit scaling preserves normalized position");
            near(actual.lambda/(mass_scale*length_scale),original_result.lambda,3e-14,"unit scaling preserves normalized material multiplier");
            near(actual.eta/(time_scale*time_scale),original_result.eta,3e-14,"unit scaling preserves normalized energy multiplier");
            near(d.last_blend,original_d.last_blend,3e-14,"dimensionless trust decisions unchanged by units");
            require(d.direction_skips==original_d.direction_skips&&d.repair_skips==original_d.repair_skips&&d.capped==original_d.capped,
                    "unit scaling preserves acceptance decisions away from numerical thresholds");
        }
    }
}

void degeneracy_audit(){
    double largest_joint_position_error=0,largest_joint_energy_row_error=0;long joint_rank_skips=0,guard_fallbacks=0;
    for(int exponent:{0,8,16,24,32,38,40,48})for(double transverse:{0.,std::ldexp(1.,-12),std::ldexp(1.,-24)}){
        Model m;const double h=1./8,at=std::ldexp(1.,-exponent);m.alpha=at*h*h;
        const V2 velocity{.25,h*transverse};const State old{V2{1,0}-velocity*h,velocity};
        Context c=context(m,old,h,0,0);c.target=physical_energy(m,{c.predictor,(c.predictor-old.x)/h})-.3;
        const Iterate start{c.predictor};const auto r=independent_rows(c,start);
        near(start.x,{1,0},0,"degeneracy fixture predictor");near(r.g,{2,transverse},0,"degeneracy fixture energy gradient");
        cs::Options o;o.schedule=cs::Schedule::Joint;cs::Diagnostics d;const Iterate actual=apply(c,start,o,0,d);
        if(d.rank_fallbacks){++joint_rank_skips;same(actual,start,0,"near-dependent joint explicitly keeps zero material baseline");}
        else{
            // Stable independent exact position oracle, avoiding subtraction of
            // lambda and g*eta when those are both huge and almost parallel.
            const long double a=at,A=1+a,gx=2,gy=transverse;
            const long double schur=gy*gy+a/A*gx*gx;
            const long double eta=-static_cast<long double>(r.H)/schur;
            const V2 expected_dx{static_cast<double>((a/A)*gx*eta),static_cast<double>(gy*eta)};
            largest_joint_position_error=std::max(largest_joint_position_error,norm(actual.x-start.x-expected_dx));
            largest_joint_energy_row_error=std::max(largest_joint_energy_row_error,std::abs(r.H+dot(r.g,actual.x-start.x)));
            require(finite(actual.x)&&std::isfinite(actual.lambda)&&std::isfinite(actual.eta),"full-rank joint stress remains finite");
        }
        o.schedule=cs::Schedule::GuardedJoint;d={};const Iterate guarded=apply(c,start,o,0,d);guard_fallbacks+=d.guarded_fallbacks;
        require(finite(guarded.x)&&std::abs(guarded.lambda)<1e5&&std::abs(guarded.eta)<1e5,"guarded joint prevents extreme inverse amplification in stress fixtures");
        o.schedule=cs::Schedule::Spectral;o.power=1;o.gain_cap=4;d={};const Iterate spectral=apply(c,start,o,0,d);
        const double output=std::hypot(spectral.lambda*std::sqrt(r.A),spectral.eta*std::sqrt(r.D));
        const double rhs=std::hypot(r.C/std::sqrt(r.A),r.H/std::sqrt(r.D));
        require(output<=4*rhs*(1+1e-13),"spectral gain remains bounded as rows approach dependence");
    }
    require(joint_rank_skips>0&&guard_fallbacks>0,"degeneracy stress exercises explicit fallback paths");
    std::cout<<"Degeneracy audit: max full-rank joint position error="<<largest_joint_position_error
             <<", max linear energy-row error="<<largest_joint_energy_row_error
             <<", joint rank fallbacks="<<joint_rank_skips<<", guarded fallbacks="<<guard_fallbacks<<'\n';
}

void hard_rod_predictor_direction_gate(){
    for(double mass:{.5,1.,5.})for(double h:{1./240,1./960})
    for(double theta:{.1,.65,1.5})for(double speed:{0.,.8,2.}){
        Model m;m.pendulum=true;m.alpha=0;m.mass=mass;m.acceleration={0,-9.81};
        const V2 x{std::sin(theta),-std::cos(theta)};const State old{x,tangent(x)*speed};
        const Context c=context(m,old,h,0,physical_energy(m,old));const Iterate start{c.predictor};
        const auto frozen=evaluate(c,start);
        near(gradient_policies::objective_gradient(c,start,frozen),{},0,"hard-rod objective gradient is exactly zero at reset predictor");
        if(speed==0){
            near(frozen.g,{},0,"resting rod has zero endpoint energy gradient at ballistic predictor");
            near(frozen.H,-.5*mass*h*h*dot(m.acceleration,m.acceleration),1e-13,"resting rod energy deficit at predictor is quadratic in time step");
        }
        for(auto schedule:schedules){
            cs::Options raw;raw.schedule=schedule;cs::Options gated=raw;gated.policy=cs::Policy::Direction;
            cs::Diagnostics dr,dg;const Iterate a=apply(c,start,raw,0,dr),b=apply(c,start,gated,0,dg);
            same(a,b,0,"direction gate is inactive for hard-rod one-predictor linearization");
            require(dg.direction_skips==0,"hard-rod predictor cannot trigger a strictly positive slope rejection");
            if(speed==0&&schedule!=cs::Schedule::Baseline)
                require(dr.core.singular==1&&dr.core.incompatible==1,"resting hard-rod energy row is explicitly incompatible at first linear order");
        }
    }
}

void deterministic_physical_stress(){
    std::uint64_t random_state=0x88bb4302c616d367ULL;
    auto uniform=[&](){random_state=random_state*6364136223846793005ULL+1442695040888963407ULL;return (random_state>>11)*0x1.0p-53;};
    double max_cap_excess=0,max_pair_roundoff=0;long updates=0,rank_cases=0;
    for(int sample=0;sample<512;++sample){
        Model m;m.mass=std::pow(10.,-3+6*uniform());m.length=std::pow(10.,-1+2*uniform());
        m.pendulum=sample%5==0;m.alpha=m.pendulum?0:std::pow(10.,-8+10*uniform());
        m.acceleration={m.length*(-10+20*uniform()),m.length*(-10+20*uniform())};
        const double theta=6.283185307179586*uniform(),radius=m.length*(.1+2*uniform()),h=std::pow(10.,-4+2*uniform());
        const State old{{radius*std::cos(theta),radius*std::sin(theta)},{m.length*(-2+4*uniform()),m.length*(-2+4*uniform())}};
        const Context c=context(m,old,h,0,physical_energy(m,old));const Iterate start{c.predictor};const auto r=independent_rows(c,start);
        const Rows numerical_rows=evaluate(c,start);
        const Pair base=baseline(r);const V2 dx0=displacement(c,r,base);const double scale=physical_scale(c);
        for(auto schedule:schedules)for(auto policy:{cs::Policy::Raw,cs::Policy::DirectionCap,cs::Policy::RepairCap}){
            cs::Options o;o.schedule=schedule;o.policy=policy;o.cap=.1;o.fraction=.5;cs::Diagnostics d;
            const Iterate actual=apply(c,start,o,0,d);++updates;rank_cases+=d.rank_fallbacks+d.core.singular;
            const V2 dx=actual.x-start.x,actual_extra=dx-dx0;
            // Use the same stored floating-point coefficients for the pure
            // assembly identity. Algebraically equivalent gravity centering
            // can round differently; physical row values have separate tests.
            const V2 paired=(numerical_rows.n*actual.lambda+numerical_rows.g*actual.eta)/m.mass;
            max_pair_roundoff=std::max(max_pair_roundoff,norm(dx-paired)/(1+norm(start.x)+norm(dx)));
            require(norm(dx-paired)<=1024*std::numeric_limits<double>::epsilon()*(1+norm(start.x)+norm(dx)),"physical stress position agrees with accepted multiplier pair");
            require(d.last_blend>=0&&d.last_blend<=1&&d.core.evaluations==1,"physical stress uses convex blending and one frozen row evaluation");
            if(schedule!=cs::Schedule::Baseline&&policy!=cs::Policy::Raw){
                const double excess=norm(actual_extra)-o.cap*scale;
                max_cap_excess=std::max(max_cap_excess,excess/(1+norm(start.x)+norm(dx0)));
                require(excess<=2e-10*(1+norm(start.x)+norm(dx0)),"physical stress cap bounds actual extra within floating-point assembly error");
                if(policy==cs::Policy::DirectionCap){
                    const double work=dot(r.q,actual_extra);
                    require(work<=2e-10*(1+norm(r.q)*(norm(start.x)+norm(dx0))),"physical stress accepted extra has nonpositive objective slope");
                }else{
                    const double correction=d.last_blend*d.last_candidate_change;
                    require(std::abs(correction)<=o.fraction*std::abs(d.last_H0)+1e-12*(1+std::abs(d.last_H0)),"physical stress repair fraction is respected");
                    require(std::abs(d.last_H0+correction)<=std::abs(d.last_H0)+1e-12*(1+std::abs(d.last_H0)),"physical stress repair improves predicted baseline residual");
                }
            }
        }
    }
    std::cout<<"Physical stress: updates="<<updates<<", rank/zero-gradient cases="<<rank_cases
             <<", max scaled actual-cap excess="<<max_cap_excess<<", max scaled pair roundoff="<<max_pair_roundoff<<'\n';
}

void refreshed_regressions(){
    namespace refreshed=refreshed_corrections;
    // Away from unresolved gradients, raw refreshed visits reproduce the old
    // per-row XPBD implementation with exactly two or three evaluations.
    const Context c=coupled_fixture();const Iterate start{c.predictor};
    for(auto schedule:{cs::Schedule::EL,cs::Schedule::LE,cs::Schedule::LEL,cs::Schedule::ELE}){
        Iterate expected=start;Counters old_count;
        if(schedule==cs::Schedule::EL){energy_update(c,expected,old_count);material_update(c,expected,old_count);}
        if(schedule==cs::Schedule::LE){material_update(c,expected,old_count);energy_update(c,expected,old_count);}
        if(schedule==cs::Schedule::LEL){material_update(c,expected,old_count);energy_update(c,expected,old_count);material_update(c,expected,old_count);}
        if(schedule==cs::Schedule::ELE){energy_update(c,expected,old_count);material_update(c,expected,old_count);energy_update(c,expected,old_count);}
        for(auto acceptance:{refreshed::Acceptance::WholeProposal,refreshed::Acceptance::AtEnergyVisit}){
            cs::Options o;o.schedule=schedule;cs::Diagnostics d;Iterate actual=start;
            require(refreshed::update(c,actual,o,0,d,acceptance),"refreshed regular update failed");
            same(actual,expected,5e-16,"raw refreshed visits preserve legacy row-by-row solver");
            require(d.core.evaluations==old_count.evaluations&&d.core.evaluations==(schedule==cs::Schedule::LEL||schedule==cs::Schedule::ELE?3:2),"refreshed schedule evaluates at every row visit");
        }
    }
    long formerly_unresolved=0;
    for(double h:{1./240,1./480,1./960}){
        Model m;const State old{{1.25,0},{}};Context rest=context(m,old,h,0,energy(m,old));
        Iterate material_only{rest.predictor};Counters count;require(material_update(rest,material_only,count),"resting material update failed");
        const Rows r=evaluate(rest,material_only);
        require(norm(r.g)<=refreshed::gradient_resolution_tolerance(rest,material_only,r),"rest material stationary gradient is below position resolution");
        formerly_unresolved+=!energy_row_zero(r);
        for(auto acceptance:{refreshed::Acceptance::WholeProposal,refreshed::Acceptance::AtEnergyVisit}){
            cs::Options o;o.schedule=cs::Schedule::LE;cs::Diagnostics d;Iterate actual{rest.predictor};
            require(refreshed::update(rest,actual,o,0,d,acceptance),"guarded refreshed rest update failed");
            same(actual,material_only,0,"unresolved energy gradient retains material result");
            require(d.core.singular==1&&d.core.incompatible==1&&d.core.scalar_solves==1&&d.core.evaluations==2,"representability guard skips division and records incompatible hard row");
        }
        // Finite energy compliance can solve through eta even at g=0.
        rest.beta=.001;cs::Options o;o.schedule=cs::Schedule::LE;cs::Diagnostics d;Iterate soft{rest.predictor};
        require(refreshed::update(rest,soft,o,0,d,refreshed::Acceptance::AtEnergyVisit),"compliant refreshed rest update failed");
        require(d.core.singular==0&&d.core.scalar_solves==2,"representability guard does not discard solvable compliant row");
        // A clearly resolved offset must not be classified as zero.
        Iterate resolved=material_only;resolved.x.x+=1e-6;
        const Rows rr=evaluate(rest,resolved);
        require(norm(rr.g)>100*refreshed::gradient_resolution_tolerance(rest,resolved,rr),"resolved gradient remains available for correction");
    }
    require(formerly_unresolved>0,"fixture reproduces old arithmetic-only tolerance failure");
    // Whole-proposal cap interpolates endpoint AND both accumulated multipliers.
    cs::Options o;o.schedule=cs::Schedule::EL;cs::Diagnostics raw_d;Iterate raw=start;
    require(refreshed::update(c,raw,o,0,raw_d,refreshed::Acceptance::WholeProposal),"raw proposal for blend test failed");
    Iterate base=start;Counters count;material_update(c,base,count);
    const V2 extra=raw.x-base.x;o.policy=cs::Policy::Cap;o.cap=.25*norm(extra)/physical_scale(c);
    cs::Diagnostics d;Iterate filtered=start;require(refreshed::update(c,filtered,o,0,d,refreshed::Acceptance::WholeProposal),"filtered whole proposal failed");
    same(filtered,{base.x+extra*.25,base.lambda+.25*(raw.lambda-base.lambda),base.eta+.25*(raw.eta-base.eta)},5e-16,"whole refreshed filter scales position and both multipliers consistently");
}
} // namespace

int main(int argc,char** argv){
    try{
        fixture_checks();
        const std::string suite=argc==2?argv[1]:"all";
        if(suite=="all"||suite=="schedules"){schedule_oracles();immediate_deferred_equivalence();old_el_parity();rank_and_guard_handling();normalized_spectral_limits();}
        if(suite=="all"||suite=="filters"){filter_properties();targeted_filter_boundaries();}
        if(suite=="all"||suite=="stress"){physical_unit_scaling();degeneracy_audit();hard_rod_predictor_direction_gate();deterministic_physical_stress();}
        if(suite=="all"||suite=="refreshed")refreshed_regressions();
        std::cout<<"PASS: correction schedule "<<suite<<" independent oracles and invariants\n";
    }catch(const std::exception& e){std::cerr<<"FAIL: "<<e.what()<<'\n';return 1;}
}
