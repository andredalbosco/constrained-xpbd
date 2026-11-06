#pragma once
#include "gradient_policies.hpp"

// Original frozen XPBD/energy equations, with alternative approximate solvers.
// All acceptance controls blend a candidate pair toward the material-only pair.
namespace correction_schedules {
using namespace energy_ordering;
enum class Schedule { Baseline, Jacobi, EL, LE, Joint, GuardedJoint, Alternating, LEL, ELE, Spectral };
enum class Policy { Raw, Direction, Cap, DirectionCap, Repair, RepairCap };
inline const char* name(Schedule s) {
    switch(s){case Schedule::Baseline:return "baseline";case Schedule::Jacobi:return "jacobi";
    case Schedule::EL:return "EL";case Schedule::LE:return "LE";case Schedule::Joint:return "joint";
    case Schedule::GuardedJoint:return "guarded_joint";case Schedule::Alternating:return "alternating";
    case Schedule::LEL:return "LEL";case Schedule::ELE:return "ELE";case Schedule::Spectral:return "spectral";}
    return "unknown";
}
inline const char* name(Policy p) {
    switch(p){case Policy::Raw:return "raw";case Policy::Direction:return "direction";case Policy::Cap:return "cap";
    case Policy::DirectionCap:return "direction_cap";case Policy::Repair:return "repair";case Policy::RepairCap:return "repair_cap";}
    return "unknown";
}
struct Options {
    Schedule schedule=Schedule::EL;
    Policy policy=Policy::Raw;
    double fraction=.5,cap=.1,power=.5,gain_cap=4;
    double rank_tol=1e-12,guard_tol=1e-3;
};
struct Diagnostics {
    Counters core;
    long spectral_solves=0,energy_attempts=0,direction_skips=0,repair_skips=0,repair_scaled=0,capped=0;
    long guarded_fallbacks=0,rank_fallbacks=0,accepted=0;
    double max_raw_extra=0,max_accepted_ratio=0,max_coupling=0,min_eigenvalue=1;
    double max_linear_energy_error=0,max_linear_material_error=0,max_prediction_error=0;
    // Last-step diagnostics for failure snapshots; values refer to frozen rows.
    double last_H0=0,last_candidate_change=0,last_blend=1,last_extra=0,last_predicted_H=0;
    double last_raw_energy_displacement=0,last_material_response_displacement=0,last_candidate_extra=0,last_actual_H0=0;
    double max_cancellation_ratio=1;
};
struct Proposal {double lambda=0,eta=0;};
inline Proposal propose(const Context& c,const Rows& r,const Options& o,long substep,Diagnostics& d) {
    const Proposal baseline{-r.rC/r.A,0};
    if(o.schedule==Schedule::Baseline){++d.core.scalar_solves;return baseline;}
    if(energy_row_zero(r)) {
        ++d.core.singular;++d.energy_attempts;++d.core.scalar_solves;
        if(std::abs(r.rE+r.B*baseline.lambda)>energy_tolerance(c,r))++d.core.incompatible;
        return baseline;
    }
    const double sqrtA=std::sqrt(r.A),sqrtD=std::sqrt(r.D);
    const double rho=std::clamp(r.B/(sqrtA*sqrtD),-1.,1.);
    const double smallest=std::max(0.,1-std::abs(rho));
    d.max_coupling=std::max(d.max_coupling,std::abs(rho));d.min_eigenvalue=std::min(d.min_eigenvalue,smallest);
    auto s=o.schedule;
    if(s==Schedule::Alternating)s=substep%2==0?Schedule::EL:Schedule::LE;
    if(s==Schedule::GuardedJoint&&smallest<o.guard_tol){s=Schedule::EL;++d.guarded_fallbacks;}
    if(s==Schedule::Joint||s==Schedule::GuardedJoint) {
        ++d.core.block_solves;++d.energy_attempts;
        if(smallest<=o.rank_tol) {
            ++d.rank_fallbacks;++d.core.rank_deficient_blocks;
            if(std::abs(r.rE+r.B*baseline.lambda)>energy_tolerance(c,r))++d.core.incompatible;
            return baseline;
        }
        // Nonnegative Schur expression avoids subtracting nearly equal D and B²/A.
        const double de=(-r.rE+r.B*r.rC/r.A)/r.schur;
        return {(-r.rC-r.B*de)/r.A,de};
    }
    if(s==Schedule::Spectral) {
        ++d.spectral_solves;++d.energy_attempts;
        const double u=r.rC/sqrtA,v=r.rE/sqrtD;
        auto filtered=[&](double eigenvalue,double rhs) {
            if(eigenvalue<=o.rank_tol) {
                ++d.core.rank_deficient_blocks;
                if(std::abs(rhs)>128*std::numeric_limits<double>::epsilon()*(1+std::abs(u)+std::abs(v)))++d.core.incompatible;
                return 0.;
            }
            double gain=std::pow(eigenvalue,-o.power);
            if(o.gain_cap>0)gain=std::min(gain,o.gain_cap);
            return gain*rhs;
        };
        const double plus=filtered(1+rho,u+v),minus=filtered(1-rho,u-v);
        return {-(plus+minus)/(2*sqrtA),-(plus-minus)/(2*sqrtD)};
    }
    Proposal z;
    auto L=[&](){z.lambda=(-r.rC-r.B*z.eta)/r.A;++d.core.scalar_solves;};
    auto E=[&](){z.eta=(-r.rE-r.B*z.lambda)/r.D;++d.core.scalar_solves;++d.energy_attempts;};
    if(s==Schedule::Jacobi){L();z.eta=-r.rE/r.D;++d.core.scalar_solves;++d.energy_attempts;}
    else if(s==Schedule::LE||s==Schedule::LEL){L();E();if(s==Schedule::LEL)L();}
    else {E();L();if(s==Schedule::ELE)E();}
    return z;
}
inline bool update(const Context& c,Iterate& z,const Options& o,long substep,Diagnostics& d) {
    const Rows r=evaluate(c,z);++d.core.evaluations;if(!r.valid)return false;
    const double mass=c.model.mass;
    const Proposal base{-r.rC/r.A,0},candidate=propose(c,r,o,substep,d);
    if(!std::isfinite(candidate.lambda)||!std::isfinite(candidate.eta))return false;
    const V2 dx0=r.n*(base.lambda/mass);
    const V2 extra=(r.n*(candidate.lambda-base.lambda)+r.g*candidate.eta)/mass;
    const double H0=r.rE+dot(r.g,dx0);
    const double change=dot(r.g,extra)+r.bt*candidate.eta;
    d.last_raw_energy_displacement=norm(r.g*(candidate.eta/mass));
    d.last_material_response_displacement=norm(r.n*((candidate.lambda-base.lambda)/mass));
    d.last_candidate_extra=norm(extra);
    const double separate=d.last_raw_energy_displacement+d.last_material_response_displacement;
    const double cancellation=separate==0?1:(norm(extra)>0?std::min(separate/norm(extra),1e300):1e300);
    d.max_cancellation_ratio=std::max(d.max_cancellation_ratio,cancellation);
    const V2 baseline_x=z.x+dx0;
    d.last_actual_H0=energy(c.model,{baseline_x,(baseline_x-c.old.x)/c.h})-c.target+r.bt*z.eta;
    double blend=1;
    const bool direction=o.policy==Policy::Direction||o.policy==Policy::DirectionCap;
    const bool repair=o.policy==Policy::Repair||o.policy==Policy::RepairCap;
    const bool cap=o.policy==Policy::Cap||o.policy==Policy::DirectionCap||o.policy==Policy::RepairCap;
    if(o.schedule!=Schedule::Baseline) {
        const V2 q=gradient_policies::objective_gradient(c,z,r);
        if(direction&&dot(q,extra)>0){blend=0;++d.direction_skips;}
        if(repair) {
            if(std::abs(H0)<=energy_tolerance(c,r)||H0*change>=0||change==0){blend=0;++d.repair_skips;}
            else if(std::abs(change)>o.fraction*std::abs(H0)) {blend=std::min(blend,o.fraction*std::abs(H0)/std::abs(change));++d.repair_scaled;}
        }
        const double physical_scale=gradient_policies::motion_scale(c),extra_size=norm(extra);
        d.max_raw_extra=std::max(d.max_raw_extra,extra_size);
        if(cap&&blend*extra_size>o.cap*physical_scale){blend=o.cap*physical_scale/extra_size;++d.capped;}
        const double accepted_size=blend*extra_size;
        const double ratio=accepted_size==0?0:(physical_scale>0?std::min(accepted_size/physical_scale,1e300):1e300);
        d.max_accepted_ratio=std::max(d.max_accepted_ratio,ratio);
        if(blend>0)++d.accepted;
    }
    const double dl=blend==1?candidate.lambda:(base.lambda+blend*(candidate.lambda-base.lambda));
    const double de=blend*candidate.eta;
    const V2 dx=(r.n*dl+r.g*de)/mass;
    d.last_H0=H0;d.last_candidate_change=change;d.last_blend=blend;d.last_extra=norm(extra)*blend;
    d.last_predicted_H=r.rE+dot(r.g,dx)+r.bt*de;
    d.max_linear_energy_error=std::max(d.max_linear_energy_error,std::abs(d.last_predicted_H));
    d.max_linear_material_error=std::max(d.max_linear_material_error,std::abs(r.rC+dot(r.n,dx)+r.at*dl));
    z.x=z.x+V2{}+dx;z.lambda+=dl;z.eta+=de;
    if(!finite(z.x)||!std::isfinite(z.lambda)||!std::isfinite(z.eta))return false;
    const double actual=energy(c.model,{z.x,(z.x-c.old.x)/c.h})-c.target+r.bt*z.eta;
    d.max_prediction_error=std::max(d.max_prediction_error,std::abs(actual-d.last_predicted_H));
    return std::isfinite(actual);
}
} // namespace correction_schedules
