#pragma once
#include "correction.hpp"

// Each visit still solves one curvature-dropped linear equation. Only the
// evaluation point is refreshed after a preceding position correction.
namespace refreshed_corrections {
using namespace energy_ordering;
namespace se=correction_schedules;
enum class Acceptance { WholeProposal, AtEnergyVisit };

// Position rounding is magnified by m/h^2 when evaluating the endpoint-energy
// gradient. After a material solve reaches an objective stationary point, this
// representability error can exceed an arithmetic-only gradient tolerance.
// Bound that uncertainty in gradient units; never perturb the denominator.
inline double gradient_resolution_tolerance(const Context& c,const Iterate& z,const Rows& r) {
    const V2 center=c.old.x+c.model.acceleration*(c.h*c.h);
    const double position_roundoff=4*std::numeric_limits<double>::epsilon()
        *(norm(z.x)+norm(center));
    return gradient_tolerance(r)+(c.model.mass/(c.h*c.h))*position_roundoff;
}

inline double filter(const Context& c,const Iterate& at,const Rows& r,
                     const se::Options& o,V2 extra,double residual,double change,
                     se::Diagnostics& d) {
    double blend=1;
    const bool direction=o.policy==se::Policy::Direction||o.policy==se::Policy::DirectionCap;
    const bool repair=o.policy==se::Policy::Repair||o.policy==se::Policy::RepairCap;
    const bool cap=o.policy==se::Policy::Cap||o.policy==se::Policy::DirectionCap||o.policy==se::Policy::RepairCap;
    if(direction&&dot(gradient_policies::objective_gradient(c,at,r),extra)>0) {
        blend=0;++d.direction_skips;
    }
    if(repair) {
        if(std::abs(residual)<=energy_tolerance(c,r)||residual*change>=0||change==0) {
            blend=0;++d.repair_skips;
        } else if(std::abs(change)>o.fraction*std::abs(residual)) {
            blend=std::min(blend,o.fraction*std::abs(residual)/std::abs(change));++d.repair_scaled;
        }
    }
    const double size=norm(extra),motion=gradient_policies::motion_scale(c);
    d.max_raw_extra=std::max(d.max_raw_extra,size);
    if(cap&&blend*size>o.cap*motion) {blend=o.cap*motion/size;++d.capped;}
    const double accepted=blend*size;
    d.max_accepted_ratio=std::max(d.max_accepted_ratio,accepted==0?0:(motion>0?accepted/motion:1e300));
    if(blend>0)++d.accepted;
    d.last_H0=residual;d.last_candidate_change=change;d.last_blend=blend;d.last_extra=accepted;
    d.last_predicted_H=residual+blend*change;
    return blend;
}

inline bool update(const Context& c,Iterate& z,const se::Options& o,long substep,
                   se::Diagnostics& d,Acceptance acceptance) {
    const Iterate start=z;
    Rows r=evaluate(c,z);++d.core.evaluations;if(!r.valid)return false;
    const Rows initial=r;
    const double mass=c.model.mass,base_dl=-r.rC/r.A;
    const Iterate baseline{z.x+r.n*(base_dl/mass),z.lambda+base_dl,z.eta};
    auto schedule=o.schedule;
    if(schedule==se::Schedule::Alternating)schedule=substep%2==0?se::Schedule::EL:se::Schedule::LE;
    const char* visits=schedule==se::Schedule::LE?"LE":schedule==se::Schedule::EL?"EL":
        schedule==se::Schedule::LEL?"LEL":schedule==se::Schedule::ELE?"ELE":nullptr;
    if(!visits)return false;
    int visit=0;
    for(const char* p=visits;*p;++p,++visit) {
        if(visit){r=evaluate(c,z);++d.core.evaluations;if(!r.valid)return false;}
        if(*p=='L') {
            const double dl=-r.rC/r.A;
            z.x+=r.n*(dl/mass);z.lambda+=dl;++d.core.scalar_solves;
        } else {
            ++d.energy_attempts;
            if(r.bt==0&&norm(r.g)<=gradient_resolution_tolerance(c,z,r)) {
                ++d.core.singular;
                if(std::abs(r.rE)>energy_tolerance(c,r))++d.core.incompatible;
                continue;
            }
            double de=-r.rE/r.D;const V2 extra=r.g*(de/mass);
            if(acceptance==Acceptance::AtEnergyVisit)
                de*=filter(c,z,r,o,extra,r.rE,dot(r.g,extra)+r.bt*de,d);
            z.x+=r.g*(de/mass);z.eta+=de;++d.core.scalar_solves;
        }
        if(!finite(z.x)||!std::isfinite(z.lambda)||!std::isfinite(z.eta))return false;
    }
    if(acceptance==Acceptance::WholeProposal) {
        const V2 extra=z.x-baseline.x;
        const double H0=initial.rE+dot(initial.g,baseline.x-start.x);
        const double change=dot(initial.g,extra)+initial.bt*(z.eta-start.eta);
        const double blend=filter(c,start,initial,o,extra,H0,change,d);
        z.x=baseline.x+extra*blend;
        z.lambda=baseline.lambda+blend*(z.lambda-baseline.lambda);
        z.eta=baseline.eta+blend*(z.eta-baseline.eta);
    }
    return finite(z.x)&&std::isfinite(z.lambda)&&std::isfinite(z.eta);
}
} // namespace refreshed_corrections
