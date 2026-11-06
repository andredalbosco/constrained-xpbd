#pragma once
#include "gradient_policies.hpp"

// One frozen EL sweep. The predictor and multipliers are reset each substep.
// These controls modify the energy increment before the material row sees it.
namespace small_steps_el {
using namespace energy_ordering;
enum class Gate { None, Raw, Effective };
struct Options {
    double relaxation=1;
    double cap=-1; // Negative disables capping; otherwise relative to old motion.
    Gate gate=Gate::Raw;
};
struct Diagnostics {
    Counters core;
    long attempts=0,accepted=0,directional_skips=0,capped=0;
    double max_raw_ratio=0,max_accepted_ratio=0,max_raw_displacement=0;
    double max_accepted_objective_work=0;
};
inline bool update(const Context& c,Iterate& z,const Options& o,Diagnostics& d) {
    const Rows r=evaluate(c,z);++d.core.evaluations;
    if(!r.valid)return false;
    const double m=c.model.mass;
    double de=0;
    ++d.attempts;
    if(energy_row_zero(r)) {
        ++d.core.singular;
        if(std::abs(r.rE)>energy_tolerance(c,r))++d.core.incompatible;
    } else {
        ++d.core.scalar_solves;
        const double raw=-r.rE/r.D;
        if(!std::isfinite(raw))return false;
        const V2 raw_dx=r.g*(raw/m);
        const V2 effect=(r.g-r.n*(r.B/r.A))*(raw/m);
        const V2 q=gradient_policies::objective_gradient(c,z,r);
        const V2 tested=o.gate==Gate::Effective?effect:raw_dx;
        const double work=dot(q,tested),scale=gradient_policies::motion_scale(c);
        const double raw_size=norm(raw_dx);
        auto ratio=[&](double size){return size==0?0:(scale>0?std::min(size/scale,1e300):1e300);};
        d.max_raw_ratio=std::max(d.max_raw_ratio,ratio(raw_size));
        d.max_raw_displacement=std::max(d.max_raw_displacement,raw_size);
        if(o.gate!=Gate::None&&work>0) {
            ++d.directional_skips;
        } else {
            de=raw*o.relaxation;
            const double proposed_size=norm(r.g*(de/m));
            if(o.cap>=0&&proposed_size>o.cap*scale) {
                de*=o.cap*scale/proposed_size;
                ++d.capped;
            }
            ++d.accepted;
            d.max_accepted_ratio=std::max(d.max_accepted_ratio,ratio(norm(r.g*(de/m))));
            const V2 accepted_dx=(o.gate==Gate::Effective?(r.g-r.n*(r.B/r.A)):r.g)*(de/m);
            d.max_accepted_objective_work=std::max(d.max_accepted_objective_work,dot(q,accepted_dx));
        }
    }
    const double dl=(-r.rC-r.B*de)/r.A;++d.core.scalar_solves;
    // Match the old frozen-EL summation order exactly when controls are disabled.
    z.x=z.x+V2{}+(r.n*dl+r.g*de)/m;z.lambda+=dl;z.eta+=de;
    return finite(z.x)&&std::isfinite(z.lambda)&&std::isfinite(z.eta);
}
} // namespace small_steps_el
