#pragma once
#include "solver.hpp"

namespace energy_ordering {
enum class Schedule { Baseline, Joint, Jacobi, LE, EL, LEL, ELE };
enum class Application { Deferred, Immediate };
inline const char* name(Schedule s) {
    switch(s) {
    case Schedule::Baseline: return "baseline";
    case Schedule::Joint: return "joint";
    case Schedule::Jacobi: return "jacobi";
    case Schedule::LE: return "gs_le";
    case Schedule::EL: return "gs_el";
    case Schedule::LEL: return "gs_lel";
    case Schedule::ELE: return "gs_ele";
    }
    return "unknown";
}

// Solve one frozen linearization. Reevaluate only on the next outer call.
// L and E are scalar coordinate solves for the TOTAL multiplier increments.
inline bool schedule_update(const Context& c,Iterate& z,Schedule schedule,
    bool retain_rx,Counters& count,int passes=1,
    Application application=Application::Deferred) {
    if(schedule==Schedule::Baseline) return material_update(c,z,count);
    if(schedule==Schedule::Joint) {
        Options options; options.method=Method::Joint; options.retain_rx=retain_rx;
        return frozen_update(c,z,options,count);
    }
    const Rows r=evaluate(c,z); ++count.evaluations;
    if(!r.valid) return false;
    const double mass=c.model.mass;
    const V2 base=retain_rx?r.R*(-1/mass):V2{};
    const double fL=-r.rC-dot(r.n,base), fE=-r.rE-dot(r.g,base);
    const bool zero=energy_row_zero(r);
    if(zero||(r.at==0&&r.bt==0&&std::abs(r.q)<=gradient_tolerance(r)))
        ++count.rank_deficient_blocks;
    double dl=0,de=0;
    const V2 original=z.x;
    if(application==Application::Immediate) z.x+=base;
    auto material=[&]() {
        const double next=(fL-r.B*de)/r.A;
        if(application==Application::Immediate) z.x+=r.n*((next-dl)/mass);
        dl=next; ++count.scalar_solves;
    };
    auto energy=[&]() {
        if(zero) {
            ++count.singular;
            if(std::abs(fE-r.B*dl)>energy_tolerance(c,r)) ++count.incompatible;
            return;
        }
        const double next=(fE-r.B*dl)/r.D;
        if(application==Application::Immediate) z.x+=r.g*((next-de)/mass);
        de=next; ++count.scalar_solves;
    };
    for(int pass=0;pass<passes;++pass) {
        if(schedule==Schedule::Jacobi) {
            // Both rows read the OLD total increments. Application order does
            // not enter either solve; only GS feeds a new value to the next row.
            const double next_dl=(fL-r.B*de)/r.A;
            double next_de=de;
            ++count.scalar_solves;
            if(zero) {
                ++count.singular;
                if(std::abs(fE-r.B*dl)>energy_tolerance(c,r)) ++count.incompatible;
            } else {
                next_de=(fE-r.B*dl)/r.D;
                ++count.scalar_solves;
            }
            if(application==Application::Immediate) {
                z.x+=r.n*((next_dl-dl)/mass);
                z.x+=r.g*((next_de-de)/mass);
            }
            dl=next_dl; de=next_de;
        } else if(schedule==Schedule::LE||schedule==Schedule::LEL) {
            material(); energy();
            if(schedule==Schedule::LEL) material();
        } else {
            energy(); material();
            if(schedule==Schedule::ELE) energy();
        }
    }
    if(application==Application::Deferred) z.x=original+base+(r.n*dl+r.g*de)/mass;
    z.lambda+=dl; z.eta+=de;
    return finite(z.x)&&std::isfinite(z.lambda)&&std::isfinite(z.eta);
}
} // namespace energy_ordering
