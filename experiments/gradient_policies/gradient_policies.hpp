#pragma once
#include "solver.hpp"

// Original endpoint-energy objective and velocity reconstruction. Only the
// approximate linearization and acceptance of the energy increment change.
namespace gradient_policies {
using namespace energy_ordering;
enum class Gradient { Fresh, LagSweep, FrozenStep, PreviousStep };
enum class Order { LE, EL };
inline const char* name(Gradient g) {
    switch(g) {
    case Gradient::Fresh:return "fresh";
    case Gradient::LagSweep:return "lag_sweep";
    case Gradient::FrozenStep:return "frozen_step";
    case Gradient::PreviousStep:return "previous_step";
    }
    return "unknown";
}
inline const char* name(Order o) {return o==Order::EL?"EL":"LE";}
struct Options {
    Order order=Order::EL;
    Gradient gradient=Gradient::Fresh;
    bool retain_rx=true;
    double gamma=-1; // Negative disables displacement-size rejection.
    bool objective_gate=false; // Reject an energy increment with positive grad(F) dot dx_E.
};
struct History {
    V2 frozen{},previous_sweep{},previous_step{};
    bool has_previous_step=false;
    void begin(const Context& c,const Iterate& z) {
        frozen=previous_sweep=evaluate(c,z).g;
    }
    void end(const Context& c,const Iterate& z) {
        previous_step=evaluate(c,z).g;has_previous_step=true;
    }
};
struct Diagnostics {
    Counters core;
    long attempts=0,accepted=0,oversized=0,stale_visits=0,directional_skips=0;
    double max_ratio=0,max_accepted_ratio=0,max_candidate=0;
    double max_accepted_objective_work=0;
};
inline double motion_scale(const Context& c) {
    V2 f=c.model.acceleration*c.model.mass;
    if(!c.model.pendulum) {
        const double radius=norm(c.old.x);
        f+=c.old.x*(-(radius-c.model.length)/(c.model.alpha*radius));
    }
    return c.h*norm(c.old.v)+c.h*c.h*norm(f)/c.model.mass;
}
inline V2 used_gradient(const Rows& r,const History& h,Gradient policy) {
    switch(policy) {
    case Gradient::Fresh:return r.g;
    case Gradient::LagSweep:return h.previous_sweep;
    case Gradient::FrozenStep:return h.frozen;
    case Gradient::PreviousStep:return h.has_previous_step?h.previous_step:r.g;
    }
    return r.g;
}
// Gradient of the original primal implicit objective. This is not the dual
// stationarity residual R: it contains neither lambda nor eta.
inline V2 objective_gradient(const Context& c,const Iterate& z,const Rows& r) {
    const V2 elastic=c.model.pendulum?V2{}:r.n*(r.C/c.model.alpha);
    return (z.x-c.predictor)*(c.model.mass/(c.h*c.h))+elastic;
}
inline bool update(const Context& c,Iterate& z,const Options& o,History& history,Diagnostics& count) {
    const Rows r=evaluate(c,z);++count.core.evaluations;
    if(!r.valid)return false;
    const V2 g=used_gradient(r,history,o.gradient);
    if(o.gradient!=Gradient::Fresh)++count.stale_visits;
    const double m=c.model.mass,B=dot(r.n,g)/m,D=dot(g,g)/m+r.bt;
    // Residuals are always CURRENT and analytical, including R's g*eta term.
    // The old gradient changes only the approximate Jacobian columns/rows.
    const V2 base=o.retain_rx?r.R*(-1/m):V2{};
    const double fL=-r.rC-dot(r.n,base),fE=-r.rE-dot(g,base);
    const bool zero=r.bt==0&&norm(g)<=gradient_tolerance(r);
    double dl=0,de=0;
    auto material=[&]() {dl=(fL-B*de)/r.A;++count.core.scalar_solves;};
    auto energy=[&]() {
        ++count.attempts;
        if(zero) {
            ++count.core.singular;
            if(std::abs(fE-B*dl)>energy_tolerance(c,r))++count.core.incompatible;
            return;
        }
        const double candidate=(fE-B*dl)/D;
        const double displacement=norm(g*(candidate/m));
        const double scale=motion_scale(c);
        // Finite diagnostic cap only; no accepted displacement is clipped.
        const double ratio=displacement==0?0:(scale>0?displacement/scale:1e300);
        count.max_ratio=std::max(count.max_ratio,std::min(ratio,1e300));
        count.max_candidate=std::max(count.max_candidate,std::min(displacement,1e300));
        const V2 dxE=g*(candidate/m);
        const double objective_work=dot(objective_gradient(c,z,r),dxE);
        if(o.objective_gate&&objective_work>0) {
            ++count.directional_skips;return; // Keep neutral/downhill directions only.
        }
        if(o.gamma>=0&&(!std::isfinite(displacement)||displacement>o.gamma*scale)) {
            ++count.oversized;return; // Leave both eta and its displacement alone.
        }
        de=candidate;++count.accepted;++count.core.scalar_solves;
        count.max_accepted_objective_work=std::max(count.max_accepted_objective_work,objective_work);
        count.max_accepted_ratio=std::max(count.max_accepted_ratio,std::min(ratio,1e300));
    };
    if(o.order==Order::EL){energy();material();}else{material();energy();}
    // Preserve the existing deferred GS summation order for the fresh control.
    z.x=z.x+base+(r.n*dl+g*de)/m;z.lambda+=dl;z.eta+=de;
    history.previous_sweep=r.g; // Snapshot before this update; reset at next step.
    return finite(z.x)&&std::isfinite(z.lambda)&&std::isfinite(z.eta);
}
} // namespace gradient_policies
