#pragma once
#include "correction.hpp"
#include "refreshed.hpp"

namespace symplectic_energy {
using namespace energy_ordering;
namespace cs=correction_schedules;
enum class Mode { XPBD, Predictor, Energy, Paired };
struct Options {
    Mode mode=Mode::Energy;
    cs::Options paired;
    double fraction=1,cap=0;
    bool reject_large=false,direction=false;
};
struct Diagnostics {
    cs::Diagnostics paired;
    long energy_solves=0,singular=0,caps=0,skips=0;
    double max_ratio=0,max_prediction_error=0,last_H=0,last_predicted_H=0,last_blend=0,last_extra=0;
};
inline double constitutive_lambda(const Model& m,V2 x,double h) {
    return -h*h*(norm(x)-m.length)/m.alpha;
}
inline bool step(const Context& c,double& stored,const Options& o,State& result,Diagnostics& d,long index=0) {
    if(c.model.pendulum)return false; // Symplectic proof here is for elastic forces.
    Iterate z{c.predictor};
    if(o.mode!=Mode::XPBD) {
        z.lambda=o.mode==Mode::Paired?stored:constitutive_lambda(c.model,c.old.x,c.h);
        z.x+=c.old.x/norm(c.old.x)*(z.lambda/c.model.mass);
    }
    if(o.mode==Mode::XPBD||o.mode==Mode::Paired) {
        auto opt=o.paired;if(o.mode==Mode::XPBD)opt.schedule=cs::Schedule::Baseline;
        if(!cs::update(c,z,opt,index,d.paired))return false;
        d.last_H=d.paired.last_H0;d.last_predicted_H=d.paired.last_predicted_H;
        d.last_blend=d.paired.last_blend;d.last_extra=d.paired.last_extra;
    } else if(o.mode==Mode::Energy) {
        const Rows r=evaluate(c,z);if(!r.valid)return false;
        const double position_bound=4*std::numeric_limits<double>::epsilon()*c.model.mass/(c.h*c.h)
            *(norm(z.x)+norm(c.old.x+c.model.acceleration*(c.h*c.h)));
        d.last_H=r.H;d.last_blend=0;d.last_extra=0;d.last_predicted_H=r.H;
        if(norm(r.g)<=gradient_tolerance(r)+position_bound)++d.singular;
        else {
            ++d.energy_solves;double de=-r.H/r.D;const V2 raw=r.g*(de/c.model.mass);
            double blend=o.fraction;const double motion=gradient_policies::motion_scale(c);
            if(o.direction&&dot(gradient_policies::objective_gradient(c,z,r),raw)>0){blend=0;++d.skips;}
            if(o.cap>0&&blend*norm(raw)>o.cap*motion) {
                if(o.reject_large){blend=0;++d.skips;}
                else {blend=o.cap*motion/norm(raw);++d.caps;}
            }
            de*=blend;z.x+=r.g*(de/c.model.mass);z.eta+=de;
            const double size=norm(raw)*blend;d.last_blend=blend;d.last_extra=size;
            if(motion>0)d.max_ratio=std::max(d.max_ratio,size/motion);
            d.last_predicted_H=r.H+dot(r.g,raw)*blend;
            const double actual=energy(c.model,{z.x,(z.x-c.old.x)/c.h})-c.target;
            d.max_prediction_error=std::max(d.max_prediction_error,std::abs(actual-d.last_predicted_H));
        }
    }
    result={z.x,(z.x-c.old.x)/c.h};
    // For predictor/E-only methods, this cache represents the force at the
    // accepted position. It is not an extra position/multiplier KKT solve.
    stored=o.mode==Mode::Paired?z.lambda:constitutive_lambda(c.model,z.x,c.h);
    return finite(result.x)&&finite(result.v)&&std::isfinite(stored);
}
} // namespace symplectic_energy
