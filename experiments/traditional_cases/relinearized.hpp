#pragma once
#include "solver.hpp"

namespace traditional {
// Each visit solves a linearized row at the current, corrected position.
// Like the frozen study, this uses the curvature-dropped, r_x-dropped update.
inline bool update_relinearized(const Context& c,Iterate& z,const Options& o,Diagnostics& stats){
    const auto& m=c.model;const int n=m.count;const Iterate start=z;
    const Rows initial=evaluate(c,start);++stats.evaluations;if(!initial.valid)return false;
    const auto baseline=material_pass(m,initial,{});
    const Vec dx0=displacement(m,initial,baseline);
    Iterate proposal=start;
    auto L=[&](){
        for(int j=0;j<n;++j){
            const Rows r=evaluate(c,proposal);++stats.evaluations;if(!r.valid)return false;
            const double dl=-r.rC[j]/r.S[j][j];
            for(int i=0;i<n;++i)proposal.x[i]+=r.J[j][i]*(dl/m.mass[i]);
            proposal.lambda[j]+=dl;++stats.scalar_solves;
        }return valid(proposal.x,n);
    };
    auto E=[&](){
        const Rows r=evaluate(c,proposal);++stats.evaluations;if(!r.valid)return false;
        // A material solve can make the true energy gradient exactly zero.
        // Include the uncertainty from representing positions in doubles.
        double tolerance=64*std::numeric_limits<double>::epsilon()*r.gradient_scale;
        for(int i=0;i<n;++i){const V2 center=c.old.x[i]+m.gravity*(c.h*c.h);
            tolerance+=4*std::numeric_limits<double>::epsilon()*m.mass[i]/(c.h*c.h)
                *(norm(proposal.x[i])+norm(center));}
        if(std::sqrt(inner(r.g,r.g,n))<=tolerance){++stats.singular;return true;}
        const double de=-r.H/r.S[n][n];
        for(int i=0;i<n;++i)proposal.x[i]+=r.g[i]*(de/m.mass[i]);
        proposal.eta+=de;++stats.scalar_solves;return valid(proposal.x,n)&&std::isfinite(proposal.eta);
    };
    bool ok=false;
    if(o.schedule==Schedule::EL)ok=E()&&L();
    else if(o.schedule==Schedule::LE)ok=L()&&E();
    else throw std::runtime_error("refreshed study supports EL or LE only");
    if(!ok)return false;
    Vec extra{};for(int i=0;i<n;++i)extra[i]=proposal.x[i]-start.x[i]-dx0[i];
    const double H0=initial.H+inner(initial.g,dx0,n),change=inner(initial.g,extra,n);
    double blend=1;
    if(o.policy==Policy::RepairCap){
        const double tol=128*std::numeric_limits<double>::epsilon()*(1+std::abs(c.target)+std::abs(initial.H));
        if(std::abs(H0)<=tol||H0*change>=0||change==0){blend=0;++stats.skips;}
        else if(std::abs(change)>o.fraction*std::abs(H0)){blend=o.fraction*std::abs(H0)/std::abs(change);++stats.repair_scaled;}
        const double size=mass_norm(m,extra),scale=motion_scale(c);
        if(blend*size>o.cap*scale){blend=o.cap*scale/size;++stats.caps;}
        if(scale>0)stats.max_ratio=std::max(stats.max_ratio,blend*size/scale);
    }else if(o.policy!=Policy::Raw)throw std::runtime_error("refreshed study supports raw or repair-cap only");
    // Blend actual reconstructed positions. A fresh gradient cannot be replaced
    // afterward by the original frozen gradient in the displacement formula.
    for(int i=0;i<n;++i)z.x[i]=start.x[i]+dx0[i]+extra[i]*blend;
    for(int j=0;j<n;++j)z.lambda[j]=start.lambda[j]+baseline[j]+blend*(proposal.lambda[j]-start.lambda[j]-baseline[j]);
    z.eta=start.eta+blend*(proposal.eta-start.eta);
    stats.last_H0=H0;stats.last_blend=blend;stats.last_extra=blend*mass_norm(m,extra);
    Vec accepted{};for(int i=0;i<n;++i)accepted[i]=z.x[i]-start.x[i];
    stats.last_predicted_H=initial.H+inner(initial.g,accepted,n);
    State endpoint;endpoint.x=z.x;for(int i=0;i<n;++i)endpoint.v[i]=(z.x[i]-c.old.x[i])/c.h;
    const double actual=energy(m,endpoint)-c.target;
    stats.max_prediction_error=std::max(stats.max_prediction_error,std::abs(actual-stats.last_predicted_H));
    return valid(z.x,n)&&std::isfinite(actual)&&std::isfinite(z.eta);
}
} // namespace traditional
