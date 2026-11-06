#pragma once
#include "solver.hpp"
namespace agenda {
struct StableReport {long attempts=0,accepted=0,rejected=0;int deepest=0;double largest_energy_error=0;};
// Accepted states satisfy a GLOBAL physical-energy band, in addition to local
// work/rod checks. Failure returns false; no violating endpoint is committed.
// This is a bound on accepted states, not a guarantee that retries terminate.
inline bool guarded_midpoint(const Model& m,const State& old,double h,double target,double scale,State& next,StableReport& diagnostics,int depth=0){
    ++diagnostics.attempts;diagnostics.deepest=std::max(diagnostics.deepest,depth);Report report;State candidate;bool valid_step=true;
    try{candidate=midpoint_step(context(m,old,h,target),4,true,report,scale);}catch(const std::exception&){valid_step=false;}
    const double energy_error=valid_step?std::abs(energy(m,candidate)-target)/scale:std::numeric_limits<double>::infinity();
    valid_step=valid_step&&valid(candidate.x,m.count)&&valid(candidate.v,m.count)&&std::isfinite(energy_error)&&energy_error<=1e-6&&report.work_defect<=1e-9;
    // Energy/work alone can accept a frozen but nonstationary state. Also
    // require the mechanical equations, using a spatial residual tolerance.
    const double spatial_scale=1+mass_norm(m,old.x);
    valid_step=valid_step&&report.rx<=1e-9*spatial_scale&&report.constitutive<=1e-9*spatial_scale;
    for(int j=0;j<m.count;++j)if(m.alpha[j]==0)valid_step=valid_step&&std::abs(norm(candidate.x[j]-(j?candidate.x[j-1]:V2{}))-m.length[j])<=1e-9*m.length[j];
    if(valid_step){next=candidate;++diagnostics.accepted;diagnostics.largest_energy_error=std::max(diagnostics.largest_energy_error,energy_error);return true;}
    ++diagnostics.rejected;if(depth==12)return false;
    State middle,result;if(!guarded_midpoint(m,old,h/2,target,scale,middle,diagnostics,depth+1))return false;
    if(!guarded_midpoint(m,middle,h/2,target,scale,result,diagnostics,depth+1))return false;
    next=result;return true;
}
} // namespace agenda
