#pragma once
#include "common.hpp"
#include "traditional_cases/relinearized.hpp"
#include <limits>

// Existing EL/LE solvers with a residual-complete extension. A frozen sweep
// freezes the linearization for that sweep, not for the entire time step.
// The bounded variant retains the historical repair test and motion cap;
// it is an accepted partial correction, not an exact constrained minimizer.
namespace compose {
namespace el_detail {
inline traditional::Model model(const Model&m){
    if(m.count>2||m.dim>2)throw std::runtime_error("EL_adapter_dimension");
    traditional::Model a;a.count=m.count;
    for(int i=0;i<m.count;++i){a.mass[i]=m.mass[i];a.length[i]=m.length[i];a.alpha[i]=m.rigid?0:1/m.stiffness[i];}
    a.gravity={m.gravity[0],m.dim==2?m.gravity[1]:0};return a;
}
inline traditional::State state(const Model&m,const State&s){traditional::State a;for(int i=0;i<m.count;++i){a.x[i]={s.x[i*m.dim],m.dim==2?s.x[i*m.dim+1]:0};a.v[i]={s.v[i*m.dim],m.dim==2?s.v[i*m.dim+1]:0};}return a;}
inline State state(const Model&m,const traditional::State&s){State a{V::Zero(dof(m)),V::Zero(dof(m))};for(int i=0;i<m.count;++i){a.x[i*m.dim]=s.x[i].x;a.v[i*m.dim]=s.v[i].x;if(m.dim==2){a.x[i*m.dim+1]=s.x[i].y;a.v[i*m.dim+1]=s.v[i].y;}}return a;}
inline bool zero_energy_gradient(const traditional::Context&c,const traditional::Iterate&z,const traditional::Rows&r){
    double tol=64*std::numeric_limits<double>::epsilon()*r.gradient_scale;
    for(int i=0;i<c.model.count;++i)tol+=4*std::numeric_limits<double>::epsilon()*c.model.mass[i]/(c.h*c.h)*(traditional::norm(z.x[i])+traditional::norm(c.old.x[i]+c.model.gravity*(c.h*c.h)));
    return std::sqrt(traditional::inner(r.g,r.g,c.model.count))<=tol;
}
inline bool retained_sweep(const traditional::Context&c,traditional::Iterate&z,const Config&cfg,traditional::Diagnostics&stats){
    using namespace traditional;const auto&m=c.model;const int n=m.count;const Iterate start=z;Rows initial=evaluate(c,start);++stats.evaluations;if(!initial.valid)return false;
    std::array<double,3> rhs{},base{},delta{};
    for(int j=0;j<=n;++j){const Vec&g=j==n?initial.g:initial.J[j];rhs[j]=-(j==n?initial.H:initial.rC[j])+weighted(m,g,initial.R);}
    // Material-only baseline at this same iterate, retaining the same R.
    for(int j=0;j<n;++j){double b=rhs[j];for(int k=0;k<n;++k)if(k!=j)b-=initial.S[j][k]*base[k];base[j]=b/initial.S[j][j];}
    Vec dx0=displacement(m,initial,base),Rinv=minv(m,initial.R);for(int i=0;i<n;++i)dx0[i]+=Rinv[i]*(-1);
    Iterate proposal=start;
    auto visit=[&](int j){
        if(cfg.relinearize){Rows r=evaluate(c,proposal);++stats.evaluations;if(!r.valid)return false;
            if(j==n&&zero_energy_gradient(c,proposal,r)){++stats.singular;return true;}
            const Vec&g=j==n?r.g:r.J[j];const double den=weighted(m,g,g)+(j==n?0:r.at[j]);
            const double d=(-(j==n?r.H:r.rC[j])+weighted(m,g,r.R))/den;
            for(int i=0;i<n;++i)proposal.x[i]+=(g[i]*d-r.R[i])/m.mass[i];
            if(j==n)proposal.eta+=d;else proposal.lambda[j]+=d;
        }else{
            if(j==n&&zero_energy_gradient(c,start,initial)){++stats.singular;return true;}
            double b=rhs[j];for(int k=0;k<=n;++k)if(k!=j)b-=initial.S[j][k]*delta[k];delta[j]=b/initial.S[j][j];
        }
        ++stats.scalar_solves;return true;
    };
    bool ok=true;if(cfg.schedule=="EL")ok=visit(n);for(int a=0;a<n&&ok;++a){int j=cfg.reverse?n-1-a:a;ok=visit(j);}if(cfg.schedule=="LE"&&ok)ok=visit(n);if(!ok)return false;
    if(!cfg.relinearize){Vec dx=displacement(m,initial,delta);for(int i=0;i<n;++i)proposal.x[i]+=dx[i]-Rinv[i];for(int j=0;j<n;++j)proposal.lambda[j]+=delta[j];proposal.eta+=delta[n];}
    Vec extra{};for(int i=0;i<n;++i)extra[i]=proposal.x[i]-start.x[i]-dx0[i];double blend=1;
    if(cfg.treatment=="bounded"){
        double H0=initial.H+inner(initial.g,dx0,n),change=inner(initial.g,extra,n);
        double tol=128*std::numeric_limits<double>::epsilon()*(1+std::abs(c.target)+std::abs(initial.H));
        if(std::abs(H0)<=tol||H0*change>=0||change==0){blend=0;++stats.skips;}
        else if(std::abs(change)>.5*std::abs(H0)){blend=.5*std::abs(H0)/std::abs(change);++stats.repair_scaled;}
        double size=mass_norm(m,extra),scale=motion_scale(c);if(blend*size>.1*scale){blend=.1*scale/size;++stats.caps;}
        if(scale>0)stats.max_ratio=std::max(stats.max_ratio,blend*size/scale);
    }
    for(int i=0;i<n;++i)z.x[i]=start.x[i]+dx0[i]+extra[i]*blend;
    for(int j=0;j<n;++j)z.lambda[j]=start.lambda[j]+base[j]+blend*(proposal.lambda[j]-start.lambda[j]-base[j]);
    z.eta=start.eta+blend*(proposal.eta-start.eta);
    return valid(z.x,n)&&std::isfinite(z.eta);
}
inline bool sweep(const traditional::Context&c,traditional::Iterate&z,const Config&cfg,traditional::Diagnostics&d){
    if(cfg.keep_rx)return retained_sweep(c,z,cfg,d);
    traditional::Options o;o.schedule=cfg.schedule=="EL"?traditional::Schedule::EL:traditional::Schedule::LE;o.policy=cfg.treatment=="bounded"?traditional::Policy::RepairCap:traditional::Policy::Raw;
    return cfg.relinearize?traditional::update_relinearized(c,z,o,d):traditional::update(c,z,o,d);
}
}
inline StepResult el_step(const Model&m,const State&old,double h,const Config&cfg,double scale){
    StepResult out;out.state=old;
    try{
        auto tm=el_detail::model(m);auto ts=el_detail::state(m,old);auto c=traditional::context(tm,ts,h,energy(m,old));traditional::Iterate z{c.predictor};traditional::Diagnostics d;
        if(m.dim==1&&c.predictor[0].x<=0){out.valid=false;out.reason="EL_distance_adapter_predictor_crossed_anchor";return out;}
        for(int it=0;it<cfg.iterations;++it){if(!el_detail::sweep(c,z,cfg,d)){out.valid=false;out.reason="EL_nonfinite_update";break;}++out.iterations;if(m.dim==1&&z.x[0].x<=0){out.valid=false;out.reason="EL_distance_adapter_crossed_anchor";break;}}
        traditional::State end;end.x=z.x;for(int i=0;i<m.count;++i)end.v[i]=(z.x[i]-ts.x[i])/h;
        out.state=el_detail::state(m,end);out.valid=out.valid&&finite(out.state);
        if(m.dim==1&&out.state.x[0]<=0){out.valid=false;out.reason="EL_distance_adapter_crossed_anchor";}
        auto r=traditional::evaluate(c,z);out.balance=std::abs(r.H)/std::max(1e-12,scale);out.rx=traditional::mass_norm(tm,traditional::minv(tm,r.R));out.visits=d.scalar_solves;out.evaluations=d.evaluations;
        out.kinematic=(out.state.v-(out.state.x-old.x)/h).norm();
        V momentum=masses(m).cwiseProduct(out.state.v-old.v)-h*force(m,out.state.x);
        if(m.rigid)for(int i=0;i<m.count;++i){V reaction=V::Zero(dof(m));for(int j=0;j<m.count;++j){reaction[i*m.dim]=r.J[j][i].x*z.lambda[j]/h;if(m.dim==2)reaction[i*m.dim+1]=r.J[j][i].y*z.lambda[j]/h;}momentum.segment(i*m.dim,m.dim)-=reaction.segment(i*m.dim,m.dim);}
        out.impulse=momentum.norm();out.work_defect=out.balance;
    }catch(const std::exception&e){out.valid=false;out.reason=e.what();}return out;
}
}
