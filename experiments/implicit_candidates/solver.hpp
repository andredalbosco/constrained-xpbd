#pragma once
#include "../research_agenda/solver.hpp"
#include <string>
#include <vector>
namespace implicit {
using namespace agenda;
enum class Family {Baseline,Energy,Work};
struct Config {
    std::string name;Family family=Family::Energy;double strength=10,work_strength=1,repay=0,cap=.1;
    int outer=1,inner=1,base_passes=1;bool EL=false,angular=false,stiffness_gate=false,local_target=false;double angular_softness=0;
};
struct StepReport {
    long visits=0,evaluations=0,backtracks=0;bool fallback=false,failed=false;double blend=1,energy=0,angular=0,ledger=0,cap_ratio=0;
    double proposal_energy=0;int base_passes=0;
};
struct AngularRow {Vec gradient{};double value=0;};
inline double angular(const Model& m,const State& s){double v=0;for(int i=0;i<m.count;++i)v+=m.mass[i]*cross(s.x[i],s.v[i]);return v;}
inline double torque(const Model& m,const Vec& x){double v=0;for(int i=0;i<m.count;++i)v+=m.mass[i]*cross(x[i],m.gravity);return v;}
inline AngularRow angular_row(const Context& c,const Vec& x){AngularRow r;const auto& m=c.model;
    const auto s=endpoint(c,x);r.value=angular(m,s)-angular(m,c.old)-.5*c.h*(torque(m,x)+torque(m,c.old.x));
    for(int i=0;i<m.count;++i)r.gradient[i]=tangent(c.old.x[i])*(m.mass[i]/c.h)+tangent(m.gravity)*(.5*c.h*m.mass[i]);return r;
}
inline double stiffness_number(const Context& c,const Vec& x){const auto r=evaluate(c,Iterate{x});double s=0;
    for(int j=0;j<c.model.count;++j)if(c.model.alpha[j]>0)s=std::max(s,c.h*c.h*weighted(c.model,r.J[j],r.J[j])/c.model.alpha[j]);return s;
}
struct Proposal {Vec x{};std::array<double,2> lambda{};Array z{},t{};};
inline Proposal material_baseline(const Context& c,int passes,StepReport& stats){Iterate z{c.predictor};Diagnostics d;for(int k=0;k<passes;++k){if(!fresh_material(c,z,d))throw std::runtime_error("invalid baseline");}
    stats.visits+=d.scalar_solves;stats.evaluations+=d.evaluations;stats.base_passes+=passes;Proposal p;p.x=z.x;p.lambda=z.lambda;return p;
}
inline Proposal energy_proposal(const Context& c,const Config& cfg,double scale,StepReport& stats){
    const auto& m=c.model;const int n=m.count,ne=n,na=n+1,nr=n+1+int(cfg.angular);Iterate current{c.predictor};double etaL=0;
    const double strength=cfg.strength/(cfg.stiffness_gate?1+stiffness_number(c,c.predictor):1);
    const double beta=scale/(c.h*c.h*strength),target=cfg.local_target?energy(m,c.old):c.target;
    for(int outer=0;outer<cfg.outer;++outer){auto r=evaluate(c,current);++stats.evaluations;if(!r.valid)throw std::runtime_error("invalid energy row");
        Dense S{};Array rhs{},delta{};std::array<Vec,N> g{};Vec R=r.R;const auto l=angular_row(c,current.x);
        for(int i=0;i<n;++i)R[i]+=l.gradient[i]*(-etaL);
        for(int j=0;j<n;++j)g[j]=r.J[j];g[ne]=r.g;if(cfg.angular)g[na]=l.gradient;
        for(int i=0;i<nr;++i){const double residual=i<n?r.rC[i]:i==ne?r.H+c.target-target+beta*current.eta:l.value;
            rhs[i]=-residual+weighted(m,g[i],R);
            for(int j=0;j<nr;++j)S[i][j]=weighted(m,g[i],g[j]);if(i<n)S[i][i]+=r.at[i];else if(i==ne)S[i][i]+=beta;}
        if(cfg.angular){const double compliance=cfg.angular_softness*S[na][na];S[na][na]+=compliance;rhs[na]-=compliance*etaL;}
        std::array<int,N> order{};int count=0;if(cfg.EL)order[count++]=ne;for(int j=0;j<n;++j)order[count++]=j;if(!cfg.EL)order[count++]=ne;if(cfg.angular)order[count++]=na;
        for(int pass=0;pass<cfg.inner;++pass)for(int q=0;q<count;++q){const int i=order[q];double v=rhs[i];for(int j=0;j<nr;++j)if(i!=j)v-=S[i][j]*delta[j];if(!(S[i][i]>0))continue;delta[i]=v/S[i][i];++stats.visits;}
        for(int i=0;i<n;++i){V2 dx=R[i]*(-1);for(int j=0;j<nr;++j)dx+=g[j][i]*delta[j];current.x[i]+=dx/m.mass[i];current.lambda[i]+=delta[i];}
        current.eta+=delta[ne];if(cfg.angular)etaL+=delta[na];
    }
    Proposal p;p.x=current.x;p.lambda=current.lambda;p.z[ne]=current.eta;if(cfg.angular)p.z[na]=etaL;return p;
}
inline Proposal work_proposal(const Context& c,const Config& cfg,double scale,const Reservoir& reservoir,StepReport& stats){
    Proposal p;p.x=c.predictor;const auto& m=c.model;const int n=m.count;WorkOptions o;o.strength=cfg.work_strength;o.repay_time=cfg.repay;
    const double beta=scale/(c.h*c.h*cfg.strength);
    for(int outer=0;outer<cfg.outer;++outer){const auto a=assemble_work(c,p.x,p.t,p.z,reservoir,o,scale);++stats.evaluations;const int nr=a.nr+int(cfg.angular);Dense S=a.S;Dense G=a.G;Array rhs{},primal=a.primal,delta{},residual=a.residual;
        const auto l=angular_row(c,p.x);
        if(cfg.angular){residual[a.nr]=l.value;for(int i=0;i<n;++i){G[a.nr][2*i]=l.gradient[i].x;G[a.nr][2*i+1]=l.gradient[i].y;primal[2*i]-=l.gradient[i].x*p.z[a.nr];primal[2*i+1]-=l.gradient[i].y*p.z[a.nr];}
            for(int j=0;j<nr;++j){double v=0;for(int k=0;k<a.nv;++k)v+=G[a.nr][k]*a.inverse[k]*G[j][k];S[a.nr][j]=S[j][a.nr]=v;}
            const double compliance=cfg.angular_softness*S[a.nr][a.nr];S[a.nr][a.nr]+=compliance;residual[a.nr]+=compliance*p.z[a.nr];}
        for(int j=n;j<a.nr;++j){residual[j]+=beta*p.z[j];S[j][j]+=beta;}
        for(int j=0;j<nr;++j){rhs[j]=-residual[j];for(int k=0;k<a.nv;++k)rhs[j]+=G[j][k]*a.inverse[k]*primal[k];}
        for(int pass=0;pass<cfg.inner;++pass)for(int i=0;i<nr;++i){double v=rhs[i];for(int j=0;j<nr;++j)if(j!=i)v-=S[i][j]*delta[j];if(!(S[i][i]>0))continue;delta[i]=v/S[i][i];++stats.visits;}
        Array dx{};for(int k=0;k<a.nv;++k){double v=-primal[k];for(int j=0;j<nr;++j)v+=G[j][k]*delta[j];dx[k]=a.inverse[k]*v;}
        for(int i=0;i<n;++i){p.x[i]+={dx[2*i],dx[2*i+1]};p.lambda[i]+=delta[i];}for(int e=0;e<a.graph.ne;++e)p.t[e]+=dx[2*n+e];for(int j=0;j<nr;++j)p.z[j]+=delta[j];
    }return p;
}
inline void project_rods(const Context& c,Proposal& p){
    for(int j=0;j<c.model.count;++j)if(c.model.alpha[j]==0){const V2 e=p.x[j]-(j?p.x[j-1]:V2{});const double length=norm(e);if(!(length>1e-14))throw std::runtime_error("collapsed rod");const V2 normal=e/length;const double w=1/c.model.mass[j]+(j?1/c.model.mass[j-1]:0);const double dl=-(length-c.model.length[j])/w;
        p.x[j]+=normal*(dl/c.model.mass[j]);if(j)p.x[j-1]+=normal*(-dl/c.model.mass[j-1]);p.lambda[j]+=dl;}
}
inline bool admissible(const Context& c,const Proposal& p,double scale){
    const auto next=endpoint(c,p.x);if(!valid(next.x,c.model.count)||!valid(next.v,c.model.count))return false;const double E=energy(c.model,next);
    if(!std::isfinite(E)||E>c.target+1e-10*scale)return false;
    for(int j=0;j<c.model.count;++j)if(c.model.alpha[j]==0&&std::abs(norm(p.x[j]-(j?p.x[j-1]:V2{}))-c.model.length[j])>1e-10*c.model.length[j])return false;return true;
}
inline void ledger(const Context& c,const Proposal& p,Reservoir& reservoir,StepReport& report,double scale){
    const auto l=layout(c.model,c.old);const auto w=work_values(c,p.x,Estimate::Lagged);const double before=energy(c.model,c.old)+total(reservoir);
    for(int i=0;i<c.model.count;++i){double b=w.dk[i].v;for(int e=0;e<l.ne;++e)if(l.edge[e].body==i)b-=w.w[e].v+p.t[e];reservoir.k[i]-=b;}
    for(int j=0;j<l.np;++j){double b=w.du[j].v;for(int e=0;e<l.ne;++e)if(l.edge[e].potential==j)b+=w.w[e].v+p.t[e];reservoir.u[j]-=b;}
    report.ledger=std::abs(energy(c.model,endpoint(c,p.x))+total(reservoir)-before)/scale;
}
inline bool advance(const Context& c,const Config& cfg,double scale,Reservoir& reservoir,State& next,StepReport& report){
    Proposal base;try{base=material_baseline(c,cfg.base_passes,report);}catch(const std::exception&){report.failed=true;return false;}
    if(cfg.family==Family::Baseline){if(!admissible(c,base,scale)){report.failed=true;return false;}next=endpoint(c,base.x);report.energy=energy(c.model,next);report.proposal_energy=report.energy;report.angular=angular_row(c,base.x).value;ledger(c,base,reservoir,report,scale);return true;}
    Proposal proposal=base;bool good=true;
    if(cfg.family!=Family::Baseline)try{proposal=cfg.family==Family::Energy?energy_proposal(c,cfg,scale,report):work_proposal(c,cfg,scale,reservoir,report);}catch(const std::exception&){good=false;}
    report.proposal_energy=good?energy(c.model,endpoint(c,proposal.x)):std::numeric_limits<double>::quiet_NaN();
    Vec extra{};for(int i=0;i<c.model.count;++i)extra[i]=proposal.x[i]-base.x[i];const double size=mass_norm(c.model,extra),motion=motion_scale(c);double blend=cfg.cap>0&&size>cfg.cap*motion?cfg.cap*motion/size:1;
    Proposal accepted=base;bool found=false;
    if(good)for(int attempt=0;attempt<20;++attempt){Proposal trial=proposal;for(int i=0;i<c.model.count;++i){trial.x[i]=base.x[i]+extra[i]*blend;trial.lambda[i]=base.lambda[i]+(proposal.lambda[i]-base.lambda[i])*blend;}
        for(int j=0;j<N;++j){trial.t[j]*=blend;trial.z[j]*=blend;}
        try{project_rods(c,trial);}catch(const std::exception&){blend*=.5;++report.backtracks;continue;}
        bool angular_ok=true;if(cfg.angular){const double original=std::abs(angular_row(c,base.x).value),candidate=std::abs(angular_row(c,trial.x).value);angular_ok=candidate<=original+1e-12*(1+std::abs(angular(c.model,c.old)));}
        Vec displacement{};for(int i=0;i<c.model.count;++i)displacement[i]=trial.x[i]-base.x[i];const bool cap_ok=cfg.cap<=0||mass_norm(c.model,displacement)<=cfg.cap*motion*(1+1e-8)+1e-15;
        if(admissible(c,trial,scale)&&angular_ok&&cap_ok){accepted=trial;found=true;break;}blend*=.5;++report.backtracks;
    }
    if(!found){if(!admissible(c,base,scale)){report.failed=true;return false;}accepted=base;blend=0;report.fallback=cfg.family!=Family::Baseline;}
    report.blend=blend;Vec actual{};for(int i=0;i<c.model.count;++i)actual[i]=accepted.x[i]-base.x[i];report.cap_ratio=motion>0?mass_norm(c.model,actual)/motion:0;
    next=endpoint(c,accepted.x);report.energy=energy(c.model,next);report.angular=angular_row(c,accepted.x).value;
    if(cfg.family==Family::Work||cfg.family==Family::Baseline)ledger(c,accepted,reservoir,report,scale);return true;
}
inline std::vector<Config> configurations(){
    std::vector<Config> out;Config b;b.name="xpbd";b.family=Family::Baseline;out.push_back(b);b.name="xpbd_4";b.base_passes=4;out.push_back(b);
    for(bool EL:{true,false})for(double strength:{.1,1.,10.,100.}){Config c;c.name=std::string(EL?"EL":"LE")+"_soft_"+(strength==.1?"01":std::to_string(int(strength)));c.EL=EL;c.strength=strength;out.push_back(c);}
    for(bool EL:{true,false}){Config c;c.name=std::string(EL?"EL":"LE")+"_stiffness_gate";c.EL=EL;c.strength=10;c.stiffness_gate=true;out.push_back(c);
        c.stiffness_gate=false;c.inner=4;c.name=std::string(EL?"EL":"LE")+"_GS4";out.push_back(c);c.inner=1;c.outer=4;c.name=std::string(EL?"EL":"LE")+"_outer4";out.push_back(c);}
    for(bool EL:{true,false})for(double strength:{1.,10.}){Config c;c.EL=EL;c.angular=true;c.strength=strength;c.name=std::string(EL?"EL":"LE")+"_angular_"+std::to_string(int(strength));out.push_back(c);}
    Config local;local.name="LE_local_target";local.local_target=true;out.push_back(local);
    for(double strength:{1.,10.,100.}){Config c;c.family=Family::Work;c.strength=strength;c.name="work_soft_"+std::to_string(int(strength));out.push_back(c);}
    Config c;c.family=Family::Work;c.name="work_angular";c.angular=true;out.push_back(c);c.name="work_repay";c.angular=false;c.repay=.25;out.push_back(c);c.name="work_repay_angular";c.angular=true;out.push_back(c);
    for(double strength:{1.,10.,100.})for(double tau:{.05,.25,1.}){if(strength==10&&tau==.25)continue;Config f=c;f.strength=strength;f.repay=tau;f.name="workL_s"+std::to_string(int(strength))+"_tau"+(tau==.05?"005":tau==.25?"025":"100");out.push_back(f);}
    for(double softness:{.1,1.,10.}){Config f=c;f.strength=1;f.repay=.05;f.angular_softness=softness;f.name="workL_soft_"+(softness==.1?std::string("01"):std::to_string(int(softness)));out.push_back(f);}
    return out;
}
} // namespace implicit
