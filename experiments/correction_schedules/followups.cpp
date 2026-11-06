#define main frozen_study_main
#include "study.cpp"
#undef main
#include "refreshed.hpp"

namespace {
struct FollowConfig {Config cfg;int refresh=0;};
struct FollowStats {
    long completed=0,frozen=0;std::string status="completed";
    double Emax=0,x2=0,v2=0,path=0,refpath=0,torque=0,Lmax=0,maxradius=0;
    se::Diagnostics diagnostic;
};
FollowStats follow_run(const Case& c,const FollowConfig& cfg,int substeps,const std::vector<State>& ref) {
    FollowStats s;State state=c.initial;const double h=1./(60*substeps),E0=energy(c.model,state),L0=angular(state,c.model.mass);
    for(size_t step=1;step<ref.size();++step) {
        const Context ctx=context(c.model,state,h,0,E0);Iterate z{ctx.predictor};
        const bool ok=cfg.refresh?refreshed_corrections::update(ctx,z,cfg.cfg.options,step-1,s.diagnostic,
            cfg.refresh==1?refreshed_corrections::Acceptance::WholeProposal:refreshed_corrections::Acceptance::AtEnergyVisit):
            se::update(ctx,z,cfg.cfg.options,step-1,s.diagnostic);
        const State next{z.x,(z.x-state.x)/h};const double E=energy(c.model,next);
        if(!ok||!finite(next.x)||!finite(next.v)||!std::isfinite(E)) {s.status="undefined";break;}
        if(norm(next.x)>1e6*c.model.length||norm(next.v)>1e6*c.model.length||std::abs(E-E0)>1e12*scale(c)
            ||std::abs(z.lambda)>1e18||std::abs(z.eta)>1e18) {s.status="safety_stop";break;}
        s.completed=step;s.Emax=std::max(s.Emax,std::abs(E-E0));
        s.x2+=dot(next.x-ref[step].x,next.x-ref[step].x);
        s.v2+=dot(next.v-ref[step].v,next.v-ref[step].v);
        s.path+=norm(next.x-state.x);s.refpath+=norm(ref[step].x-ref[step-1].x);
        s.frozen+=norm(next.x-state.x)<1e-10*c.model.length;
        auto torque=[&](V2 x){return c.model.mass*(x.x*c.model.acceleration.y-x.y*c.model.acceleration.x);};
        s.torque+=.5*h*(torque(state.x)+torque(next.x));
        s.Lmax=std::max(s.Lmax,std::abs(angular(next,c.model.mass)-L0-s.torque));
        s.maxradius=std::max(s.maxradius,norm(next.x));state=next;
    }
    return s;
}
std::vector<FollowConfig> follow_configs(const std::string& mode) {
    std::vector<FollowConfig> out;
    for(const auto& c:configurations())
        if(c.name=="baseline"||c.name=="EL_repair_cap"||c.name=="LE_repair_cap"||c.name=="joint_repair_cap")out.push_back({c,0});
    if(mode=="refresh") {
        for(auto s:{se::Schedule::EL,se::Schedule::LE,se::Schedule::LEL,se::Schedule::ELE,se::Schedule::Alternating})
        for(auto p:{se::Policy::Raw,se::Policy::Direction,se::Policy::DirectionCap,se::Policy::RepairCap})
        for(int refresh:{0,1,2}) {
            se::Options o;o.schedule=s;o.policy=p;
            const std::string label=std::string(se::name(s))+"_"+se::name(p)+(refresh==0?"_frozen":refresh==1?"_refresh_whole":"_refresh_local");
            out.push_back({{label,o},refresh});
        }
    } else if(mode=="parameters") {
        for(auto s:{se::Schedule::EL,se::Schedule::LE,se::Schedule::Joint})
        for(double fraction:{.1,.25,.5,.75,1.})for(double cap:{.001,.01,.1,1.}) {
            se::Options o;o.schedule=s;o.policy=se::Policy::RepairCap;o.fraction=fraction;o.cap=cap;
            std::ostringstream label;label<<se::name(s)<<"_f"<<fraction<<"_c"<<cap;
            out.push_back({{label.str(),o},0});
        }
    } else if(mode=="long") {
        for(const auto& c:configurations())if(c.name=="LE_direction_cap"||c.name=="spectral_repair_cap")out.push_back({c,0});
    } else throw std::runtime_error("mode must be refresh, parameters, or long");
    return out;
}
}
int main(int argc,char** argv) {
    try {
        if(argc<3)throw std::runtime_error("usage: followups mode output.csv [duration]");
        const std::string mode=argv[1];const auto configs=follow_configs(mode);
        const double duration=argc>3?std::stod(argv[3]):mode=="long"?100:10;
        if(!(duration>0&&duration<=1000))throw std::runtime_error("duration out of range");
        std::ofstream out(argv[2]);if(!out)throw std::runtime_error("cannot open output");out<<std::setprecision(17);
        out<<"case,method,schedule,policy,refresh,fraction,cap,substeps,h,duration,completed,requested,status,max_energy_percent,rms_position,rms_velocity,path_ratio,max_torque_error,max_radius,frozen_steps,evaluations,scalar_solves,block_solves,spectral_solves,singular,incompatible,direction_skips,repair_skips,repair_scaled,capped,max_extra_ratio\n";
        long count=0,stops=0;
        for(const auto& c:cases())for(int substeps:{4,8,16}) {
            const double h=1./(60*substeps);const int steps=std::lround(duration/h);
            const auto ref=reference(c,h,steps,16),fine=reference(c,h,steps,32);
            double referror=0;for(int i=0;i<=steps;++i)referror=std::max(referror,norm(ref[i].x-fine[i].x));
            if(referror>1e-6)throw std::runtime_error("reference convergence failed");
            for(const auto& cfg:configs) {
                const auto s=follow_run(c,cfg,substeps,ref);const auto& d=s.diagnostic;const auto& o=cfg.cfg.options;
                out<<c.name<<','<<cfg.cfg.name<<','<<se::name(o.schedule)<<','<<se::name(o.policy)<<','<<cfg.refresh<<','<<o.fraction<<','<<o.cap<<','<<substeps<<','<<h<<','<<duration<<','<<s.completed<<','<<steps<<','<<s.status<<','<<100*s.Emax/scale(c)<<','<<std::sqrt(s.x2/std::max(1L,s.completed))<<','<<std::sqrt(s.v2/std::max(1L,s.completed))<<','<<s.path/s.refpath<<','<<s.Lmax<<','<<s.maxradius<<','<<s.frozen<<','<<d.core.evaluations<<','<<d.core.scalar_solves<<','<<d.core.block_solves<<','<<d.spectral_solves<<','<<d.core.singular<<','<<d.core.incompatible<<','<<d.direction_skips<<','<<d.repair_skips<<','<<d.repair_scaled<<','<<d.capped<<','<<d.max_accepted_ratio<<'\n';
                ++count;stops+=s.completed<steps;
            }
        }
        std::cout<<mode<<": "<<count<<" runs, "<<stops<<" stops\n";return 0;
    } catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}
