#define main correction_study_main
#include "study.cpp"
#undef main
#include "symplectic_energy.hpp"
namespace sy=symplectic_energy;
namespace {
struct EnergyConfig {std::string name;sy::Options options;};
std::vector<EnergyConfig> energy_configs(bool longer) {
    std::vector<EnergyConfig> v;sy::Options o;
    o.mode=sy::Mode::XPBD;v.push_back({"xpbd",o});
    o.mode=sy::Mode::Predictor;v.push_back({"symplectic",o});
    o.mode=sy::Mode::Energy;v.push_back({"E_raw",o});
    o.fraction=.5;v.push_back({"E_half",o});
    o.fraction=.1;v.push_back({"E_tenth",o});
    o.fraction=1;o.cap=.1;v.push_back({"E_cap",o});
    o.fraction=.5;v.push_back({"E_half_cap",o});
    o.cap=.01;v.push_back({"E_half_cap_small",o});
    o.fraction=1;o.cap=.1;o.reject_large=true;v.push_back({"E_skip",o});
    o.reject_large=false;o.direction=true;v.push_back({"E_direction_cap",o});
    for(auto s:{se::Schedule::EL,se::Schedule::LE,se::Schedule::Joint})for(auto p:{se::Policy::Raw,se::Policy::RepairCap}) {
        o={};o.mode=sy::Mode::Paired;o.paired.schedule=s;o.paired.policy=p;
        v.push_back({std::string("paired_")+se::name(s)+"_"+se::name(p),o});
    }
    if(longer) {
        const std::set<std::string> selected{"xpbd","symplectic","E_half_cap","E_half_cap_small","E_skip","E_direction_cap","paired_EL_repair_cap","paired_LE_repair_cap"};
        v.erase(std::remove_if(v.begin(),v.end(),[&](const EnergyConfig& c){return !selected.count(c.name);}),v.end());
    }
    return v;
}
void sanity() {
    Model m;State s{{1.25,0},{}};const double h=1./480;
    const auto c=context(m,s,h,0,energy(m,s));double lambda=sy::constitutive_lambda(m,s.x,h);
    sy::Options o;o.mode=sy::Mode::Predictor;sy::Diagnostics d;State next;
    if(!sy::step(c,lambda,o,next,d))throw std::runtime_error("predictor failed");
    const V2 expected_v=s.v+acceleration(m,s.x)*h;
    if(norm(next.x-(s.x+expected_v*h))>1e-14||norm(next.v-expected_v)>1e-12)throw std::runtime_error("not the symplectic baseline");
    o.mode=sy::Mode::Energy;lambda=sy::constitutive_lambda(m,s.x,h);d={};
    if(!sy::step(c,lambda,o,next,d))throw std::runtime_error("raw energy step failed");
    const double a=h*h/m.alpha/m.mass,expected=.5*(1-a)*.25;
    if(std::abs(next.x.x-(1+expected))>1e-7)throw std::runtime_error("raw first-step analytic oracle failed");
    o.fraction=.5;o.cap=.1;lambda=sy::constitutive_lambda(m,s.x,h);d={};
    if(!sy::step(c,lambda,o,next,d)||d.max_ratio>.1+1e-10)throw std::runtime_error("energy cap failed");
    if(std::abs(lambda-sy::constitutive_lambda(m,next.x,h))>1e-14)throw std::runtime_error("force cache not refreshed");
    // An equilibrium has a zero gradient and must not produce a spurious kick.
    s={{1,0},{}};const auto eq=context(m,s,h,0,0);lambda=0;d={};
    if(!sy::step(eq,lambda,o,next,d)||norm(next.x-s.x)>1e-14||d.singular!=1)throw std::runtime_error("equilibrium guard failed");
    // Central-force predictor conserves angular momentum without correction.
    s={{1.1,0},{0,.8}};const double L=angular(s,1);o.mode=sy::Mode::Predictor;
    for(int j=0;j<4800;++j){const auto ctx=context(m,s,h,0,energy(m,s));if(!sy::step(ctx,lambda,o,next,d))throw std::runtime_error("orbit baseline failed");s=next;}
    if(std::abs(angular(s,1)-L)>1e-10)throw std::runtime_error("symplectic central-force angular invariant failed");
}
void write_trace(std::ostream& out,const Case& c,const EnergyConfig& cfg,double h,int step,State s,State r,sy::Diagnostics& d,double torque) {
    out<<c.name<<','<<cfg.name<<','<<step<<','<<step*h<<','<<s.x.x<<','<<s.x.y<<','<<s.v.x<<','<<s.v.y<<','<<energy(c.model,s)<<','<<r.x.x<<','<<r.x.y<<','<<energy(c.model,r)<<','<<angular(s,c.model.mass)<<','<<angular(r,c.model.mass)<<','<<torque<<','<<d.last_blend<<','<<d.last_extra<<','<<d.last_H<<','<<d.last_predicted_H<<'\n';
}
}
int main(int argc,char** argv) {
    try {
        sanity();if(argc>1&&std::string(argv[1])=="--test"){std::cout<<"PASS symplectic baseline, angular invariant, first energy step, cap, cache, equilibrium\n";return 0;}
        const std::filesystem::path dir=argc>1?argv[1]:"symplectic_energy_results";
        const double duration=argc>2?std::stod(argv[2]):10;const bool longer=duration>10;
        const auto configs=energy_configs(longer);std::filesystem::create_directories(dir);
        std::ofstream summary(dir/"summary.csv"),trace(dir/"traces.csv"),checks(dir/"reference_checks.csv");
        summary<<std::setprecision(17);trace<<std::setprecision(17);checks<<std::setprecision(17);
        summary<<"case,method,substeps,h,duration,completed,requested,status,Escale,max_energy_percent,final_energy_percent,rms_position,rms_velocity,reference_excursion,path_ratio,max_torque_error,final_L,initial_L,energy_solves,singular,caps,skips,max_correction_ratio,max_prediction_error,peak_time,peak_energy_percent\n";
        trace<<"case,method,step,time,x,y,vx,vy,energy,ref_x,ref_y,ref_energy,L,ref_L,torque,blend,extra,H,predicted_H\n";
        checks<<"case,h,duration,max_position_difference,max_energy_difference\n";
        auto scenes=cases();scenes.erase(std::remove_if(scenes.begin(),scenes.end(),[](const Case& c){return c.model.pendulum;}),scenes.end());
        Model free;scenes.push_back({"circular_orbit",free,{{1.1,0},{0,std::sqrt(11.)}}});
        int runs=0,stops=0;
        for(const auto& c:scenes)for(int substeps:{4,8,16}) {
            const double h=1./(60*substeps);const int count=std::lround(duration/h);
            const auto ref=reference(c,h,count,16),fine=reference(c,h,count,32);
            double refdiff=0,eref=0;for(int j=0;j<=count;++j){refdiff=std::max(refdiff,norm(ref[j].x-fine[j].x));eref=std::max(eref,std::abs(energy(c.model,ref[j])-energy(c.model,c.initial)));}
            if(refdiff>1e-6||eref>1e-6*scale(c))throw std::runtime_error("reference convergence failed");
            checks<<c.name<<','<<h<<','<<duration<<','<<refdiff<<','<<eref<<'\n';
            for(const auto& cfg:configs) {
                State s=c.initial;double stored=sy::constitutive_lambda(c.model,s.x,h);sy::Diagnostics d;
                const double E0=energy(c.model,s),L0=angular(s,c.model.mass),es=scale(c);
                double emax=0,x2=0,v2=0,exc=0,path=0,refpath=0,torque=0,tmax=0,peakt=0,peake=0;
                int done=0;std::string status="completed";
                const bool save=!longer&&substeps==8&&(cfg.name=="symplectic"||cfg.name=="E_half_cap"||cfg.name=="E_skip"||cfg.name=="E_direction_cap"||cfg.name=="paired_LE_repair_cap"||cfg.name=="E_raw");
                if(save)write_trace(trace,c,cfg,h,0,s,ref[0],d,torque);
                for(int j=1;j<=count;++j) {
                    const Context ctx=context(c.model,s,h,0,E0);State next;
                    const bool ok=sy::step(ctx,stored,cfg.options,next,d,j-1);const double E=energy(c.model,next);
                    if(!ok||!std::isfinite(E)){status="undefined";break;}
                    if(norm(next.x)>1e6||norm(next.v)>1e6||std::abs(E-E0)>1e12*es||std::abs(stored)>1e18){status="safety_stop";break;}
                    done=j;const double error=std::abs(E-E0);
                    if(error>emax){emax=error;peakt=j*h;peake=100*(E-E0)/es;}
                    x2+=dot(next.x-ref[j].x,next.x-ref[j].x);v2+=dot(next.v-ref[j].v,next.v-ref[j].v);exc+=dot(ref[j].x-c.initial.x,ref[j].x-c.initial.x);
                    path+=norm(next.x-s.x);refpath+=norm(ref[j].x-ref[j-1].x);
                    auto force_torque=[&](V2 x){return c.model.mass*(x.x*c.model.acceleration.y-x.y*c.model.acceleration.x);};
                    torque+=.5*h*(force_torque(s.x)+force_torque(next.x));tmax=std::max(tmax,std::abs(angular(next,c.model.mass)-L0-torque));
                    if(save)write_trace(trace,c,cfg,h,j,next,ref[j],d,torque);s=next;
                }
                if((cfg.options.mode==sy::Mode::Predictor||cfg.options.mode==sy::Mode::XPBD)&&done!=count)throw std::runtime_error("baseline run failed");
                summary<<c.name<<','<<cfg.name<<','<<substeps<<','<<h<<','<<duration<<','<<done<<','<<count<<','<<status<<','<<es<<','<<100*emax/es<<','<<100*(energy(c.model,s)-E0)/es<<','<<std::sqrt(x2/std::max(1,done))<<','<<std::sqrt(v2/std::max(1,done))<<','<<std::sqrt(exc/std::max(1,done))<<','<<(refpath>0?path/refpath:0)<<','<<tmax<<','<<angular(s,c.model.mass)<<','<<L0<<','<<d.energy_solves+d.paired.energy_attempts<<','<<d.singular+d.paired.core.singular<<','<<d.caps+d.paired.capped<<','<<d.skips+d.paired.direction_skips+d.paired.repair_skips<<','<<std::max(d.max_ratio,d.paired.max_accepted_ratio)<<','<<std::max(d.max_prediction_error,d.paired.max_prediction_error)<<','<<peakt<<','<<peake<<'\n';
                ++runs;stops+=done!=count;
            }
        }
        std::cout<<runs<<" runs of "<<duration<<" seconds; "<<stops<<" explicit stops\n";return 0;
    }catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}
