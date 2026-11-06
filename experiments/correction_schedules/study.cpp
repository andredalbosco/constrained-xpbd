#include "correction.hpp"
#include <sstream>
#include <set>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <vector>

using namespace energy_ordering;
namespace se=correction_schedules;
namespace {
struct Case { std::string name; Model model; State initial; };
std::vector<Case> cases() {
    Model m; m.acceleration={0,-9.81};
    const double equilibrium=m.length+m.mass*9.81*m.alpha;
    auto initial=[](double radius,double theta,double radial,double tangential) {
        const V2 n{std::sin(theta),-std::cos(theta)};
        return State{n*radius,n*radial+tangent(n)*tangential};
    };
    auto list=std::vector<Case>{{"vertical_oscillation",m,initial(equilibrium+.2,0,0,0)},
        {"small_swing",m,initial(equilibrium,.15,0,0)},
        {"swing_rest",m,initial(1.25,.65,0,0)},
        {"swing_tangent",m,initial(equilibrium,.65,0,1.2)},
        {"coupled_motion",m,initial(1.35,.65,-.7,.9)}};
    Model free=m;free.acceleration={};
    list.push_back({"spring_release",free,{{1.25,0},{}}});
    list.push_back({"spring_orbit",free,{{1.1,0},{0,.8}}});
    Model rod=m;rod.pendulum=true;rod.alpha=0;
    list.push_back({"pendulum_rest",rod,initial(1,.65,0,0)});
    list.push_back({"pendulum_kick",rod,initial(1,.4,0,1.2)});
    return list;
}
double scale(const Case& c) {
    // Minimum of spring plus gravity occurs directly below the anchor,
    // at radius L + m*g/k. This is a shifted potential, not a rigid pendulum.
    const double g=norm(c.model.acceleration),m=c.model.mass;
    const double minimum=-m*g*c.model.length-(c.model.pendulum?0:.5*m*m*g*g*c.model.alpha);
    return energy(c.model,c.initial)-minimum;
}
V2 acceleration(const Model& m,V2 x) {
    return m.acceleration+x*(-(norm(x)-m.length)/(m.alpha*m.mass*norm(x)));
}
State reference_step(const Model& m,State s,double h) {
    if(m.pendulum) {
        const double theta=std::atan2(s.x.x,-s.x.y),omega=dot(tangent(s.x),s.v)/(m.length*m.length);
        auto f=[&](double a){return -norm(m.acceleration)*std::sin(a)/m.length;};
        const double k1a=omega,k1w=f(theta),k2a=omega+h*k1w/2,k2w=f(theta+h*k1a/2);
        const double k3a=omega+h*k2w/2,k3w=f(theta+h*k2a/2),k4a=omega+h*k3w,k4w=f(theta+h*k3a);
        const double a=theta+h*(k1a+2*k2a+2*k3a+k4a)/6,w=omega+h*(k1w+2*k2w+2*k3w+k4w)/6;
        const V2 x{m.length*std::sin(a),-m.length*std::cos(a)};
        return {x,tangent(x)*w};
    }
    const V2 k1x=s.v, k1v=acceleration(m,s.x);
    const V2 k2x=s.v+k1v*(h/2), k2v=acceleration(m,s.x+k1x*(h/2));
    const V2 k3x=s.v+k2v*(h/2), k3v=acceleration(m,s.x+k2x*(h/2));
    const V2 k4x=s.v+k3v*h, k4v=acceleration(m,s.x+k3x*h);
    return {s.x+(k1x+k2x*2+k3x*2+k4x)*(h/6),
        s.v+(k1v+k2v*2+k3v*2+k4v)*(h/6)};
}
std::vector<State> reference(const Case& c,double h,int steps,int subdivisions) {
    std::vector<State> states{c.initial};
    State s=c.initial;
    for(int i=0;i<steps;++i) {
        for(int j=0;j<subdivisions;++j) s=reference_step(c.model,s,h/subdivisions);
        states.push_back(s);
    }
    return states;
}
struct Config {std::string name;se::Options options;};
std::vector<Config> configurations() {
    std::vector<Config> list;
    se::Options baseline;baseline.schedule=se::Schedule::Baseline;list.push_back({"baseline",baseline});
    const std::vector<se::Schedule> schedules{se::Schedule::Jacobi,se::Schedule::EL,se::Schedule::LE,se::Schedule::Joint,
        se::Schedule::GuardedJoint,se::Schedule::Alternating,se::Schedule::LEL,se::Schedule::ELE,se::Schedule::Spectral};
    const std::vector<se::Policy> policies{se::Policy::Raw,se::Policy::Direction,se::Policy::Cap,se::Policy::DirectionCap,se::Policy::Repair,se::Policy::RepairCap};
    for(auto schedule:schedules)for(auto policy:policies) {
        se::Options o;o.schedule=schedule;o.policy=policy;
        list.push_back({std::string(se::name(schedule))+"_"+se::name(policy),o});
        if(schedule==se::Schedule::Spectral) {
            o.power=1;list.push_back({"spectral_p1_"+std::string(se::name(policy)),o});
            o.power=.5;o.gain_cap=0;list.push_back({"spectral_unbounded_"+std::string(se::name(policy)),o});
        }
    }
    return list;
}
struct Stats {
    int completed=0,reversals=0,reference_reversals=0,frozen_steps=0,peak_step=0;
    std::string status="completed";
    double energy_error=0,rms_position=0,rms_velocity=0,path=0,ref_path=0;
    double max_rC=0,max_rx=0,final_energy=0,attempted_energy=0;
    double torque_impulse=0,ref_torque_impulse=0,max_torque_error=0,max_ref_torque_error=0;
    double min_length=1e300,max_length=0,reference_excursion=0;
    double peak_energy=0,peak_predicted_H=0,peak_H0=0,peak_blend=1,peak_extra=0;
    se::Diagnostics diagnostic;
};
double angular(State s,double mass){return mass*(s.x.x*s.v.y-s.x.y*s.v.x);}
double reversal_coordinate(const Case& c,State s) {
    return c.name=="vertical_oscillation"?s.v.y:dot(tangent(s.x),s.v)/norm(s.x);
}
void record_reversal(double value,int& last,int& count) {
    const int sign=value>1e-8?1:(value< -1e-8?-1:0);
    if(sign){if(last&&last!=sign)++count;last=sign;}
}
void prefix(std::ostream& out,const Case& c,const Config& cfg,int substeps,double duration) {
    const auto& o=cfg.options;
    out<<c.name<<','<<cfg.name<<','<<se::name(o.schedule)<<','<<se::name(o.policy)<<','<<o.fraction<<','<<o.cap<<','<<o.power<<','<<o.gain_cap<<','
        <<substeps<<','<<1./(60*substeps)<<",1,"<<duration<<','<<scale(c)<<','<<energy(c.model,c.initial)<<','<<c.model.mass<<','<<c.model.length<<','<<c.model.alpha<<','<<norm(c.model.acceleration)<<','<<c.model.pendulum;
}
void trace_row(std::ostream& out,const Case& c,const Config& cfg,int substeps,double duration,int step,State state,State ref,const Rows& r,const Stats& s) {
    prefix(out,c,cfg,substeps,duration);const auto& d=s.diagnostic;
    out<<','<<step<<','<<step/(60.*substeps)<<','<<state.x.x<<','<<state.x.y<<','<<state.v.x<<','<<state.v.y<<','
       <<energy(c.model,state)<<','<<ref.x.x<<','<<ref.x.y<<','<<energy(c.model,ref)<<','<<r.rC<<','
       <<d.energy_attempts<<','<<d.accepted<<','<<d.direction_skips<<','<<d.repair_skips<<','<<d.repair_scaled<<','<<d.capped<<','<<d.core.singular<<','
       <<angular(state,c.model.mass)<<','<<s.torque_impulse<<','<<angular(ref,c.model.mass)<<','<<s.ref_torque_impulse<<','
       <<d.last_H0<<','<<d.last_candidate_change<<','<<d.last_blend<<','<<d.last_predicted_H<<','<<d.last_extra<<','<<d.last_raw_energy_displacement<<','<<d.last_material_response_displacement<<','<<d.last_candidate_extra<<','<<d.last_actual_H0<<'\n';
}
Stats run(const Case& c,const Config& cfg,int substeps,double duration,const std::vector<State>& ref,std::ostream* traces) {
    Stats s;State state=c.initial;const double h=1./(60*substeps);
    const double E0=energy(c.model,state),L0=angular(state,c.model.mass),es=scale(c);
    s.min_length=s.max_length=norm(state.x);int sign=0,refsign=0;
    record_reversal(reversal_coordinate(c,state),sign,s.reversals);
    record_reversal(reversal_coordinate(c,state),refsign,s.reference_reversals);
    if(traces)trace_row(*traces,c,cfg,substeps,duration,0,state,ref[0],Rows{},s);
    for(int step=1;step<static_cast<int>(ref.size());++step) {
        // Fixed hard energy target; reset both multipliers each small substep.
        const Context ctx=context(c.model,state,h,0,E0);Iterate z{ctx.predictor};
        const bool ok=se::update(ctx,z,cfg.options,step-1,s.diagnostic);
        const State next{z.x,(z.x-state.x)/h};s.attempted_energy=energy(c.model,next);
        if(!ok||!finite(next.x)||!finite(next.v)||!std::isfinite(s.attempted_energy)) {s.status="aborted_undefined";break;}
        if(norm(next.x)>1e6*c.model.length||norm(next.v)>1e6*c.model.length||std::abs(s.attempted_energy-E0)>1e12*es
            ||std::abs(z.lambda)>1e18||std::abs(z.eta)>1e18) {s.status="aborted_safety_limit";break;}
        const Rows r=evaluate(ctx,z);if(!r.valid){s.status="aborted_undefined";break;}
        s.completed=step;
        const double error=std::abs(energy(c.model,next)-E0);
        if(error>s.energy_error){s.energy_error=error;s.peak_step=step;s.peak_energy=energy(c.model,next);
            s.peak_predicted_H=s.diagnostic.last_predicted_H;s.peak_H0=s.diagnostic.last_H0;s.peak_blend=s.diagnostic.last_blend;s.peak_extra=s.diagnostic.last_extra;}
        s.rms_position+=dot(next.x-ref[step].x,next.x-ref[step].x);
        s.rms_velocity+=dot(next.v-ref[step].v,next.v-ref[step].v);
        s.reference_excursion+=dot(ref[step].x-c.initial.x,ref[step].x-c.initial.x);
        s.path+=norm(next.x-state.x);s.ref_path+=norm(ref[step].x-ref[step-1].x);
        s.frozen_steps+=norm(next.x-state.x)<1e-10*c.model.length;
        s.max_rC=std::max(s.max_rC,std::abs(r.rC));s.max_rx=std::max(s.max_rx,norm(r.R)/(h*h));
        auto torque=[&](V2 x){return c.model.mass*(x.x*c.model.acceleration.y-x.y*c.model.acceleration.x);};
        s.torque_impulse+=.5*h*(torque(state.x)+torque(next.x));
        s.ref_torque_impulse+=.5*h*(torque(ref[step-1].x)+torque(ref[step].x));
        s.max_torque_error=std::max(s.max_torque_error,std::abs(angular(next,c.model.mass)-L0-s.torque_impulse));
        s.max_ref_torque_error=std::max(s.max_ref_torque_error,std::abs(angular(ref[step],c.model.mass)-L0-s.ref_torque_impulse));
        s.min_length=std::min(s.min_length,norm(next.x));s.max_length=std::max(s.max_length,norm(next.x));
        record_reversal(reversal_coordinate(c,next),sign,s.reversals);
        record_reversal(reversal_coordinate(c,ref[step]),refsign,s.reference_reversals);
        if(traces)trace_row(*traces,c,cfg,substeps,duration,step,next,ref[step],r,s);
        state=next;
    }
    if(s.completed){s.rms_position=std::sqrt(s.rms_position/s.completed);s.rms_velocity=std::sqrt(s.rms_velocity/s.completed);s.reference_excursion=std::sqrt(s.reference_excursion/s.completed);}
    s.final_energy=energy(c.model,state);
    if(s.status=="completed"&&s.diagnostic.core.incompatible)s.status="completed_with_incompatible_updates";
    return s;
}
} // namespace
int main(int argc,char** argv) {
    try {
        const std::filesystem::path dir=argc>1?argv[1]:"results/correction_schedules";
        // Optional comma-separated methods: run only these and save full 480Hz traces.
        std::set<std::string> selected;
        if(argc>2&&std::string(argv[2])!="all"){std::stringstream ss(argv[2]);std::string name;while(std::getline(ss,name,','))selected.insert(name);}
        const double duration=argc>3?std::stod(argv[3]):10;
        if(!(duration>0&&duration<=1000))throw std::runtime_error("duration must be in (0,1000]");
        auto configs=configurations();
        if(!selected.empty()) {
            for(const auto& name:selected)if(std::none_of(configs.begin(),configs.end(),[&](const Config& c){return c.name==name;}))throw std::runtime_error("unknown method: "+name);
            configs.erase(std::remove_if(configs.begin(),configs.end(),[&](const Config& c){return !selected.count(c.name);}),configs.end());
        }
        std::filesystem::create_directories(dir);
        std::ofstream summary(dir/"summary.csv"),traces(dir/"traces.csv"),refs(dir/"reference_checks.csv");
        if(!summary||!traces||!refs)throw std::runtime_error("cannot open outputs");
        summary<<std::setprecision(17);traces<<std::setprecision(17);refs<<std::setprecision(17);
        const char* pre="case,method,schedule,policy,fraction,cap,power,gain_cap,substeps,h,linearizations,duration,Escale,target_energy,mass,length,alpha,gravity,pendulum";
        summary<<pre<<",completed_steps,requested_steps,status,max_abs_energy_error,final_energy,rms_position_error,rms_velocity_error,path_length,reference_path_length,path_ratio,velocity_reversals,reference_reversals,frozen_steps,max_material_residual,max_rx,energy_attempts,accepted,direction_skips,repair_skips,repair_scaled,capped,singular_updates,incompatible_updates,guarded_fallbacks,rank_fallbacks,max_raw_extra,max_accepted_ratio,max_coupling,min_eigenvalue,max_linear_energy_error,max_linear_material_error,max_prediction_error,evaluations,scalar_solves,block_solves,spectral_solves,attempted_energy,max_torque_balance_error,reference_torque_balance_error,min_length,max_length,reference_rms_excursion,peak_step,peak_energy,peak_predicted_H,peak_H0,peak_blend,peak_extra,max_cancellation_ratio\n";
        traces<<pre<<",step,time,x,y,vx,vy,energy,ref_x,ref_y,ref_energy,rC,energy_attempts,accepted,direction_skips,repair_skips,repair_scaled,capped,singular_updates,angular_momentum,torque_impulse,ref_angular_momentum,ref_torque_impulse,H0,candidate_change,blend,predicted_H,extra_displacement,raw_energy_displacement,material_response_displacement,candidate_extra_displacement,actual_H0\n";
        refs<<"case,substeps,h,max_position_difference,max_energy_error\n";
        int runs=0,aborted=0;
        for(const auto& c:cases())for(int substeps:{4,8,16}) {
            const double h=1./(60*substeps);const int steps=static_cast<int>(std::lround(duration/h));
            const auto ref=reference(c,h,steps,16),fine=reference(c,h,steps,32);
            double pos_error=0,E_error=0;
            for(int i=0;i<=steps;++i){pos_error=std::max(pos_error,norm(ref[i].x-fine[i].x));E_error=std::max(E_error,std::abs(energy(c.model,ref[i])-energy(c.model,c.initial)));}
            if(pos_error>1e-7||E_error>1e-7*scale(c))throw std::runtime_error("reference check failed");
            refs<<c.name<<','<<substeps<<','<<h<<','<<pos_error<<','<<E_error<<'\n';
            for(const auto& cfg:configs) {
                const Stats s=run(c,cfg,substeps,duration,ref,!selected.empty()&&substeps==8?&traces:nullptr);const auto& d=s.diagnostic;
                prefix(summary,c,cfg,substeps,duration);
                summary<<','<<s.completed<<','<<steps<<','<<s.status<<','<<s.energy_error<<','<<s.final_energy<<','<<s.rms_position<<','<<s.rms_velocity<<','
                    <<s.path<<','<<s.ref_path<<','<<(s.ref_path>0?s.path/s.ref_path:0)<<','<<s.reversals<<','<<s.reference_reversals<<','<<s.frozen_steps<<','
                    <<s.max_rC<<','<<s.max_rx<<','<<d.energy_attempts<<','<<d.accepted<<','<<d.direction_skips<<','<<d.repair_skips<<','<<d.repair_scaled<<','<<d.capped<<','<<d.core.singular<<','<<d.core.incompatible<<','
                    <<d.guarded_fallbacks<<','<<d.rank_fallbacks<<','<<d.max_raw_extra<<','<<d.max_accepted_ratio<<','<<d.max_coupling<<','<<d.min_eigenvalue<<','
                    <<d.max_linear_energy_error<<','<<d.max_linear_material_error<<','<<d.max_prediction_error<<','<<d.core.evaluations<<','<<d.core.scalar_solves<<','<<d.core.block_solves<<','<<d.spectral_solves<<','<<s.attempted_energy<<','
                    <<s.max_torque_error<<','<<s.max_ref_torque_error<<','<<s.min_length<<','<<s.max_length<<','<<s.reference_excursion<<','<<s.peak_step<<','<<s.peak_energy<<','<<s.peak_predicted_H<<','<<s.peak_H0<<','<<s.peak_blend<<','<<s.peak_extra<<','<<d.max_cancellation_ratio<<'\n';
                ++runs;aborted+=s.completed<steps;
            }
        }
        std::cout<<runs<<" schedule/filter runs; "<<aborted<<" guarded stops; output "<<dir<<'\n';
        return 0;
    }catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}
