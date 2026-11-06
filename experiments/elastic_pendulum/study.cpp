#include "gradient_policies.hpp"
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <vector>

using namespace energy_ordering;
using namespace gradient_policies;
namespace {
struct Case { std::string name; Model model; State initial; };
std::vector<Case> cases() {
    Model m; m.acceleration={0,-9.81};
    const double equilibrium=m.length+m.mass*9.81*m.alpha;
    auto initial=[](double radius,double theta,double radial,double tangential) {
        const V2 n{std::sin(theta),-std::cos(theta)};
        return State{n*radius,n*radial+tangent(n)*tangential};
    };
    return {{"vertical_oscillation",m,initial(equilibrium+.2,0,0,0)},
        {"small_swing",m,initial(equilibrium,.15,0,0)},
        {"swing_rest",m,initial(1.25,.65,0,0)},
        {"swing_tangent",m,initial(equilibrium,.65,0,1.2)},
        {"coupled_motion",m,initial(1.35,.65,-.7,.9)}};
}
double scale(const Case& c) {
    // Minimum of spring plus gravity occurs directly below the anchor,
    // at radius L + m*g/k. This is a shifted potential, not a rigid pendulum.
    const double g=norm(c.model.acceleration),m=c.model.mass;
    const double minimum=-m*g*c.model.length-.5*m*m*g*g*c.model.alpha;
    return energy(c.model,c.initial)-minimum;
}
V2 acceleration(const Model& m,V2 x) {
    return m.acceleration+x*(-(norm(x)-m.length)/(m.alpha*m.mass*norm(x)));
}
State reference_step(const Model& m,State s,double h) {
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
struct Stats {
    int completed=0,reversals=0,reference_reversals=0,frozen_steps=0;
    std::string status="completed";
    double energy_error=0,rms_position=0,rms_velocity=0,path=0,ref_path=0;
    double max_rC=0,max_rx=0,max_angular_change=0,final_energy=0,attempted_energy=0;
    double torque_impulse=0,ref_torque_impulse=0,max_torque_error=0,max_ref_torque_error=0;
    double min_length=1e300,max_length=0;
    Diagnostics diagnostic;
};
double angular(State s,double mass){return mass*(s.x.x*s.v.y-s.x.y*s.v.x);}
double reversal_coordinate(const Case& c,State s) {
    return c.name=="vertical_oscillation" ? s.v.y : dot(tangent(s.x),s.v)/norm(s.x);
}
void record_reversal(double value,int& last,int& count) {
    const int sign=value>1e-8?1:(value< -1e-8?-1:0);
    if(sign){if(last&&last!=sign)++count;last=sign;}
}
void trace_row(std::ostream& out,const Case& c,const gradient_policies::Options& o,bool baseline,double h,int sweeps,
    int step,State state,State ref,const Rows& r,const Stats& s) {
    out<<c.name<<','<<(baseline?"baseline":name(o.order))<<','<<(baseline?"none":name(o.gradient))
       <<','<<(baseline?"baseline":(o.objective_gate?"directional":(o.gamma<0?"unguarded":"size_reject")))<<','<<o.gamma<<','<<o.retain_rx<<','<<h<<','<<sweeps<<','<<step<<','<<step*h<<','
       <<scale(c)<<','<<energy(c.model,c.initial)<<','<<state.x.x<<','<<state.x.y<<','<<state.v.x<<','<<state.v.y<<','
       <<energy(c.model,state)<<','<<ref.x.x<<','<<ref.x.y<<','<<energy(c.model,ref)<<','<<r.rC<<','<<norm(r.R)/(h*h)<<','
       <<s.diagnostic.attempts<<','<<s.diagnostic.accepted<<','<<s.diagnostic.oversized<<','
       <<s.diagnostic.core.singular<<','<<s.diagnostic.core.incompatible<<','
       <<angular(state,c.model.mass)<<','<<s.torque_impulse<<','<<angular(ref,c.model.mass)<<','<<s.ref_torque_impulse<<','<<s.diagnostic.directional_skips<<'\n';
}
Stats run(const Case& c,const gradient_policies::Options& o,bool baseline,double h,int sweeps,const std::vector<State>& ref,std::ostream* traces) {
    Stats s;State state=c.initial;History history;
    const double E0=energy(c.model,state),L0=angular(state,c.model.mass),es=scale(c);
    s.min_length=s.max_length=norm(state.x);
    int sign=0,refsign=0;
    record_reversal(reversal_coordinate(c,state),sign,s.reversals);
    record_reversal(reversal_coordinate(c,state),refsign,s.reference_reversals);
    if(traces)trace_row(*traces,c,o,baseline,h,sweeps,0,state,ref[0],Rows{},s);
    for(int step=1;step<static_cast<int>(ref.size());++step) {
        const Context ctx=context(c.model,state,h,0,E0);
        Iterate z{ctx.predictor};history.begin(ctx,z);bool ok=true;
        for(int outer=0;outer<sweeps;++outer) {
            ok=baseline?material_update(ctx,z,s.diagnostic.core):update(ctx,z,o,history,s.diagnostic);
            const State trial{z.x,(z.x-state.x)/h};
            s.attempted_energy=energy(c.model,trial);
            if(!ok||!finite(trial.x)||!finite(trial.v)||!std::isfinite(s.attempted_energy)) {
                s.status="aborted_undefined";ok=false;break;
            }
            if(norm(trial.x)>1e6*c.model.length||norm(trial.v)>1e6*c.model.length
                ||std::abs(s.attempted_energy-E0)>1e12*es
                ||std::abs(z.lambda)>1e18||std::abs(z.eta)>1e18) {
                s.status="aborted_safety_limit";ok=false;break;
            }
        }
        if(!ok)break;
        const Rows r=evaluate(ctx,z);
        if(!r.valid){s.status="aborted_undefined";break;}
        const State next{z.x,(z.x-state.x)/h};history.end(ctx,z);
        s.completed=step;s.energy_error=std::max(s.energy_error,std::abs(energy(c.model,next)-E0));
        s.rms_position+=dot(next.x-ref[step].x,next.x-ref[step].x);
        s.rms_velocity+=dot(next.v-ref[step].v,next.v-ref[step].v);
        s.path+=norm(next.x-state.x);s.ref_path+=norm(ref[step].x-ref[step-1].x);
        s.frozen_steps+=norm(next.x-state.x)<1e-10*c.model.length;
        s.max_rC=std::max(s.max_rC,std::abs(r.rC));s.max_rx=std::max(s.max_rx,norm(r.R)/(h*h));
        s.max_angular_change=std::max(s.max_angular_change,std::abs(angular(next,c.model.mass)-L0));
        auto torque=[&](V2 x){return c.model.mass*(x.x*c.model.acceleration.y-x.y*c.model.acceleration.x);};
        s.torque_impulse+=.5*h*(torque(state.x)+torque(next.x));
        s.ref_torque_impulse+=.5*h*(torque(ref[step-1].x)+torque(ref[step].x));
        s.max_torque_error=std::max(s.max_torque_error,std::abs(angular(next,c.model.mass)-L0-s.torque_impulse));
        s.max_ref_torque_error=std::max(s.max_ref_torque_error,std::abs(angular(ref[step],c.model.mass)-L0-s.ref_torque_impulse));
        s.min_length=std::min(s.min_length,norm(next.x));s.max_length=std::max(s.max_length,norm(next.x));
        record_reversal(reversal_coordinate(c,next),sign,s.reversals);
        record_reversal(reversal_coordinate(c,ref[step]),refsign,s.reference_reversals);
        if(traces)trace_row(*traces,c,o,baseline,h,sweeps,step,next,ref[step],r,s);
        state=next;
    }
    if(s.completed){s.rms_position=std::sqrt(s.rms_position/s.completed);s.rms_velocity=std::sqrt(s.rms_velocity/s.completed);}
    s.final_energy=energy(c.model,state);
    if(s.status=="completed"&&s.diagnostic.core.incompatible)s.status="completed_with_incompatible_updates";
    return s;
}
} // namespace
int main(int argc,char** argv) {
    try {
        const std::filesystem::path dir=argc>1?argv[1]:"results/elastic_pendulum";
        const std::string requested=argc>2?argv[2]:"reject";
        if(requested!="reject"&&requested!="directional"&&requested!="unguarded")
            throw std::runtime_error("energy policy must be reject, directional, or unguarded");
        std::filesystem::create_directories(dir);
        std::ofstream summary(dir/"summary.csv"),traces(dir/"traces.csv"),refs(dir/"reference_checks.csv");
        if(!summary||!traces||!refs)throw std::runtime_error("cannot open outputs");
        summary<<std::setprecision(17);traces<<std::setprecision(17);refs<<std::setprecision(17);
        summary<<"case,order,gradient,energy_policy,gamma,retain_rx,h,sweeps,Escale,target_energy,completed_steps,requested_steps,status,"
            "max_abs_energy_error,final_energy,rms_position_error,rms_velocity_error,path_length,reference_path_length,path_ratio,"
            "velocity_reversals,reference_reversals,frozen_steps,max_material_residual,max_rx,max_angular_change,"
            "energy_attempts,energy_accepted,oversized_skips,singular_updates,incompatible_updates,"
            "max_candidate_ratio,max_accepted_ratio,max_candidate_displacement,evaluations,scalar_solves,attempted_energy,max_torque_balance_error,reference_torque_balance_error,min_length,max_length,directional_skips,max_accepted_objective_work\n";
        traces<<"case,order,gradient,energy_policy,gamma,retain_rx,h,sweeps,step,time,Escale,target_energy,x,y,vx,vy,energy,"
            "ref_x,ref_y,ref_energy,rC,rx,energy_attempts,energy_accepted,oversized_skips,singular_updates,incompatible_updates,angular_momentum,torque_impulse,ref_angular_momentum,ref_torque_impulse,directional_skips\n";
        refs<<"case,h,max_position_difference,max_energy_error\n";
        int runs=0,aborted=0;
        for(const auto& c:cases())for(double h:{1./60,1./120,1./240}) {
            const int steps=static_cast<int>(std::lround(10/h));
            const auto ref=reference(c,h,steps,32),fine=reference(c,h,steps,64);
            double pos_error=0,E_error=0;
            for(int i=0;i<=steps;++i){pos_error=std::max(pos_error,norm(ref[i].x-fine[i].x));
                E_error=std::max(E_error,std::abs(energy(c.model,ref[i])-energy(c.model,c.initial)));}
            if(pos_error>1e-7||E_error>1e-7*scale(c))throw std::runtime_error("reference check failed");
            refs<<c.name<<','<<h<<','<<pos_error<<','<<E_error<<'\n';
            for(int sweeps:{1,5,20}) {
                auto execute=[&](gradient_policies::Options o,bool baseline=false) {
                    o.objective_gate=requested=="directional";
                    const Stats s=run(c,o,baseline,h,sweeps,ref,(h==1./60||h==1./240)&&sweeps==5?&traces:nullptr);
                    const auto& d=s.diagnostic;
                    summary<<c.name<<','<<(baseline?"baseline":name(o.order))<<','<<(baseline?"none":name(o.gradient))
                        <<','<<(baseline?"baseline":(o.objective_gate?"directional":(o.gamma<0?"unguarded":"size_reject")))<<','<<o.gamma<<','<<o.retain_rx<<','<<h<<','<<sweeps<<','<<scale(c)<<','<<energy(c.model,c.initial)<<','
                        <<s.completed<<','<<steps<<','<<s.status<<','<<s.energy_error<<','<<s.final_energy<<','<<s.rms_position<<','<<s.rms_velocity<<','
                        <<s.path<<','<<s.ref_path<<','<<(s.ref_path>0?s.path/s.ref_path:0)<<','<<s.reversals<<','<<s.reference_reversals<<','
                        <<s.frozen_steps<<','<<s.max_rC<<','<<s.max_rx<<','<<s.max_angular_change<<','<<d.attempts<<','<<d.accepted<<','<<d.oversized<<','
                        <<d.core.singular<<','<<d.core.incompatible<<','<<d.max_ratio<<','<<d.max_accepted_ratio<<','<<d.max_candidate<<','
                        <<d.core.evaluations<<','<<d.core.scalar_solves<<','<<s.attempted_energy<<','<<s.max_torque_error<<','<<s.max_ref_torque_error<<','<<s.min_length<<','<<s.max_length<<','<<d.directional_skips<<','<<d.max_accepted_objective_work<<'\n';
                    ++runs;aborted+=s.completed<steps;
                };
                execute({Order::EL,Gradient::Fresh,false,-1},true);
                const std::vector<double> gammas=requested=="reject"
                    ?std::vector<double>{-1.,.25,1.,4.,16.}:std::vector<double>{-1.};
                for(double gamma:gammas)execute({Order::EL,Gradient::Fresh,false,gamma});
            }
        }
        std::cout<<runs<<" runs; "<<aborted<<" guarded stops; output "<<dir<<'\n';
        return 0;
    } catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}
