#include "frozen_schedules.hpp"
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <vector>

using namespace energy_ordering;
namespace {
struct Case { std::string name; Model model; State initial; };
std::vector<Case> cases() {
    const double theta=.65;
    Model spring, pendulum;
    pendulum.pendulum=true; pendulum.alpha=0; pendulum.acceleration={0,-9.81};
    const V2 x{std::sin(theta),-std::cos(theta)}, t{std::cos(theta),std::sin(theta)};
    return {{"spring_rest",spring,{{1.25,0},{0,0}}},
        {"spring_radial",spring,{{1.25,0},{-.7,0}}},
        {"spring_tangent",spring,{{1.25,0},{0,1}}},
        {"pendulum_rest",pendulum,{x,{0,0}}},
        {"pendulum_tangent",pendulum,{x,t*1.2}}};
}
double scale(const Case& c) {
    // Energy above the minimum potential, independent of the gravity zero.
    return energy(c.model,c.initial)+c.model.mass*norm(c.model.acceleration)*c.model.length;
}
V2 acceleration(const Model& m,V2 x) {
    return x*(-(norm(x)-m.length)/(m.alpha*m.mass*norm(x)));
}
State reference_step(const Model& m,State s,double h) {
    if (m.pendulum) {
        const double theta=std::atan2(s.x.x,-s.x.y);
        const double omega=dot(s.v,tangent(s.x))/(m.length*m.length);
        auto f=[&](double th) {return -norm(m.acceleration)/m.length*std::sin(th);};
        const double a1=f(theta),w2=omega+h*a1/2,a2=f(theta+h*omega/2);
        const double w3=omega+h*a2/2,a3=f(theta+h*w2/2);
        const double w4=omega+h*a3,a4=f(theta+h*w3);
        const double th=theta+h*(omega+2*w2+2*w3+w4)/6;
        const double om=omega+h*(a1+2*a2+2*a3+a4)/6;
        V2 x{m.length*std::sin(th),-m.length*std::cos(th)};
        return {x,tangent(x)*om};
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
struct Stats {
    int completed=0,first_incompatible=0,infeasible_steps=0;
    std::string status="completed";
    double max_energy_error=0,max_rC=0,max_rx=0,rms_position=0,path_length=0;
    double application_abs=0,application_rel=0,attempted_energy=0,attempted_radius=0;
    Counters counts;
};
void trace_header(std::ostream& out) {
    out<<"case,method,application,retain_rx,h,sweeps,softness,beta,Escale,target_energy,step,time,x,y,vx,vy,energy,"
        "ref_x,ref_y,ref_energy,rC,rE,rx,singular_updates,incompatible_updates\n";
}
void trace_row(std::ostream& out,const Case& c,Schedule method,Application application,bool retain,double h,int sweeps,double soft,
    int step,State s,State ref,const Rows& r,const Stats& stats) {
    out<<c.name<<','<<name(method)<<','<<(application==Application::Deferred?"deferred":"immediate")<<','<<retain<<','<<h<<','<<sweeps<<','<<soft<<','<<soft*scale(c)<<','
        <<scale(c)<<','<<energy(c.model,c.initial)<<','<<step<<','<<step*h<<','<<s.x.x<<','<<s.x.y<<','
        <<s.v.x<<','<<s.v.y<<','<<energy(c.model,s)<<','<<ref.x.x<<','<<ref.x.y<<','<<energy(c.model,ref)<<','
        <<r.rC<<','<<r.rE<<','<<norm(r.R)/(h*h)<<','<<stats.counts.singular<<','<<stats.counts.incompatible<<'\n';
}
Stats run(const Case& c,Schedule method,Application application,bool retain,double h,int sweeps,double soft,
    const std::vector<State>& ref,std::ostream* trace) {
    Stats stats; State state=c.initial;
    const double E0=energy(c.model,state),es=scale(c);
    if(trace) trace_row(*trace,c,method,application,retain,h,sweeps,soft,0,state,ref[0],Rows{},stats);
    const int steps=static_cast<int>(ref.size())-1;
    for(int step=1;step<=steps;++step) {
        const Context ctx=context(c.model,state,h,soft*es,E0);
        if(c.model.pendulum&&soft==0&&method!=Schedule::Baseline&&!pendulum_feasible(ctx))
            ++stats.infeasible_steps;
        Iterate z{ctx.predictor}; bool ok=true;
        for(int outer=0;outer<sweeps;++outer) {
            const Iterate before=z;
            ok=schedule_update(ctx,z,method,retain,stats.counts,1,application);
            // Compare immediate/deferred application on the SAME frozen input.
            // Their complete unstable trajectories need not agree in floating point.
            if(method!=Schedule::Joint&&method!=Schedule::Baseline) {
                Iterate immediate=before; Counters unused;
                const Application alternative=application==Application::Deferred?Application::Immediate:Application::Deferred;
                const bool immediate_ok=schedule_update(ctx,immediate,method,retain,unused,1,alternative);
                if(ok&&immediate_ok) {
                    const double error=norm(immediate.x-z.x);
                    stats.application_abs=std::max(stats.application_abs,error);
                    stats.application_rel=std::max(stats.application_rel,error/(1+norm(z.x)+norm(before.x)));
                    if(immediate.lambda!=z.lambda||immediate.eta!=z.eta)
                        throw std::runtime_error("application mode changed frozen multiplier solve");
                }
            }
            if(stats.counts.incompatible&&!stats.first_incompatible) stats.first_incompatible=step;
            const State trial{z.x,(z.x-state.x)/h};
            stats.attempted_energy=energy(c.model,trial); stats.attempted_radius=norm(z.x);
            if(!ok||!finite(trial.x)||!std::isfinite(stats.attempted_energy)) {
                stats.status="aborted_undefined";ok=false;break;
            }
            if(norm(trial.x)>1e6*c.model.length||norm(trial.v)>1e6*c.model.length
                ||std::abs(stats.attempted_energy-E0)>1e12*es
                ||std::abs(z.lambda)>1e18||std::abs(z.eta)>1e18) {
                stats.status="aborted_safety_limit";ok=false;break;
            }
        }
        if(!ok) break;
        const Rows r=evaluate(ctx,z);
        if(!r.valid) {stats.status="aborted_undefined";break;}
        const State next{z.x,(z.x-state.x)/h};
        stats.completed=step;
        stats.max_energy_error=std::max(stats.max_energy_error,std::abs(energy(c.model,next)-E0));
        stats.max_rC=std::max(stats.max_rC,std::abs(r.rC));
        stats.max_rx=std::max(stats.max_rx,norm(r.R)/(h*h));
        stats.rms_position+=dot(next.x-ref[step].x,next.x-ref[step].x);
        stats.path_length+=norm(next.x-state.x);
        if(trace) trace_row(*trace,c,method,application,retain,h,sweeps,soft,step,next,ref[step],r,stats);
        state=next;
    }
    if(stats.completed) stats.rms_position=std::sqrt(stats.rms_position/stats.completed);
    if(stats.status=="completed"&&stats.counts.incompatible) stats.status="completed_with_incompatible_updates";
    return stats;
}
} // namespace
int main(int argc,char** argv) {
    try {
        const std::filesystem::path dir=argc>1?argv[1]:"results/frozen_schedules_jacobi";
        std::filesystem::create_directories(dir);
        std::ofstream traces(dir/"traces.csv"),summary(dir/"summary.csv"),refs(dir/"reference_checks.csv");
        if(!traces||!summary||!refs) throw std::runtime_error("cannot open result files");
        traces<<std::setprecision(17);summary<<std::setprecision(17);refs<<std::setprecision(17);
        trace_header(traces);
        summary<<"case,method,application,retain_rx,h,sweeps,softness,Escale,target_energy,completed_steps,requested_steps,status,"
            "max_abs_energy_error,max_abs_rC,max_rx,rms_position_error,path_length,evaluations,scalar_solves,block_solves,"
            "singular_updates,incompatible_updates,rank_deficient_blocks,first_incompatible_step,infeasible_steps,"
            "application_max_absolute_difference,application_max_relative_difference,attempted_energy,attempted_radius\n";
        refs<<"case,h,max_position_difference,max_energy_error\n";
        int runs=0,aborted=0;
        for(const auto& c:cases()) for(double h:{1./240,1./60}) {
            const int steps=static_cast<int>(std::lround(10/h));
            const auto ref=reference(c,h,steps,32),fine=reference(c,h,steps,64);
            double position_error=0,energy_error=0;
            for(int i=0;i<=steps;++i) {
                position_error=std::max(position_error,norm(ref[i].x-fine[i].x));
                energy_error=std::max(energy_error,std::abs(energy(c.model,ref[i])-energy(c.model,c.initial)));
            }
            if(position_error>1e-7||energy_error>1e-7*scale(c)) throw std::runtime_error("reference check failed");
            refs<<c.name<<','<<h<<','<<position_error<<','<<energy_error<<'\n';
            for(int sweeps:{1,5,20}) for(Schedule method:{Schedule::Baseline,Schedule::Joint,Schedule::Jacobi,Schedule::LE,
                Schedule::EL,Schedule::LEL,Schedule::ELE}) for(double soft:{0.,.01}) for(bool retain:{false,true})
                for(Application application:{Application::Deferred,Application::Immediate}) {
                if(method==Schedule::Baseline&&(soft!=0||retain)) continue;
                // Jacobi computes both corrections before applying their sum.
                // A second application mode would only change summation order.
                if(application==Application::Immediate&&(method==Schedule::Baseline||method==Schedule::Joint
                    ||method==Schedule::Jacobi)) continue;
                const bool record=h==1./60&&sweeps==5;
                const Stats s=run(c,method,application,retain,h,sweeps,soft,ref,record?&traces:nullptr);
                summary<<c.name<<','<<name(method)<<','<<(application==Application::Deferred?"deferred":"immediate")<<','<<retain<<','<<h<<','<<sweeps<<','<<soft<<','<<scale(c)<<','
                    <<energy(c.model,c.initial)<<','<<s.completed<<','<<steps<<','<<s.status<<','<<s.max_energy_error<<','
                    <<s.max_rC<<','<<s.max_rx<<','<<s.rms_position<<','<<s.path_length<<','<<s.counts.evaluations<<','
                    <<s.counts.scalar_solves<<','<<s.counts.block_solves<<','<<s.counts.singular<<','<<s.counts.incompatible<<','
                    <<s.counts.rank_deficient_blocks<<','<<s.first_incompatible<<','<<s.infeasible_steps<<','
                    <<s.application_abs<<','<<s.application_rel<<','<<s.attempted_energy<<','<<s.attempted_radius<<'\n';
                ++runs;aborted+=s.completed<steps;
            }
        }
        std::cout<<runs<<" runs; "<<aborted<<" guarded stops. Results: "<<dir<<'\n';
        return 0;
    } catch(const std::exception& e) {std::cerr<<e.what()<<'\n';return 1;}
}
