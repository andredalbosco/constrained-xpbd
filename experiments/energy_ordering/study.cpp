#include "solver.hpp"

#include <array>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
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
struct Metrics {
    std::string status="completed";
    int completed=0, first_singular=0, first_incompatible=0, first_infeasible=0, infeasible=0;
    double max_dE=0,max_rC=0,max_rE=0,max_rx=0,max_speed=0,max_displacement=0;
    double path_length=0,rms_position=0,late_kinetic=0,max_dual=0;
    double final_dE=0,final_rC=0,final_rE=0,final_rx=0,final_kinetic=0;
    double attempted_norm=0, attempted_energy=0;
    long late_samples=0,rows_closed=0,kkt_closed=0;
    Counters count;
};
void trace_header(std::ostream& out) {
    out<<"case,method,retain_rx,inner,h,sweeps,softness,beta,Escale,step,time,x,y,vx,vy,energy_error,rC,rE,rx,kinetic,lambda,eta,singular_updates,incompatible_updates,ref_x,ref_y,ref_kinetic\n";
}
void trace_row(std::ostream& out,const Case& c,const Options& o,double h,double softness,double beta,
    int step,State state,const Rows& r,const Iterate& z,const Metrics& stats,State ref) {
    out<<c.name<<','<<name(o.method)<<','<<o.retain_rx<<','<<o.inner<<','<<h<<','<<o.sweeps<<','
        <<softness<<','<<beta<<','<<scale(c)<<','<<step<<','<<step*h<<','<<state.x.x<<','<<state.x.y<<','
        <<state.v.x<<','<<state.v.y<<','<<r.H<<','<<r.rC<<','<<r.rE<<','<<norm(r.R)/(h*h)<<','
        <<.5*c.model.mass*dot(state.v,state.v)<<','<<z.lambda<<','<<z.eta<<','<<stats.count.singular<<','
        <<stats.count.incompatible<<','<<ref.x.x<<','<<ref.x.y<<','<<.5*c.model.mass*dot(ref.v,ref.v)<<'\n';
}
Metrics run(const Case& c,const Options& o,double h,double softness,int steps,
            const std::vector<State>& ref,std::ostream* trace) {
    Metrics stats;
    const double es=scale(c), E0=energy(c.model,c.initial), beta=softness*es;
    State s=c.initial;
    if(trace) {
        Rows r; Iterate z{s.x,0,0};
        trace_row(*trace,c,o,h,softness,beta,0,s,r,z,stats,ref[0]);
    }
    for(int step=1;step<=steps;++step) {
        const Context ctx=context(c.model,s,h,beta,E0);
        if(c.model.pendulum && beta==0 && o.method!=Method::Baseline && !pendulum_feasible(ctx)) {
            ++stats.infeasible;
            if(!stats.first_infeasible) stats.first_infeasible=step;
        }
        Iterate z{ctx.predictor,0,0};
        bool ok=true;
        for(int k=0;k<o.sweeps;++k) {
            ok=sweep(ctx,z,o,stats.count);
            if(stats.count.singular && !stats.first_singular) stats.first_singular=step;
            if(stats.count.incompatible && !stats.first_incompatible) stats.first_incompatible=step;
            // Termination guard, never clipping, damping, or accepting a modified step.
            const State trial{z.x,(z.x-s.x)/h};
            stats.attempted_norm=norm(z.x);
            stats.attempted_energy=energy(c.model,trial)-E0;
            if(!ok || !finite(z.x) || !std::isfinite(stats.attempted_energy)) {
                stats.status="aborted_nonfinite_or_undefined"; ok=false; break;
            }
            if(norm(z.x)>1e6*c.model.length || norm(trial.v)>1e6*c.model.length
                || std::abs(stats.attempted_energy)>1e12*es
                || std::abs(z.lambda)>1e18 || std::abs(z.eta)>1e18) {
                stats.status="aborted_safety_limit"; ok=false; break;
            }
        }
        if(!ok) break;
        const Rows r=evaluate(ctx,z);
        if(!r.valid) { stats.status="aborted_nonfinite_or_undefined"; break; }
        const State next{z.x,(z.x-s.x)/h};
        stats.completed=step;
        stats.final_dE=r.H; stats.final_rC=r.rC; stats.final_rE=r.rE;
        stats.final_rx=norm(r.R)/(h*h);
        stats.final_kinetic=.5*c.model.mass*dot(next.v,next.v);
        stats.max_dE=std::max(stats.max_dE,std::abs(r.H));
        stats.max_rC=std::max(stats.max_rC,std::abs(r.rC));
        stats.max_rE=std::max(stats.max_rE,std::abs(r.rE));
        stats.max_rx=std::max(stats.max_rx,stats.final_rx);
        stats.max_speed=std::max(stats.max_speed,norm(next.v));
        stats.max_displacement=std::max(stats.max_displacement,norm(next.x-c.initial.x));
        stats.max_dual=std::max({stats.max_dual,std::abs(z.lambda),std::abs(z.eta)});
        stats.path_length+=norm(next.x-s.x);
        stats.rms_position+=dot(next.x-ref[step].x,next.x-ref[step].x);
        if(step>steps*3/4) { stats.late_kinetic+=stats.final_kinetic; ++stats.late_samples; }
        const bool closed=std::abs(r.rC)<1e-8*c.model.length
            && (o.method==Method::Baseline||std::abs(r.rE)<1e-8*es);
        stats.rows_closed+=closed;
        stats.kkt_closed+=closed && stats.final_rx<1e-8*es/c.model.length;
        if(trace) trace_row(*trace,c,o,h,softness,beta,step,next,r,z,stats,ref[step]);
        s=next;
    }
    if(stats.completed) stats.rms_position=std::sqrt(stats.rms_position/stats.completed);
    if(stats.late_samples) stats.late_kinetic/=stats.late_samples;
    if(stats.status=="completed"&&stats.count.incompatible) stats.status="completed_with_incompatible_updates";
    else if(stats.status=="completed"&&stats.count.singular) stats.status="completed_with_redundant_updates";
    return stats;
}
void summary_header(std::ostream& out) {
    out<<"case,method,retain_rx,inner,h,sweeps,softness,beta,Escale,requested_steps,completed_steps,status,"
        "max_abs_energy_error,final_energy_error,max_abs_rC,max_abs_rE,max_rx,final_rC,final_rE,final_rx,"
        "max_speed,max_displacement,path_length,rms_position_error,late_mean_kinetic,final_kinetic,max_abs_dual,"
        "rows_closed_steps,kkt_closed_steps,evaluations,scalar_solves,block_solves,singular_updates,"
        "incompatible_updates,first_singular_step,first_incompatible_step,infeasible_steps,first_infeasible_step,"
        "attempted_norm,attempted_energy_error,rank_deficient_blocks,late_samples\n";
}
void summary_row(std::ostream& out,const Case& c,const Options& o,double h,double soft,int steps,const Metrics& s) {
    out<<c.name<<','<<name(o.method)<<','<<o.retain_rx<<','<<o.inner<<','<<h<<','<<o.sweeps<<','
        <<soft<<','<<soft*scale(c)<<','<<scale(c)<<','<<steps<<','<<s.completed<<','<<s.status<<','
        <<s.max_dE<<','<<s.final_dE<<','<<s.max_rC<<','<<s.max_rE<<','<<s.max_rx<<','
        <<s.final_rC<<','<<s.final_rE<<','<<s.final_rx<<','<<s.max_speed<<','<<s.max_displacement<<','
        <<s.path_length<<','<<s.rms_position<<','<<s.late_kinetic<<','<<s.final_kinetic<<','<<s.max_dual<<','
        <<s.rows_closed<<','<<s.kkt_closed<<','<<s.count.evaluations<<','<<s.count.scalar_solves<<','
        <<s.count.block_solves<<','<<s.count.singular<<','<<s.count.incompatible<<','<<s.first_singular<<','
        <<s.first_incompatible<<','<<s.infeasible<<','<<s.first_infeasible<<','<<s.attempted_norm<<','
        <<s.attempted_energy<<','<<s.count.rank_deficient_blocks<<','<<s.late_samples<<'\n';
}
} // namespace
int main(int argc,char** argv) {
    try {
        const std::filesystem::path dir=argc>1?argv[1]:"results/energy_ordering";
        std::filesystem::create_directories(dir);
        std::ofstream summary(dir/"summary.csv"),traces(dir/"traces.csv"),refs(dir/"reference_checks.csv");
        summary<<std::setprecision(17); traces<<std::setprecision(17); refs<<std::setprecision(17);
        if(!summary||!traces||!refs) throw std::runtime_error("cannot open outputs");
        summary_header(summary); trace_header(traces);
        refs<<"case,h,max_refinement_position_error,max_reference_energy_error\n";
        int runs=0,aborted=0;
        for(const Case& c:cases()) for(double h:{1./240,1./60}) {
            const int steps=static_cast<int>(std::lround(10/h));
            const auto ref=reference(c,h,steps,32), finer=reference(c,h,steps,64);
            double refinement=0,energy_drift=0;
            for(int i=0;i<=steps;++i) {
                refinement=std::max(refinement,norm(ref[i].x-finer[i].x));
                energy_drift=std::max(energy_drift,std::abs(energy(c.model,ref[i])-energy(c.model,c.initial)));
            }
            refs<<c.name<<','<<h<<','<<refinement<<','<<energy_drift<<'\n';
            if(refinement>1e-7 || energy_drift>1e-7*scale(c))
                throw std::runtime_error("reference refinement/energy check failed");
            for(int sweeps:{1,5,20}) for(Method method:{Method::Baseline,Method::Joint,Method::Immediate,
                    Method::Frozen,Method::Reverse,Method::Symmetric}) {
                for(double soft:{0.,.01,1.}) {
                    if(method==Method::Baseline && soft!=0) continue;
                    for(bool retain:{false,true}) {
                        if(retain && method!=Method::Joint && method!=Method::Frozen) continue;
                        for(int inner:{1,5,20}) {
                            if(method!=Method::Frozen && inner!=1) continue;
                            Options opt{method,sweeps,inner,retain};
                            // Full timestep traces for a representative coarse-step slice.
                            const bool record=h==1./60&&sweeps==5&&inner==1
                                &&(!retain||method==Method::Joint);
                            const auto m=run(c,opt,h,soft,steps,ref,record?&traces:nullptr);
                            summary_row(summary,c,opt,h,soft,steps,m);
                            ++runs; aborted+=m.completed<steps;
                        }
                    }
                }
            }
        }
        std::cout<<runs<<" runs; "<<aborted<<" stopped by explicit numerical guards; "
                 <<"reference refinement checks passed. Outputs: "<<dir<<'\n';
        return 0;
    } catch(const std::exception& e) {std::cerr<<e.what()<<'\n';return 1;}
}
