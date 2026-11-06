#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

struct V {
    double x=0,y=0;
    V operator+(V b) const { return {x+b.x,y+b.y}; }
    V operator-(V b) const { return {x-b.x,y-b.y}; }
    V operator*(double a) const { return {a*x,a*y}; }
    V operator/(double a) const { return {x/a,y/a}; }
    double dot(V b) const { return x*b.x+y*b.y; }
    double norm2() const { return dot(*this); }
    double norm() const { return std::hypot(x,y); }
    bool finite() const { return std::isfinite(x)&&std::isfinite(y); }
};
struct State { V x,v; };
struct Rows { V j,g; double c,H,rC,a,b,d,schur,condition,schur_ratio; };
struct Stats {
    double max_condition=1,min_schur_ratio=1,max_proposed=0,max_accepted=0;
    double max_lambda=0,max_eta=0,max_x=0,max_speed=0,max_H=0,max_material=0;
    double max_stationarity=0;
    int rejected=0,singular=0,negative_raw_determinants=0,backtracks=0;
    bool finite_failure=false;
};
double potential(V x,double k) { const double c=x.norm()-1; return .5*k*c*c; }
double energy(State s,double k) { return .5*s.v.norm2()+potential(s.x,k); }

Rows rows(V x,const State& old,double lambda,double h,double k,double target_energy) {
    if(!x.finite() || x.norm()<1e-14) throw std::runtime_error("nonfinite_or_collapsed_position");
    const V j=x/x.norm(), dx=x-old.x;
    const double c=x.norm()-1, alpha=1/(k*h*h);
    const double H=.5*dx.norm2()/(h*h)+potential(x,k)-target_energy;
    // Always the analytical endpoint H gradient, regardless of multiplier consistency.
    const V g=dx/(h*h)+j*(k*c);
    const double a=1+alpha, parallel=j.dot(g);
    const V perpendicular=g-j*parallel;
    const double d=g.norm2();
    // Positive Schur complement; no cancellation of a*d - b*b.
    const double schur=perpendicular.norm2()+(alpha/a)*parallel*parallel;
    const double ratio=d>0?schur/d:0;
    const double rho=std::sqrt(std::max(0.0,1-ratio));
    const double condition=ratio>0?(1+rho)*(1+rho)/ratio:std::numeric_limits<double>::infinity();
    if(!std::isfinite(H)||!g.finite()||!std::isfinite(c+alpha*lambda))
        throw std::runtime_error("nonfinite_energy_or_row");
    return {j,g,c,H,c+alpha*lambda,a,parallel,d,schur,condition,ratio};
}

void record(Stats& s,const Rows& r) {
    if(r.d>0) {
        s.max_condition=std::max(s.max_condition,r.condition);
        s.min_schur_ratio=std::min(s.min_schur_ratio,r.schur_ratio);
        if(r.a*r.d-r.b*r.b<=0) ++s.negative_raw_determinants;
    }
}
bool tiny_gradient(const Rows& r,V x,const State& old,double h,double k) {
    const double scale=1+(x-old.x).norm()/(h*h)+std::abs(k*r.c);
    return r.g.norm()<=1e-12*scale;
}

struct Trace {
    std::ofstream* stream=nullptr;
    std::string scenario,method;
    double h=0,k=0;
    void write(int iteration,const char* phase,V x,const State& old,double lambda,double eta,
               const Rows& r,const Stats& s) {
        if(!stream) return;
        *stream<<scenario<<','<<method<<','<<k<<','<<h<<','<<iteration<<','<<phase<<','
               <<x.x<<','<<x.y<<','<<(x.x-old.x.x)/h<<','<<(x.y-old.x.y)/h<<','
               <<r.H<<','<<r.rC<<','<<lambda<<','<<eta<<','<<r.g.norm()<<','
               <<(r.d>0?r.condition:0)<<','<<s.max_proposed<<','<<s.max_accepted<<','
               <<s.rejected<<','<<s.singular<<'\n';
    }
};

State advance(State old,double h,double k,int iterations,const std::string& method,
              Stats& stats,Trace trace,double& final_H,double& final_material,double target_energy) {
    V x=old.x+old.v*h;
    double lambda=0,eta=0;
    const bool joint=method=="joint"||method=="joint_backtracking";
    const bool sequential=method=="sequential"||method=="sequential_backtracking"||method=="sequential_raw";
    const bool backtracking=method=="joint_backtracking"||method=="sequential_backtracking";
    const bool raw=method=="sequential_raw";
    const double energy_tolerance=1e-12*(1+std::abs(target_energy));
    for(int iteration=0;iteration<iterations;++iteration) {
        Rows r=rows(x,old,lambda,h,k,target_energy);
        record(stats,r);
        trace.write(iteration,"start",x,old,lambda,eta,r,stats);
        if(std::abs(r.rC)<1e-13 && (method=="xpbd"||std::abs(r.H)<energy_tolerance)) break;
        if(joint) {
            if(tiny_gradient(r,x,old,h,k)) {
                if(std::abs(r.H)>energy_tolerance) ++stats.singular;
                ++stats.rejected;
                continue;
            }
            const double d_eta=(-r.H+r.b*r.rC/r.a)/r.schur;
            const double d_lambda=(-r.rC-r.b*d_eta)/r.a;
            const V dx=r.j*d_lambda+r.g*d_eta;
            stats.max_proposed=std::max(stats.max_proposed,dx.norm());
            bool accepted=false;
            double scale=1;
            const double merit=r.rC*r.rC/r.a+r.H*r.H/r.d;
            for(int attempt=0;attempt<(backtracking?40:1);++attempt,scale*=.5) {
                const V candidate=x+dx*scale;
                const double next_lambda=lambda+scale*d_lambda;
                Rows next=rows(candidate,old,next_lambda,h,k,target_energy);
                const double next_merit=next.rC*next.rC/r.a+next.H*next.H/r.d;
                if(!backtracking || next_merit<1e-26 || next_merit<=merit*(1-1e-4*scale)) {
                    x=candidate;lambda=next_lambda;eta+=scale*d_eta;
                    stats.max_accepted=std::max(stats.max_accepted,scale*dx.norm());
                    accepted=true;break;
                }
                ++stats.backtracks;
            }
            if(!accepted) ++stats.rejected;
            r=rows(x,old,lambda,h,k,target_energy);
            trace.write(iteration,"joint",x,old,lambda,eta,r,stats);
        } else {
            // Ordinary XPBD material update, relinearized independently of the energy row.
            const double d_lambda=-r.rC/r.a;
            const V dx=r.j*d_lambda;
            stats.max_proposed=std::max(stats.max_proposed,dx.norm());
            x=x+dx;lambda+=d_lambda;
            stats.max_accepted=std::max(stats.max_accepted,dx.norm());
            r=rows(x,old,lambda,h,k,target_energy);
            trace.write(iteration,"material",x,old,lambda,eta,r,stats);
            if(sequential && std::abs(r.H)>energy_tolerance) {
                record(stats,r);
                if((!raw&&tiny_gradient(r,x,old,h,k)) || r.d==0) {
                    ++stats.singular;++stats.rejected;
                } else {
                    const double d_eta=-r.H/r.d;
                    const V dx_energy=r.g*d_eta;
                    stats.max_proposed=std::max(stats.max_proposed,dx_energy.norm());
                    bool accepted=false;
                    double scale=1;
                    for(int attempt=0;attempt<(backtracking?40:1);++attempt,scale*=.5) {
                        const V candidate=x+dx_energy*scale;
                        const Rows next=rows(candidate,old,lambda,h,k,target_energy);
                        if(!backtracking || std::abs(next.H)<energy_tolerance
                           || std::abs(next.H)<=std::abs(r.H)*(1-1e-4*scale)) {
                            x=candidate;eta+=scale*d_eta;
                            stats.max_accepted=std::max(stats.max_accepted,scale*dx_energy.norm());
                            accepted=true;break;
                        }
                        ++stats.backtracks;
                    }
                    if(!accepted) ++stats.rejected;
                }
                r=rows(x,old,lambda,h,k,target_energy);
                trace.write(iteration,"energy",x,old,lambda,eta,r,stats);
            }
        }
        stats.max_lambda=std::max(stats.max_lambda,std::abs(lambda));
        stats.max_eta=std::max(stats.max_eta,std::abs(eta));
        if(!x.finite()||!std::isfinite(lambda)||!std::isfinite(eta)) {
            stats.finite_failure=true;
            throw std::runtime_error("nonfinite_iterate");
        }
    }
    const Rows r=rows(x,old,lambda,h,k,target_energy);
    final_H=r.H;final_material=r.rC;
    stats.max_H=std::max(stats.max_H,std::abs(r.H));
    stats.max_material=std::max(stats.max_material,std::abs(r.rC));
    const V stationarity=x-(old.x+old.v*h)-r.j*lambda-r.g*eta;
    if(!stationarity.finite()) {
        stats.finite_failure=true;
        throw std::runtime_error("nonfinite_stationarity_residual");
    }
    stats.max_stationarity=std::max(stats.max_stationarity,stationarity.norm());
    const State next{x,(x-old.x)/h};
    stats.max_x=std::max(stats.max_x,next.x.norm());
    stats.max_speed=std::max(stats.max_speed,next.v.norm());
    if(!next.x.finite()||!next.v.finite()||!std::isfinite(energy(next,k))) {
        stats.finite_failure=true;
        throw std::runtime_error("nonfinite_accepted_state");
    }
    return next;
}

State initial(const std::string& scenario) {
    if(scenario=="turning") return {{1.25,0},{0,0}};
    if(scenario=="rotating") return {{1.25,0},{0,1}};
    return {{1.25,0},{-.7,1}};
}

void run_case(const std::string& scenario,const std::string& method,double k,double h,
              int iterations,double duration,const std::string& target_policy,std::ofstream& out) {
    State state=initial(scenario);
    const State start=state;
    const double e0=energy(state,k);
    Stats stats;
    stats.max_x=state.x.norm();stats.max_speed=state.v.norm();
    double max_error=0,max_motion=0,final_H=0,final_material=0;
    int completed=0,unresolved=0;
    std::string status="completed";
    const int steps=static_cast<int>(std::lround(duration/h));
    try {
        for(int step=0;step<steps;++step) {
            const double previous_energy=energy(state,k);
            const double target_energy=target_policy=="initial"?e0:previous_energy;
            state=advance(state,h,k,iterations,method,stats,{},final_H,final_material,target_energy);
            ++completed;
            if(method!="xpbd" && (std::abs(final_H)>1e-10*(1+std::abs(target_energy))
                                  ||std::abs(final_material)>1e-10)) ++unresolved;
            max_error=std::max(max_error,std::abs(energy(state,k)-e0));
            max_motion=std::max(max_motion,(state.x-start.x).norm());
            if(state.x.norm()>1e100 || state.v.norm()>1e100) {
                status="magnitude_guard";break;
            }
        }
    } catch(const std::exception& error) {
        status=error.what();
        if(status.find("nonfinite")!=std::string::npos) stats.finite_failure=true;
    }
    if(status=="completed"&&scenario=="turning"&&max_motion<1e-9) status="frozen_turning_state";
    if(status=="completed"&&stats.singular>0) status="completed_with_singular_energy_updates";
    if(status=="completed") status=unresolved>0?"completed_unresolved":"completed_rows_converged";
    if(status.rfind("completed",0)==0 && (stats.max_x>10||stats.max_speed>100)) status="runaway_accepted";
    out<<scenario<<','<<method<<','<<target_policy<<','<<k<<','<<h<<','<<iterations<<','<<completed<<','<<completed*h<<','<<status<<','
       <<e0<<','<<energy(state,k)-e0<<','<<max_error<<','<<max_error/e0<<','<<max_motion<<','
       <<stats.max_H<<','<<stats.max_material<<','<<stats.max_lambda<<','<<stats.max_eta<<','
       <<stats.max_condition<<','<<stats.min_schur_ratio<<','<<stats.max_x<<','<<stats.max_speed<<','
       <<stats.max_proposed<<','<<stats.max_accepted<<','<<stats.rejected<<','<<stats.singular<<','
       <<stats.negative_raw_determinants<<','<<stats.backtracks<<','<<stats.finite_failure<<','<<unresolved<<','
       <<stats.max_stationarity<<'\n';
}

void mandatory_checks() {
    const double k=100,h=1./60;
    const State old=initial("turning");
    Stats baseline;
    double H=0,rC=0;
    const State next=advance(old,h,k,20,"xpbd",baseline,{},H,rC,energy(old,k));
    const double expected_q=.25/(1+k*h*h);
    if(std::abs(next.x.x-(1+expected_q))>1e-14 || std::abs(next.v.x-(expected_q-.25)/h)>1e-12)
        throw std::runtime_error("mandatory analytical BE recurrence check failed");
    const Rows after_material=rows(next.x,old,baseline.max_lambda*-1,h,k,energy(old,k));
    if(after_material.g.norm()>1e-10 || !(after_material.H<-.08))
        throw std::runtime_error("mandatory sequential zero-gradient check failed");
    Stats coupled;
    const State frozen=advance(old,h,k,20,"joint",coupled,{},H,rC,energy(old,k));
    if((frozen.x-old.x).norm()>1e-14 || frozen.v.norm()>1e-12
       ||std::abs(coupled.max_lambda-h*h*k*.25)>1e-14||std::abs(coupled.max_eta-h*h)>1e-14
       ||std::abs(H)>1e-12||std::abs(rC)>1e-12)
        throw std::runtime_error("mandatory joint frozen-branch check failed");
    State state=old;
    double q=.25,v=0;
    for(int i=0;i<600;++i) {
        const double before=energy(state,k);
        state=advance(state,h,k,1,"xpbd",baseline,{},H,rC,before);
        const double next_q=(q+h*v)/(1+k*h*h),next_v=(next_q-q)/h;
        q=next_q;v=next_v;
        if(!state.x.finite()||!state.v.finite()||std::abs(state.x.x-(1+q))>1e-12
           ||std::abs(state.v.x-v)>1e-10||state.x.norm()>1.251||state.v.norm()>2.501
           ||energy(state,k)>before+1e-10)
            throw std::runtime_error("mandatory bounded baseline trajectory check failed");
    }
    std::cout<<"Mandatory checks passed: BE recurrence/bounds, joint frozen branch/multipliers, sequential zero gradient.\n";
}

void selected_trajectories(double duration) {
    std::ofstream out("coupling_selected_trajectories.csv");
    out.exceptions(std::ios::badbit|std::ios::failbit);
    out<<std::setprecision(17);
    out<<"method,energy_target,iterations,time,energy,relative_energy_change,x,y,vx,vy,H,rC,condition\n";
    const double k=1000,h=1./60;
    for(const std::string method:{"xpbd","joint","sequential"})
        for(const std::string target_policy:{"initial","previous_actual"})
            for(int iterations:{1,5}) {
                State state=initial("rotating");
                const double e0=energy(state,k);
                Stats stats;
                for(int step=0;step<static_cast<int>(std::lround(duration/h));++step) {
                    double H=0,rC=0;
                    const double target=target_policy=="initial"?e0:energy(state,k);
                    state=advance(state,h,k,iterations,method,stats,{},H,rC,target);
                    out<<method<<','<<target_policy<<','<<iterations<<','<<(step+1)*h<<','
                       <<energy(state,k)<<','<<energy(state,k)/e0-1<<','<<state.x.x<<','<<state.x.y<<','
                       <<state.v.x<<','<<state.v.y<<','<<H<<','<<rC<<','<<stats.max_condition<<'\n';
                }
            }
    out.close();
}

int main(int argc,char**argv) {
    try {
        const double duration=argc>1?std::stod(argv[1]):10;
        if(!(duration>0)||!std::isfinite(duration)) throw std::runtime_error("invalid duration");
        mandatory_checks();
        std::ofstream out("coupling_results.csv"),traces("coupling_traces.csv");
        out.exceptions(std::ios::badbit|std::ios::failbit);
        traces.exceptions(std::ios::badbit|std::ios::failbit);
        out<<std::setprecision(17);traces<<std::setprecision(17);
        out<<"scenario,method,energy_target,k,h,iterations,completed_steps,duration,status,initial_energy,final_energy_change,"
             "max_energy_error,max_relative_energy_error,max_motion,max_step_energy_residual,max_material_residual,"
             "max_abs_lambda,max_abs_eta,max_scaled_condition,min_schur_ratio,max_position,max_speed,max_proposed_update,"
             "max_accepted_update,rejected_updates,singular_energy_updates,nonpositive_raw_determinants,backtracks,"
             "finite_failure,unresolved_steps,max_stationarity_residual\n";
        traces<<"scenario,method,k,h,iteration,phase,x,y,vx,vy,H,rC,lambda,eta,energy_gradient_norm,"
                "scaled_condition,max_proposed_update,max_accepted_update,rejected,singular\n";
        const std::vector<std::string> methods={"xpbd","joint","sequential","joint_backtracking",
                                               "sequential_backtracking","sequential_raw"};
        for(const auto& scenario:{"turning","rotating","mixed"})
            for(double k:{100.,1000.})
                for(double h:{1./240,1./60,1./15})
                    for(const auto& method:methods)
                        for(int iterations:{1,5,20,80})
                            for(const std::string target_policy:{"initial","previous_actual"})
                                run_case(scenario,method,k,h,iterations,duration,target_policy,out);
        for(const auto& scenario:{"turning","rotating","mixed"})
            for(const auto& method:methods)
                for(int setting:{0,1,2}) {
                    const double h=setting==1?1./15:1./60,k=setting?1000:100;
                    State state=initial(scenario);
                    Stats stats;
                    double H=0,rC=0;
                    try { advance(state,h,k,80,method,stats,{&traces,scenario,method,h,k},H,rC,energy(state,k)); }
                    catch(const std::exception& error) {
                        std::cerr<<"trace "<<scenario<<' '<<method<<" k="<<k<<" h="<<h<<": "<<error.what()<<'\n';
                    }
                }
        selected_trajectories(duration);
        out.close();traces.close();
        std::cout<<"Wrote 864 duration-study rows and one-step traces. Outcomes are recorded, not forced to pass.\n";
    } catch(const std::exception& error) {
        std::cerr<<error.what()<<'\n';return 1;
    }
}
