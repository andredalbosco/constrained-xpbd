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
    V operator+(V b) const{return{x+b.x,y+b.y};}
    V operator-(V b) const{return{x-b.x,y-b.y};}
    V operator*(double a) const{return{x*a,y*a};}
    V operator/(double a) const{return{x/a,y/a};}
    double dot(V b)const{return x*b.x+y*b.y;}
    double cross(V b)const{return x*b.y-y*b.x;}
    double norm2()const{return dot(*this);}
    double norm()const{return std::hypot(x,y);}
    bool finite()const{return std::isfinite(x)&&std::isfinite(y);}
};
struct State{V x,v;};
double energy(State s,V gravity){return .5*s.v.norm2()-gravity.dot(s.x);}
State initial(const std::string& name){
    const double theta=name=="equilibrium"?0:.65;
    const double speed=name=="gravity_tangent"?1.2:name=="free_rotation"?2:0;
    return{{std::sin(theta),-std::cos(theta)},
           {speed*std::cos(theta),speed*std::sin(theta)}};
}
struct Problem{
    State old; V g,center,pred; double h,scale,energy_offset=0;bool fixed_target=false;
    Problem(State s,V gravity,double dt,double target=std::numeric_limits<double>::quiet_NaN()):old(s),g(gravity),h(dt){
        center=old.x+g*(h*h); pred=center+old.v*h;
        scale=std::max({1.0,g.norm(),std::abs(energy(old,g)),.5*old.v.norm2()});
        fixed_target=std::isfinite(target);
        if(fixed_target){energy_offset=energy(old,g)-target;scale=std::max(scale,std::abs(target));}
    }
    double c(V x)const{return x.norm()-1;}
    // Complete-square form of the requested physical energy constraint. This
    // avoids manufacturing a tiny nonzero derivative at a resting predictor.
    double H(V x)const{
        return .5*((x-center).norm2()/(h*h)-old.v.norm2()-h*h*g.norm2())+energy_offset;
    }
    double physical_H(V x)const{
        const V dx=x-old.x;
        return .5*(dx.norm2()/(h*h)-old.v.norm2())-g.dot(dx)+energy_offset;
    }
    V grad(V x)const{return(x-center)/(h*h);}
    bool feasible()const{
        const double radius_squared=h*h*old.v.norm2()+h*h*h*h*g.norm2()-2*h*h*energy_offset;
        if(radius_squared<0)return false;
        const double r=std::sqrt(radius_squared);
        const double d=center.norm();
        return d<=1+r+1e-12 && d+std::min(1.0,r)>=std::max(1.0,r)-1e-12;
    }
};
struct Step{
    State next; std::string status="ok";
    double c=0,H=0,rx=0,max_condition=0,max_update=0,max_c=0,max_H=0,max_x=0;
    double min_sine=1,physical_form_difference=0,lambda=0,eta=0;
    double min_gradient_norm=std::numeric_limits<double>::infinity();
    int near_parallel=0,zero_energy_gradient=0,exact_singular=0,iterations=0,backtracks=0;
    bool converged=false,feasible=true;
};
const double infinity=std::numeric_limits<double>::infinity();

bool converged(const Problem& p,V x){return std::abs(p.c(x))<1e-12 && std::abs(p.H(x))<1e-11*p.scale;}

void sample(const Problem& p,V x,Step& out){
    out.c=p.c(x);out.H=p.physical_H(x);
    out.max_c=std::max(out.max_c,std::abs(out.c));
    out.max_H=std::max(out.max_H,std::abs(out.H));
    out.max_x=std::max(out.max_x,x.norm());
    out.physical_form_difference=std::max(out.physical_form_difference,std::abs(p.H(x)-out.H));
    if(!x.finite() || !std::isfinite(out.H)){out.status="nonfinite";return;}
    if(x.norm()==0){out.status="zero_radius";return;}
    const V J=x/x.norm(),G=p.grad(x);
    const double gn=G.norm();
    out.min_gradient_norm=std::min(out.min_gradient_norm,gn);
    if(gn==0){++out.zero_energy_gradient;out.max_condition=infinity;return;}
    const double sine=std::min(1.0,std::abs(J.cross(G))/gn);
    const double cosine=std::min(1.0,std::abs(J.dot(G))/gn);
    out.min_sine=std::min(out.min_sine,sine);
    if(sine<1e-6)++out.near_parallel;
    const double condition=sine==0?infinity:(1+cosine)*(1+cosine)/(sine*sine);
    out.max_condition=std::max(out.max_condition,condition);
}

void trace(std::ofstream* stream,const std::string& scene,const std::string& method,
           const Problem& p,int iteration,const std::string& phase,V x,const Step& out){
    if(!stream)return;
    const V G=p.grad(x);
    *stream<<(p.fixed_target?"fixed_initial":"previous_actual")<<','<<scene<<','<<method<<','<<p.h<<','<<iteration<<','<<phase<<','<<out.status<<','
        <<x.x<<','<<x.y<<','<<p.c(x)<<','<<p.H(x)<<','<<p.physical_H(x)<<','
        <<G.norm()<<','<<out.min_sine<<','<<out.max_condition<<','<<out.max_update<<','
        <<out.lambda<<','<<out.eta<<'\n';
}

Step advance(State old,V gravity,double h,int limit,const std::string& method,
             const std::string& scene,std::ofstream* traces=nullptr,
             double target=std::numeric_limits<double>::quiet_NaN()){
    Problem p(old,gravity,h,target); Step out;out.feasible=p.feasible(); V x=p.pred;
    sample(p,x,out);trace(traces,scene,method,p,-1,"prediction",x,out);
    const auto material=[&](){
        if(x.norm()==0){out.status="zero_radius";return;}
        const V J=x/x.norm();const double dl=-p.c(x);
        x=x+J*dl;out.lambda+=dl;out.max_update=std::max(out.max_update,std::abs(dl));
    };
    if(method=="xpbd"){
        material();out.iterations=1;sample(p,x,out);
        trace(traces,scene,method,p,0,"material",x,out);
    }else{
        if(method=="joint_seed" || method=="joint_seed_backtrack"){
            material();sample(p,x,out);trace(traces,scene,method,p,-1,"material_seed",x,out);
        }
        for(int iteration=0;iteration<limit && out.status=="ok";++iteration){
            if(converged(p,x))break;
            ++out.iterations;
            if(method=="sequential"){
                material();sample(p,x,out);trace(traces,scene,method,p,iteration,"material",x,out);
                if(out.status!="ok")break;
                if(std::abs(p.H(x))>=1e-11*p.scale){
                    const V G=p.grad(x);const double denominator=G.norm2();
                    if(denominator==0){++out.exact_singular;out.status="zero_energy_gradient";break;}
                    const double de=-p.H(x)/denominator;
                    const V dx=G*de;out.eta+=de;x=x+dx;
                    out.max_update=std::max(out.max_update,dx.norm());
                }
                sample(p,x,out);trace(traces,scene,method,p,iteration,"energy",x,out);
            }else{
                const V J=x/x.norm(),G=p.grad(x),tangent{-J.y,J.x};
                const double cross=J.cross(G),parallel=J.dot(G);
                // Orthogonal-component Schur: |G_perp|^2. No d-b*b cancellation.
                const double schur=cross*cross;
                const double rhs=-p.H(x)+parallel*p.c(x);
                if(schur==0){++out.exact_singular;out.status=G.norm()==0?"zero_energy_gradient":"rank_one_joint";break;}
                const double de=rhs/schur,dl=-p.c(x)-parallel*de;
                // Algebraically J*dl+G*de, evaluated without cancellation of huge parallel components.
                const V dx=J*(-p.c(x))+tangent*(rhs/cross);
                out.max_update=std::max(out.max_update,dx.norm());
                if(!dx.finite() || !std::isfinite(de) || !std::isfinite(dl)){
                    out.status="nonfinite_update";break;
                }
                double step_scale=1;
                if(method=="joint_seed_backtrack"){
                    const double merit=p.c(x)*p.c(x)+std::pow(p.H(x)/p.scale,2);
                    bool accepted=false;
                    for(int attempt=0;attempt<30;++attempt,step_scale*=.5){
                        const V candidate=x+dx*step_scale;
                        const double candidate_merit=p.c(candidate)*p.c(candidate)+std::pow(p.H(candidate)/p.scale,2);
                        if(std::isfinite(candidate_merit) && candidate_merit<merit){accepted=true;break;}
                        ++out.backtracks;
                    }
                    if(!accepted){out.status="line_search_rejected";break;}
                }
                x=x+dx*step_scale;out.lambda+=dl*step_scale;out.eta+=de*step_scale;
                sample(p,x,out);trace(traces,scene,method,p,iteration,"joint",x,out);
            }
        }
    }
    out.next={x,(x-old.x)/h};out.converged=converged(p,x);
    if(x.norm()>0)out.rx=(x-p.pred-(x/x.norm())*out.lambda-p.grad(x)*out.eta).norm();
    sample(p,x,out);
    if(!std::isfinite(out.lambda) || !std::isfinite(out.eta) || !std::isfinite(out.rx))
        out.status="nonfinite_dual_or_rx";
    trace(traces,scene,method,p,out.iterations,"final",x,out);
    return out;
}

void validate_baseline(){
    for(const std::string scene:{"gravity_rest","equilibrium","free_rotation"}){
        const State old=initial(scene);const V g=scene=="free_rotation"?V{}:V{0,-9.81};
        const double h=1.0/60;const Step result=advance(old,g,h,20,"xpbd",scene);
        const V predicted=old.x+old.v*h+g*(h*h),expected=predicted/predicted.norm();
        if((result.next.x-expected).norm()>1e-14 || std::abs(result.next.x.norm()-1)>1e-14)
            throw std::runtime_error("baseline XPBD projection mismatch");
        if((result.next.v-(expected-old.x)/h).norm()>1e-12)
            throw std::runtime_error("baseline XPBD velocity mismatch");
        if(scene=="equilibrium" && result.next.v.norm()>1e-14)
            throw std::runtime_error("baseline equilibrium moved");
        if(scene=="free_rotation"){
            const double exact_energy=(1-1/std::sqrt(1+4*h*h))/(h*h);
            if(std::abs(energy(result.next,g)-exact_energy)>1e-11)
                throw std::runtime_error("baseline free-rotation energy mismatch");
        }
    }
    for(double h:{1.0/240,1.0/60,1.0/15,1.0/5}){
        const V gravity{0,-9.81};const State rest=initial("gravity_rest");
        const Problem p(rest,gravity,h);
        const Problem shifted(rest,gravity,h,energy(rest,gravity)-1);
        const Problem impossible(rest,gravity,h,energy(rest,gravity)-1000);
        if((shifted.grad(p.pred)-p.grad(p.pred)).norm()!=0 ||
           std::abs(shifted.H(p.pred)-p.H(p.pred)-1)>1e-12 || impossible.feasible())
            throw std::runtime_error("fixed-target shift or negative-radius feasibility check failed");
        if(p.grad(p.pred).norm()!=0 || !(p.H(p.pred)<0) ||
           std::abs(p.H(p.pred)+.5*h*h*gravity.norm2())>1e-13)
            throw std::runtime_error("resting predictor gradient/energy identity failed");
        const V projected=p.pred/p.pred.norm(),J=projected/projected.norm(),G=p.grad(projected);
        if(std::abs(J.cross(G))>1e-10*std::max(1.0,G.norm()) || !(p.H(projected)<0))
            throw std::runtime_error("projected resting gradients are not parallel");
        if(advance(rest,gravity,h,80,"joint_raw","gravity_rest").status!="zero_energy_gradient")
            throw std::runtime_error("raw resting joint solve did not expose its required singularity");
        if(advance(rest,gravity,h,80,"sequential","gravity_rest").converged)
            throw std::runtime_error("resting sequential test unexpectedly lost its unresolved row");
        for(const std::string scene:{"gravity_rest","gravity_tangent","equilibrium","free_rotation"}){
            State model=initial(scene);const V g=scene=="free_rotation"?V{}:gravity;
            for(int frame=0;frame<static_cast<int>(std::lround(10/h));++frame){
                const Step step=advance(model,g,h,1,"xpbd",scene);
                model=step.next;
                if(step.status!="ok" || !model.x.finite() || !model.v.finite() || std::abs(model.x.norm()-1)>1e-13)
                    throw std::runtime_error("baseline ten-second length/finite-state validation failed");
            }
        }
    }
}

int main(){
    try{
        validate_baseline();
        std::ofstream summary("summary.csv"),traces("one_step_traces.csv");
        if(!summary || !traces)throw std::runtime_error("cannot open study outputs");
        summary<<std::setprecision(17);traces<<std::setprecision(17);std::cout<<std::setprecision(10);
        summary<<"energy_target,scene,method,h,iteration_limit,accepted_steps,status,unconverged_steps,infeasible_inputs,"
            "final_energy_change,max_accepted_energy_error,max_accepted_C,max_accepted_H,max_accepted_rx,"
            "max_accepted_radial_velocity,max_accepted_angular_change,max_scaled_condition,min_angle_sine,"
            "near_parallel_samples,zero_gradient_samples,exact_singular_events,max_attempted_update,"
            "max_attempted_C,max_attempted_H,max_attempted_radius,last_candidate_C,last_candidate_H,"
            "backtracks,max_energy_form_difference,completion_class,max_accepted_speed,max_accepted_radius,"
            "rx_unclosed_steps,minimum_energy_gradient_norm\n";
        traces<<"energy_target,scene,method,h,iteration,phase,status,x,y,C,H_complete_square,H_physical,energy_gradient_norm,"
            "minimum_angle_sine,maximum_scaled_condition,maximum_attempted_update,lambda,eta\n";
        for(const std::string target_mode:{"previous_actual","fixed_initial"})
        for(const std::string scene:{"gravity_rest","gravity_tangent","equilibrium","free_rotation"})
            for(const std::string method:{"xpbd","joint_raw","joint_seed","sequential","joint_seed_backtrack"})
                for(double h:{1.0/240,1.0/60,1.0/15,1.0/5})for(int limit:{1,5,20,80}){
                    State model=initial(scene);const V gravity=scene=="free_rotation"?V{}:V{0,-9.81};
                    const double e0=energy(model,gravity),l0=model.x.cross(model.v);
                    const double target=target_mode=="fixed_initial"?e0:std::numeric_limits<double>::quiet_NaN();
                    double max_e=0,max_c=0,max_H=0,max_rx=0,max_radial=0,max_l=0,max_cond=0,min_sine=1;
                    double max_update=0,max_attempt_c=0,max_attempt_H=0,max_radius=0,form_difference=0;
                    double max_speed=0,max_accepted_radius=0;
                    double min_gradient_norm=infinity;
                    int accepted=0,incomplete=0,infeasible=0,parallel=0,zeros=0,singular=0,backtracks=0;
                    int unclosed_rx=0;
                    Step report;
                    for(int k=0;k<static_cast<int>(std::lround(10/h));++k){
                        report=advance(model,gravity,h,limit,method,scene,k==0&&limit==80?&traces:nullptr,target);
                        infeasible+=!report.feasible;parallel+=report.near_parallel;zeros+=report.zero_energy_gradient;
                        singular+=report.exact_singular;backtracks+=report.backtracks;
                        max_cond=std::max(max_cond,report.max_condition);min_sine=std::min(min_sine,report.min_sine);
                        min_gradient_norm=std::min(min_gradient_norm,report.min_gradient_norm);
                        max_update=std::max(max_update,report.max_update);
                        max_attempt_c=std::max(max_attempt_c,report.max_c);
                        max_attempt_H=std::max(max_attempt_H,report.max_H);
                        max_radius=std::max(max_radius,report.max_x);
                        form_difference=std::max(form_difference,report.physical_form_difference);
                        if(report.status!="ok"){
                            if(k>0 && limit==80)
                                advance(model,gravity,h,limit,method,scene+"_failure_at_"+std::to_string(k),&traces,target);
                            break;
                        }
                        model=report.next;++accepted;incomplete+=!report.converged;
                        max_e=std::max(max_e,std::abs(energy(model,gravity)-e0));
                        max_c=std::max(max_c,std::abs(report.c));max_H=std::max(max_H,std::abs(report.H));
                        max_rx=std::max(max_rx,report.rx);
                        max_speed=std::max(max_speed,model.v.norm());
                        max_accepted_radius=std::max(max_accepted_radius,model.x.norm());
                        unclosed_rx+=report.rx>1e-11;
                        max_radial=std::max(max_radial,std::abs(model.x.dot(model.v)/model.x.norm()));
                        max_l=std::max(max_l,std::abs(model.x.cross(model.v)-l0));
                    }
                    std::string classification;
                    if(max_accepted_radius>10 || max_speed>100)classification="runaway_accepted";
                    else if(max_update>10 || max_radius>10)classification="runaway_trial";
                    else if(report.status=="zero_energy_gradient" || report.status=="rank_one_joint")classification="singular";
                    else if(report.status!="ok")classification="rejected_or_nonfinite";
                    else if(incomplete)classification="bounded_constraints_nonconverged";
                    else if(unclosed_rx)classification="bounded_constraints_converged_rx_unclosed";
                    else classification="bounded_constraints_and_rx_converged";
                    summary<<target_mode<<','<<scene<<','<<method<<','<<h<<','<<limit<<','<<accepted<<','<<report.status<<','
                        <<incomplete<<','<<infeasible<<','<<energy(model,gravity)-e0<<','<<max_e<<','<<max_c<<','
                        <<max_H<<','<<max_rx<<','<<max_radial<<','<<max_l<<','<<max_cond<<','<<min_sine<<','
                        <<parallel<<','<<zeros<<','<<singular<<','<<max_update<<','<<max_attempt_c<<','
                        <<max_attempt_H<<','<<max_radius<<','<<report.c<<','<<report.H<<','<<backtracks<<','
                        <<form_difference<<','<<classification<<','<<max_speed<<','<<max_accepted_radius<<','
                        <<unclosed_rx<<','<<min_gradient_norm<<'\n';
                    if(limit==80)std::cout<<target_mode<<' '<<scene<<' '<<method<<" h="<<h<<" accepted="<<accepted
                        <<" status="<<report.status<<" maxE="<<max_e<<" maxC="<<max_c
                        <<" nearparallel="<<parallel<<" maxupdate="<<max_update<<'\n';
                }
        summary.flush();traces.flush();if(!summary||!traces)throw std::runtime_error("failed writing outputs");
        std::cout<<"Baseline checks passed; failures of experimental variants are recorded, not regularized.\n";
    }catch(const std::exception&error){std::cerr<<error.what()<<'\n';return 1;}
}
