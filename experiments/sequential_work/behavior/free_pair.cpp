#include "xpbd.hpp"
#include <Eigen/Geometry>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>

using namespace cxpbd;

Model rotating_pair() {
    Model model;
    Particle a,b;
    a.mass=2; b.mass=1;
    const double total=a.mass+b.mass;
    const Vec3 center(.3,.4,.2),relative(1.25,0,0),translation(.25,-.1,.15),velocity(0,1,.2);
    a.x=center+(b.mass/total)*relative;
    b.x=center-(a.mass/total)*relative;
    a.v=translation+(b.mass/total)*velocity;
    b.v=translation-(a.mass/total)*velocity;
    model.particles={a,b};
    Spring spring;
    spring.a=0; spring.b=1; spring.rest_length=1; spring.compliance=.01;
    model.springs={spring};
    return model;
}

struct Measures {
    double K=0,U=0,R=0,S=0;
    Vec3 P=Vec3::Zero(),L=Vec3::Zero(),Lcm=Vec3::Zero();
    double E() const { return K+U; }
};

Measures measure(const Model& model) {
    Measures out;
    Vec3 center=Vec3::Zero();
    double mass=0;
    for(const auto& p:model.particles) {
        mass+=p.mass; center+=p.mass*p.x;
        out.K+=.5*p.mass*p.v.squaredNorm(); out.R+=p.reservoir;
        out.P+=p.mass*p.v; out.L+=p.x.cross(p.mass*p.v);
    }
    center/=mass;
    out.Lcm=out.L-center.cross(out.P);
    for(const auto& s:model.springs) {
        const double c=(model.particles[s.a].x-model.particles[s.b].x).norm()-s.rest_length;
        out.U+=.5*c*c/s.compliance; out.S+=s.reservoir;
    }
    return out;
}

struct Bias {
    double total=0,quadratic=0,geometry=0,history=0,error=0;
};

Bias bias(const Model& old,const Model& next) {
    const auto& s=old.springs[0];
    const Vec3 rold=old.particles[0].x-old.particles[1].x;
    const Vec3 rnew=next.particles[0].x-next.particles[1].x;
    const Vec3 n=rold/rold.norm(),dr=rnew-rold;
    const double cold=rold.norm()-s.rest_length,cnew=rnew.norm()-s.rest_length,dc=cnew-cold;
    const double analytic=-cold/s.compliance;
    const double used=old.previous_dt>0?s.previous_lambda/(old.previous_dt*old.previous_dt):analytic;
    Bias out;
    out.total=.5*(cnew*cnew-cold*cold)/s.compliance+used*n.dot(dr);
    out.quadratic=.5*dc*dc/s.compliance;
    out.geometry=cold/s.compliance*(dc-n.dot(dr));
    out.history=(used-analytic)*n.dot(dr);
    out.error=out.total-out.quadratic-out.geometry-out.history;
    return out;
}

void vector_csv(std::ofstream& out,const Vec3& v) { out<<','<<v.x()<<','<<v.y()<<','<<v.z(); }

void write_row(std::ofstream& out,const std::string& method,double h,double compliance,int index,
               const Model& model,const Measures& start,const Diagnostics& d,const Bias& b) {
    const auto m=measure(model);
    out<<method<<','<<h<<','<<compliance<<','<<index<<','<<index*h<<','<<m.E()<<','<<m.K<<','<<m.U
       <<','<<m.E()-start.E()<<','<<m.R<<','<<m.S<<','<<m.E()+m.R+m.S-start.E();
    vector_csv(out,m.P); vector_csv(out,m.L); vector_csv(out,m.Lcm);
    out<<','<<(m.P-start.P).norm()<<','<<(m.L-start.L).norm()<<','<<(m.Lcm-start.Lcm).norm()
       <<','<<model.particles[0].reservoir<<','<<model.particles[1].reservoir
       <<','<<b.total<<','<<b.quadratic<<','<<b.geometry<<','<<b.history<<','<<b.error
       <<','<<d.max_material_residual<<','<<d.max_feedback_residual<<','<<d.position_residual_norm
       <<','<<d.energy_balance_error<<','<<d.backtracks<<','<<d.rejected_updates;
    for(const auto& p:model.particles) { vector_csv(out,p.x); vector_csv(out,p.v); }
    out<<'\n';
}

void require_finite(const Model& model,const Diagnostics& d) {
    for(const auto& p:model.particles)
        if(!p.x.allFinite()||!p.v.allFinite()||!std::isfinite(p.reservoir))
            throw std::runtime_error("nonfinite particle");
    for(const auto& s:model.springs)
        if(!std::isfinite(s.lambda)||!std::isfinite(s.previous_lambda)||!std::isfinite(s.reservoir))
            throw std::runtime_error("nonfinite spring");
    const double values[]={d.mechanical_energy,d.reservoir_energy,d.energy_balance_error,
        d.max_balance_residual,d.max_feedback_residual,d.max_material_residual,d.position_residual_norm};
    for(double value:values) if(!std::isfinite(value)) throw std::runtime_error("nonfinite diagnostic");
}

int main() {
    try {
        std::ofstream trajectory("trajectories.csv"),summary("summary.csv");
        trajectory.exceptions(std::ios::badbit|std::ios::failbit);
        summary.exceptions(std::ios::badbit|std::ios::failbit);
        trajectory<<std::setprecision(17); summary<<std::setprecision(17);
        trajectory<<"method,h,beta_gamma,step,time,E,K,U,delta_E,R,S,extended_error,Px,Py,Pz,Lx,Ly,Lz,"
            "Lcm_x,Lcm_y,Lcm_z,P_error,L_error,Lcm_error,R0,R1,d,d_quadratic,d_geometry,d_history,"
            "d_identity_error,material_residual,feedback_residual,position_residual,step_account_error,"
            "backtracks,rejections,x0,y0,z0,vx0,vy0,vz0,x1,y1,z1,vx1,vy1,vz1\n";
        summary<<"method,h,beta_gamma,sweeps,duration,status,completed_steps,E0,Efinal,final_delta_E,"
            "max_E_error,P0_norm,L0_norm,Lcm0_norm,max_P_error,final_P_error,max_L_error,final_L_error,"
            "max_Lcm_error,final_Lcm_error,Rfinal,Sfinal,max_abs_extended_error,max_material_residual,"
            "max_feedback_residual,max_position_residual,final_material_residual,final_feedback_residual,"
            "final_position_residual,sum_d_quadratic,sum_d_geometry,sum_d_history,max_d_identity_error,"
            "positive_d_steps,negative_d_steps,backtracks,rejections\n";
        for(double h:{1./240,1./60}) for(double compliance:{0.,1.,1000.,1e6}) {
            const std::string method=compliance==0?"xpbd":"sequential";
            Model model=rotating_pair();
            const auto start=measure(model);
            Settings settings;
            settings.dt=h;settings.iterations=20;
            settings.mode=compliance==0?Mode::XPBD:Mode::Active;
            settings.active_order=ActiveOrder::Sequential;
            settings.kinetic_reservoir_compliance=compliance==0?1000:compliance;
            settings.elastic_reservoir_compliance=compliance==0?1000:compliance;
            double maxe=0,maxp=0,maxl=0,maxlcm=0,maxextended=0,maxmaterial=0,maxfeedback=0,maxposition=0;
            double sumquadratic=0,sumgeometry=0,sumhistory=0,maxidentity=0;
            int positive=0,negative=0,completed=0;
            long long backtracks=0,rejections=0;
            std::string status="completed";
            Diagnostics last;
            write_row(trajectory,method,h,compliance,0,model,start,last,{});
            try {
                for(int index=1;index<=static_cast<int>(std::lround(10/h));++index) {
                    Model candidate=model;
                    const auto d=step(candidate,settings);
                    require_finite(candidate,d);
                    const auto m=measure(candidate);
                    const auto b=bias(model,candidate);
                    if(!std::isfinite(b.error)) throw std::runtime_error("nonfinite bias diagnostic");
                    maxe=std::max(maxe,std::abs(m.E()-start.E()));
                    maxp=std::max(maxp,(m.P-start.P).norm());
                    maxl=std::max(maxl,(m.L-start.L).norm());
                    maxlcm=std::max(maxlcm,(m.Lcm-start.Lcm).norm());
                    maxextended=std::max(maxextended,std::abs(m.E()+m.R+m.S-start.E()));
                    maxmaterial=std::max(maxmaterial,d.max_material_residual);
                    maxfeedback=std::max(maxfeedback,d.max_feedback_residual);
                    maxposition=std::max(maxposition,d.position_residual_norm);
                    sumquadratic+=b.quadratic;sumgeometry+=b.geometry;sumhistory+=b.history;
                    maxidentity=std::max(maxidentity,std::abs(b.error));
                    positive+=b.total>1e-14;negative+=b.total<-1e-14;
                    backtracks+=d.backtracks;rejections+=d.rejected_updates;
                    model=std::move(candidate);last=d;++completed;
                    write_row(trajectory,method,h,compliance,index,model,start,d,b);
                }
            } catch(const std::exception& e) { status=e.what(); }
            const auto final=measure(model);
            if(method=="xpbd"&&maxp>1e-9) throw std::runtime_error("baseline momentum check failed");
            summary<<method<<','<<h<<','<<compliance<<",20,10,"<<status<<','<<completed
                <<','<<start.E()<<','<<final.E()<<','<<final.E()-start.E()<<','<<maxe
                <<','<<start.P.norm()<<','<<start.L.norm()<<','<<start.Lcm.norm()
                <<','<<maxp<<','<<(final.P-start.P).norm()<<','<<maxl<<','<<(final.L-start.L).norm()
                <<','<<maxlcm<<','<<(final.Lcm-start.Lcm).norm()<<','<<final.R<<','<<final.S
                <<','<<maxextended<<','<<maxmaterial<<','<<maxfeedback<<','<<maxposition
                <<','<<last.max_material_residual<<','<<last.max_feedback_residual<<','<<last.position_residual_norm
                <<','<<sumquadratic<<','<<sumgeometry<<','<<sumhistory<<','<<maxidentity
                <<','<<positive<<','<<negative<<','<<backtracks<<','<<rejections<<'\n';
            std::cout<<method<<" h="<<h<<" beta=gamma="<<compliance<<" status="<<status
                <<" Efinal="<<final.E()<<" maxP="<<maxp<<" maxL="<<maxl<<" S="<<final.S<<'\n';
        }
        trajectory.close();summary.close();
    } catch(const std::exception& e) { std::cerr<<e.what()<<'\n';return 1; }
}
