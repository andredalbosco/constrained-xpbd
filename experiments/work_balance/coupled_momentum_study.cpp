// Reuses the curvature-dropped midpoint solver; no second solver implementation.
#include "midpoint_gs.hpp"
#include <fstream>
#include <iomanip>
#include <iostream>

Vec3 momentum(const Model& model) {
    Vec3 value=Vec3::Zero();
    for(const auto& p:model.particles) value+=p.mass*p.v;
    return value;
}

Vec3 angular_momentum(const Model& model) {
    Vec3 value=Vec3::Zero();
    for(const auto& p:model.particles) value+=p.x.cross(p.mass*p.v);
    return value;
}

Model free_triangle(bool deforming) {
    Model model;
    Particle a,b,c;
    a.mass=1.0; b.mass=2.0; c.mass=1.5;
    a.x=Vec3(-.7,-.3,0); b.x=Vec3(.9,-.1,.1); c.x=Vec3(.1,1,.2);
    model.particles={a,b,c};
    Vec3 center=Vec3::Zero();
    double mass=0;
    for(const auto& p:model.particles) { center+=p.mass*p.x; mass+=p.mass; }
    center/=mass;
    const Vec3 translation(.4,-.2,.3),rotation(.3,-.5,1.1);
    for(auto& p:model.particles) p.v=translation+rotation.cross(p.x-center);
    if(deforming) {
        model.particles[0].v+=Vec3(.2,.05,-.1);
        model.particles[1].v-=Vec3(.1,.025,-.05);
    }
    add_spring(model,0,1,.01);
    add_spring(model,1,2,.015);
    add_spring(model,2,0,.008);
    model.springs[0].rest_length*=.9;
    model.springs[1].rest_length*=1.05;
    model.springs[2].rest_length*=.95;
    return model;
}

int main() {
    try {
        std::ofstream output("coupled_momentum_summary.csv");
        if(!output) throw std::runtime_error("cannot open coupled momentum output");
        output<<std::setprecision(17)
            <<"initial_motion,method,h,outer_limit,inner_sweeps,final_energy_change,max_energy_error,"
              "final_P_error,max_P_error,final_L_error,max_L_error,max_step_P_error,max_step_L_error,"
              "max_scaled_rx,max_constitutive,max_kinetic_work_defect,unconverged_steps,backtracks\n";
        for(bool deforming:{false,true}) for(const std::string method:{"drop_rx","retain_rx","retain_rx_symmetric"})
            for(double h:{1.0/240,1.0/60,1.0/15}) for(int limit:{1,5,20,80}) for(int inner:{1,10}) {
                Model model=free_triangle(deforming);
                const double initial_energy=mechanical_energy(model);
                const Vec3 initial_p=momentum(model),initial_l=angular_momentum(model);
                double max_e=0,max_p=0,max_l=0,max_step_p=0,max_step_l=0,max_rx=0,max_c=0,max_work=0;
                int incomplete=0,backtracks=0;
                const int steps=static_cast<int>(std::lround(10/h));
                for(int frame=0;frame<steps;++frame) {
                    const Model before=model;
                    const StepReport report=advance(model,h,limit,inner,method!="drop_rx",method=="retain_rx_symmetric");
                    const double ep=momentum(model).operator-(momentum(before)).norm();
                    const double el=angular_momentum(model).operator-(angular_momentum(before)).norm();
                    max_e=std::max(max_e,std::abs(mechanical_energy(model)-initial_energy));
                    max_p=std::max(max_p,(momentum(model)-initial_p).norm());
                    max_l=std::max(max_l,(angular_momentum(model)-initial_l).norm());
                    max_step_p=std::max(max_step_p,ep); max_step_l=std::max(max_step_l,el);
                    max_rx=std::max(max_rx,report.rx); max_c=std::max(max_c,report.constitutive);
                    max_work=std::max(max_work,std::abs(report.kinetic_work_defect));
                    incomplete+=(report.rx>1e-11 || report.constitutive>1e-11);
                    backtracks+=report.backtracks;
                    double linear_factor=0,angular_factor=0;
                    for(std::size_t i=0;i<model.particles.size();++i) {
                        const auto& p=model.particles[i];
                        const Vec3 middle=.5*(p.x+before.particles[i].x);
                        linear_factor+=p.mass*p.mass;
                        angular_factor+=p.mass*p.mass*middle.squaredNorm();
                        if(!p.x.allFinite() || !p.v.allFinite()) throw std::runtime_error("nonfinite coupled state");
                    }
                    // With central midpoint discrete gradients, momentum defects are
                    // 4/h times the net force/torque of rx. report.rx is ||M^-1 rx||.
                    const double p_bound=(4/h)*std::sqrt(linear_factor)*report.rx;
                    const double l_bound=(4/h)*std::sqrt(angular_factor)*report.rx;
                    if(ep>p_bound+2e-10 || el>l_bound+2e-10)
                        throw std::runtime_error("momentum defect exceeded position-residual bound");
                }
                if(max_p>1e-7) throw std::runtime_error("pair corrections lost linear momentum");
                if(method!="drop_rx" && limit==80 && inner==10 && h<.02 &&
                   (max_e>1e-6 || max_l>1e-6 || incomplete!=0))
                    throw std::runtime_error("resolved coupled run failed E/L convergence");
                output<<(deforming?"rigid_plus_deformation":"rigid")<<','<<method<<','<<h<<','<<limit<<','<<inner<<','
                    <<mechanical_energy(model)-initial_energy<<','<<max_e<<','
                    <<(momentum(model)-initial_p).norm()<<','<<max_p<<','
                    <<(angular_momentum(model)-initial_l).norm()<<','<<max_l<<','<<max_step_p<<','<<max_step_l<<','
                    <<max_rx<<','<<max_c<<','<<max_work<<','<<incomplete<<','<<backtracks<<'\n';
                if(limit==80 && inner==10) std::cout<<(deforming?"deforming":"rigid")<<' '<<method<<" h="<<h
                    <<" maxE="<<max_e<<" maxP="<<max_p<<" maxL="<<max_l<<" rx="<<max_rx
                    <<" incomplete="<<incomplete<<'\n';
            }
        output.flush();
        if(!output) throw std::runtime_error("failed writing coupled momentum output");
    } catch(const std::exception& error) {
        std::cerr<<error.what()<<'\n'; return 1;
    }
}
