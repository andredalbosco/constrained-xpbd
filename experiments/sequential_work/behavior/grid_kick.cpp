#include "scenes.hpp"
#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>

using namespace cxpbd;

namespace {
void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}
double kinetic(const Model& model) {
    double value=0;
    for (const auto& p:model.particles) value+=.5*p.mass*p.v.squaredNorm();
    return value;
}
double elastic(const Model& model) {
    double value=0;
    for (const auto& s:model.springs) {
        const double c=(model.particles[s.a].x-model.particles[s.b].x).norm()-s.rest_length;
        if(s.compliance>0) value+=.5*c*c/s.compliance;
    }
    return value;
}
double gravity_potential(const Model& model) {
    double value=0;
    for (const auto& p:model.particles) value-=p.mass*model.gravity.dot(p.x);
    return value;
}
double particle_account(const Model& model) {
    double value=0;
    for (const auto& p:model.particles) value+=p.reservoir;
    return value;
}
double spring_account(const Model& model) {
    double value=0;
    for (const auto& s:model.springs) value+=s.reservoir;
    return value;
}
void check(const Model& model,const Model& initial,const Diagnostics& d,bool planar) {
    for(std::size_t i=0;i<model.particles.size();++i) {
        const auto& p=model.particles[i];
        require(p.x.allFinite()&&p.v.allFinite()&&std::isfinite(p.reservoir),"nonfinite particle state");
        if(planar) require(p.x.z()==0&&p.v.z()==0,"grid left its plane");
        if(p.pinned) require((p.x-initial.particles[i].x).squaredNorm()==0&&p.v.squaredNorm()==0,"pin moved");
    }
    for(const auto& s:model.springs)
        require(std::isfinite(s.lambda)&&std::isfinite(s.previous_lambda)&&std::isfinite(s.reservoir),"nonfinite spring state");
    for(double value:{d.mechanical_energy,d.reservoir_energy,d.reservoir_norm,d.energy_balance_error,
        d.max_balance_residual,d.max_feedback_residual,d.max_material_residual,d.position_residual_norm,
        kinetic(model),elastic(model),gravity_potential(model)})
        require(std::isfinite(value),"nonfinite diagnostic");
}
void snapshot(std::ofstream& out,const Model& model,int step_index,double time) {
    for(std::size_t i=0;i<model.particles.size();++i) {
        const auto& p=model.particles[i];
        out<<step_index<<','<<time<<','<<i<<','<<p.x.x()<<','<<p.x.y()<<','<<p.x.z()<<','
           <<p.v.x()<<','<<p.v.y()<<','<<p.v.z()<<','<<p.pinned<<','<<p.reservoir<<'\n';
    }
}
}

int main(int argc,char** argv) {
    try {
        require(argc==4||argc==5,"usage: grid_kick xpbd|sequential beta output_directory [cloth]");
        const std::string scene=argc==5?argv[4]:"grid_kick";
        require(scene=="grid_kick"||scene=="cloth","invalid scene");
        const std::string method=argv[1];
        require(method=="xpbd"||method=="sequential","invalid method");
        const double beta=std::stod(argv[2]);
        require(std::isfinite(beta)&&beta>0,"invalid beta");
        const std::filesystem::path directory=argv[3];
        std::filesystem::create_directories(directory);
        Model model=scene=="grid_kick"?make_grid():make_cloth();
        const double pi=std::acos(-1.);
        for(std::size_t i=0;i<model.particles.size();++i) {
            auto& p=model.particles[i];
            if(scene=="grid_kick"&&!p.pinned) p.v.y()=-2*std::sin(pi*static_cast<double>(i%6)/5)*std::sin(pi*static_cast<double>(i/6)/5);
        }
        const Model initial=model;
        const double initial_energy=kinetic(model)+elastic(model)+gravity_potential(model);
        const double initial_kinetic=kinetic(model);
        require(std::abs(initial_kinetic-(scene=="grid_kick"?1.25:0))<1e-13,"initial kinetic energy check failed");
        require(std::abs(initial_energy-(scene=="grid_kick"?32.1515:59.59575))<1e-12,"initial total energy check failed");
        Settings settings;
        settings.dt=1./240;settings.iterations=20;
        settings.mode=method=="xpbd"?Mode::XPBD:Mode::Active;
        settings.active_order=ActiveOrder::Sequential;
        settings.kinetic_reservoir_compliance=beta;
        settings.elastic_reservoir_compliance=beta;
        std::ofstream trajectory(directory/"trajectory.csv"),snapshots(directory/"snapshots.csv"),summary(directory/"summary.csv");
        for(auto* stream:{&trajectory,&snapshots,&summary}) {
            stream->exceptions(std::ios::badbit|std::ios::failbit);
            *stream<<std::setprecision(17);
        }
        trajectory<<"step,time,mechanical_energy,kinetic_energy,elastic_energy,gravity_potential,particle_reservoir,spring_reservoir,reservoir_norm,extended_energy_drift,relative_stretch,material_residual,feedback_residual,position_residual,middle_x,middle_y,middle_vy,bottom_x,bottom_y,bottom_vy,max_speed,backtracks,rejected_updates\n";
        snapshots<<"step,time,particle,x,y,z,vx,vy,vz,pinned,reservoir\n";
        Diagnostics d;
        double max_energy_error=0,max_kinetic=kinetic(model),max_stretch=0,max_extended=0,max_material=0,max_feedback=0,max_position=0,max_speed=0;
        double tail_kinetic_sum=0,tail_kinetic_max=0,tail_middle_min=1e100,tail_middle_max=-1e100;
        double tail_bottom_min=1e100,tail_bottom_max=-1e100;
        long long backtracks=0,rejected=0;int tail_count=0,completed=0;
        std::string status="completed";
        try {
            for(int i=0;i<=2400;++i) {
                if(i>0) { d=step(model,settings);check(model,initial,d,scene=="grid_kick");++completed; }
                const double t=i*settings.dt,k=kinetic(model),u=elastic(model),vg=gravity_potential(model),e=k+u+vg;
                const double r=particle_account(model),s=spring_account(model),extended=e+r+s-initial_energy;
                const double stretch=maximum_relative_stretch(model);
                double speed=0;for(const auto& p:model.particles) speed=std::max(speed,p.v.norm());
                max_energy_error=std::max(max_energy_error,std::abs(e-initial_energy));
                max_kinetic=std::max(max_kinetic,k);max_stretch=std::max(max_stretch,stretch);
                max_extended=std::max(max_extended,std::abs(extended));max_speed=std::max(max_speed,speed);
                max_material=std::max(max_material,d.max_material_residual);
                max_feedback=std::max(max_feedback,d.max_feedback_residual);
                max_position=std::max(max_position,d.position_residual_norm);
                backtracks+=d.backtracks;rejected+=d.rejected_updates;
                const auto& middle=model.particles[scene=="grid_kick"?20:40];
                const auto& bottom=model.particles[scene=="grid_kick"?32:76];
                if(t>=8) {
                    tail_kinetic_sum+=k;tail_kinetic_max=std::max(tail_kinetic_max,k);++tail_count;
                    tail_middle_min=std::min(tail_middle_min,middle.x.y());tail_middle_max=std::max(tail_middle_max,middle.x.y());
                    tail_bottom_min=std::min(tail_bottom_min,bottom.x.y());tail_bottom_max=std::max(tail_bottom_max,bottom.x.y());
                }
                trajectory<<i<<','<<t<<','<<e<<','<<k<<','<<u<<','<<vg<<','<<r<<','<<s<<','<<d.reservoir_norm<<','
                  <<extended<<','<<stretch<<','<<d.max_material_residual<<','<<d.max_feedback_residual<<','<<d.position_residual_norm<<','
                  <<middle.x.x()<<','<<middle.x.y()<<','<<middle.v.y()<<','<<bottom.x.x()<<','<<bottom.x.y()<<','<<bottom.v.y()<<','
                  <<speed<<','<<d.backtracks<<','<<d.rejected_updates<<'\n';
                if(i%60==0) snapshot(snapshots,model,i,t);
            }
        } catch(const std::exception& error) { status=error.what(); }
        summary<<"scene,order,dt,beta,gamma,iterations,status,completed_steps,duration,initial_energy,initial_kinetic,final_energy,final_energy_change,final_kinetic_energy,max_kinetic_energy,max_energy_error,max_relative_stretch,max_extended_energy_drift,max_material_residual,max_feedback_residual,max_position_residual,final_particle_reservoir,final_spring_reservoir,final_reservoir_norm,max_speed,mean_kinetic_last_2s,max_kinetic_last_2s,middle_y_range_last_2s,bottom_y_range_last_2s,backtracks,rejected_updates,pins_checked_each_step,plane_checked_each_step,finite_checked_each_step\n";
        const double final_energy=kinetic(model)+elastic(model)+gravity_potential(model);
        summary<<scene<<','<<method<<','<<settings.dt<<','<<beta<<','<<beta<<",20,"<<status<<','<<completed<<','<<completed*settings.dt<<','
          <<initial_energy<<','<<initial_kinetic<<','<<final_energy<<','<<final_energy-initial_energy<<','<<kinetic(model)<<','<<max_kinetic<<','<<max_energy_error<<','
          <<max_stretch<<','<<max_extended<<','<<max_material<<','<<max_feedback<<','<<max_position<<','<<particle_account(model)<<','<<spring_account(model)<<','
          <<d.reservoir_norm<<','<<max_speed<<','<<(tail_count?tail_kinetic_sum/tail_count:0)<<','<<tail_kinetic_max<<','
          <<(tail_count?tail_middle_max-tail_middle_min:0)<<','<<(tail_count?tail_bottom_max-tail_bottom_min:0)<<','<<backtracks<<','<<rejected<<",1,"<<(scene=="grid_kick")<<",1\n";
        trajectory.close();snapshots.close();summary.close();
        std::cout<<method<<" beta="<<beta<<" "<<status<<" final dE="<<final_energy-initial_energy<<" final K="<<kinetic(model)<<'\n';
        return status=="completed"?0:1;
    } catch(const std::exception& error) { std::cerr<<error.what()<<'\n';return 1; }
}
