#include "scenes.hpp"
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

using namespace cxpbd;
struct Case { Mode mode; ActiveOrder order; double compliance; std::string name; };

int main(int argc,char** argv) {
    try {
        if(argc!=2)throw std::runtime_error("usage: trajectory_probe output_directory");
        const std::filesystem::path output=argv[1];std::filesystem::create_directories(output);
        std::vector<Case> cases={{Mode::XPBD,ActiveOrder::Paired,1000,"xpbd_paired_b1000"}};
        for(double c:{1.,10.,100.,1000.,1e6}) {
            const std::string label=c==1e6?"1e+06":std::to_string(static_cast<int>(c));
            cases.push_back({Mode::Active,ActiveOrder::Sequential,c,"active_sequential_b"+label});
        }
        for(double c:{1.,1000.})
            cases.push_back({Mode::Active,ActiveOrder::PairedLocal,c,"active_paired-local_b"+std::to_string(static_cast<int>(c))});
        double max_analytic_x_error=0,max_analytic_v_error=0;
        for(const std::string scene:{"spring","pendulum"})for(int hz:{240,60})for(const auto& test:cases) {
            Settings s;s.dt=1./hz;s.iterations=20;s.mode=test.mode;s.active_order=test.order;
            s.kinetic_reservoir_compliance=s.elastic_reservoir_compliance=test.compliance;
            Model model=make_scene(scene);const double initial_energy=mechanical_energy(model);
            const std::string name=scene+"_h1-"+std::to_string(hz)+"_"+test.name;
            std::ofstream csv(output/(name+".csv"));if(!csv)throw std::runtime_error("cannot open "+name);
            csv<<std::setprecision(17)
               <<"step,t,x,y,z,vx,vy,vz,C,angle,speed,radial_velocity,kinetic,elastic,gravity,mechanical,energy_above_minimum,particle_reservoir,spring_reservoir,feedback_residual,material_residual,position_residual,extended_drift,backtracks,rejected_updates,analytic_be_x,analytic_be_v,continuous_x,continuous_v\n";
            Diagnostics d;
            for(int k=0;k<=10*hz;++k) {
                if(k)d=step(model,s);
                const auto& p=model.particles.at(1);const auto& spring=model.springs.at(0);
                double kinetic=0,elastic=0,gravity=0,particle_r=0,spring_r=0;
                for(const auto& a:model.particles){kinetic+=.5*a.mass*a.v.squaredNorm();gravity-=a.mass*model.gravity.dot(a.x);particle_r+=a.reservoir;}
                for(const auto& a:model.springs){elastic+=spring_energy(model,a);spring_r+=a.reservoir;}
                const double mechanical=kinetic+elastic+gravity,t=k*s.dt;
                const Vec3 relative=p.x-model.particles[0].x;
                const double C=relative.norm()-spring.rest_length;
                const double angle=std::atan2(relative.x(),-relative.y());
                const double be_amplitude=.25*std::pow(1+100*s.dt*s.dt,-.5*k);
                const double phase=k*std::atan(10*s.dt);
                const double bx=1+be_amplitude*std::cos(phase),bv=-10*be_amplitude*std::sin(phase);
                if(scene=="spring"&&test.mode==Mode::XPBD) {
                    max_analytic_x_error=std::max(max_analytic_x_error,std::abs(p.x.x()-bx));
                    max_analytic_v_error=std::max(max_analytic_v_error,std::abs(p.v.x()-bv));
                }
                if(!p.x.allFinite()||!p.v.allFinite()||!std::isfinite(mechanical))throw std::runtime_error("nonfinite "+name);
                csv<<k<<','<<t<<','<<p.x.x()<<','<<p.x.y()<<','<<p.x.z()<<','<<p.v.x()<<','<<p.v.y()<<','<<p.v.z()<<','
                   <<C<<','<<angle<<','<<p.v.norm()<<','<<relative.dot(p.v)/relative.norm()<<','<<kinetic<<','<<elastic<<','<<gravity<<','
                   <<mechanical<<','<<(mechanical+(scene=="pendulum"?9.81:0))<<','<<particle_r<<','<<spring_r<<','
                   <<d.max_feedback_residual<<','<<d.max_material_residual<<','<<d.position_residual_norm<<','
                   <<mechanical+particle_r+spring_r-initial_energy<<','<<d.backtracks<<','<<d.rejected_updates<<','
                   <<(scene=="spring"?bx:0)<<','<<(scene=="spring"?bv:0)<<','
                   <<(scene=="spring"?1+.25*std::cos(10*t):0)<<','<<(scene=="spring"?-2.5*std::sin(10*t):0)<<'\n';
            }
            csv.flush();if(!csv)throw std::runtime_error("failed writing "+name);
            std::cout<<name<<" complete\n";
        }
        if(max_analytic_x_error>1e-11||max_analytic_v_error>1e-10)
            throw std::runtime_error("XPBD spring disagrees with analytic backward Euler");
        std::cout<<std::setprecision(17)<<"analytic_be_max_x_error="<<max_analytic_x_error
                 <<", analytic_be_max_v_error="<<max_analytic_v_error<<'\n';
    } catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}
