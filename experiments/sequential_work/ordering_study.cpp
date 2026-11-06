#include "scenes.hpp"
#include <Eigen/Geometry>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

using namespace cxpbd;

namespace {

void require(bool value,const char* message) {
    if(!value) throw std::runtime_error(message);
}

Model rotating_pair() {
    Model model;
    Particle a,b;
    a.mass=2;b.mass=1;
    const double total=a.mass+b.mass;
    const Vec3 center(.3,.4,.2),relative(1.25,0,0),translation(.25,-.1,.15),velocity(0,1,.2);
    a.x=center+(b.mass/total)*relative;
    b.x=center-(a.mass/total)*relative;
    a.v=translation+(b.mass/total)*velocity;
    b.v=translation-(a.mass/total)*velocity;
    model.particles={a,b};
    add_spring(model,0,1,.01);
    model.springs[0].rest_length=1;
    return model;
}

Model spring_triangle() {
    Model model;
    Particle a,b,c;
    a.mass=1;b.mass=2;c.mass=.7;
    a.x=Vec3(-.7,0,0);b.x=Vec3(.8,0,0);c.x=Vec3(.15,1.1,.25);
    model.particles={a,b,c};
    add_spring(model,0,1,.01);
    add_spring(model,1,2,.02);
    add_spring(model,2,0,.015);
    model.particles[1].x+=Vec3(.15,-.05,.1);
    Vec3 center=Vec3::Zero();
    double mass=0;
    for(const auto& p:model.particles) { center+=p.mass*p.x;mass+=p.mass; }
    center/=mass;
    const Vec3 omega(.2,-.3,.9),translation(.2,-.1,.1);
    for(auto& p:model.particles) p.v=translation+omega.cross(p.x-center);
    model.particles[2].v+=Vec3(-.15,.05,.1);
    return model;
}

Model make_test_scene(const std::string& scene) {
    if(scene=="spring") return make_spring();
    if(scene=="pendulum") return make_pendulum();
    if(scene=="free_pair") return rotating_pair();
    if(scene=="free_triangle") return spring_triangle();
    throw std::runtime_error("unknown ordering-study scene");
}

double independent_potential(const Model& model,const Spring& spring) {
    if(spring.compliance==0) return 0;
    const double c=(model.particles[spring.a].x-model.particles[spring.b].x).norm()-spring.rest_length;
    return c*c/(2*spring.compliance);
}
double independent_energy(const Model& model) {
    double value=0;
    for(const auto& p:model.particles) value+=.5*p.mass*p.v.squaredNorm()-p.mass*model.gravity.dot(p.x);
    for(const auto& s:model.springs) value+=independent_potential(model,s);
    return value;
}
double particle_reservoir(const Model& model) {
    double value=0;
    for(const auto& p:model.particles) value+=p.reservoir;
    return value;
}
double spring_reservoir(const Model& model) {
    double value=0;
    for(const auto& s:model.springs) value+=s.reservoir;
    return value;
}
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

double local_account_residual(const Model& old,const Model& next) {
    std::vector<Vec3> force(old.particles.size(),Vec3::Zero());
    double residual=0;
    for(std::size_t j=0;j<old.springs.size();++j) {
        const auto& s=old.springs[j];
        const Vec3 edge=old.particles[s.a].x-old.particles[s.b].x;
        const double length=edge.norm();
        Vec3 f=Vec3::Zero();
        if(length>1e-14*s.rest_length) {
            const double magnitude=old.previous_dt>0
                ?s.previous_lambda/(old.previous_dt*old.previous_dt)
                :(s.compliance>0?-(length-s.rest_length)/s.compliance:0);
            f=magnitude*edge/length;
        }
        force[s.a]+=f;force[s.b]-=f;
        const Vec3 dx=(next.particles[s.a].x-old.particles[s.a].x)
                    -(next.particles[s.b].x-old.particles[s.b].x);
        const double d=independent_potential(next,next.springs[j])-independent_potential(old,s)+f.dot(dx);
        residual=std::max(residual,std::abs(d+next.springs[j].reservoir-s.reservoir));
    }
    for(std::size_t i=0;i<old.particles.size();++i) {
        const auto& p=old.particles[i];
        const auto& q=next.particles[i];
        const double b=.5*p.mass*(q.v.squaredNorm()-p.v.squaredNorm())
            -(force[i]+p.mass*old.gravity).dot(q.x-p.x);
        residual=std::max(residual,std::abs(b+q.reservoir-p.reservoir));
    }
    return residual;
}

void require_finite(const Model& model,const Diagnostics& d) {
    for(const auto& p:model.particles)
        require(p.x.allFinite()&&p.v.allFinite()&&std::isfinite(p.reservoir),"nonfinite particle");
    for(const auto& s:model.springs)
        require(std::isfinite(s.lambda)&&std::isfinite(s.previous_lambda)&&std::isfinite(s.reservoir),"nonfinite spring");
    const double values[]={d.mechanical_energy,d.reservoir_energy,d.reservoir_norm,
        d.energy_balance_error,d.max_balance_residual,d.max_feedback_residual,
        d.max_material_residual,d.position_residual_norm,independent_energy(model),
        particle_reservoir(model),spring_reservoir(model)};
    for(double value:values) require(std::isfinite(value),"nonfinite diagnostic");
}

Settings settings_for(const std::string& method,double h,int iterations,double beta,double gamma) {
    Settings settings;
    settings.dt=h;settings.iterations=iterations;
    settings.kinetic_reservoir_compliance=beta;
    settings.elastic_reservoir_compliance=gamma;
    settings.mode=method=="xpbd"?Mode::XPBD:Mode::Active;
    if(method=="paired") settings.active_order=ActiveOrder::Paired;
    if(method=="paired_local") settings.active_order=ActiveOrder::PairedLocal;
    if(method=="sequential") settings.active_order=ActiveOrder::Sequential;
    return settings;
}

void baseline_checks() {
    Model spring=make_spring();
    Settings settings=settings_for("xpbd",1./60,1,1,1);
    double q=.25,v=0;
    for(int i=0;i<600;++i) {
        const double next_q=(q+settings.dt*v)/(1+100*settings.dt*settings.dt);
        const double next_v=(next_q-q)/settings.dt;
        const auto d=step(spring,settings);
        require_finite(spring,d);
        require(std::abs(spring.particles[1].x.x()-(1+next_q))<1e-12,"XPBD analytical position mismatch");
        require(std::abs(spring.particles[1].v.x()-next_v)<1e-10,"XPBD analytical velocity mismatch");
        require(spring.particles[0].x.squaredNorm()==0,"stationary spring pin moved");
        q=next_q;v=next_v;
    }
    for(const std::string scene:{"free_pair","free_triangle"}) {
        Model model=make_test_scene(scene);
        const Vec3 p0=momentum(model);
        for(int i=0;i<240;++i) {
            const auto d=step(model,settings);
            require_finite(model,d);
            require((momentum(model)-p0).norm()<1e-9,"ordinary XPBD linear momentum mismatch");
        }
    }
    std::cout<<"Mandatory checks passed: scalar XPBD recurrence, stationary pin, free-system linear momentum.\n";
}

struct Summary {
    double max_energy_error=0,max_extended_error=0,max_step_ledger=0,max_local_ledger=0;
    double max_material=0,max_feedback=0,max_position_residual=0,max_reservoir_norm=0;
    double max_momentum=0,max_angular=0,max_hard_error=0,max_stretch=0,max_pin_error=0;
    double max_x=0,max_speed=0,max_motion=0;
    long long backtracks=0,rejected=0;
    int completed=0;
    bool finite_failure=false;
    std::string status="completed";
};

void run(const std::string& scene,const std::string& method,double h,int iterations,
         double beta,double gamma,double duration,std::ofstream& output) {
    Model model=make_test_scene(scene);
    const Model initial=model;
    const Settings settings=settings_for(method,h,iterations,beta,gamma);
    const bool isolated=scene=="free_pair"||scene=="free_triangle";
    const double e0=independent_energy(model);
    const double augmented0=e0+particle_reservoir(model)+spring_reservoir(model);
    const Vec3 p0=momentum(model),l0=angular_momentum(model);
    Summary summary;
    try {
        for(int step_index=0;step_index<static_cast<int>(std::lround(duration/h));++step_index) {
            // Commit only an accepted finite call; errors retain the last accepted snapshot.
            Model candidate=model;
            const Diagnostics d=step(candidate,settings);
            require_finite(candidate,d);
            const double energy=independent_energy(candidate);
            const double augmented=energy+particle_reservoir(candidate)+spring_reservoir(candidate);
            const double local=local_account_residual(model,candidate);
            require(std::isfinite(local),"nonfinite independent ledger");
            summary.max_energy_error=std::max(summary.max_energy_error,std::abs(energy-e0));
            summary.max_extended_error=std::max(summary.max_extended_error,std::abs(augmented-augmented0));
            summary.max_step_ledger=std::max(summary.max_step_ledger,std::abs(d.energy_balance_error));
            summary.max_local_ledger=std::max(summary.max_local_ledger,local);
            summary.max_material=std::max(summary.max_material,d.max_material_residual);
            summary.max_feedback=std::max(summary.max_feedback,d.max_feedback_residual);
            summary.max_position_residual=std::max(summary.max_position_residual,d.position_residual_norm);
            summary.max_reservoir_norm=std::max(summary.max_reservoir_norm,d.reservoir_norm);
            summary.max_momentum=std::max(summary.max_momentum,(momentum(candidate)-p0).norm());
            summary.max_angular=std::max(summary.max_angular,(angular_momentum(candidate)-l0).norm());
            summary.max_stretch=std::max(summary.max_stretch,maximum_relative_stretch(candidate));
            for(const auto& s:candidate.springs)
                if(s.compliance==0)
                    summary.max_hard_error=std::max(summary.max_hard_error,
                        std::abs((candidate.particles[s.a].x-candidate.particles[s.b].x).norm()-s.rest_length));
            for(std::size_t i=0;i<candidate.particles.size();++i) {
                const auto& p=candidate.particles[i];
                summary.max_x=std::max(summary.max_x,p.x.norm());
                summary.max_speed=std::max(summary.max_speed,p.v.norm());
                const double motion=(p.x-initial.particles[i].x).norm();
                summary.max_motion=std::max(summary.max_motion,motion);
                if(p.pinned) summary.max_pin_error=std::max(summary.max_pin_error,motion);
            }
            summary.backtracks+=d.backtracks;summary.rejected+=d.rejected_updates;
            model=std::move(candidate);++summary.completed;
            if(summary.max_x>1e50||summary.max_speed>1e50) {
                summary.status="magnitude_guard";break;
            }
        }
    } catch(const std::exception& e) {
        summary.status=e.what();
        summary.finite_failure=summary.status.find("nonfinite")!=std::string::npos
                             ||summary.status.find("non-finite")!=std::string::npos;
    }
    if(summary.status=="completed"&&(summary.max_x>10||summary.max_speed>100))
        summary.status="large_accepted_excursion";
    else if(summary.status=="completed"&&summary.max_motion<1e-9)
        summary.status="frozen";
    output<<scene<<','<<method<<','<<h<<','<<iterations<<','<<beta<<','<<gamma<<','<<isolated<<','
          <<summary.completed<<','<<summary.completed*h<<','<<summary.status<<','<<summary.finite_failure<<','
          <<e0<<','<<independent_energy(model)-e0<<','<<summary.max_energy_error<<','
          <<summary.max_energy_error/std::max(1.,std::abs(e0))<<','
          <<particle_reservoir(model)<<','<<spring_reservoir(model)<<','<<summary.max_reservoir_norm<<','
          <<summary.max_extended_error<<','<<summary.max_step_ledger<<','<<summary.max_local_ledger<<','
          <<summary.max_material<<','<<summary.max_feedback<<','<<summary.max_position_residual<<','
          <<summary.max_momentum<<','<<summary.max_angular<<','<<summary.max_hard_error<<','
          <<summary.max_stretch<<','<<summary.max_pin_error<<','<<summary.max_x<<','<<summary.max_speed<<','
          <<summary.max_motion<<','<<summary.backtracks<<','<<summary.rejected<<'\n';
}

} // namespace

int main(int argc,char** argv) {
    try {
        const double duration=argc>1?std::stod(argv[1]):10;
        require(std::isfinite(duration)&&duration>0,"invalid duration");
        baseline_checks();
        std::ofstream output("ordering_results.csv");
        output.exceptions(std::ios::failbit|std::ios::badbit);
        output<<std::setprecision(17);
        output<<"scene,method,h,iterations,beta,gamma,isolated_system,completed_steps,duration,status,finite_failure,"
                "initial_energy,final_energy_change,max_energy_error,max_relative_energy_error,"
                "final_particle_reservoir,final_spring_reservoir,max_reservoir_norm,max_extended_energy_error,"
                "max_step_ledger_error,max_independent_local_ledger_error,max_material_residual,"
                "max_feedback_residual,max_position_residual,max_momentum_change,max_angular_momentum_change,"
                "max_hard_length_error,max_relative_stretch,max_pin_error,max_position,max_speed,max_motion,"
                "backtracks,rejected_updates\n";
        const std::vector<std::pair<double,double>> compliances={{1,1},{10,10},{1000,1000},
            {1,1000},{1000,1},{1e6,1e6}};
        for(const std::string scene:{"spring","pendulum","free_pair","free_triangle"})
            for(double h:{1./240,1./60,1./15})
                for(int iterations:{1,5,20}) {
                    run(scene,"xpbd",h,iterations,1000,1000,duration,output);
                    for(const std::string method:{"paired","paired_local","sequential"})
                        for(const auto& pair:compliances)
                            run(scene,method,h,iterations,pair.first,pair.second,duration,output);
                }
        output.close();
        std::cout<<"Wrote 684 ordering cases. Completion/ledger closure do not imply physical conservation.\n";
    } catch(const std::exception& e) {
        std::cerr<<e.what()<<'\n';return 1;
    }
}
