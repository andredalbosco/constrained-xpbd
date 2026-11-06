#include "midpoint_gs.hpp"

Model two_springs() {
    Model m;
    Particle a, b, c;
    a.x = Vec3(-1.0,0.0,0.0); a.pinned = true;
    b.x = Vec3(0.8,0.5,0.0); b.pinned = true;
    c.x = Vec3(0.0,0.8,0.0); c.v = Vec3(0.7,-0.2,0.0);
    m.particles = {a,b,c};
    add_spring(m,0,2,0.001); add_spring(m,1,2,0.001);
    m.springs[0].rest_length *= 0.9; m.springs[1].rest_length *= 0.9;
    return m;
}

int main() {
    try {
        std::cout << std::setprecision(15)
            << "scene,method,h,outer_limit,inner_sweeps,steps,final_energy_change,max_energy_error,"
               "max_scaled_rx,max_constitutive,max_kinetic_work_defect,max_identity_error,"
               "max_constraint_error,max_endpoint_radial_velocity,unconverged_steps,max_outer_used,backtracks\n";
        for (const std::string scene : {"spring","pendulum","two_springs"}) {
            for (const std::string method : {"drop_rx","retain_rx","retain_rx_symmetric"}) {
                for (double h : {1.0/240.0,1.0/60.0,1.0/15.0}) {
                  for (int linear_sweeps : {1,10}) {
                    for (int limit : {1,5,20,80}) {
                        Model model = scene == "two_springs" ? two_springs() : make_scene(scene);
                        const double initial = mechanical_energy(model);
                        const int steps = static_cast<int>(std::lround(10.0/h));
                        double max_e=0,max_rx=0,max_c=0,max_w=0,max_identity=0,max_length=0,max_radial=0;
                        int incomplete=0,used=0,backtracks=0;
                        for (int k=0;k<steps;++k) {
                            const StepReport r=advance(model,h,limit,linear_sweeps,method!="drop_rx",method=="retain_rx_symmetric");
                            const double energy_error=std::abs(mechanical_energy(model)-initial);
                            if (!std::isfinite(energy_error)) throw std::runtime_error("nonfinite GS state");
                            max_e=std::max(max_e,energy_error); max_rx=std::max(max_rx,r.rx);
                            max_c=std::max(max_c,r.constitutive); max_w=std::max(max_w,std::abs(r.kinetic_work_defect));
                            max_identity=std::max(max_identity,r.identity_error);
                            incomplete += (r.rx>1e-11 || r.constitutive>1e-11);
                            used=std::max(used,r.iterations); backtracks+=r.backtracks;
                            if(scene=="pendulum") {
                                const Vec3 p=model.particles[1].x-model.particles[0].x;
                                max_length=std::max(max_length,std::abs(p.norm()-1.0));
                                max_radial=std::max(max_radial,std::abs(p.normalized().dot(model.particles[1].v)));
                            }
                        }
                        if (max_identity>1e-8) throw std::runtime_error("residual-to-kinetic-work identity failed");
                        if(method!="drop_rx" && limit==80 && h<0.02 && (max_e>1e-6 || incomplete>0))
                            throw std::runtime_error("well-resolved curvature-dropped GS did not conserve energy");
                        std::cout<<scene<<','<<method<<','<<h<<','<<limit<<','<<linear_sweeps<<','<<steps<<','
                            <<mechanical_energy(model)-initial<<','<<max_e<<','<<max_rx<<','<<max_c<<','<<max_w<<','
                            <<max_identity<<','<<max_length<<','<<max_radial<<','<<incomplete<<','<<used<<','<<backtracks<<'\n';
                    }
                  }
                }
            }
        }
        return 0;
    } catch(const std::exception& e) { std::cerr<<e.what()<<'\n'; return 1; }
}
