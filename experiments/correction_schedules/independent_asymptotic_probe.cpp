#define main parent_study_main
#include "study.cpp"
#undef main
int main() {
    std::cout<<std::setprecision(17);
    std::cout<<"case,method,substeps,h,rms_position_error,max_energy_percent,max_torque_error,final_energy,path_ratio,cap_count,max_extra_ratio,status\n";
    for(const auto& c:cases()) {
        if(c.name!="spring_orbit" && c.name!="coupled_motion")continue;
        for(int substeps:{32,64,128,256}) {
            const double h=1./(60*substeps); const int steps=600*substeps;
            const auto ref=reference(c,h,steps,8);
            for(const auto& cfg:configurations()) {
                if(cfg.name!="baseline" && cfg.name!="LE_repair_cap")continue;
                const Stats s=run(c,cfg,substeps,10,ref,nullptr);
                std::cout<<c.name<<','<<cfg.name<<','<<substeps<<','<<h<<','<<s.rms_position<<','<<100*s.energy_error/scale(c)<<','<<s.max_torque_error<<','<<s.final_energy<<','<<s.path/s.ref_path<<','<<s.diagnostic.capped<<','<<s.diagnostic.max_accepted_ratio<<','<<s.status<<'\n';
            }
        }
    }
}
