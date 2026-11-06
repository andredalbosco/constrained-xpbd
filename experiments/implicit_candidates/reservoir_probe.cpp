#include "solver.hpp"
#include "../research_agenda/cases.hpp"
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
using namespace implicit;
int main(int argc,char** argv){try{const std::filesystem::path dir=argc>1?argv[1]:"results";std::filesystem::create_directories(dir);std::ofstream out(dir/"reservoir_probe.csv");out<<std::setprecision(17)<<"case,method,hz,duration,completed,steps,max_local_L1,first_10s_max,last_10s_max,max_cumulative_ledger_error\n";
    for(const auto& test:cases()){
    for(const auto& cfg:configurations())if(cfg.name=="workL_s1_tau005"||cfg.name=="workL_s10_tau005"||cfg.name=="workL_soft_01"){
        State state=test.initial;Reservoir reservoir;const double E0=energy(test.model,state),scale=E0-minimum_energy(test.model);double largest=0,first=0,last=0,ledger_error=0;int completed=0;constexpr int hz=480,steps=60*hz;
        for(int k=1;k<=steps;++k){State next;StepReport r;if(!advance(context(test.model,state,1./hz,E0),cfg,scale,reservoir,next,r))break;
            double local=0;for(double v:reservoir.k)local+=std::abs(v);for(double v:reservoir.u)local+=std::abs(v);local/=scale;largest=std::max(largest,local);if(k<=10*hz)first=std::max(first,local);if(k>50*hz)last=std::max(last,local);
            ledger_error=std::max(ledger_error,std::abs(energy(test.model,next)+total(reservoir)-E0)/scale);state=next;completed=k;
        }
        out<<test.name<<','<<cfg.name<<",480,60,"<<completed<<','<<steps<<','<<largest<<','<<first<<','<<last<<','<<ledger_error<<'\n';
    }
    }
    std::cout<<"Local reservoir and cumulative ledger probes complete\n";
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
