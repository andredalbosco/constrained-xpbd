#include "dual_gs.hpp"
#include "oracle.hpp"
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>

using compose::Case;
using Clock=std::chrono::steady_clock;
static double elapsed(Clock::time_point start){return std::chrono::duration<double,std::milli>(Clock::now()-start).count();}
static double full_residual(const dual_gs::Residuals&r){return std::max({r.material,r.kinetic,r.storage,r.geometry,r.tangency,r.stationarity});}
static double state_gap(const compose::Model&m,const compose::State&a,const compose::State&b,double E){
    return std::sqrt((compose::masses(m).array()*(a.v-b.v).array().square()).sum()/E+(a.x-b.x).squaredNorm()/m.count);
}

struct Selection {
    dual_gs::Result result;
    int selected_prefix=0,largest_successful_prefix=0,next_prefix=-1,search_runs=0;
    long long next_evaluations=-1;
    double search_ms=0;
    std::string boundary,next_failure;
};

// Each candidate starts from the identical f0 minimum. Counts are monotone along
// deterministic prefixes until the first guard failure. Never expand past it.
static Selection select_prefix(const dual_gs::Problem&p,const std::string&method,
                               const std::string&schedule,int ceiling){
    Selection selected;const int cap=2*ceiling;std::map<int,dual_gs::Result> cache;
    const auto start=Clock::now();
    auto run=[&](int visits)->const dual_gs::Result&{
        auto found=cache.find(visits);if(found!=cache.end())return found->second;
        ++selected.search_runs;
        return cache.emplace(visits,dual_gs::solve(p,method,schedule,cap+1,1,visits)).first->second;
    };
    auto admissible=[&](const dual_gs::Result&r){return r.failure.empty()&&r.row_evaluations<=ceiling;};
    int low=0,high=1;run(0);
    while(admissible(run(high))){low=high;if(high==cap)break;high=std::min(cap,high*2);}
    if(low<high)while(high-low>1){const int mid=low+(high-low)/2;if(admissible(run(mid)))low=mid;else high=mid;}
    selected.largest_successful_prefix=low;selected.selected_prefix=low;selected.result=run(low);
    if(low==cap){selected.boundary="visit_cap";}
    else{
        const auto& next=run(low+1);selected.next_prefix=low+1;selected.next_evaluations=next.row_evaluations;selected.next_failure=next.failure;
        if(!next.failure.empty()&&next.row_evaluations<=ceiling){
            selected.boundary="first_guard_failure";selected.selected_prefix=low+1;selected.result=next;
        }else{
            selected.boundary="evaluation_ceiling";
            if(next.row_evaluations<=ceiling)throw std::runtime_error("uncertified_prefix_boundary");
        }
    }
    if(selected.result.row_evaluations>ceiling)throw std::runtime_error("selected_prefix_exceeds_ceiling");
    // Independently check cached prefix count monotonicity up to the first
    // observed failure. Failed prefixes may saturate and are never expanded.
    long long prior=-1;
    for(const auto& pair:cache){if(pair.second.row_evaluations<prior)throw std::runtime_error("nonmonotone_prefix_counts");prior=pair.second.row_evaluations;if(!pair.second.failure.empty())break;}
    selected.search_ms=elapsed(start);return selected;
}

int main(int argc,char**argv){try{
    const std::string directory=argc>1?argv[1]:"results";std::filesystem::create_directories(directory);
    std::ofstream out(directory+"/fair_cost.csv"),oracles(directory+"/fair_cost_oracles.csv");
    if(!out||!oracles)throw std::runtime_error("cannot_open_cost_outputs");out<<std::setprecision(17);oracles<<std::setprecision(17);
    out<<"case,N,h,method,schedule,evaluation_ceiling,selected_prefix,largest_successful_prefix,actual_visits,row_evaluations,coordinate_writes,restoration_passes,completed_sweeps,next_prefix,next_prefix_evaluations,boundary,next_prefix_failure,setup_ms,solve_ms,search_runs,search_ms,material,kinetic,storage,geometry,tangency,stationarity,energy,full_residual,objective,oracle_converged,oracle_gap,failure\n";
    oracles<<"case,N,h,converged,primal,stationarity,objective,curvature,regular,local_minimum_certified,evaluations,solve_ms,reason\n";
    const std::vector<std::string> methods={"refreshed_drop","refreshed_local_rx","frozen_full_rx"};
    const std::vector<std::string> schedules={"material_storage_kinetic","kinetic_storage_material","storage_kinetic_material","edge","edge_storage_first","edge_reverse"};
    const std::vector<Case> cases={compose::make_case("chain_mixed"),hour_local::chain_case(3),hour_local::chain_case(8)};
    for(const auto&c:cases)for(double h:{.02,.08}){
        double E=std::max(1e-8,compose::energy(c.model,c.initial)-compose::minimum_energy(c.model));auto start=Clock::now();
        dual_gs::Problem p(c.model,c.initial,h,E);double setup_ms=elapsed(start);auto oracle=dual_gs::solve_oracle(c.model,c.initial,h,E);
        oracles<<c.name<<','<<c.model.count<<','<<h<<','<<oracle.converged<<','<<oracle.primal<<','<<oracle.stationarity<<','<<oracle.objective<<','<<oracle.curvature<<','<<oracle.regular<<','<<oracle.local_minimum_certified<<','<<oracle.evaluations<<','<<oracle.milliseconds<<','<<oracle.reason<<'\n';
        for(const auto&method:methods)for(const auto&schedule:schedules)for(int ceiling:{100,300,1000,3000}){
            auto s=select_prefix(p,method,schedule,ceiling);const auto&r=s.result;
            const double gap=oracle.converged?state_gap(c.model,p.state(r.z),oracle.state,E):std::numeric_limits<double>::quiet_NaN();
            double objective=std::numeric_limits<double>::quiet_NaN();try{objective=p.objective(r.z);}catch(const std::exception&){}
            out<<c.name<<','<<c.model.count<<','<<h<<','<<method<<','<<schedule<<','<<ceiling<<','<<s.selected_prefix<<','<<s.largest_successful_prefix<<','<<r.visits<<','<<r.row_evaluations<<','<<r.coordinate_writes<<','<<r.stationarity_passes<<','<<r.completed_sweeps<<','<<s.next_prefix<<','<<s.next_evaluations<<','<<s.boundary<<','<<s.next_failure<<','<<setup_ms<<','<<r.milliseconds<<','<<s.search_runs<<','<<s.search_ms<<','<<r.material<<','<<r.kinetic<<','<<r.storage<<','<<r.geometry<<','<<r.tangency<<','<<r.stationarity<<','<<r.energy<<','<<full_residual(r)<<','<<objective<<','<<oracle.converged<<','<<gap<<','<<r.failure<<'\n';
        }
        std::cerr<<"fair_cost "<<c.name<<" h="<<h<<" oracle="<<oracle.converged<<'\n';out.flush();oracles.flush();
    }
    return 0;
}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}
