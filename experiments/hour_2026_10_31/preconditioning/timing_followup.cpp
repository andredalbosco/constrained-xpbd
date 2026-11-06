#include "solver.hpp"
#include <fstream>
#include <iomanip>
#include <iostream>
struct Variant{std::string label,metric,forcing="fixed",linear="transformed";int chunk=0;bool step_cache=false;};
int main(){try{
    std::ofstream out("timing_followup.csv");out<<std::setprecision(12)<<"N,regime,variant,budget,repeat,completed,converged,max_residual,max_energy,products,residuals,outer,metric_updates,setup_us,cpu_ms\n";
    std::vector<Variant> variants={{"diagonal","diagonal"},{"vertex_right","vertex_right"},{"vertex_rows","vertex_rows"},{"vertex_chunk2","vertex_rows","fixed","transformed",2},{"vertex_ew2","vertex_rows","ew2","physical_stop"},{"edge_right","edge_right"},{"edge_step_cache","edge_right","fixed","transformed",0,true},{"curvature_graph","curvature_graph"},{"curvature_step_cache","curvature_graph","fixed","transformed",0,true},{"curvature_ew2","curvature_graph","ew2","physical_stop"},{"direct_block","direct_block"}};
    for(int repeat=0;repeat<7;++repeat){std::rotate(variants.begin(),variants.begin()+1,variants.end());for(int N:{32,128})for(std::string regime:{"mass100","both100"}){Problem p=problem(N,regime,"mixed",100);double h=.04;
        for(auto v:variants)for(int budget:(v.metric=="direct_block"?std::vector<int>{0}:std::vector<int>{12,48})){
            double start=cpu_us();Metric metric(p,h,v.metric);metric.cache_per_step=v.step_cache;if(v.step_cache)metric.refresh_period=1000000;Stats stats;stats.setup_us=metric.initial_setup_us;double cpu=cpu_us()-start;State s=p.initial;int done=0,converged=0;double maxr=0,maxe=0;
            for(int k=0;k<10;++k){start=cpu_us();Step r=solve(p,s,h,metric,budget,predictor(p,s,h,"ballistic"),v.metric=="direct_block",v.chunk,v.linear,v.forcing);cpu+=cpu_us()-start;stats+=r.stats;if(!r.finite)break;s=r.state;++done;converged+=r.converged;maxr=std::max(maxr,r.residual);maxe=std::max(maxe,std::abs(energy(p.m,s)-energy(p.m,p.initial))/p.escale);}
            out<<N<<','<<regime<<','<<v.label<<','<<budget<<','<<repeat<<','<<done<<','<<converged<<','<<maxr<<','<<maxe<<','<<stats.products<<','<<stats.residuals<<','<<stats.outer<<','<<stats.metric_updates<<','<<stats.setup_us<<','<<cpu/1000<<'\n';
        }out.flush();std::cout<<"followup timing repeat="<<repeat<<" N="<<N<<" "<<regime<<std::endl;
    }}return 0;
}catch(const std::exception&e){std::cerr<<e.what()<<std::endl;return 1;}}
