#include "solver.hpp"
#include <fstream>
#include <iomanip>
#include <iostream>
int main(int argc,char**argv){try{
    bool timing=argc>1;std::ofstream out(timing?"cache_timing.csv":"cache.csv"),test(timing?"cache_timing_identity.csv":"cache_identity.csv");out<<std::setprecision(12)<<"N,regime,metric,budget,cache,repeat,completed,converged,max_residual,max_energy,products,residuals,metric_updates,setup_us,cpu_ms\n";test<<std::setprecision(12)<<"N,metric,cache,converged,root_difference\n";
    for(int repeat=0;repeat<(timing?5:1);++repeat)for(int N:{32,128})for(std::string regime:{"mass100","both100"}){Problem p=problem(N,regime,"mixed",100);double h=.04;
        std::vector<std::string> names={"vertex_rows","vertex_right","edge_right","curvature_graph"};std::rotate(names.begin(),names.begin()+repeat%names.size(),names.end());
        for(auto name:names)for(int budget:{12,48})for(std::string cache:{"fresh","every2","step","trajectory"}){double start=cpu_us();Metric metric(p,h,name);if(cache=="every2")metric.refresh_period=2;if(cache=="step"||cache=="trajectory")metric.refresh_period=1000000;metric.cache_per_step=cache=="step";
            Stats stats;stats.setup_us=metric.initial_setup_us;State s=p.initial;int done=0,converged=0;double maxr=0,maxe=0;
            for(int k=0;k<10;++k){Step r=solve(p,s,h,metric,budget,predictor(p,s,h,"ballistic"));stats+=r.stats;if(!r.finite)break;s=r.state;++done;converged+=r.converged;maxr=std::max(maxr,r.residual);maxe=std::max(maxe,std::abs(energy(p.m,s)-energy(p.m,p.initial))/p.escale);}
            out<<N<<','<<regime<<','<<name<<','<<budget<<','<<cache<<','<<repeat<<','<<done<<','<<converged<<','<<maxr<<','<<maxe<<','<<stats.products<<','<<stats.residuals<<','<<stats.metric_updates<<','<<stats.setup_us<<','<<(cpu_us()-start)/1000<<'\n';
        }out.flush();std::cout<<"cache repeat="<<repeat<<" N="<<N<<" "<<regime<<std::endl;
    }
    if(!timing)for(int N:{32,128}){Problem p=problem(N,"mass100","mixed",100);double h=.04;V guess=predictor(p,p.initial,h,"ballistic");Metric direct(p,h,"direct_block");Step ref=solve(p,p.initial,h,direct,0,guess,true);
        for(std::string name:{"vertex_rows","vertex_right","edge_right","curvature_graph"})for(std::string cache:{"every2","step","trajectory"}){Metric metric(p,h,name);metric.refresh_period=cache=="every2"?2:1000000;metric.cache_per_step=cache=="step";Step r=solve(p,p.initial,h,metric,0,guess);double diff=p.msqrt.cwiseProduct(r.state.x-ref.state.x).norm()+h*p.msqrt.cwiseProduct(r.state.v-ref.state.v).norm();test<<N<<','<<name<<','<<cache<<','<<r.converged<<','<<diff<<'\n';if(!r.converged||diff>1e-8)throw std::runtime_error("cached preconditioner root mismatch");}
    }return 0;
}catch(const std::exception&e){std::cerr<<e.what()<<std::endl;return 1;}}
