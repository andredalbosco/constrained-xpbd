#include "solver.hpp"
#include <fstream>
#include <iomanip>
#include <iostream>
int main(){try{
    std::ofstream out("forcing.csv"),tests("forcing_identity.csv");out<<std::setprecision(12)<<"N,regime,contrast,h,metric,budget,policy,requested,completed,converged,max_residual,max_energy,products,residuals,outer,rejects,setup_us,cpu_ms\n";tests<<std::setprecision(12)<<"N,metric,converged,root_difference\n";
    for(int N:{8,32,128})for(std::string regime:{"uniform","mass100","stiff100","both100"})for(double contrast:(regime=="uniform"?std::vector<double>{1}:std::vector<double>{100,1000}))for(double h:{.01,.04}){
        Problem p=problem(N,regime,"mixed",contrast);int steps=int(std::round(.4/h));
        for(std::string name:{"diagonal","graph","vertex_rows","vertex_right","edge_right","curvature_graph"})for(int budget:{12,48})for(std::string policy:{"fixed_transformed","fixed_physical","loose_physical","ew2_physical"}){
            Metric metric(p,h,name);Stats stats;State s=p.initial;int done=0,converged=0;double maxr=0,maxe=0,cpu=cpu_us();std::string forcing=policy=="ew2_physical"?"ew2":policy=="loose_physical"?"loose":"fixed",linear=policy=="fixed_transformed"?"transformed":"physical_stop";
            for(int k=0;k<steps;++k){Step r=solve(p,s,h,metric,budget,predictor(p,s,h,"ballistic"),false,0,linear,forcing);stats+=r.stats;if(!r.finite)break;s=r.state;++done;converged+=r.converged;maxr=std::max(maxr,r.residual);maxe=std::max(maxe,std::abs(energy(p.m,s)-energy(p.m,p.initial))/p.escale);}
            cpu=cpu_us()-cpu;out<<N<<','<<regime<<','<<contrast<<','<<h<<','<<name<<','<<budget<<','<<policy<<','<<steps<<','<<done<<','<<converged<<','<<maxr<<','<<maxe<<','<<stats.products<<','<<stats.residuals<<','<<stats.outer<<','<<stats.rejects<<','<<stats.setup_us<<','<<cpu/1000<<'\n';
        }out.flush();std::cout<<"forcing N="<<N<<" "<<regime<<" contrast="<<contrast<<" h="<<h<<std::endl;
    }
    for(int N:{8,32,128}){Problem p=problem(N,"mass100","mixed",100);double h=.04;V guess=predictor(p,p.initial,h,"ballistic");Metric direct(p,h,"direct_block");Step ref=solve(p,p.initial,h,direct,0,guess,true);
        for(std::string name:{"vertex_rows","vertex_right","edge_right","curvature_graph"}){Metric metric(p,h,name);Step r=solve(p,p.initial,h,metric,0,guess,false,0,"physical_stop","ew2");double diff=p.msqrt.cwiseProduct(r.state.x-ref.state.x).norm()+h*p.msqrt.cwiseProduct(r.state.v-ref.state.v).norm();tests<<N<<','<<name<<','<<r.converged<<','<<diff<<'\n';if(!r.converged||diff>1e-8)throw std::runtime_error("forcing root mismatch");}
    }return 0;
}catch(const std::exception&e){std::cerr<<e.what()<<std::endl;return 1;}}
