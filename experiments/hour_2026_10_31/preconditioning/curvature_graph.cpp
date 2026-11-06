#include "solver.hpp"
#include <fstream>
#include <iomanip>
#include <iostream>
int main(){try{
    std::ofstream out("curvature_graph.csv"),tests("curvature_graph_identity.csv");out<<std::setprecision(12)<<"N,regime,contrast,h,budget,completed,converged,max_residual,max_energy,products,residuals,setup_us,cpu_ms\n";tests<<std::setprecision(12)<<"N,regime,converged,root_difference,positive_eigenvalue,inverse_error\n";
    for(int N:{8,32,128})for(std::string regime:{"uniform","mass100","stiff100","both100"})for(double contrast:(regime=="uniform"?std::vector<double>{1}:std::vector<double>{100,1000}))for(double h:{.01,.04}){
        Problem p=problem(N,regime,"mixed",contrast);int steps=int(std::round(.4/h));
        for(int budget:{12,48}){Metric metric(p,h,"curvature_graph");Stats stats;State s=p.initial;int done=0,converged=0;double maxr=0,maxe=0,cpu=cpu_us();
            for(int k=0;k<steps;++k){Step r=solve(p,s,h,metric,budget,predictor(p,s,h,"ballistic"));stats+=r.stats;if(!r.finite)break;s=r.state;++done;converged+=r.converged;maxr=std::max(maxr,r.residual);maxe=std::max(maxe,std::abs(energy(p.m,s)-energy(p.m,p.initial))/p.escale);}
            cpu=cpu_us()-cpu;out<<N<<','<<regime<<','<<contrast<<','<<h<<','<<budget<<','<<done<<','<<converged<<','<<maxr<<','<<maxe<<','<<stats.products<<','<<stats.residuals<<','<<stats.setup_us<<','<<cpu/1000<<'\n';
        }
        if(contrast==100&&h==.01){Metric metric(p,h,"curvature_graph"),direct(p,h,"direct_block");V guess=predictor(p,p.initial,h,"ballistic");Stats stats;Eval e=evaluate(p,p.initial,guess,h,true,stats);metric.update(e.J,stats);Mat L(2*N,2*N),Q=L;for(int k=0;k<2*N;++k){V unit=V::Zero(2*N);unit[k]=1;L.col(k)=metric.left(unit);Q.col(k)=metric.right(unit);}double mineig=Eigen::SelfAdjointEigenSolver<Mat>(Q*Q.transpose()).eigenvalues().minCoeff(),inverse_error=(L-Q.transpose()).norm();Step ref=solve(p,p.initial,h,direct,0,guess,true),r=solve(p,p.initial,h,metric,0,guess);double diff=p.msqrt.cwiseProduct(r.state.x-ref.state.x).norm()+h*p.msqrt.cwiseProduct(r.state.v-ref.state.v).norm();tests<<N<<','<<regime<<','<<r.converged<<','<<diff<<','<<mineig<<','<<inverse_error<<'\n';if(!r.converged||diff>1e-8||mineig<=0||inverse_error>1e-9)throw std::runtime_error("curvature graph identity");}
        out.flush();std::cout<<"curvature graph N="<<N<<" "<<regime<<" contrast="<<contrast<<" h="<<h<<std::endl;
    }return 0;
}catch(const std::exception&e){std::cerr<<e.what()<<std::endl;return 1;}}
