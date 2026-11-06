#include "solver.hpp"
#include <fstream>
#include <iomanip>
#include <iostream>

int main(){try{
    Eigen::setNbThreads(1);std::ofstream out("krylov_ablation.csv"),linear("linear_probe.csv"),identities("krylov_identity.csv");
    out<<std::setprecision(12)<<"N,regime,contrast,metric,budget,chunk,policy,completed,converged,max_residual,max_energy,products,residuals,outer,rejects,setup_us,linear_us,cpu_ms\n";
    linear<<std::setprecision(12)<<"N,regime,contrast,metric,budget,policy,physical_relative,transformed_relative,direction_relative,products,cpu_us\n";
    identities<<std::setprecision(12)<<"N,metric,policy,converged,root_difference\n";
    struct Case {int N;std::string regime;double contrast;};
    std::vector<Case> cases={{32,"mass100",100},{128,"mass100",100},{32,"stiff100",100},{32,"both100",100},{128,"both100",100},{128,"both100",1000}};
    for(auto c:cases){Problem p=problem(c.N,c.regime,"mixed",c.contrast);double h=.04;
        for(std::string name:{"diagonal","graph","vertex_rows","vertex_right","edge_rows","edge_right"}){
            Metric metric(p,h,name);Stats st;V guess=predictor(p,p.initial,h,"ballistic");Eval e=evaluate(p,p.initial,guess,h,true,st);metric.update(e.J,st);bool ok;V exact=e.J.solve(-e.r,ok);
            for(int budget:{2,4,6,8,12,24,48})for(std::string policy:{"transformed","physical_stop","physical_min"}){Stats stats;double start=cpu_us();V dx=gmres(e.J,e.r,metric,budget,1e-3,stats,policy);double cpu=cpu_us()-start;V remaining=e.r+e.J.apply(dx);linear<<c.N<<','<<c.regime<<','<<c.contrast<<','<<name<<','<<budget<<','<<policy<<','<<p.misqrt.cwiseProduct(remaining).norm()/p.misqrt.cwiseProduct(e.r).norm()<<','<<metric.left(remaining).norm()/metric.left(e.r).norm()<<','<<(dx-exact).norm()/std::max(1e-30,exact.norm())<<','<<stats.products<<','<<cpu<<'\n';}
            for(int budget:{12,48})for(int chunk:(budget==12?std::vector<int>{2,4,6,12}:std::vector<int>{4,8,12,24,48}))for(std::string policy:{"transformed","physical_stop","physical_min"}){
                Metric m(p,h,name);Stats stats;stats.setup_us=m.initial_setup_us;State s=p.initial;int done=0,converged=0;double maxr=0,maxe=0,start=cpu_us();
                for(int step=0;step<10;++step){Step r=solve(p,s,h,m,budget,predictor(p,s,h,"ballistic"),false,chunk,policy);stats+=r.stats;if(!r.finite)break;s=r.state;++done;converged+=r.converged;maxr=std::max(maxr,r.residual);maxe=std::max(maxe,std::abs(energy(p.m,s)-energy(p.m,p.initial))/p.escale);}
                double cpu=cpu_us()-start;out<<c.N<<','<<c.regime<<','<<c.contrast<<','<<name<<','<<budget<<','<<chunk<<','<<policy<<','<<done<<','<<converged<<','<<maxr<<','<<maxe<<','<<stats.products<<','<<stats.residuals<<','<<stats.outer<<','<<stats.rejects<<','<<stats.setup_us<<','<<stats.linear_us<<','<<cpu/1000<<'\n';
            }
        }out.flush();linear.flush();std::cout<<"ablation N="<<c.N<<" "<<c.regime<<" contrast="<<c.contrast<<std::endl;
    }
    for(int N:{8,32}){Problem p=problem(N,"mass100","mixed",100);double h=.04;V guess=predictor(p,p.initial,h,"ballistic");Metric direct(p,h,"direct_block");Step reference=solve(p,p.initial,h,direct,0,guess,true);
        for(std::string name:{"vertex_rows","vertex_right","edge_rows","edge_right"})for(std::string policy:{"transformed","physical_stop","physical_min"}){Metric metric(p,h,name);Step r=solve(p,p.initial,h,metric,0,guess,false,0,policy);double diff=p.msqrt.cwiseProduct(r.state.x-reference.state.x).norm()+h*p.msqrt.cwiseProduct(r.state.v-reference.state.v).norm();identities<<N<<','<<name<<','<<policy<<','<<r.converged<<','<<diff<<'\n';if(!r.converged||diff>1e-8)throw std::runtime_error("Krylov policy root mismatch");}
    }return 0;
}catch(const std::exception&e){std::cerr<<e.what()<<std::endl;return 1;}}
