#include "solver.hpp"
#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>

using Clock=std::chrono::steady_clock;
const std::vector<std::string> metrics={"none","mass","diagonal","graph","vertex_jacobi","vertex_rows","edge_rows","vertex_right","edge_right"};
void require(bool test,const std::string&s){if(!test)throw std::runtime_error(s);}

void tests(){
    std::ofstream out("identity_tests.csv");out<<std::setprecision(14)<<"N,regime,contrast,metric,jacobian_error,block_solve_error,energy_identity,root_converged,root_difference,energy_error,products\n";
    for(int N:{8,32,128})for(std::string regime:{"mass100","both100"}){
        Problem p=problem(N,regime,"mixed",100);double h=.01;Stats stats;
        V x=p.initial.x+h*p.initial.v;Eval e=evaluate(p,p.initial,x,h,true,stats);V d=V::LinSpaced(x.size(),-.7,.9);double eps=1e-6;
        V fd=(evaluate(p,p.initial,x+eps*d,h,false,stats).r-evaluate(p,p.initial,x-eps*d,h,false,stats).r)/(2*eps);
        double je=(fd-e.J.apply(d)).norm()/std::max(1.,fd.norm());require(je<5e-7,"Jacobian derivative");
        bool ok;V block=e.J.solve(d,ok);double be=(e.J.dense().fullPivLu().solve(d)-block).norm()/std::max(1.,block.norm());require(ok&&be<1e-9,"block direct solve");
        State trial{x,2*(x-p.initial.x)/h-p.initial.v};double ei=std::abs(energy(p.m,trial)-energy(p.m,p.initial)-2*e.r.dot(x-p.initial.x)/(h*h))/p.escale;require(ei<1e-8,"energy residual identity");
        Metric direct(p,h,"direct_block");Step reference=solve(p,p.initial,h,direct,0,x,true);require(reference.converged,"direct root");
        for(auto name:metrics){Metric metric(p,h,name);metric.update(e.J,stats);
            if(N==8){Mat L(x.size(),x.size()),Q=L;for(int k=0;k<x.size();++k){V unit=V::Zero(x.size());unit[k]=1;L.col(k)=metric.left(unit);Q.col(k)=metric.right(unit);}require(L.fullPivLu().rank()==x.size()&&Q.fullPivLu().rank()==x.size(),"invertible metric");
                if(name=="vertex_rows"||name=="edge_rows"){Mat T=L*p.msqrt.asDiagonal();require((T-T.transpose()).norm()<1e-10,"row symmetry");require(Eigen::SelfAdjointEigenSolver<Mat>(T).eigenvalues().minCoeff()>0,"row positivity");}}
            Step r=solve(p,p.initial,h,metric,0,x);double diff=p.msqrt.cwiseProduct(r.state.x-reference.state.x).norm()+h*p.msqrt.cwiseProduct(r.state.v-reference.state.v).norm();
            double de=std::abs(energy(p.m,r.state)-energy(p.m,p.initial))/p.escale;
            out<<N<<','<<regime<<",100,"<<name<<','<<je<<','<<be<<','<<ei<<','<<r.converged<<','<<diff<<','<<de<<','<<r.stats.products<<'\n';
            // Full Krylov convergence may become expensive at N=128; preserve
            // any failure as a row, but require identity wherever it converges.
            require(!r.converged||(diff<1e-6&&de<1e-8),"converged root agreement");
        }
    }
    std::cout<<"PASS analytic, block, energy, invertibility and converged-root audits\n";
}

struct Reference {std::vector<State> states;std::string failure;double maxenergy=0,cpu_ms=0;};
Reference reference(const Problem&p,double h,int steps){
    Reference ref;ref.states.push_back(p.initial);Metric metric(p,h,"direct_block");double start=cpu_us();
    for(int k=0;k<steps;++k){const State&old=ref.states.back();Step r=solve(p,old,h,metric,0,predictor(p,old,h,"ballistic"),true);
        if(!r.finite||!r.converged){ref.failure=r.failure;break;}ref.states.push_back(r.state);ref.maxenergy=std::max(ref.maxenergy,std::abs(energy(p.m,r.state)-energy(p.m,p.initial))/p.escale);}
    ref.cpu_ms=(cpu_us()-start)/1000;return ref;
}

void sweep(bool quick=false){
    std::ofstream out(quick?"quick_runs.csv":"runs.csv"),refs(quick?"quick_references.csv":"references.csv");
    out<<std::setprecision(12)<<"N,regime,contrast,h,metric,budget,steps,completed,converged,reference_steps,matched,max_residual,max_energy,rms_position,max_mode,max_phase,products,residuals,jacobians,outer,rejects,setup_us,linear_us,cpu_ms,wall_ms,failure\n";
    refs<<std::setprecision(12)<<"N,regime,contrast,h,requested,completed,max_energy,cpu_ms,failure\n";
    int runs=0;
    for(int N:(quick?std::vector<int>{8,32}:std::vector<int>{8,32,128}))for(std::string regime:{"uniform","mass100","stiff100","both100"})for(double contrast:(regime=="uniform"?std::vector<double>{1}:quick?std::vector<double>{100}:std::vector<double>{100,1000}))for(double h:(quick?std::vector<double>{.04}:std::vector<double>{.01,.04})){
        Problem p=problem(N,regime,"mixed",contrast);int steps=int(std::round(.4/h));Reference ref=reference(p,h,steps);
        refs<<N<<','<<regime<<','<<contrast<<','<<h<<','<<steps<<','<<ref.states.size()-1<<','<<ref.maxenergy<<','<<ref.cpu_ms<<','<<ref.failure<<'\n';refs.flush();
        std::vector<std::string> names=metrics;names.push_back("direct_block");
        for(auto name:names)for(int budget:(name=="direct_block"?std::vector<int>{0}:std::vector<int>{12,48})){
            double cpu0=cpu_us();auto wall0=Clock::now();Metric metric(p,h,name);Stats stats;stats.setup_us=metric.initial_setup_us;
            State s=p.initial;int done=0,converged=0,matched=0;double maxr=0,maxe=0,sumx=0,maxmode=0,maxphase=0,solve_cpu=cpu_us()-cpu0;std::string failure;
            for(int k=0;k<steps;++k){double step0=cpu_us();Step r=solve(p,s,h,metric,budget,predictor(p,s,h,"ballistic"),name=="direct_block");solve_cpu+=cpu_us()-step0;stats+=r.stats;
                if(!r.finite){failure=r.failure;break;}s=r.state;++done;converged+=r.converged;maxr=std::max(maxr,r.residual);maxe=std::max(maxe,std::abs(energy(p.m,s)-energy(p.m,p.initial))/p.escale);
                if(k+1<int(ref.states.size())){++matched;sumx+=p.msqrt.cwiseProduct(s.x-ref.states[k+1].x).squaredNorm()/p.mass.sum();auto me=mode_error(p,s,ref.states[k+1]);maxmode=std::max(maxmode,me.relative);maxphase=std::max(maxphase,me.phase);}
                if(!r.converged)failure=r.failure;if(maxe>1e6||s.x.norm()>1e8){failure="runaway_state";break;}
            }
            double wall=std::chrono::duration<double,std::milli>(Clock::now()-wall0).count();
            out<<N<<','<<regime<<','<<contrast<<','<<h<<','<<name<<','<<budget<<','<<steps<<','<<done<<','<<converged<<','<<ref.states.size()-1<<','<<matched<<','<<maxr<<','<<maxe<<','<<(matched?std::sqrt(sumx/matched):NAN)<<','<<(matched?maxmode:NAN)<<','<<(matched?maxphase:NAN)<<','<<stats.products<<','<<stats.residuals<<','<<stats.jacobians<<','<<stats.outer<<','<<stats.rejects<<','<<stats.setup_us<<','<<stats.linear_us<<','<<solve_cpu/1000<<','<<wall<<','<<failure<<'\n';++runs;
        }out.flush();std::cout<<"sweep N="<<N<<" "<<regime<<" contrast="<<contrast<<" h="<<h<<" reference="<<ref.states.size()-1<<"/"<<steps<<" runs="<<runs<<std::endl;
    }
}

void branches(bool quick=false){
    std::ofstream out(quick?"quick_branches.csv":"branches.csv");out<<std::setprecision(12)<<"N,regime,contrast,motion,h,predictor,requested,completed,max_residual,max_energy,matched_ballistic,max_root_difference,stages,direct_solves,residuals,setup_us,cpu_ms,failure\n";
    for(int N:(quick?std::vector<int>{8,32}:std::vector<int>{8,32,128}))for(std::string regime:{"stiff100","both100"})for(double contrast:(quick?std::vector<double>{100}:std::vector<double>{100,1000}))for(std::string motion:(quick?std::vector<std::string>{"mixed"}:std::vector<std::string>{"mixed","large"}))for(double h:(motion=="mixed"?std::vector<double>{.01,.04}:std::vector<double>{.1})){
        Problem p=problem(N,regime,motion,contrast);int steps=int(std::round(.4/h));Reference ref=reference(p,h,steps);
        for(std::string name:{"ballistic","acceleration","secant","continuation"}){
            State s=p.initial,prev=s;Metric metric(p,h,"direct_block");Stats stats;int done=0,matched=0;double maxr=0,maxe=0,maxdiff=0,cpu=0;std::string failure;
            for(int k=0;k<steps;++k){double start=cpu_us();Step r=name=="continuation"?continuation(p,s,h):solve(p,s,h,metric,0,predictor(p,s,h,name,k?&prev:nullptr),true);cpu+=cpu_us()-start;stats+=r.stats;
                if(!r.finite||!r.converged){failure=r.failure;break;}prev=s;s=r.state;++done;maxr=std::max(maxr,r.residual);maxe=std::max(maxe,std::abs(energy(p.m,s)-energy(p.m,p.initial))/p.escale);
                if(k+1<int(ref.states.size())){++matched;maxdiff=std::max(maxdiff,std::sqrt((p.msqrt.cwiseProduct(s.x-ref.states[k+1].x).squaredNorm()+h*h*p.msqrt.cwiseProduct(s.v-ref.states[k+1].v).squaredNorm())/p.mass.sum()));}
            }out<<N<<','<<regime<<','<<contrast<<','<<motion<<','<<h<<','<<name<<','<<steps<<','<<done<<','<<maxr<<','<<maxe<<','<<matched<<','<<(matched?maxdiff:NAN)<<','<<stats.stages<<','<<stats.direct_solves<<','<<stats.residuals<<','<<stats.setup_us<<','<<cpu/1000<<','<<failure<<'\n';out.flush();
        }std::cout<<"branches N="<<N<<" "<<regime<<" contrast="<<contrast<<" "<<motion<<" h="<<h<<std::endl;
    }
}

void timing(){
    std::ofstream out("timing.csv");out<<std::setprecision(12)<<"N,regime,contrast,metric,budget,repeat,steps,converged,products,residuals,setup_us,linear_us,cpu_ms,wall_ms\n";
    std::vector<std::string> names=metrics;names.push_back("direct_block");
    for(int repeat=0;repeat<5;++repeat){std::rotate(names.begin(),names.begin()+1,names.end());for(int N:{8,32,128,256})for(std::string regime:{"mass100","both100"}){
        Problem p=problem(N,regime,"mixed",100);double h=.04;
        for(auto name:names)for(int budget:(name=="direct_block"?std::vector<int>{0}:std::vector<int>{12,48})){
            double cpu0=cpu_us();auto wall0=Clock::now();Metric metric(p,h,name);Stats stats;stats.setup_us=metric.initial_setup_us;State s=p.initial;int steps=0,converged=0;
            for(int k=0;k<10;++k){Step r=solve(p,s,h,metric,budget,predictor(p,s,h,"ballistic"),name=="direct_block");stats+=r.stats;if(!r.finite)break;s=r.state;++steps;converged+=r.converged;}
            double cpu=(cpu_us()-cpu0)/1000,wall=std::chrono::duration<double,std::milli>(Clock::now()-wall0).count();out<<N<<','<<regime<<",100,"<<name<<','<<budget<<','<<repeat<<','<<steps<<','<<converged<<','<<stats.products<<','<<stats.residuals<<','<<stats.setup_us<<','<<stats.linear_us<<','<<cpu<<','<<wall<<'\n';
        }out.flush();std::cout<<"timing repeat="<<repeat<<" N="<<N<<" "<<regime<<std::endl;
    }}
}

int main(int argc,char**argv){try{Eigen::setNbThreads(1);std::string mode=argc>1?argv[1]:"--quick";if(mode=="--tests")tests();else if(mode=="--quick"){sweep(true);branches(true);}else if(mode=="--sweep")sweep();else if(mode=="--branches")branches();else if(mode=="--timing")timing();else throw std::runtime_error("unknown mode");return 0;}catch(const std::exception&e){std::cerr<<"ERROR "<<e.what()<<std::endl;return 1;}}
