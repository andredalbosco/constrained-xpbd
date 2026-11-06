#include "solvers.hpp"
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <iostream>
#include <thread>
#include <mutex>
#include <atomic>

using namespace local_work;
namespace fs=std::filesystem;
struct Job{std::string id,name,target,kin,method,order;double rho=10;int budget=128;};
int main(int argc,char**argv){try{
    fs::path out=argc>1?argv[1]:"results";fs::create_directories(out/"local_traces");std::vector<Job>jobs;
    for(std::string name:{"scalar_release","scalar_mixed","chain_mixed"})for(std::string target:{"live_dg","implicit_midpoint_dg"})for(std::string kin:{"soft_endpoint","soft_midpoint"}){
        auto add=[&](std::string method,std::string order,double rho,int budget){jobs.push_back({"local_"+std::to_string(jobs.size()),name,target,kin,method,order,rho,budget});};
        add("joint","joint",0,128);for(int budget:{128,512})for(std::string order:{"forward","reverse","symmetric"})add("local_AL",order,10,budget);if(name=="scalar_mixed")add("local_AL","reverse",100,512);
    }
    std::ofstream summary(out/"local_solver_dynamics.csv");summary<<std::setprecision(16)<<"id,case,target,kinematics,method,rho,order,budget,h,duration,status,completed,steps,nonconverged_steps,max_primal,max_stationarity,max_energy_error,rms_position,path_ratio,final_position_error,evaluations,compute_ms\n";
    std::atomic<int>next{0};std::mutex mutex;std::vector<std::thread>threads;
    for(int thread=0;thread<4;++thread)threads.emplace_back([&](){for(;;){int index=next.fetch_add(1);if(index>=int(jobs.size()))break;const auto&job=jobs[index];auto c=make_case(job.name);const auto&m=c.model;double h=1./60,duration=.5,scale=std::max(1e-8,energy(m,c.initial)-minimum_energy(m)),E0=energy(m,c.initial);auto cfg=configuration(job.target,job.kin);State s=c.initial,reference=s;int completed=0,bad=0,evaluations=0;double maxP=0,maxS=0,maxE=0,error=0,path=0,refpath=0,ms=0;std::string status="completed";
        std::ofstream trace(out/"local_traces"/(job.id+".csv"));trace<<std::setprecision(16)<<"time,energy,reference_energy";for(int k=0;k<dof(m);++k)trace<<",x"<<k<<",v"<<k<<",ref_x"<<k<<",ref_v"<<k;trace<<'\n';auto row=[&](double t){trace<<t<<','<<energy(m,s)<<','<<energy(m,reference);for(int k=0;k<dof(m);++k)trace<<','<<s.x[k]<<','<<s.v[k]<<','<<reference.x[k]<<','<<reference.v[k];trace<<'\n';};row(0);
        for(int step=0;step<30;++step){try{ExplicitProblem p(m,s,h,cfg,scale,job.target);Solve solution=job.method=="joint"?joint_solve(p,job.budget):local_al(p,job.budget,job.rho,false,job.order);State current=p.state(solution.z),refnext=reference_step(m,reference,h,64);evaluations+=solution.evaluations;ms+=solution.milliseconds;if(!solution.converged)++bad;maxP=std::max(maxP,solution.primal);maxS=std::max(maxS,solution.dual_residual);if(!finite(current)){status="nonfinite";break;}double ee=std::abs(energy(m,current)-E0)/scale;if(ee>1e4||current.x.norm()>1e4||current.v.norm()>1e5){status="runaway";break;}maxE=std::max(maxE,ee);path+=(current.x-s.x).norm();refpath+=(refnext.x-reference.x).norm();error+=(current.x-refnext.x).squaredNorm()/m.count;s=std::move(current);reference=std::move(refnext);++completed;row(completed*h);}catch(const std::exception&e){status=e.what();break;}}
        if(status=="completed"&&bad)status="completed_inexact";std::ostringstream line;line<<std::setprecision(16)<<job.id<<','<<job.name<<','<<job.target<<','<<job.kin<<','<<job.method<<','<<job.rho<<','<<job.order<<','<<job.budget<<','<<h<<','<<duration<<','<<status<<','<<completed<<",30,"<<bad<<','<<maxP<<','<<maxS<<','<<maxE<<','<<(completed?std::sqrt(error/completed):0)<<','<<(refpath?path/refpath:0)<<','<<(s.x-reference.x).norm()<<','<<evaluations<<','<<ms<<'\n';std::lock_guard<std::mutex>lock(mutex);summary<<line.str();summary.flush();if(index%20==0)std::cout<<"AL rollout "<<index<<'/'<<jobs.size()<<" done\n";
    }});
    for(auto&t:threads)t.join();std::cout<<jobs.size()<<" AL/joint rollouts complete\n";return 0;
}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}
