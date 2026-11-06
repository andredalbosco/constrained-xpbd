#include "cross_problem.hpp"
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


struct Job{std::string id,name,target,objective,method;};
int main(int argc,char**argv){try{
    fs::path out=argc>1?argv[1]:"results";fs::create_directories(out/"cross_traces");
    std::ofstream validation(out/"objective_work_cross_validation.csv");validation<<std::setprecision(16)<<"case,target,gradient_error,hessian_error,unchanged_constraints,passed\n";
    for(std::string name:{"scalar_release","radial_orbit","chain_mixed"})for(std::string target:{"frozen_power","implicit_midpoint_dg"}){
        auto c=make_case(name);double E=std::max(1e-8,energy(c.model,c.initial)-minimum_energy(c.model));auto cfg=configuration(target,"soft_midpoint");
        CrossProblem original(c.model,c.initial,1./60,cfg,E,target,false),quarter(c.model,c.initial,1./60,cfg,E,target,true);
        V z=quarter.initial();for(int i=0;i<z.size();++i)z[i]+=.017*std::sin(i+1.);auto e=quarter.evaluate(z),o=original.evaluate(z);double ge=0,he=0,eps=1e-5;
        for(int i=0;i<z.size();++i){V a=z,b=z;a[i]+=eps;b[i]-=eps;auto ea=quarter.evaluate(a),eb=quarter.evaluate(b);ge=std::max(ge,std::abs((ea.objective.a-eb.objective.a)/(2*eps)-e.objective.g[i]));he=std::max(he,((ea.objective.g-eb.objective.g)/(2*eps)-e.objective.H.col(i)).cwiseAbs().maxCoeff());}
        bool same=(residuals(e)-residuals(o)).norm()==0&&(jacobian(e)-jacobian(o)).norm()==0;bool pass=ge<1e-6&&he<1e-6&&same;validation<<name<<','<<target<<','<<ge<<','<<he<<','<<same<<','<<pass<<'\n';if(!pass)throw std::runtime_error("cross derivative validation failed");
    }
    std::vector<Job>jobs;
    for(std::string name:{"scalar_release","radial_orbit","chain_mixed"})for(std::string target:{"frozen_power","implicit_midpoint_dg"})for(std::string objective:{"original","quarter"})for(std::string method:{"joint","local_AL"})jobs.push_back({"cross_"+std::to_string(jobs.size()),name,target,objective,method});
    std::ofstream summary(out/"objective_work_cross.csv");summary<<std::setprecision(16)<<"id,case,target,objective,kinematics,beta,gamma,method,rho,order,budget,h,duration,status,completed,steps,nonconverged_steps,max_primal,max_stationarity,max_energy_error,rms_position,path_ratio,first_displacement_ratio,final_position_error,min_converged_reduced_curvature,evaluations,compute_ms\n";
    std::atomic<int>next{0};std::mutex mutex;std::vector<std::thread>threads;
    for(int worker=0;worker<4;++worker)threads.emplace_back([&](){for(;;){int index=next.fetch_add(1);if(index>=int(jobs.size()))break;const auto&job=jobs[index];auto c=make_case(job.name);const auto&m=c.model;double h=1./60,duration=1,E0=energy(m,c.initial),scale=std::max(1e-8,E0-minimum_energy(m));auto cfg=configuration(job.target,"soft_midpoint");State s=c.initial,reference=s;int completed=0,bad=0,evaluations=0,budget=job.method=="joint"?128:512;double maxP=0,maxS=0,maxE=0,error=0,path=0,refpath=0,firstRatio=0,minCurvature=1e100,ms=0;std::string status="completed";
        std::ofstream trace(out/"cross_traces"/(job.id+".csv"));trace<<std::setprecision(16)<<"time,energy,reference_energy";for(int k=0;k<dof(m);++k)trace<<",x"<<k<<",v"<<k<<",ref_x"<<k<<",ref_v"<<k;trace<<'\n';auto row=[&](double t){trace<<t<<','<<energy(m,s)<<','<<energy(m,reference);for(int k=0;k<dof(m);++k)trace<<','<<s.x[k]<<','<<s.v[k]<<','<<reference.x[k]<<','<<reference.v[k];trace<<'\n';};row(0);
        for(int step=0;step<60;++step){try{CrossProblem p(m,s,h,cfg,scale,job.target,job.objective=="quarter");Solve solution=job.method=="joint"?joint_solve(p,budget):local_al(p,budget,10,false,"forward");State current=p.state(solution.z),refnext=reference_step(m,reference,h,64);evaluations+=solution.evaluations;ms+=solution.milliseconds;if(!solution.converged)++bad;else minCurvature=std::min(minCurvature,reduced_curvature(p,solution));maxP=std::max(maxP,solution.primal);maxS=std::max(maxS,solution.dual_residual);if(!finite(current)){status="nonfinite";break;}double ee=std::abs(energy(m,current)-E0)/scale;if(ee>1e4||current.x.norm()>1e4||current.v.norm()>1e5){status="runaway";break;}maxE=std::max(maxE,ee);path+=(current.x-s.x).norm();refpath+=(refnext.x-reference.x).norm();if(!step)firstRatio=(current.x-s.x).norm()/std::max(1e-30,(refnext.x-reference.x).norm());error+=(current.x-refnext.x).squaredNorm()/m.count;s=std::move(current);reference=std::move(refnext);++completed;row(completed*h);}catch(const std::exception&e){status=e.what();break;}}
        if(status=="completed"&&bad)status="completed_inexact";if(minCurvature==1e100)minCurvature=0;
        std::ostringstream line;line<<std::setprecision(16)<<job.id<<','<<job.name<<','<<job.target<<','<<job.objective<<",soft_midpoint,1,1,"<<job.method<<','<<(job.method=="joint"?0:10)<<','<<(job.method=="joint"?"joint":"forward")<<','<<budget<<','<<h<<','<<duration<<','<<status<<','<<completed<<",60,"<<bad<<','<<maxP<<','<<maxS<<','<<maxE<<','<<(completed?std::sqrt(error/completed):0)<<','<<(refpath?path/refpath:0)<<','<<firstRatio<<','<<(s.x-reference.x).norm()<<','<<minCurvature<<','<<evaluations<<','<<ms<<'\n';std::lock_guard<std::mutex>lock(mutex);summary<<line.str();summary.flush();std::cout<<job.id<<' '<<job.objective<<' '<<job.method<<' '<<job.name<<' '<<job.target<<' '<<status<<'\n';
    }});
    for(auto&t:threads)t.join();std::cout<<jobs.size()<<" cross-objective rollouts complete\n";return 0;
}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}
