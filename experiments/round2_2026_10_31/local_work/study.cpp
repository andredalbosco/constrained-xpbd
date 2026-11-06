#include "solvers.hpp"
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <mutex>
#include <atomic>
#include <thread>

using namespace local_work;
namespace fs=std::filesystem;
template<class F>void parallel(int count,F fn){std::atomic<int>next{0};std::vector<std::thread>threads;for(int t=0;t<4;++t)threads.emplace_back([&](){for(;;){int j=next.fetch_add(1);if(j>=count)break;fn(j);}});for(auto&t:threads)t.join();}
double escale(const Case&c){return std::max(1e-8,energy(c.model,c.initial)-minimum_energy(c.model));}
void validate(const fs::path&out){
    std::ofstream f(out/"validation.csv");f<<std::setprecision(16)<<"case,target,objective_gradient_error,constraint_gradient_error,kinetic_position_derivative,storage_velocity_derivative,force_time_derivative_error,rod_taylor_work,pass\n";bool all=true;
    for(auto name:{"scalar_mixed","radial_orbit","rod_kick","elastic_mixed","chain_mixed"})for(auto target:{"live_dg","taylor2","predicted_dg","implicit_dg","implicit_midpoint_dg"}){
        auto c=make_case(name);auto cfg=configuration(target);ExplicitProblem p(c.model,c.initial,1./120,cfg,escale(c),target);V z=p.initial();for(int i=0;i<z.size();++i)z[i]+=.003*std::sin(i+1.);auto e=p.evaluate(z);double ge=0,ce=0,kx=0,uv=0;
        for(int k=0;k<p.p.n;++k){double eps=1e-6;V a=z,b=z;a[k]+=eps;b[k]-=eps;auto ea=p.evaluate(a),eb=p.evaluate(b);ge=std::max(ge,std::abs((ea.objective.a-eb.objective.a)/(2*eps)-e.objective.g[k]));for(int j=0;j<int(e.rows.size());++j)ce=std::max(ce,std::abs((ea.rows[j].value.a-eb.rows[j].value.a)/(2*eps)-e.rows[j].value.g[k]));}
        for(auto&r:e.rows){if(r.kind==independent_detail::Kinetic)kx=std::max(kx,r.value.g.head(p.p.d).norm());if(r.kind==independent_detail::Storage||r.kind==independent_detail::GravityStorage)uv=std::max(uv,r.value.g.segment(p.p.d,p.p.d).norm());}
        auto forces=components(c.model,c.initial),plus=components(c.model,reference_step(c.model,c.initial,1e-6,4)),minus=components(c.model,reference_step(c.model,c.initial,-1e-6,4));double fd=0;for(int j=0;j<int(forces.entries.size());++j)fd=std::max(fd,((plus.entries[j].f-minus.entries[j].f)/(2e-6)-forces.entries[j].fdot).norm());double rod=c.model.rigid?work_prediction(c.model,c.initial,1./60,"taylor2")[0]:0;
        bool pass=ge<1e-6&&ce<1e-6&&kx<1e-10&&uv<1e-10&&fd<1e-5&&std::abs(rod)<1e-12;all=all&&pass;f<<name<<','<<target<<','<<ge<<','<<ce<<','<<kx<<','<<uv<<','<<fd<<','<<rod<<','<<pass<<'\n';
    }
    if(!all)throw std::runtime_error("validation_failed");std::cout<<"derivative/locality validation passed\n";
}
void predictions(const fs::path&out){
    std::ofstream f(out/"prediction_accuracy.csv"),r(out/"prediction_reference_checks.csv");f<<std::setprecision(16);r<<std::setprecision(16);
    f<<"case,scenario,kfactor,mass_ratio,amplitude,phase,h,target,incidence,body,potential,estimate,reference_work,absolute_error,normalized_error,energy_scale\n";
    r<<"case,scenario,phase,h,max_work_refinement,max_normalized_work_refinement,max_storage_defect,max_kinetic_defect\n";
    struct Scenario{std::string name;double k,m,a;};std::vector<Scenario>scenarios={{"nominal",1,1,1},{"fast_rotation",1,1,4},{"stiff",100,1,1},{"heavy",1,100,1},{"stiff_heavy",10,10,1}};
    for(auto name:{"scalar_release","scalar_mixed","radial_orbit","rod_kick","elastic_mixed","chain_mixed"})for(auto sc:scenarios){auto c=make_case(name,sc.k,sc.m,sc.a);double scale=escale(c);for(double phase:{0.,.08,.23}){State old=phase?reference_step(c.model,c.initial,phase,512):c.initial;for(double h:{1./30,1./120,1./480}){
        auto coarse=integrate_work(c.model,old,h,128),ref=integrate_work(c.model,old,h,256);auto layout=components(c.model,old);V du=potential_parts(c.model,ref.end.x)-potential_parts(c.model,old.x),storage=du,kin=V::Zero(c.model.count);
        for(int i=0;i<c.model.count;++i)kin[i]=.5*c.model.mass[i]*(ref.end.v.segment(i*c.model.dim,c.model.dim).squaredNorm()-old.v.segment(i*c.model.dim,c.model.dim).squaredNorm());
        for(int j=0;j<int(layout.entries.size());++j){storage[layout.entries[j].potential]+=ref.work[j];kin[layout.entries[j].body]-=ref.work[j];}
        double error=(coarse.work-ref.work).lpNorm<Eigen::Infinity>();r<<name<<','<<sc.name<<','<<phase<<','<<h<<','<<error<<','<<error/scale<<','<<storage.lpNorm<Eigen::Infinity>()/scale<<','<<kin.lpNorm<Eigen::Infinity>()/scale<<'\n';
        for(const auto&law:prediction_names()){V pred=work_prediction(c.model,old,h,law);for(int j=0;j<pred.size();++j){double error=std::abs(pred[j]-ref.work[j]);f<<name<<','<<sc.name<<','<<sc.k<<','<<sc.m<<','<<sc.a<<','<<phase<<','<<h<<','<<law<<','<<j<<','<<layout.entries[j].body<<','<<layout.entries[j].potential<<','<<pred[j]<<','<<ref.work[j]<<','<<error<<','<<error/scale<<','<<scale<<'\n';}}
    }}std::cout<<"prediction "<<name<<' '<<sc.name<<" done\n";}
}
struct Fixture{std::string id,target,kin;Case c;State old;double h=1./60;};
void al_comparison(const fs::path&out){
    std::vector<Fixture>fixtures;for(std::string name:{"scalar_release","scalar_mixed","chain_mixed"})for(std::string target:{"live_dg","taylor2","predicted_dg","implicit_midpoint_dg"})for(std::string kin:{"soft_endpoint","soft_midpoint"}){auto c=make_case(name);fixtures.push_back({name+"_"+target+"_"+kin,target,kin,c,c.initial});}
    std::ofstream f(out/"local_solver_comparison.csv");f<<std::setprecision(16)<<"fixture,case,target,kinematics,method,rho_initial,adaptive,order,budget,evaluations,block_visits,dual_sweeps,rho_final,objective,primal_residual,stationarity,converged,reason,joint_reference_converged,joint_reference_objective,objective_gap,state_distance_to_joint,regularizations,compute_ms,reduced_min_eigenvalue,joint_reduced_min_eigenvalue\n";std::mutex mutex;
    parallel(fixtures.size(),[&](int index){const auto&fixture=fixtures[index];auto cfg=configuration(fixture.target,fixture.kin);ExplicitProblem p(fixture.c.model,fixture.old,fixture.h,cfg,escale(fixture.c),fixture.target);auto ref=joint_solve(p,4096);std::ostringstream rows;rows<<std::setprecision(16);
        auto row=[&](const Solve&s,std::string method,double rho,bool adaptive,std::string order,int budget){double gap=ref.converged?s.objective-ref.objective:std::numeric_limits<double>::quiet_NaN();rows<<fixture.id<<','<<fixture.c.name<<','<<fixture.target<<','<<fixture.kin<<','<<method<<','<<rho<<','<<adaptive<<','<<order<<','<<budget<<','<<s.evaluations<<','<<s.visits<<','<<s.sweeps<<','<<s.rho<<','<<s.objective<<','<<s.primal<<','<<s.dual_residual<<','<<s.converged<<','<<s.reason<<','<<ref.converged<<','<<ref.objective<<','<<gap<<','<<(s.z-ref.z).norm()<<','<<s.regularizations<<','<<s.milliseconds<<','<<reduced_curvature(p,s)<<','<<reduced_curvature(p,ref)<<'\n';};
        row(ref,"joint_reference",0,false,"joint",4096);
        for(int budget:{32,128,512}){row(joint_solve(p,budget),"joint_budget",0,false,"joint",budget);for(double rho:{.1,1.,10.,100.})for(bool adaptive:{false,true})for(std::string order:{"forward","reverse","symmetric"})row(local_al(p,budget,rho,adaptive,order),"local_AL",rho,adaptive,order,budget);}
        std::lock_guard<std::mutex>lock(mutex);f<<rows.str();f.flush();std::cout<<"local AL "<<fixture.id<<" reference="<<ref.converged<<" done\n";
    });
}
struct DynamicsJob{Case c;std::string target,kin;double h,duration;};
void dynamics(const fs::path&out){
    std::vector<DynamicsJob>jobs;std::vector<std::string>laws=prediction_names();laws.push_back("live_start");laws.push_back("live_dg");for(std::string name:{"scalar_release","radial_orbit","rod_kick","elastic_mixed","chain_mixed"})for(double k:{1.,10.}){if(name=="rod_kick"&&k!=1)continue;for(double h:{1./60,1./240})for(std::string kin:{"soft_endpoint","soft_midpoint"})for(const auto&target:laws){if(target=="projected_dg"&&name!="rod_kick")continue;if(name=="rod_kick"&&(target=="implicit_dg"||target=="implicit_midpoint_dg"))continue;jobs.push_back({make_case(name,k),target,kin,h,1.});}}
    std::ofstream f(out/"frozen_target_dynamics.csv");f<<std::setprecision(16)<<"case,kfactor,target,kinematics,tangency,h,duration,status,completed,steps,nonconverged_steps,max_primal_residual,max_stationarity,max_energy_error,rms_position,path_ratio,first_displacement_ratio,final_position_error,evaluations,compute_ms\n";std::mutex mutex;
    parallel(jobs.size(),[&](int index){const auto&job=jobs[index];const auto&m=job.c.model;auto cfg=configuration(job.target,job.kin);cfg.tangency=m.rigid;State s=job.c.initial,ref=s;double E0=energy(m,s),scale=escale(job.c),maxP=0,maxS=0,maxE=0,sumerror=0,path=0,refpath=0,first=0,ms=0;int steps=std::lround(job.duration/job.h),completed=0,bad=0,eval=0;std::string status="completed";
        for(int k=0;k<steps;++k){try{ExplicitProblem p(m,s,job.h,cfg,scale,job.target);auto solved=joint_solve(p,128);auto next=p.state(solved.z),refnext=reference_step(m,ref,job.h,64);eval+=solved.evaluations;ms+=solved.milliseconds;if(!solved.converged)++bad;maxP=std::max(maxP,solved.primal);maxS=std::max(maxS,solved.dual_residual);if(!finite(next)){status="nonfinite";break;}double ee=std::abs(energy(m,next)-E0)/scale;if(ee>1e4||next.x.norm()>1e4||next.v.norm()>1e5){status="runaway";break;}maxE=std::max(maxE,ee);double dx=(next.x-s.x).norm(),rdx=(refnext.x-ref.x).norm();if(k==0)first=rdx>1e-14?dx/rdx:0;path+=dx;refpath+=rdx;sumerror+=(next.x-refnext.x).squaredNorm()/m.count;s=std::move(next);ref=std::move(refnext);++completed;}catch(const std::exception&e){status=e.what();break;}}
        std::ostringstream row;row<<std::setprecision(16)<<job.c.name<<','<<m.stiffness[0]/100<<','<<job.target<<','<<job.kin<<','<<cfg.tangency<<','<<job.h<<','<<job.duration<<','<<status<<','<<completed<<','<<steps<<','<<bad<<','<<maxP<<','<<maxS<<','<<maxE<<','<<(completed?std::sqrt(sumerror/completed):0)<<','<<(refpath?path/refpath:0)<<','<<first<<','<<(s.x-ref.x).norm()<<','<<eval<<','<<ms<<'\n';std::lock_guard<std::mutex>lock(mutex);f<<row.str();f.flush();if(index%20==0)std::cout<<"dynamics job "<<index<<'/'<<jobs.size()<<" done\n";
    });
}
void reference_checks(const fs::path&out){std::ofstream f(out/"dynamics_reference_checks.csv");f<<std::setprecision(16)<<"case,kfactor,h,position_refinement,velocity_refinement,normalized_energy_error,pass\n";bool all=true;
    for(std::string name:{"scalar_release","radial_orbit","rod_kick","elastic_mixed","chain_mixed"})for(double k:{1.,10.}){if(name=="rod_kick"&&k!=1)continue;auto c=make_case(name,k);for(double h:{1./60,1./240}){int steps=std::lround(1/h);auto a=reference_step(c.model,c.initial,1.,steps*64),b=reference_step(c.model,c.initial,1.,steps*128);double dx=(a.x-b.x).norm(),dv=(a.v-b.v).norm(),de=std::abs(energy(c.model,b)-energy(c.model,c.initial))/escale(c);bool pass=dx<1e-7&&dv<1e-6&&de<1e-7;all=all&&pass;f<<name<<','<<k<<','<<h<<','<<dx<<','<<dv<<','<<de<<','<<pass<<'\n';}}
    if(!all)throw std::runtime_error("dynamics_reference_refinement_failed");std::cout<<"dynamics references refined successfully\n";
}
int main(int argc,char**argv){try{std::string mode=argc>1?argv[1]:"all";fs::path out=argc>2?argv[2]:"results";fs::create_directories(out);if(mode=="all"||mode=="validate"){validate(out);reference_checks(out);}if(mode=="all"||mode=="predictions")predictions(out);if(mode=="all"||mode=="al")al_comparison(out);if(mode=="all"||mode=="dynamics")dynamics(out);return 0;}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}
