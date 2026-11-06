#include "dual_gs.hpp"
#include "oracle.hpp"
#include "problem.hpp"
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>

using namespace compose;
using Clock=std::chrono::steady_clock;
static const std::vector<std::string> methods={"refreshed_drop","refreshed_local_rx","frozen_full_rx"};
static const std::vector<std::string> schedules={"material_storage_kinetic","kinetic_storage_material","storage_kinetic_material","edge","edge_storage_first","edge_reverse"};
static double elapsed(Clock::time_point a){return std::chrono::duration<double,std::milli>(Clock::now()-a).count();}
static double scale(const Case& c){return std::max(1e-8,energy(c.model,c.initial)-minimum_energy(c.model));}
static double fullres(const dual_gs::Residuals&r){return std::max({r.material,r.kinetic,r.storage,r.geometry,r.tangency,r.stationarity});}
static std::vector<Case> primary(){std::vector<Case> out;for(const auto&name:case_names())out.push_back(make_case(name));out.push_back(hour_local::chain_case(3));return out;}
static double state_gap(const Model&m,const State&a,const State&b,double E){return std::sqrt((masses(m).array()*((a.v-b.v).array().square())).sum()/E+(a.x-b.x).squaredNorm()/m.count);}
static std::ofstream output(const std::string&dir,const std::string&file){std::ofstream out(dir+"/"+file);if(!out)throw std::runtime_error("cannot_write_"+file);out<<std::setprecision(17);return out;}

static void write_row(std::ostream&out,const Case&c,double h,const std::string&method,const std::string&schedule,double omega,int budget,const dual_gs::Problem&p,const dual_gs::Result&r,double setup,const State&ref,double refgap,const dual_gs::OracleResult&o){
    State s=p.state(r.z);double objective=p.objective(r.z);double rootgap=o.converged?state_gap(c.model,s,o.state,p.E):std::numeric_limits<double>::quiet_NaN();
    out<<c.name<<','<<c.model.count<<','<<h<<','<<method<<','<<schedule<<','<<omega<<','<<budget<<','<<r.visits<<','<<r.completed_sweeps<<','<<r.row_evaluations<<','<<r.coordinate_writes<<','<<r.stationarity_passes<<','<<setup<<','<<r.milliseconds<<','<<r.material<<','<<r.kinetic<<','<<r.storage<<','<<r.geometry<<','<<r.tangency<<','<<r.stationarity<<','<<r.energy<<','<<fullres(r)<<','<<objective<<','<<rootgap<<','<<state_gap(c.model,s,ref,p.E)<<','<<refgap<<','<<r.converged<<','<<r.failure<<'\n';
}

static void one_steps(const std::string&dir,bool controls){
    auto out=output(dir,controls?"controls.csv":"one_step.csv");
    out<<"case,N,h,method,schedule,omega,budget,visits,sweeps,row_evaluations,coordinate_writes,restoration_passes,setup_ms,solve_ms,material,kinetic,storage,geometry,tangency,stationarity,energy,full_residual,objective,oracle_gap,physical_gap,reference_gap,converged,failure\n";
    auto oracleout=output(dir,controls?"control_oracles.csv":"oracles.csv");
    oracleout<<"case,N,h,converged,primal,stationarity,curvature,rank,objective,physical_gap,reference_gap\n";
    auto cases=primary();std::vector<double> steps={.02,.08,.2};std::vector<int> budgets={1,2,4,8,16,32,128};std::vector<double> omegas={1};
    if(controls){cases.clear();for(const auto&name:{"scalar_release","radial_orbit","rod_release","elastic_mixed","chain_mixed"}){auto c=make_case(name,10);c.name+="_stiff10";cases.push_back(c);}auto mass=make_case("chain_mixed",1,10);mass.name+="_mass10";cases.push_back(mass);for(const auto&name:{"scalar_release","rod_release"}){auto c=make_case(name,1,1,.001);c.name+="_nearrest";cases.push_back(c);}cases.push_back(hour_local::chain_case(3));steps={.08};budgets={8,32,128,512};omegas={1,.5};}
    for(const auto&c:cases)for(double h:steps){double E=scale(c);auto start=Clock::now();dual_gs::Problem p(c.model,c.initial,h,E);double setup=elapsed(start);int q=p.rows.size();auto ref=reference_step(c.model,c.initial,h,128),fine=reference_step(c.model,c.initial,h,256);double refgap=state_gap(c.model,ref,fine,E);auto oracle=dual_gs::solve_oracle(c.model,c.initial,h,E);
        oracleout<<c.name<<','<<c.model.count<<','<<h<<','<<oracle.converged<<','<<oracle.primal<<','<<oracle.stationarity<<','<<oracle.curvature<<','<<oracle.jacobian_rank<<','<<oracle.objective<<','<<state_gap(c.model,oracle.state,fine,E)<<','<<refgap<<'\n';
        for(const auto&method:methods)for(const auto&schedule:schedules)for(double omega:omegas)for(int b:budgets){auto r=dual_gs::solve(p,method,schedule,b*q+1,omega,b*q);write_row(out,c,h,method,schedule,omega,b*q,p,r,setup,fine,refgap,oracle);}
        std::cerr<<"one_step "<<c.name<<" h="<<h<<" oracle="<<oracle.converged<<'\n';
    }
}

static void rollouts(const std::string&dir,bool longer){
    auto out=output(dir,longer?"long_rollouts.csv":"short_rollouts.csv");
    out<<"case,N,h,method,schedule,omega,budget_multiplier,target_time,accepted_time,steps,failed,reason,max_energy,max_step_energy,max_material,max_kinetic,max_storage,max_geometry,max_tangency,max_stationarity,rms_position,rms_velocity,final_position,final_velocity,reference_gap,visits,row_evaluations,coordinate_writes,setup_ms,solve_ms\n";
    auto trajectory=output(dir,longer?"long_trajectories.csv":"short_trajectories.csv");
    trajectory<<"case,h,method,schedule,omega,budget_multiplier,t,energy,stationarity,position_error,velocity_error,x0,v0\n";
    auto referenceout=output(dir,longer?"long_references.csv":"short_references.csv");referenceout<<"case,h,t,x0,v0\n";
    auto cases=primary();if(longer){cases.clear();for(const auto&name:{"scalar_mixed","radial_orbit","rod_kick","elastic_mixed","chain_mixed"})cases.push_back(make_case(name));cases.push_back(hour_local::chain_case(3));}
    std::vector<std::string> orderlist=longer?std::vector<std::string>{"material_storage_kinetic","storage_kinetic_material","edge_reverse"}:schedules;
    std::vector<int> budgets=longer?std::vector<int>{8,32}:std::vector<int>{4,16};
    std::vector<double> hs=longer?std::vector<double>{.02,.08}:std::vector<double>{.02,.08};double duration=longer?5:1;
    for(const auto&c:cases)for(double h:hs){double E=scale(c);int steps=std::lround(duration/h);std::vector<State> refs{c.initial},fines{c.initial};double cert=0;
        for(int i=0;i<steps;++i){refs.push_back(reference_step(c.model,refs.back(),h,64));fines.push_back(reference_step(c.model,fines.back(),h,128));cert=std::max(cert,state_gap(c.model,refs.back(),fines.back(),E));}
        for(int i=0;i<=steps;++i)referenceout<<c.name<<','<<h<<','<<i*h<<','<<fines[i].x[0]<<','<<fines[i].v[0]<<'\n';
        auto rollout_methods=methods;rollout_methods.push_back("dense_oracle");
        for(const auto&method:rollout_methods)for(const auto&schedule:orderlist)for(int budget:budgets){if(method=="dense_oracle"&&(schedule!=orderlist.front()||budget!=budgets.front()))continue;int reported_budget=method=="dense_oracle"?0:budget;State state=c.initial;int accepted=0;std::string reason;bool failed=false;double maxenergy=0,stepenergy=0,px=0,pv=0,fx=0,fv=0,setupms=0,solvems=0;dual_gs::Residuals maxr;long long visits=0,evals=0,writes=0;
            for(int step=0;step<steps;++step){try{auto start=Clock::now();dual_gs::Problem p(c.model,state,h,E);setupms+=elapsed(start);dual_gs::Result r;if(method=="dense_oracle"){auto oracle=dual_gs::solve_oracle(c.model,state,h,E);r.z=oracle.z.size()?oracle.z:p.initial();r.y=oracle.y_lifted;r.milliseconds=oracle.milliseconds;if(!oracle.converged)r.failure="oracle_"+oracle.reason;else {static_cast<dual_gs::Residuals&>(r)=p.residuals(r.z,r.y);r.converged=true;}}else r=dual_gs::solve(p,method,schedule,budget*int(p.rows.size())+1,1,budget*int(p.rows.size()));solvems+=r.milliseconds;visits+=r.visits;evals+=r.row_evaluations;writes+=r.coordinate_writes;State next=p.state(r.z);
                if(!r.failure.empty()){failed=true;reason=r.failure;break;}if(!finite(next)){failed=true;reason="nonfinite_state";break;}double drift=std::abs(energy(c.model,next)-energy(c.model,c.initial))/E;
                if(drift>1e4||next.x.norm()>1e5||next.v.norm()>1e6){failed=true;reason="trajectory_runaway_guard";break;}
                state=next;++accepted;maxenergy=std::max(maxenergy,drift);stepenergy=std::max(stepenergy,std::abs(r.energy));maxr.material=std::max(maxr.material,r.material);maxr.kinetic=std::max(maxr.kinetic,r.kinetic);maxr.storage=std::max(maxr.storage,r.storage);maxr.geometry=std::max(maxr.geometry,r.geometry);maxr.tangency=std::max(maxr.tangency,r.tangency);maxr.stationarity=std::max(maxr.stationarity,r.stationarity);
                fx=(state.x-fines[step+1].x).norm()/std::sqrt(c.model.count);fv=std::sqrt((masses(c.model).array()*(state.v-fines[step+1].v).array().square()).sum()/E);px+=fx*fx;pv+=fv*fv;
                trajectory<<c.name<<','<<h<<','<<method<<','<<schedule<<",1,"<<reported_budget<<','<<(step+1)*h<<','<<drift<<','<<r.stationarity<<','<<fx<<','<<fv<<','<<state.x[0]<<','<<state.v[0]<<'\n';
            }catch(const std::exception&e){failed=true;reason=e.what();break;}}
            out<<c.name<<','<<c.model.count<<','<<h<<','<<method<<','<<schedule<<",1,"<<reported_budget<<','<<steps*h<<','<<accepted*h<<','<<accepted<<','<<failed<<','<<reason<<','<<maxenergy<<','<<stepenergy<<','<<maxr.material<<','<<maxr.kinetic<<','<<maxr.storage<<','<<maxr.geometry<<','<<maxr.tangency<<','<<maxr.stationarity<<','<<std::sqrt(px/std::max(1,accepted))<<','<<std::sqrt(pv/std::max(1,accepted))<<','<<fx<<','<<fv<<','<<cert<<','<<visits<<','<<evals<<','<<writes<<','<<setupms<<','<<solvems<<'\n';
        }
        std::cerr<<"rollouts "<<c.name<<" h="<<h<<" duration="<<duration<<'\n';
    }
}

static void scaling(const std::string&dir){auto out=output(dir,"scaling.csv");out<<"N,h,method,schedule,budget_multiplier,visits,row_evaluations,coordinate_writes,restoration_passes,setup_ms,solve_ms,full_residual,energy,failure\n";
    for(int N:{3,8,32}){auto c=hour_local::chain_case(N);double E=scale(c);for(double h:{.02,.08}){auto begin=Clock::now();dual_gs::Problem p(c.model,c.initial,h,E);double setup=elapsed(begin);for(const auto&method:methods)for(const auto&schedule:schedules)for(int budget:{1,4,16,64,256}){auto r=dual_gs::solve(p,method,schedule,budget*int(p.rows.size())+1,1,budget*int(p.rows.size()));out<<N<<','<<h<<','<<method<<','<<schedule<<','<<budget<<','<<r.visits<<','<<r.row_evaluations<<','<<r.coordinate_writes<<','<<r.stationarity_passes<<','<<setup<<','<<r.milliseconds<<','<<fullres(r)<<','<<r.energy<<','<<r.failure<<'\n';}}}
}

static void traces(const std::string&dir){auto out=output(dir,"visit_trace.csv");out<<"case,method,schedule,visit,row,kind,index,dual_value,x_middle,v_middle,material,kinetic,storage,stationarity,energy\n";auto c=hour_local::chain_case(3);dual_gs::Problem p(c.model,c.initial,.08,scale(c));for(const auto&method:methods)for(const auto&schedule:schedules){auto order=p.schedule_indices(schedule);for(int visit=1;visit<=int(order.size())*3;++visit){auto r=dual_gs::solve(p,method,schedule,100,1,visit);int row=order[(visit-1)%order.size()];auto state=p.state(r.z);out<<c.name<<','<<method<<','<<schedule<<','<<visit<<','<<row<<','<<dual_gs::kind_name(p.rows[row].kind)<<','<<p.rows[row].index<<','<<r.y[row]<<','<<state.x[2]<<','<<state.v[2]<<','<<r.material<<','<<r.kinetic<<','<<r.storage<<','<<r.stationarity<<','<<r.energy<<'\n';}}}

int main(int argc,char**argv){try{std::string mode=argc>1?argv[1]:"one",dir=argc>2?argv[2]:"results";std::filesystem::create_directories(dir);if(mode=="one")one_steps(dir,false);else if(mode=="controls")one_steps(dir,true);else if(mode=="short")rollouts(dir,false);else if(mode=="long")rollouts(dir,true);else if(mode=="scaling")scaling(dir);else if(mode=="trace")traces(dir);else throw std::runtime_error("unknown_mode");return 0;}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}
