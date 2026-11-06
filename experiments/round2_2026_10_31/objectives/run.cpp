#include "objectives.hpp"
#include <fstream>
#include <iostream>
#include <iomanip>
#include <filesystem>
#include <chrono>
#include <map>
#include <omp.h>
using namespace objectives2;
namespace fs=std::filesystem;
inline double scalar_value(const Objective&p,double x,double v){double m=p.m.mass[0],k=p.m.stiffness[0],l=p.m.length[0],u=(x-p.old.x[0])/p.h,vn=p.old.v[0],a=p.aold[0],U=.5*k*(x-l)*(x-l),r=v+vn-2*u;
 if(p.variant.objective=="projection")return (.5*m*std::pow((x-p.target.x[0])/p.h,2)+.5*m*std::pow(v-p.target.v[0],2))/p.E;
 if(p.variant.objective.find("ls_")==0){double f=-k*(.5*(x+p.old.x[0])-l),rk=u-.5*(v+vn),rp=v-vn-p.h*f/m;return .5*m*(rk*rk+rp*rp)/p.E;}
 if(p.variant.objective=="original"){if(p.variant.kin=="soft_endpoint")r=v-u;return (.5*m*(u-vn)*(u-vn)+U+.5*m*r*r/p.cfg.beta)/p.E;}
 double aa=p.variant.objective=="quarter"?.25:(p.variant.objective=="half_half"?.5:0),b=p.variant.objective=="quarter"?.25:(p.variant.objective=="half_half"?.5:1);
 if(p.variant.objective=="midpotential")U=.5*k*std::pow(.5*(x+p.old.x[0])-l,2);return (.5*m*std::pow(u-vn-aa*p.h*a,2)+b*U+.5*m*r*r/p.cfg.beta)/p.E;
}
void validate(){double ge=0,he=0,identity=0;int checks=0;
 for(auto name:{"scalar_mixed","radial_orbit","rod_kick","chain_mixed"})for(auto&q:variants()){
  auto c=make_case(name);double E=std::max(1e-8,energy(c.model,c.initial)-minimum_energy(c.model));Objective p(c.model,c.initial,.025,q,E);V z=p.initial();for(int i=0;i<z.size();++i)z[i]+=.01*std::sin(i+1);auto e=p.evaluate(z);double eps=1e-5;
  for(int i=0;i<z.size();++i){V zp=z,zm=z;zp[i]+=eps;zm[i]-=eps;auto a=p.evaluate(zp),b=p.evaluate(zm);double g=(a.objective.a-b.objective.a)/(2*eps);V H=(a.objective.g-b.objective.g)/(2*eps);ge=std::max(ge,std::abs(g-e.objective.g[i])/(1+std::abs(g)));he=std::max(he,(H-e.objective.H.col(i)).norm()/(1+H.norm()));for(int j=0;j<int(e.rows.size());++j){double ng=(a.rows[j].value.a-b.rows[j].value.a)/(2*eps);V nh=(a.rows[j].value.g-b.rows[j].value.g)/(2*eps);ge=std::max(ge,std::abs(ng-e.rows[j].value.g[i])/(1+std::abs(ng)));he=std::max(he,(nh-e.rows[j].value.H.col(i)).norm()/(1+nh.norm()));}++checks;}
  if(c.model.dim==1)identity=std::max(identity,std::abs(scalar_value(p,e.x[0].a,e.v[0].a)-e.objective.a));
 }
 std::ofstream f("validation.csv");f<<"derivative_columns,max_relative_gradient,max_relative_hessian,scalar_value_AD_discrepancy\n"<<checks<<','<<std::setprecision(16)<<ge<<','<<he<<','<<identity<<'\n';if(ge>2e-5||he>2e-5||identity>1e-10)throw std::runtime_error("derivative_validation_failed");
}
void local_and_landscape(){std::ofstream out("local.csv"),dense("landscape.csv");out<<std::setprecision(17);dense<<std::setprecision(17);
 out<<"case,variant,hz,beta,energy_hard,angular_hard,converged,finite,position_error,velocity_error,energy_error,angular_error,acceleration_ratio,balance,stationarity,objective\n";
 dense<<"case,variant,hz,beta,min_x,min_v,min_objective,single_objective_gap,multistart_objective_gap,multistart_converged,reference_x,reference_v,min_position_error,min_velocity_error\n";
 for(auto name:{"scalar_release","scalar_mixed","radial_orbit","rod_release","elastic_mixed","chain_mixed"})for(int hz:{20,40,80,160,320})for(auto&q:variants()){
  auto c=make_case(name);double h=1./hz,E=std::max(1e-8,energy(c.model,c.initial)-minimum_energy(c.model));Objective p(c.model,c.initial,h,q,E);auto r=solve(p,p.initial(q.objective=="original"),64);auto ref=reference_step(c.model,c.initial,h,256);double ratio=std::string(name)=="scalar_release"?r.state.v[0]/ref.v[0]:0;
  out<<name<<','<<q.id<<','<<hz<<','<<p.cfg.beta<<','<<q.energy<<','<<q.angular<<','<<r.converged<<','<<r.finite<<','<<(r.state.x-ref.x).norm()<<','<<(r.state.v-ref.v).norm()<<','<<std::abs(energy(c.model,r.state)-energy(c.model,c.initial))/E<<','<<std::abs(angular(c.model,r.state)-angular(c.model,c.initial))<<','<<ratio<<','<<r.balance<<','<<r.stationarity<<','<<r.objective<<'\n';
  if(c.model.dim!=1||!q.energy)continue;
  double En=energy(c.model,c.initial),A=std::sqrt(2*En/c.model.stiffness[0]),B=std::sqrt(2*En/c.model.mass[0]),tau=2*std::acos(-1.);auto state=[&](double t){State s{V(1),V(1)};s.x[0]=c.model.length[0]+A*std::cos(t);s.v[0]=B*std::sin(t);return s;};auto f=[&](double t){auto s=state(t);return scalar_value(p,s.x[0],s.v[0]);};int n=65536,best=0;double value=f(0);for(int k=1;k<n;++k){double v=f(tau*k/n);if(v<value){value=v;best=k;}}double lo=tau*(best-1)/n,hi=tau*(best+1)/n;for(int k=0;k<80;++k){double a=lo+(hi-lo)/3,b=hi-(hi-lo)/3;if(f(a)<f(b))hi=b;else lo=a;}double t=.5*(lo+hi);State minimum=state(t);double minobj=f(t);
  std::vector<V> seeds={p.initial(false),p.initial(true),p.encode(c.initial)};State flipped=p.target;flipped.v=-flipped.v;seeds.push_back(p.encode(flipped));double bestobj=std::numeric_limits<double>::infinity();int ok=0;for(auto&z:seeds){auto rr=solve(p,z,96);if(rr.converged){++ok;bestobj=std::min(bestobj,rr.objective);}}
  dense<<name<<','<<q.id<<','<<hz<<','<<p.cfg.beta<<','<<minimum.x[0]<<','<<minimum.v[0]<<','<<minobj<<','<<r.objective-minobj<<','<<bestobj-minobj<<','<<ok<<','<<ref.x[0]<<','<<ref.v[0]<<','<<std::abs(minimum.x[0]-ref.x[0])<<','<<std::abs(minimum.v[0]-ref.v[0])<<'\n';
 }
}
struct Job{std::string name;Variant q;double duration;int hz=120;};
void trajectories(){std::vector<Job>jobs;
 // All objectives receive the common 10-second comparison. Long runs retain
 // bad controls as well as each distinct repair mechanism.
 for(auto name:case_names())for(auto&q:variants())jobs.push_back({name,q,10});
 std::vector<std::string>long_ids={"original_endpoint_b1","original_midpoint_b1","original_midpoint_h2","original_midpoint_h4","quarter_b1","midpotential_b1","projection_xpbd_E","projection_verlet_E","projection_verlet_EL","projection_midpoint_EL","ls_midpoint_E","ls_dg_E","ls_dg_EL"};
 for(auto name:case_names())for(auto&q:variants())if(std::find(long_ids.begin(),long_ids.end(),q.id)!=long_ids.end())jobs.push_back({name,q,60});
 std::ofstream out("summary.csv");out<<std::setprecision(16)<<"case,variant,hz,duration,beta_initial,status,completed,steps,reference_valid,reference_dx,reference_dv,max_energy_error,max_angular_balance,rms_position,rms_velocity,path_ratio,max_rod,max_tangency,max_balance,max_stationarity,unconverged_steps,compute_ms\n";fs::create_directories("traces");
 // References are cached once per physical case and duration, not per method.
 struct Ref{std::vector<State>s;double dx=0,dv=0;bool valid=true;};std::map<std::string,Ref>refs;
 for(auto&job:jobs){std::string key=job.name+std::to_string(int(job.duration));if(refs.count(key))continue;auto c=make_case(job.name);Ref ref;State fine=c.initial,coarse=fine;ref.s.push_back(fine);for(int k=0;k<int(job.duration*job.hz);++k){fine=reference_step(c.model,fine,1./job.hz,64);coarse=reference_step(c.model,coarse,1./job.hz,32);ref.dx=std::max(ref.dx,(fine.x-coarse.x).norm());ref.dv=std::max(ref.dv,(fine.v-coarse.v).norm());ref.s.push_back(fine);}ref.valid=ref.dx<1e-6&&ref.dv<1e-5;refs.emplace(key,std::move(ref));}
 int finished=0;omp_set_num_threads(8);
 #pragma omp parallel for schedule(dynamic,1)
 for(int ji=0;ji<int(jobs.size());++ji){auto&job=jobs[ji];auto c=make_case(job.name);auto&m=c.model;State s=c.initial;double h=1./job.hz,E0=energy(m,s),E=std::max(1e-8,E0-minimum_energy(m)),L0=angular(m,s),torques=0;const auto&ref=refs.at(job.name+std::to_string(int(job.duration)));double maxE=0,maxL=0,sx=0,sv=0,path=0,rpath=0,rod=0,tan=0,balance=0,station=0,ms=0;int complete=0,unconv=0;std::string status="complete";bool save=job.duration==10;
  std::ofstream trace;if(save){trace.open("traces/"+job.name+"__"+job.q.id+".csv");trace<<std::setprecision(16)<<"time,energy,angular";for(int i=0;i<dof(m);++i)trace<<",x"<<i<<",v"<<i<<",ref_x"<<i<<",ref_v"<<i;trace<<'\n';}
  for(int k=1;k<int(ref.s.size());++k){auto start=std::chrono::steady_clock::now();Solve r;try{r=step(m,s,h,job.q,E,24);}catch(const std::exception&e){r.finite=false;r.reason=e.what();}ms+=std::chrono::duration<double,std::milli>(std::chrono::steady_clock::now()-start).count();if(!r.finite){status="solver_failure";break;}double err=std::abs(energy(m,r.state)-E0)/E;if(!std::isfinite(err)||err>1e5||r.state.x.norm()>1e5||r.state.v.norm()>1e5){status="runaway";break;}++complete;if(!r.converged)++unconv;maxE=std::max(maxE,err);torques+=.5*h*(torque(m,s)+torque(m,r.state));maxL=std::max(maxL,std::abs(angular(m,r.state)-L0-torques));sx+=(r.state.x-ref.s[k].x).squaredNorm()/m.count;sv+=(r.state.v-ref.s[k].v).squaredNorm()/m.count;path+=(r.state.x-s.x).norm();rpath+=(ref.s[k].x-ref.s[k-1].x).norm();balance=std::max(balance,r.balance);station=std::max(station,r.stationarity);if(m.rigid){rod=std::max(rod,std::abs(r.state.x.norm()-m.length[0]));tan=std::max(tan,std::abs(r.state.x.dot(r.state.v)/r.state.x.norm()));}s=r.state;
   if(save&&k%2==0){trace<<k*h<<','<<energy(m,s)<<','<<angular(m,s);for(int i=0;i<dof(m);++i)trace<<','<<s.x[i]<<','<<s.v[i]<<','<<ref.s[k].x[i]<<','<<ref.s[k].v[i];trace<<'\n';}
  }
  #pragma omp critical
  {out<<job.name<<','<<job.q.id<<','<<job.hz<<','<<job.duration<<','<<beta_for(job.q,m,h)<<','<<status<<','<<complete<<','<<int(job.duration*job.hz)<<','<<ref.valid<<','<<ref.dx<<','<<ref.dv<<','<<maxE<<','<<maxL<<','<<(complete?std::sqrt(sx/complete):0)<<','<<(complete?std::sqrt(sv/complete):0)<<','<<(rpath?path/rpath:0)<<','<<rod<<','<<tan<<','<<balance<<','<<station<<','<<unconv<<','<<ms<<'\n';out.flush();if(++finished%20==0)std::cout<<"trajectories "<<finished<<'/'<<jobs.size()<<std::endl;}
 }
}
int main(int argc,char**argv){try{std::string mode=argc>1?argv[1]:"all";if(mode=="--help"){std::cout<<"objective_study [--validate-only|--local-only|--trajectories-only]\nNo option runs all three phases. CSVs and traces are written in the working directory.\n";return 0;}if(mode!="all"&&mode!="--validate-only"&&mode!="--local-only"&&mode!="--trajectories-only")throw std::runtime_error("unknown_mode");if(mode=="all"||mode=="--validate-only"){validate();std::cout<<"derivatives validated"<<std::endl;}if(mode=="all"||mode=="--local-only"){local_and_landscape();std::cout<<"local and dense scalar studies complete"<<std::endl;}if(mode=="all"||mode=="--trajectories-only")trajectories();}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}
