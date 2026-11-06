#include "common.hpp"
#include "independent.hpp"
#include "mechanics.hpp"
#include "el.hpp"
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <mutex>
#include <omp.h>
using namespace compose;
namespace fs=std::filesystem;
struct Job {Config c;std::string tier,case_name;int hz=120;double duration=10,kfactor=1,mratio=1,perturb=1;bool trace=false;};
std::vector<std::string> split(const std::string&s){std::vector<std::string> a;std::stringstream f(s);std::string t;while(std::getline(f,t,',')){if(!t.empty()&&t.back()=='\r')t.pop_back();a.push_back(t);}return a;}
std::vector<Job> read_jobs(const std::string&path){std::ifstream f(path);if(!f)throw std::runtime_error("cannot_open_manifest");std::string line;std::getline(f,line);auto h=split(line);std::map<std::string,int> cols;for(int i=0;i<int(h.size());++i)cols[h[i]]=i;std::vector<Job> jobs;
    while(std::getline(f,line)){if(line.empty())continue;auto a=split(line);auto s=[&](std::string k){return a.at(cols.at(k));};auto d=[&](std::string k){return std::stod(s(k));};auto n=[&](std::string k){return std::stoi(s(k));};Job j;auto&c=j.c;c.id=s("id");j.tier=s("tier");j.case_name=s("case");c.family=s("family");c.work=s("work");c.treatment=s("treatment");c.kinematics=s("kinematics");c.schedule=s("schedule");c.curvature=s("curvature");c.beta=d("beta");c.gamma=d("gamma");c.iterations=n("iterations");c.tangency=n("tangency");c.keep_rx=n("keep_rx");c.relinearize=n("relinearize");c.reverse=n("reverse");j.hz=n("hz");j.duration=d("duration");j.kfactor=d("stiffness_factor");j.mratio=d("mass_ratio");j.perturb=d("perturbation");j.trace=n("trace");if(j.hz<=0||j.duration<=0||c.iterations<1)throw std::runtime_error("invalid_parameters");jobs.push_back(j);}return jobs;
}
std::string key(const Job&j){std::ostringstream o;o<<std::setprecision(17)<<j.case_name<<'_'<<j.hz<<'_'<<j.duration<<'_'<<j.kfactor<<'_'<<j.mratio<<'_'<<j.perturb;return o.str();}
struct Reference {std::vector<State> state;double position=0,velocity=0,energy_error=0;bool certified=false;int subdivisions=0;};
Reference reference(const Case&c,int hz,int steps,double scale){Reference r;const double h=1./hz;double omega=std::sqrt(4*c.model.stiffness.maxCoeff()/c.model.mass.minCoeff());r.subdivisions=std::max(16,int(std::ceil(h*omega/.02)));State fine=c.initial,coarse=c.initial;r.state.push_back(fine);double E=energy(c.model,fine);
    for(int k=0;k<steps;++k){fine=reference_step(c.model,fine,h,2*r.subdivisions);coarse=reference_step(c.model,coarse,h,r.subdivisions);r.position=std::max(r.position,(fine.x-coarse.x).norm());r.velocity=std::max(r.velocity,(fine.v-coarse.v).norm());r.energy_error=std::max(r.energy_error,std::abs(energy(c.model,fine)-E)/scale);r.state.push_back(fine);}r.certified=r.position<1e-6&&r.velocity<1e-5&&r.energy_error<1e-7;return r;}
StepResult advance(const Model&m,const State&s,double h,const Config&c,double scale){if(c.family=="independent")return independent_step(m,s,h,c,scale);if(c.family=="EL")return el_step(m,s,h,c,scale);return mechanics_step(m,s,h,c,scale);}
struct Stats {std::string status="completed",reason="none";int completed=0;long visits=0,evals=0,iterations=0;double maxE=0,finalE=0,rmsx=0,rmsv=0,path=0,refpath=0,rod=0,tangent=0,angular_error=0,balance=0,rx=0,kinematic=0,impulse=0,work=0,compute_ms=0,firstcross=-1,refcross=-1,maxspeed=0,maxacc=0,objective_steps=0;};
void trace_header(std::ostream&o,int n){o<<"time,energy,reference_energy";for(int i=0;i<n;++i)o<<",x"<<i<<",v"<<i<<",ref_x"<<i<<",ref_v"<<i;o<<"\n";}
void trace_row(std::ostream&o,const Model&m,double t,const State&s,const State&r){o<<t<<','<<energy(m,s)<<','<<energy(m,r);for(int i=0;i<dof(m);++i)o<<','<<s.x[i]<<','<<s.v[i]<<','<<r.x[i]<<','<<r.v[i];o<<'\n';}
Stats run(const Job&j,const Case&c,const Reference&ref,const fs::path&out){Stats t;State s=c.initial;const Model&m=c.model;const double h=1./j.hz,E0=energy(m,s),scale=std::max(1e-10,E0-minimum_energy(m));double L0=angular(m,s),torque_integral=0,sumsqx=0,sumsqv=0;std::ofstream trace;
    if(j.trace){trace.open(out/"traces"/(j.c.id+".csv"));trace<<std::setprecision(16);trace_header(trace,dof(m));trace_row(trace,m,0,s,ref.state[0]);}
    for(int k=1;k<int(ref.state.size());++k){StepResult z;auto start=std::chrono::steady_clock::now();try{z=advance(m,s,h,j.c,scale);}catch(const std::exception&e){z.valid=false;z.reason=e.what();}t.compute_ms+=std::chrono::duration<double,std::milli>(std::chrono::steady_clock::now()-start).count();
        t.visits+=z.visits;t.evals+=z.evaluations;t.iterations+=z.iterations;if(!z.valid){t.status="solver_failed";t.reason=z.reason;break;}const State&next=z.state;
        if(next.x.size()!=dof(m)||next.v.size()!=dof(m)||!finite(next)){t.status="nonfinite";t.reason="invalid_state";break;}double E=energy(m,next),err=std::abs(E-E0)/scale;
        if(err>1e6||next.x.norm()>1e6||next.v.norm()>1e6){t.status="runaway";t.reason="safety_limit";t.maxE=std::max(t.maxE,err);break;}
        ++t.completed;t.maxE=std::max(t.maxE,err);t.finalE=(E-E0)/scale;sumsqx+=(next.x-ref.state[k].x).squaredNorm();sumsqv+=(next.v-ref.state[k].v).squaredNorm();t.path+=(next.x-s.x).norm();t.refpath+=(ref.state[k].x-ref.state[k-1].x).norm();
        t.maxspeed=std::max(t.maxspeed,next.v.norm());t.maxacc=std::max(t.maxacc,(next.v-s.v).norm()/h);
        if(m.rigid){V e=edge(m,next.x,0);t.rod=std::max(t.rod,std::abs(e.norm()-m.length[0]));t.tangent=std::max(t.tangent,std::abs(e.dot(next.v)/e.norm()));}
        torque_integral+=.5*h*(torque(m,s)+torque(m,next));t.angular_error=std::max(t.angular_error,std::abs(angular(m,next)-L0-torque_integral));
        t.balance=std::max(t.balance,std::abs(z.balance));t.rx=std::max(t.rx,std::abs(z.rx));t.kinematic=std::max(t.kinematic,std::abs(z.kinematic));t.impulse=std::max(t.impulse,std::abs(z.impulse));t.work=std::max(t.work,std::abs(z.work_defect));
        if(m.dim==1){double a=s.x[0]-m.length[0],b=next.x[0]-m.length[0];if(t.firstcross<0&&a>0&&b<=0)t.firstcross=(k-1+a/(a-b))*h;double ar=ref.state[k-1].x[0]-m.length[0],br=ref.state[k].x[0]-m.length[0];if(t.refcross<0&&ar>0&&br<=0)t.refcross=(k-1+ar/(ar-br))*h;}
        s=next;if(j.trace&&(k%std::max(1,j.hz/60)==0||k+1==int(ref.state.size())))trace_row(trace,m,k*h,s,ref.state[k]);
    }
    if(t.completed){t.rmsx=std::sqrt(sumsqx/(t.completed*m.count));t.rmsv=std::sqrt(sumsqv/(t.completed*m.count));}return t;
}
std::string summary(const Job&j,const Case&c,const Reference&r,const Stats&t){std::ostringstream o;o<<std::setprecision(16);auto cfg=j.c;std::string reason=t.reason;std::replace(reason.begin(),reason.end(),',',';');std::replace(reason.begin(),reason.end(),'\n',' ');
    o<<cfg.id<<','<<j.tier<<','<<j.case_name<<','<<cfg.family<<','<<cfg.work<<','<<cfg.treatment<<','<<cfg.kinematics<<','<<cfg.schedule<<','<<cfg.curvature<<','<<cfg.beta<<','<<cfg.gamma<<','<<cfg.iterations<<','<<cfg.tangency<<','<<cfg.keep_rx<<','<<cfg.relinearize<<','<<cfg.reverse<<','<<j.hz<<','<<j.duration<<','<<j.kfactor<<','<<j.mratio<<','<<j.perturb<<','<<t.status<<','<<reason<<','<<t.completed<<','<<int(std::lround(j.duration*j.hz))<<','<<std::max(1e-10,energy(c.model,c.initial)-minimum_energy(c.model))<<','<<t.maxE<<','<<t.finalE<<','<<t.rmsx<<','<<t.rmsv<<','<<(t.refpath?t.path/t.refpath:0)<<','<<t.rod<<','<<t.tangent<<','<<t.angular_error<<','<<t.balance<<','<<t.rx<<','<<t.kinematic<<','<<t.impulse<<','<<t.work<<','<<t.visits<<','<<t.evals<<','<<t.iterations<<','<<t.compute_ms<<','<<t.firstcross<<','<<t.refcross<<','<<t.maxspeed<<','<<t.maxacc<<','<<r.certified<<','<<r.position<<','<<r.velocity<<','<<r.energy_error<<','<<j.trace;return o.str();}
int main(int argc,char**argv){try{std::string manifest,out,filter;int threads=4,limit=0;for(int i=1;i<argc;++i){std::string a=argv[i];auto next=[&](){if(++i>=argc)throw std::runtime_error("missing_argument");return std::string(argv[i]);};if(a=="--manifest")manifest=next();else if(a=="--out")out=next();else if(a=="--threads")threads=std::stoi(next());else if(a=="--limit")limit=std::stoi(next());else if(a=="--case")filter=next();else if(a=="--help"){std::cout<<"study --manifest jobs.csv --out DIRECTORY [--threads N] [--limit N] [--case NAME]\nAll physical/solver parameters are explicit in the CSV manifest. Deterministic; no random initialization.\n";return 0;}else throw std::runtime_error("unknown_argument_"+a);}if(manifest.empty()||out.empty()||threads<1)throw std::runtime_error("need_manifest_out_positive_threads");
    fs::create_directories(fs::path(out)/"traces");auto jobs=read_jobs(manifest);std::map<std::string,std::vector<Job>> groups;int accepted=0;for(auto&j:jobs){if(!filter.empty()&&j.case_name!=filter)continue;if(limit&&accepted>=limit)break;groups[key(j)].push_back(j);++accepted;}
    std::ofstream results(fs::path(out)/"summary.csv"),refs(fs::path(out)/"reference_checks.csv");results<<"id,tier,case,family,work,treatment,kinematics,schedule,curvature,beta,gamma,iterations,tangency,keep_rx,relinearize,reverse,hz,duration,stiffness_factor,mass_ratio,perturbation,status,reason,completed,steps,energy_scale,max_energy_rel,final_energy_rel,rms_position,rms_velocity,path_ratio,max_rod_error,max_tangency,max_angular_balance,max_balance,max_rx,max_kinematic,max_impulse,max_work_defect,visits,evaluations,total_iterations,compute_ms,first_crossing,reference_crossing,max_speed,max_acceleration,reference_valid,reference_dx,reference_dv,reference_energy,trace\n";
    refs<<"key,subdivisions,fine_subdivisions,max_dx,max_dv,max_energy_rel,certified\n";omp_set_num_threads(threads);int done=0;std::mutex output_mutex;
    for(const auto&entry:groups){const auto&k=entry.first;const auto&group=entry.second;auto j=group[0];auto c=make_case(j.case_name,j.kfactor,j.mratio,j.perturb);double scale=std::max(1e-10,energy(c.model,c.initial)-minimum_energy(c.model));auto ref=reference(c,j.hz,int(std::lround(j.duration*j.hz)),scale);refs<<std::setprecision(16)<<k<<','<<ref.subdivisions<<','<<2*ref.subdivisions<<','<<ref.position<<','<<ref.velocity<<','<<ref.energy_error<<','<<ref.certified<<'\n';refs.flush();
        #pragma omp parallel for schedule(dynamic,1)
        for(int i=0;i<int(group.size());++i){auto t=run(group[i],c,ref,out);auto row=summary(group[i],c,ref,t);std::lock_guard<std::mutex> lock(output_mutex);results<<row<<'\n';results.flush();++done;if(done%25==0||done==accepted)std::cout<<"completed "<<done<<'/'<<accepted<<" last="<<group[i].c.id<<" status="<<t.status<<std::endl;}
    }std::cout<<"finished "<<done<<" configurations"<<std::endl;
}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}
