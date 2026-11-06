#include "solver.hpp"
#include "../research_agenda/cases.hpp"
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <set>
using namespace implicit;
const std::set<std::string> selected{"xpbd","xpbd_4","EL_soft_10","LE_soft_10","EL_angular_10","LE_angular_10","work_soft_10","work_angular","work_repay","work_repay_angular","workL_s1_tau005","workL_s1_tau025","workL_s1_tau100","workL_s10_tau005","workL_s10_tau100","workL_s100_tau005","workL_s100_tau025","workL_s100_tau100","workL_soft_01","workL_soft_1","workL_soft_10"};
struct Stats {int completed=0,fallbacks=0;long visits=0,eval=0,backtracks=0;std::string status="completed";double energy=0,rms=0,maxL=0,localL=0,rod=0,path=0,refpath=0,impulse=0,ledger=0,blend=0,ratio=0,reservoir=0;};
void trace(std::ostream& out,const Case& c,const Config& cfg,int hz,int k,const State& s,const State& ref,const Reservoir& reservoir,const StepReport& r,const Stats& stats){
    out<<c.name<<','<<cfg.name<<','<<hz<<','<<k<<','<<double(k)/hz<<','<<energy(c.model,s)<<','<<total(reservoir)<<','<<angular(c.model,s)<<','<<stats.impulse<<','<<r.blend<<','<<r.fallback;
    for(int i=0;i<2;++i)out<<','<<s.x[i].x<<','<<s.x[i].y<<','<<s.v[i].x<<','<<s.v[i].y<<','<<ref.x[i].x<<','<<ref.x[i].y<<','<<ref.v[i].x<<','<<ref.v[i].y;out<<'\n';
}
Stats run(const Case& c,const Config& cfg,int hz,const std::vector<State>& ref,std::ostream* traces,int stride=1){
    const auto& m=c.model;State s=c.initial;Reservoir reservoir;Stats st;const double h=1./hz,E0=energy(m,s),scale=E0-minimum_energy(m),L0=angular(m,s);StepReport last;if(traces)trace(*traces,c,cfg,hz,0,s,ref[0],reservoir,last,st);
    for(int k=1;k<int(ref.size());++k){State next;StepReport r;const auto ctx=context(m,s,h,E0);
        if(!advance(ctx,cfg,scale,reservoir,next,r)){st.status="cannot_advance_safely";break;}
        st.completed=k;st.energy=std::max(st.energy,std::abs(energy(m,next)-E0)/scale);st.localL=std::max(st.localL,std::abs(r.angular));
        st.impulse+=h*.5*(torque(m,s.x)+torque(m,next.x));st.maxL=std::max(st.maxL,std::abs(angular(m,next)-L0-st.impulse));
        for(int i=0;i<m.count;++i){const V2 e=next.x[i]-ref[k].x[i];st.rms+=dot(e,e);st.path+=norm(next.x[i]-s.x[i]);st.refpath+=norm(ref[k].x[i]-ref[k-1].x[i]);if(m.alpha[i]==0)st.rod=std::max(st.rod,std::abs(norm(next.x[i]-(i?next.x[i-1]:V2{}))-m.length[i]));}
        st.fallbacks+=r.fallback;st.visits+=r.visits;st.eval+=r.evaluations;st.backtracks+=r.backtracks;st.ledger=std::max(st.ledger,r.ledger);st.blend+=r.blend;st.ratio=std::max(st.ratio,r.cap_ratio);st.reservoir=std::max(st.reservoir,std::abs(total(reservoir))/scale);
        s=next;last=r;if(traces&&(k%stride==0||k+1==int(ref.size())))trace(*traces,c,cfg,hz,k,s,ref[k],reservoir,r,st);
    }
    if(traces&&st.status!="completed"&&st.completed%stride)trace(*traces,c,cfg,hz,st.completed,s,ref[st.completed],reservoir,last,st);
    if(st.completed){st.rms=std::sqrt(st.rms/(st.completed*m.count));st.blend/=st.completed;}return st;
}
void stress(const std::filesystem::path& dir,double duration){
    std::ofstream out(dir/"stability.csv");out<<std::setprecision(17)<<"case,method,hz,stiffness_factor,light_mass_ratio,perturbation,duration,completed,steps,status,max_energy_gain,fallbacks,backtracks,mean_blend\n";
    int count=0,failed=0;
    for(const auto& original:cases()){
    if(original.name!="spring_release"&&original.name!="spring_orbit"&&original.name!="pendulum_rest"&&original.name!="double_large")continue;
    for(int hz:{240,960})for(double stiffness:{1.,100.,10000.})for(double ratio:{1.,100.})for(double perturbation:{0.,1e-8})for(const auto& cfg:configurations())if(selected.count(cfg.name)){
        auto c=original;for(int i=0;i<c.model.count;++i)if(c.model.alpha[i]>0)c.model.alpha[i]/=stiffness;c.model.mass[c.model.count-1]/=ratio;
        const double ca=std::cos(perturbation),sa=std::sin(perturbation);for(int i=0;i<c.model.count;++i){const auto x=c.initial.x[i],v=c.initial.v[i];c.initial.x[i]={ca*x.x-sa*x.y,sa*x.x+ca*x.y};c.initial.v[i]={ca*v.x-sa*v.y,sa*v.x+ca*v.y};}
        State s=c.initial;Reservoir reservoir;const double E0=energy(c.model,s),scale=E0-minimum_energy(c.model);int completed=0,fallback=0;long backtracks=0;double gain=0,blend=0;std::string status="completed";
        const int steps=std::lround(hz*duration);
        for(int k=1;k<=steps;++k){State next;StepReport r;if(!advance(context(c.model,s,1./hz,E0),cfg,scale,reservoir,next,r)){status="cannot_advance_safely";break;}
            gain=std::max(gain,(energy(c.model,next)-E0)/scale);fallback+=r.fallback;backtracks+=r.backtracks;blend+=r.blend;s=next;completed=k;}
        if(status!="completed")++failed;++count;out<<c.name<<','<<cfg.name<<','<<hz<<','<<stiffness<<','<<ratio<<','<<perturbation<<','<<duration<<','<<completed<<','<<steps<<','<<status<<','<<gain<<','<<fallback<<','<<backtracks<<','<<(completed?blend/completed:0)<<'\n';
    }
    }
    std::cout<<"Stability: "<<count<<" runs, "<<failed<<" unable to advance safely\n";
}
int main(int argc,char** argv){try{
    const std::filesystem::path dir=argc>1?argv[1]:"results";const double duration=argc>2?std::stod(argv[2]):10;const bool long_run=argc>3&&std::string(argv[3])=="--long";std::filesystem::create_directories(dir);
    if(argc>3&&std::string(argv[3])=="--stress"){stress(dir,duration);return 0;}
    std::ofstream summary(dir/"summary.csv"),traces(dir/"traces.csv"),meta(dir/"cases.csv"),checks(dir/"reference_checks.csv");for(auto* out:{&summary,&traces,&meta,&checks}){if(!*out)throw std::runtime_error("output open failed");*out<<std::setprecision(17);}
    summary<<"case,method,hz,duration,count,energy_scale,initial_energy,status,completed,steps,max_energy_rel,rms_position,max_angular_balance,max_local_angular,max_rod_error,path_ratio,fallbacks,backtracks,mean_blend,max_extra_ratio,scalar_visits,evaluations,max_ledger_error,max_reservoir_rel,reference_valid,trace_stride\n";
    traces<<"case,method,hz,step,time,energy,reservoir,angular,torque_impulse,blend,fallback,x1,y1,vx1,vy1,ref_x1,ref_y1,ref_vx1,ref_vy1,x2,y2,vx2,vy2,ref_x2,ref_y2,ref_vx2,ref_vy2\n";
    meta<<"case,count,m1,m2,l1,l2,alpha1,alpha2,gx,gy,initial_energy,energy_scale\n";checks<<"case,hz,max_position_difference,max_velocity_difference,max_energy_error,reference_valid\n";
    for(const auto& c:cases()){const auto& m=c.model;const double E0=energy(m,c.initial),scale=E0-minimum_energy(m);meta<<c.name<<','<<m.count<<','<<m.mass[0]<<','<<m.mass[1]<<','<<m.length[0]<<','<<m.length[1]<<','<<m.alpha[0]<<','<<m.alpha[1]<<','<<m.gravity.x<<','<<m.gravity.y<<','<<E0<<','<<scale<<'\n';
        for(int hz:long_run?std::vector<int>{480}:std::vector<int>{240,480,960}){const int steps=std::lround(hz*duration);const auto ref=reference(c,1./hz,steps,16),fine=reference(c,1./hz,steps,32);double dx=0,dv=0,de=0;
            for(int k=0;k<=steps;++k){for(int i=0;i<m.count;++i){dx=std::max(dx,norm(ref[k].x[i]-fine[k].x[i]));dv=std::max(dv,norm(ref[k].v[i]-fine[k].v[i]));}de=std::max(de,std::abs(energy(m,fine[k])-E0)/scale);}
            const bool reference_valid=dx<=1e-6&&dv<=1e-5&&de<=1e-7;checks<<c.name<<','<<hz<<','<<dx<<','<<dv<<','<<de<<','<<reference_valid<<'\n';if(!reference_valid&&!long_run)throw std::runtime_error("reference refinement failed");
            for(const auto& cfg:configurations()){if(long_run&&!selected.count(cfg.name))continue;const auto st=run(c,cfg,hz,fine,hz==480&&selected.count(cfg.name)?&traces:nullptr,long_run?8:1);
                summary<<c.name<<','<<cfg.name<<','<<hz<<','<<duration<<','<<m.count<<','<<scale<<','<<E0<<','<<st.status<<','<<st.completed<<','<<steps<<','<<st.energy<<','<<st.rms<<','<<st.maxL<<','<<st.localL<<','<<st.rod<<','<<(st.refpath?st.path/st.refpath:0)<<','<<st.fallbacks<<','<<st.backtracks<<','<<st.blend<<','<<st.ratio<<','<<st.visits<<','<<st.eval<<','<<st.ledger<<','<<st.reservoir<<','<<reference_valid<<','<<(long_run?8:1)<<'\n';}
        }std::cout<<c.name<<" complete\n";
    }
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
