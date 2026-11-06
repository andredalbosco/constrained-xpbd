#include "cases.hpp"
#include "../traditional_cases/relinearized.hpp"
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <set>
using namespace agenda;
enum class Family {XPBD,Ledger,Legacy,Endpoint,Work,Midpoint};
struct Config {std::string name;Family family=Family::XPBD;bool EL=true,keep_rx=false;int dual=0,iterations=1;double cap=.1;WorkOptions work;};
std::vector<Config> configurations(){
    std::vector<Config> out;
    auto add=[&](std::string name,Family family,bool EL=true){Config c;c.name=name;c.family=family;c.EL=EL;out.push_back(c);};
    add("xpbd",Family::XPBD);add("xpbd_ledger",Family::Ledger);add("EL_bounded",Family::Legacy);add("LE_bounded",Family::Legacy,false);
    for(bool EL:{true,false}){
        for(double cap:{0.,.1}){Config c;c.name=std::string(EL?"EL":"LE")+"_rx_"+(cap?"bounded":"raw");c.family=Family::Endpoint;c.EL=EL;c.keep_rx=true;c.cap=cap;out.push_back(c);}
        for(int dual:{0,1,2}){Config c;c.name=std::string(EL?"EL":"LE")+(dual==0?"_refresh_bounded":dual==1?"_dual_complete":"_dual_approx");c.family=Family::Endpoint;c.EL=EL;c.dual=dual;out.push_back(c);}
    }
    auto work=[&](std::string name,Estimate estimate,double strength,double cap,int sweeps,double repay=0,bool kfirst=false){Config c;c.name=name;c.family=Family::Work;c.work={estimate,strength,cap,repay,sweeps,kfirst};out.push_back(c);};
    work("work_lagged_weak",Estimate::Lagged,.1,.1,1);
    work("work_lagged",Estimate::Lagged,1,.1,1);
    work("work_lagged_strong",Estimate::Lagged,10,.1,1);
    work("work_lagged_raw",Estimate::Lagged,1,0,1);
    work("work_column",Estimate::Column,1,.1,1);
    work("work_row",Estimate::Row,1,.1,1);
    work("work_secant",Estimate::Secant,1,.1,1);
    work("work_lagged_GS4",Estimate::Lagged,1,.1,4);
    work("work_lagged_GS16",Estimate::Lagged,1,.1,16);
    work("work_kinetic_first",Estimate::Lagged,1,.1,1,0,true);
    work("reservoir_fast",Estimate::Lagged,1,.1,1,.25);
    work("reservoir_slow",Estimate::Lagged,1,.1,1,1);
    for(int iterations:{1,4,8}){Config c;c.name="midpoint_"+std::to_string(iterations);c.family=Family::Midpoint;c.keep_rx=true;c.iterations=iterations;out.push_back(c);}
    Config drop;drop.name="midpoint_drop_rx";drop.family=Family::Midpoint;drop.iterations=4;out.push_back(drop);return out;
}
double angular(const Model& m,const State& s){double a=0;for(int i=0;i<m.count;++i)a+=m.mass[i]*cross(s.x[i],s.v[i]);return a;}
double torque(const Model& m,const State& s){double a=0;for(int i=0;i<m.count;++i)a+=m.mass[i]*cross(s.x[i],m.gravity);return a;}
struct Stats {int completed=0;std::string status="completed";double max_energy=0,rms=0,max_angular=0,max_rod=0,path=0,refpath=0,impulse=0,max_reservoir=0,max_local_reservoir=0,max_ledger=0,max_cumulative_ledger=0,max_work=0,final_energy=0,attempted_energy=0,milliseconds=0;Report r;};
State step(const Context& c,const Config& cfg,double scale,Reservoir& reservoir,Report& r){
    if(cfg.family==Family::Work)return work_step(c,cfg.work,scale,reservoir,r);
    if(cfg.family==Family::Midpoint)return midpoint_step(c,cfg.iterations,cfg.keep_rx,r,scale);
    if(cfg.family==Family::Endpoint)return energy_step(c,cfg.EL,cfg.keep_rx,cfg.dual,cfg.cap,r);
    if(cfg.family==Family::XPBD||cfg.family==Family::Ledger){const auto z=baseline(c);r.lambda=z.lambda;const auto next=endpoint(c,z.x);
        if(cfg.family==Family::Ledger){const auto graph=layout(c.model,c.old);const auto v=work_values(c,z.x,Estimate::Lagged);const double before=energy(c.model,c.old)+total(reservoir);
            for(int i=0;i<c.model.count;++i){double b=v.dk[i].v;for(int e=0;e<graph.ne;++e)if(graph.edge[e].body==i)b-=v.w[e].v;reservoir.k[i]-=b;r.work_defect=std::max(r.work_defect,std::abs(b)/scale);}
            for(int p=0;p<graph.np;++p){double b=v.du[p].v;for(int e=0;e<graph.ne;++e)if(graph.edge[e].potential==p)b+=v.w[e].v;reservoir.u[p]-=b;r.work_defect=std::max(r.work_defect,std::abs(b)/scale);}
            r.ledger_error=std::abs(energy(c.model,next)+total(reservoir)-before)/scale;}
        return next;}
    Options o;o.schedule=cfg.EL?Schedule::EL:Schedule::LE;o.policy=Policy::RepairCap;Diagnostics d;Iterate z{c.predictor};if(!update(c,z,o,d))throw std::runtime_error("legacy failed");r.lambda=z.lambda;r.nu=z.eta;r.visits=d.scalar_solves;r.evaluations=d.evaluations;r.caps=d.caps;return endpoint(c,z.x);
}
void trace(std::ostream& out,const Case& c,const Config& cfg,int hz,int k,const State& s,const State& ref,const Reservoir& reservoir,const Stats& st){
    out<<c.name<<','<<cfg.name<<','<<hz<<','<<k<<','<<double(k)/hz<<','<<energy(c.model,s)<<','<<energy(c.model,ref)<<','<<total(reservoir)<<','<<angular(c.model,s)<<','<<st.impulse;
    for(int i=0;i<2;++i)out<<','<<s.x[i].x<<','<<s.x[i].y<<','<<s.v[i].x<<','<<s.v[i].y<<','<<ref.x[i].x<<','<<ref.x[i].y<<','<<ref.v[i].x<<','<<ref.v[i].y;
    for(double v:reservoir.k)out<<','<<v;for(double v:reservoir.u)out<<','<<v;
    out<<','<<st.r.work_defect<<','<<st.r.ledger_error<<'\n';
}
Stats run(const Case& c,const Config& cfg,int hz,const std::vector<State>& ref,std::ostream* traces,int trace_stride=1){
    Stats st;State s=c.initial;Reservoir reservoir;const auto& m=c.model;const double h=1./hz,E0=energy(m,s),scale=E0-minimum_energy(m),L0=angular(m,s);if(traces)trace(*traces,c,cfg,hz,0,s,ref[0],reservoir,st);
    const auto start=std::chrono::steady_clock::now();
    for(int k=1;k<int(ref.size());++k){Report report;State next;Reservoir pending=reservoir;const auto cxt=context(m,s,h,E0);
        try{next=step(cxt,cfg,scale,pending,report);}catch(const std::exception&){st.status="aborted_invalid";break;}
        st.attempted_energy=energy(m,next);
        if(!valid(next.x,m.count)||!valid(next.v,m.count)||!std::isfinite(st.attempted_energy)){st.status="aborted_nonfinite";break;}
        if(mass_norm(m,next.x)>1e4||mass_norm(m,next.v)>1e5||std::abs(st.attempted_energy-E0)>1e8*scale){st.status="aborted_safety_limit";break;}
        st.completed=k;st.max_energy=std::max(st.max_energy,std::abs(st.attempted_energy-E0)/scale);
        for(int i=0;i<m.count;++i){const auto e=next.x[i]-ref[k].x[i];st.rms+=dot(e,e);st.path+=norm(next.x[i]-s.x[i]);st.refpath+=norm(ref[k].x[i]-ref[k-1].x[i]);
            const double C=norm(next.x[i]-(i?next.x[i-1]:V2{}))-m.length[i];if(m.alpha[i]==0)st.max_rod=std::max(st.max_rod,std::abs(C));}
        st.impulse+=h*.5*(torque(m,s)+torque(m,next));st.max_angular=std::max(st.max_angular,std::abs(angular(m,next)-L0-st.impulse));
        st.max_reservoir=std::max(st.max_reservoir,std::abs(total(pending))/scale);st.max_ledger=std::max(st.max_ledger,report.ledger_error);st.max_work=std::max(st.max_work,report.work_defect);
        double local=0;for(double v:pending.k)local+=std::abs(v);for(double v:pending.u)local+=std::abs(v);st.max_local_reservoir=std::max(st.max_local_reservoir,local/scale);
        if(cfg.family==Family::Work||cfg.family==Family::Ledger)st.max_cumulative_ledger=std::max(st.max_cumulative_ledger,std::abs(st.attempted_energy+total(pending)-E0)/scale);
        st.r.visits+=report.visits;st.r.evaluations+=report.evaluations;st.r.caps+=report.caps;st.r.singular+=report.singular;st.r.work_defect=report.work_defect;st.r.ledger_error=report.ledger_error;
        st.r.max_extra_ratio=std::max(st.r.max_extra_ratio,report.max_extra_ratio);st.r.kinematic=std::max(st.r.kinematic,report.kinematic);st.r.rx=std::max(st.r.rx,report.rx);
        reservoir=pending;s=next;if(traces&&(k%trace_stride==0||k+1==int(ref.size())))trace(*traces,c,cfg,hz,k,s,ref[k],reservoir,st);
    }
    st.milliseconds=std::chrono::duration<double,std::milli>(std::chrono::steady_clock::now()-start).count();if(st.completed)st.rms=std::sqrt(st.rms/(st.completed*m.count));st.final_energy=energy(m,s);return st;
}
int main(int argc,char** argv){try{
    const std::filesystem::path dir=argc>1?argv[1]:"results";const double duration=argc>2?std::stod(argv[2]):10;const bool long_run=argc>3&&std::string(argv[3])=="--long";
    std::filesystem::create_directories(dir);std::ofstream summary(dir/"summary.csv"),traces(dir/"traces.csv"),meta(dir/"cases.csv"),checks(dir/"reference_checks.csv");
    for(auto* o:{&summary,&traces,&meta,&checks}){if(!*o)throw std::runtime_error("output open failed");*o<<std::setprecision(17);}
    summary<<"case,method,hz,duration,count,energy_scale,initial_energy,status,completed,steps,max_energy_rel,final_energy_rel,rms_position,max_angular_balance,max_rod_error,path_ratio,max_reservoir_rel,max_local_reservoir_rel,max_ledger_error,max_cumulative_ledger_error,max_work_defect,attempted_energy,visits,evaluations,caps,singular,max_extra_ratio,max_rx,max_kinematic_error,wall_ms,reference_valid,trace_stride\n";
    traces<<"case,method,hz,step,time,energy,reference_energy,reservoir,angular,torque_impulse,x1,y1,vx1,vy1,ref_x1,ref_y1,ref_vx1,ref_vy1,x2,y2,vx2,vy2,ref_x2,ref_y2,ref_vx2,ref_vy2,RK1,RK2,RU1,RU2,RU3,RU4,work_defect,ledger_error\n";
    meta<<"case,count,m1,m2,l1,l2,alpha1,alpha2,gx,gy,initial_energy,energy_scale\n";checks<<"case,hz,max_position_difference,max_velocity_difference,max_energy_rel,trajectory_valid\n";
    const std::set<std::string> selected{"xpbd","xpbd_ledger","EL_bounded","LE_bounded","EL_rx_bounded","LE_rx_bounded","work_lagged","work_column","work_row","work_secant","work_lagged_GS16","reservoir_fast","reservoir_slow","midpoint_1","midpoint_4","midpoint_drop_rx"};
    for(const auto& c:cases()){const auto& m=c.model;const double E0=energy(m,c.initial),scale=E0-minimum_energy(m);meta<<c.name<<','<<m.count<<','<<m.mass[0]<<','<<m.mass[1]<<','<<m.length[0]<<','<<m.length[1]<<','<<m.alpha[0]<<','<<m.alpha[1]<<','<<m.gravity.x<<','<<m.gravity.y<<','<<E0<<','<<scale<<'\n';
        for(int hz:long_run?std::vector<int>{480}:std::vector<int>{240,480,960}){
            const int steps=std::lround(duration*hz);const auto ref=reference(c,1./hz,steps,16),fine=reference(c,1./hz,steps,32);double dx=0,dv=0,de=0;
            for(int k=0;k<=steps;++k){for(int i=0;i<m.count;++i){dx=std::max(dx,norm(ref[k].x[i]-fine[k].x[i]));dv=std::max(dv,norm(ref[k].v[i]-fine[k].v[i]));}de=std::max(de,std::abs(energy(m,fine[k])-E0)/scale);}
            const bool reference_valid=dx<=1e-6&&dv<=1e-5&&de<=1e-7;checks<<c.name<<','<<hz<<','<<dx<<','<<dv<<','<<de<<','<<reference_valid<<'\n';if(!reference_valid&&!long_run)throw std::runtime_error("reference refinement failed");
            if(!reference_valid)std::cout<<c.name<<": long-horizon trajectory reference not certified; RMS diagnostic only"<<std::endl;
            for(const auto& cfg:configurations()){if(long_run&&!selected.count(cfg.name))continue;const auto st=run(c,cfg,hz,fine,hz==480&&selected.count(cfg.name)?&traces:nullptr,long_run?8:1);
                summary<<c.name<<','<<cfg.name<<','<<hz<<','<<duration<<','<<m.count<<','<<scale<<','<<E0<<','<<st.status<<','<<st.completed<<','<<steps<<','<<st.max_energy<<','<<(st.final_energy-E0)/scale<<','<<st.rms<<','<<st.max_angular<<','<<st.max_rod<<','<<(st.refpath?st.path/st.refpath:0)<<','<<st.max_reservoir<<','<<st.max_local_reservoir<<','<<st.max_ledger<<','<<st.max_cumulative_ledger<<','<<st.max_work<<','<<st.attempted_energy<<','<<st.r.visits<<','<<st.r.evaluations<<','<<st.r.caps<<','<<st.r.singular<<','<<st.r.max_extra_ratio<<','<<st.r.rx<<','<<st.r.kinematic<<','<<st.milliseconds<<','<<reference_valid<<','<<(long_run?8:1)<<'\n';
            }
        }std::cout<<c.name<<" complete"<<std::endl;
    }
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
