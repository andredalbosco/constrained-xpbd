#include "solver.hpp"
#include "relinearized.hpp"
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <set>
using namespace traditional;
struct Case {std::string name;Model model;State initial;};
V2 direction(double angle){return {std::sin(angle),-std::cos(angle)};}
std::vector<Case> cases(){
    Model spring,rod,chain;rod.alpha[0]=0;rod.gravity={0,-9.81};chain.count=2;chain.gravity=rod.gravity;
    State release,orbit,rest,kick;release.x[0]={1.25,0};orbit.x[0]={1.1,0};orbit.v[0]={0,.8};
    rest.x[0]=direction(.65);kick.x[0]=direction(.4);kick.v[0]=tangent(kick.x[0])*1.2;
    auto double_state=[](double a,double b,double dr1,double dr2,double vr1,double vt1,double vr2,double vt2){
        State s;const auto n1=direction(a),n2=direction(b);
        s.x[0]=n1*(1.1962+dr1);s.x[1]=s.x[0]+n2*(1.0981+dr2);
        s.v[0]=n1*vr1+tangent(n1)*vt1;s.v[1]=s.v[0]+n2*vr2+tangent(n2)*vt2;return s;};
    return {{"spring_release",spring,release},{"spring_orbit",spring,orbit},{"pendulum_rest",rod,rest},{"pendulum_kick",rod,kick},
        {"double_small",chain,double_state(.15,-.1,0,0,0,0,0,0)},
        {"double_large",chain,double_state(.65,-.45,.1,.2,0,0,0,0)},
        {"double_mixed",chain,double_state(.45,-.3,.05,.1,-.2,.5,.3,.8)}};
}
State reference_step(const Model& m,State s,double h){
    if(m.count==1&&m.alpha[0]==0){
        const double l=m.length[0],a=std::atan2(s.x[0].x,-s.x[0].y),w=dot(tangent(s.x[0]),s.v[0])/(l*l);
        auto f=[&](double t){return -norm(m.gravity)*std::sin(t)/l;};
        const double a1=w,w1=f(a),a2=w+h*w1/2,w2=f(a+h*a1/2),a3=w+h*w2/2,w3=f(a+h*a2/2),a4=w+h*w3,w4=f(a+h*a3);
        const double an=a+h*(a1+2*a2+2*a3+a4)/6,wn=w+h*(w1+2*w2+2*w3+w4)/6;
        s.x[0]=direction(an)*l;s.v[0]=tangent(s.x[0])*wn;return s;
    }
    const Vec a1=acceleration(m,s.x);State s2=s,s3=s,s4=s;
    for(int i=0;i<m.count;++i){s2.x[i]+=s.v[i]*(h/2);s2.v[i]+=a1[i]*(h/2);}const Vec a2=acceleration(m,s2.x);
    for(int i=0;i<m.count;++i){s3.x[i]+=s2.v[i]*(h/2);s3.v[i]+=a2[i]*(h/2);}const Vec a3=acceleration(m,s3.x);
    for(int i=0;i<m.count;++i){s4.x[i]+=s3.v[i]*h;s4.v[i]+=a3[i]*h;}const Vec a4=acceleration(m,s4.x);
    for(int i=0;i<m.count;++i){s.x[i]+=(s.v[i]+s2.v[i]*2+s3.v[i]*2+s4.v[i])*(h/6);s.v[i]+=(a1[i]+a2[i]*2+a3[i]*2+a4[i])*(h/6);}return s;
}
std::vector<State> reference(const Case& c,double h,int steps,int subdivisions){
    std::vector<State> out;out.reserve(steps+1);out.push_back(c.initial);State s=c.initial;
    for(int j=0;j<steps;++j){for(int k=0;k<subdivisions;++k)s=reference_step(c.model,s,h/subdivisions);out.push_back(s);}return out;
}
struct Config{std::string name;Options options;bool refresh=false;};
std::vector<Config> configurations(){
    Options fresh;fresh.fresh=true;std::vector<Config> out{{"xpbd",fresh},{"frozen_xpbd",{}}};
    for(auto s:std::vector<std::pair<std::string,Schedule>>{{"EL",Schedule::EL},{"LE",Schedule::LE},{"jacobi",Schedule::Jacobi},{"joint",Schedule::Joint}})
        for(auto p:std::vector<std::pair<std::string,Policy>>{{"raw",Policy::Raw},{"direction",Policy::Direction},{"repair_cap",Policy::RepairCap}}){Options o;o.schedule=s.second;o.policy=p.second;out.push_back({s.first+"_"+p.first,o});}
    for(auto p:std::vector<std::pair<std::string,Policy>>{{"EL_raw_gate",Policy::RawEnergyGate},{"EL_raw_gate_cap",Policy::RawEnergyGateCap}}){Options o;o.schedule=Schedule::EL;o.policy=p.second;out.push_back({p.first,o});}
    for(auto s:std::vector<std::pair<std::string,Schedule>>{{"LEL_repair_cap",Schedule::LEL},{"ELE_repair_cap",Schedule::ELE}}){Options o;o.schedule=s.second;o.policy=Policy::RepairCap;out.push_back({s.first,o});}return out;
}
double angular(const Model& m,const State& s){double L=0;for(int i=0;i<m.count;++i)L+=m.mass[i]*cross(s.x[i],s.v[i]);return L;}
double torque(const Model& m,const State& s){double t=0;for(int i=0;i<m.count;++i)t+=m.mass[i]*cross(s.x[i],m.gravity);return t;}
struct Stats{
    std::string status="completed";int completed=0;double max_energy=0,rms_position=0,rms_velocity=0,reference_excursion=0,path=0,ref_path=0,
    max_angular=0,max_reference_angular=0,max_constraint=0,max_rod=0,final_energy=0,attempted_energy=0,impulse=0,ref_impulse=0;Diagnostics d;
};
void trace(std::ostream& out,const Case& c,const Config& cfg,int hz,int step,const State& s,const State& ref,const Iterate& z,const Stats& stats){
    out<<c.name<<','<<cfg.name<<','<<hz<<','<<step<<','<<double(step)/hz<<','<<energy(c.model,s)<<','<<energy(c.model,ref)<<','
       <<angular(c.model,s)<<','<<stats.impulse<<','<<angular(c.model,ref)<<','<<stats.ref_impulse<<','<<stats.d.last_blend;
    for(int i=0;i<2;++i)out<<','<<s.x[i].x<<','<<s.x[i].y<<','<<s.v[i].x<<','<<s.v[i].y<<','<<ref.x[i].x<<','<<ref.x[i].y<<','<<ref.v[i].x<<','<<ref.v[i].y;
    out<<','<<z.lambda[0]<<','<<z.lambda[1]<<','<<z.eta<<'\n';
}
Stats run(const Case& c,const Config& cfg,int hz,const std::vector<State>& ref,std::ostream* traces){
    Stats st;State s=c.initial;const auto& m=c.model;const double h=1./hz,E0=energy(m,s),es=E0-minimum_energy(m),L0=angular(m,s);
    if(traces)trace(*traces,c,cfg,hz,0,s,ref[0],Iterate{},st);
    for(int k=1;k<int(ref.size());++k){
        const Context ctx=context(m,s,h,E0);Iterate z{ctx.predictor};const bool ok=cfg.refresh?update_relinearized(ctx,z,cfg.options,st.d):update(ctx,z,cfg.options,st.d);
        State next;next.x=z.x;for(int i=0;i<m.count;++i)next.v[i]=(next.x[i]-s.x[i])/h;
        st.attempted_energy=energy(m,next);
        if(!ok||!valid(next.x,m.count)||!valid(next.v,m.count)||!std::isfinite(st.attempted_energy)){st.status="aborted_undefined";break;}
        if(mass_norm(m,next.x)>1e6||mass_norm(m,next.v)>1e6||std::abs(st.attempted_energy-E0)>1e12*es){st.status="aborted_safety_limit";break;}
        const Rows rows=evaluate(ctx,z);if(!rows.valid){st.status="aborted_undefined";break;}
        st.completed=k;st.max_energy=std::max(st.max_energy,std::abs(st.attempted_energy-E0)/es);
        for(int i=0;i<m.count;++i){const auto ex=next.x[i]-ref[k].x[i],ev=next.v[i]-ref[k].v[i],exc=ref[k].x[i]-c.initial.x[i];
            st.rms_position+=dot(ex,ex);st.rms_velocity+=dot(ev,ev);st.reference_excursion+=dot(exc,exc);
            st.path+=norm(next.x[i]-s.x[i]);st.ref_path+=norm(ref[k].x[i]-ref[k-1].x[i]);
            st.max_constraint=std::max(st.max_constraint,std::abs(rows.rC[i]));if(m.alpha[i]==0)st.max_rod=std::max(st.max_rod,std::abs(rows.C[i]));}
        st.impulse+=h*.5*(torque(m,s)+torque(m,next));st.ref_impulse+=h*.5*(torque(m,ref[k-1])+torque(m,ref[k]));
        st.max_angular=std::max(st.max_angular,std::abs(angular(m,next)-L0-st.impulse));
        st.max_reference_angular=std::max(st.max_reference_angular,std::abs(angular(m,ref[k])-L0-st.ref_impulse));
        if(traces)trace(*traces,c,cfg,hz,k,next,ref[k],z,st);s=next;
    }
    if(st.completed){st.rms_position=std::sqrt(st.rms_position/(st.completed*m.count));st.rms_velocity=std::sqrt(st.rms_velocity/(st.completed*m.count));st.reference_excursion=std::sqrt(st.reference_excursion/(st.completed*m.count));}
    st.final_energy=energy(m,s);return st;
}
int main(int argc,char** argv){try{
    const std::filesystem::path dir=argc>1?argv[1]:"results";const double duration=argc>2?std::stod(argv[2]):10;
    const bool el_le_focus=argc>3&&std::string(argv[3])=="--el-le-focus";
    const bool el_focus=el_le_focus||(argc>3&&std::string(argv[3])=="--el-focus");
    const std::set<std::string> el_methods=el_le_focus?
        std::set<std::string>{"xpbd","EL_raw","LE_raw","EL_repair_cap","LE_repair_cap","EL_refresh_raw","LE_refresh_raw","EL_refresh_repair_cap","LE_refresh_repair_cap"}:
        std::set<std::string>{"xpbd","EL_raw","EL_direction","EL_raw_gate","EL_raw_gate_cap","EL_repair_cap"};
    auto configs=configurations();
    if(el_le_focus)for(auto schedule:{Schedule::EL,Schedule::LE})for(auto policy:{Policy::Raw,Policy::RepairCap}){
        Options o;o.schedule=schedule;o.policy=policy;
        configs.push_back({std::string(schedule==Schedule::EL?"EL":"LE")+"_refresh_"+(policy==Policy::Raw?"raw":"repair_cap"),o,true});
    }
    if(!(duration>0&&duration<=1000))throw std::runtime_error("duration outside (0,1000]");std::filesystem::create_directories(dir);
    std::ofstream summary(dir/"summary.csv"),traces(dir/"traces.csv"),checks(dir/"reference_checks.csv"),meta(dir/"cases.csv");
    for(auto* out:{&summary,&traces,&checks,&meta}){if(!*out)throw std::runtime_error("cannot open output");*out<<std::setprecision(17);}
    summary<<"case,method,hz,h,duration,count,energy_scale,initial_energy,status,completed,steps,max_energy_rel,final_energy_rel,rms_position,rms_velocity,reference_excursion,path_ratio,max_angular_balance,max_reference_angular_balance,max_constraint_residual,max_rod_error,final_energy,attempted_energy,evaluations,scalar_solves,blocks,rank_fallbacks,singular,skips,caps,repair_scaled,max_extra_ratio\n";
    traces<<"case,method,hz,step,time,energy,reference_energy,angular,torque_impulse,reference_angular,reference_torque_impulse,blend,x1,y1,vx1,vy1,ref_x1,ref_y1,ref_vx1,ref_vy1,x2,y2,vx2,vy2,ref_x2,ref_y2,ref_vx2,ref_vy2,lambda1,lambda2,eta\n";
    checks<<"case,hz,subdivisions,refined_subdivisions,max_position_difference,max_velocity_difference,max_energy_rel\n";
    meta<<"case,count,m1,m2,l1,l2,alpha1,alpha2,gx,gy,x1,y1,vx1,vy1,x2,y2,vx2,vy2,initial_energy,energy_scale\n";
    const std::set<std::string> selected{"xpbd","frozen_xpbd","EL_raw_gate","EL_raw_gate_cap","EL_repair_cap","LE_repair_cap","joint_repair_cap"};
    for(const auto& c:cases()){
        if(el_focus&&c.name!="spring_release"&&c.name!="pendulum_rest"&&c.name!="spring_orbit"&&c.name!="pendulum_kick")continue;
        const auto& m=c.model;const double E0=energy(m,c.initial),es=E0-minimum_energy(m);
        meta<<c.name<<','<<m.count<<','<<m.mass[0]<<','<<m.mass[1]<<','<<m.length[0]<<','<<m.length[1]<<','<<m.alpha[0]<<','<<m.alpha[1]<<','<<m.gravity.x<<','<<m.gravity.y;
        for(int i=0;i<2;++i)meta<<','<<c.initial.x[i].x<<','<<c.initial.x[i].y<<','<<c.initial.v[i].x<<','<<c.initial.v[i].y;meta<<','<<E0<<','<<es<<'\n';
        for(int hz:{240,480,960}){
            const int steps=std::lround(duration*hz);const auto ref=reference(c,1./hz,steps,16),fine=reference(c,1./hz,steps,32);
            double dx=0,dv=0,de=0;for(int k=0;k<=steps;++k){for(int i=0;i<m.count;++i){dx=std::max(dx,norm(ref[k].x[i]-fine[k].x[i]));dv=std::max(dv,norm(ref[k].v[i]-fine[k].v[i]));}de=std::max(de,std::abs(energy(m,fine[k])-E0)/es);}
            checks<<c.name<<','<<hz<<",16,32,"<<dx<<','<<dv<<','<<de<<'\n';
            if(dx>1e-6||dv>1e-5||de>1e-7)throw std::runtime_error("reference refinement failed: "+c.name);
            for(const auto& cfg:configs){
                if(el_focus&&!el_methods.count(cfg.name))continue;
                const auto s=run(c,cfg,hz,fine,hz==480&&(el_focus||selected.count(cfg.name))?&traces:nullptr);const auto& d=s.d;
                summary<<c.name<<','<<cfg.name<<','<<hz<<','<<1./hz<<','<<duration<<','<<m.count<<','<<es<<','<<E0<<','<<s.status<<','<<s.completed<<','<<steps<<','<<s.max_energy<<','<<(s.final_energy-E0)/es<<','<<s.rms_position<<','<<s.rms_velocity<<','<<s.reference_excursion<<','<<(s.ref_path?s.path/s.ref_path:0)<<','<<s.max_angular<<','<<s.max_reference_angular<<','<<s.max_constraint<<','<<s.max_rod<<','<<s.final_energy<<','<<s.attempted_energy<<','<<d.evaluations<<','<<d.scalar_solves<<','<<d.blocks<<','<<d.rank_fallbacks<<','<<d.singular<<','<<d.skips<<','<<d.caps<<','<<d.repair_scaled<<','<<d.max_ratio<<'\n';
            }
        }
        std::cout<<c.name<<" complete"<<std::endl;
    }
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}return 0;}
