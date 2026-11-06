#include "frozen_schedules.hpp"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

using namespace energy_ordering;
namespace {
enum class Relax { WholeSweep, EnergyOnly };
struct Options {
    Schedule order;
    Relax relax;
    double omega;
    bool half_seed,retain;
    int sweeps;
};
bool update(const Context& c,Iterate& z,const Options& o,Counters& counts) {
    if(o.relax==Relax::WholeSweep) {
        const Iterate before=z;
        if(!schedule_update(c,z,o.order,o.retain,counts))return false;
        if(o.omega!=1) {
            z.x=before.x+(z.x-before.x)*o.omega;
            z.lambda=before.lambda+(z.lambda-before.lambda)*o.omega;
            z.eta=before.eta+(z.eta-before.eta)*o.omega;
        }
        return finite(z.x)&&std::isfinite(z.lambda)&&std::isfinite(z.eta);
    }
    const Rows r=evaluate(c,z);++counts.evaluations;
    if(!r.valid)return false;
    const double m=c.model.mass;
    const V2 base=o.retain?r.R*(-1/m):V2{};
    const double fL=-r.rC-dot(r.n,base),fE=-r.rE-dot(r.g,base);
    double dl=0,de=0;
    auto material=[&]() {dl=(fL-r.B*de)/r.A;++counts.scalar_solves;};
    auto energy=[&]() {
        if(energy_row_zero(r)) {
            ++counts.singular;
            if(std::abs(fE-r.B*dl)>energy_tolerance(c,r))++counts.incompatible;
        } else {
            const double solved=(fE-r.B*dl)/r.D;
            de+=o.omega*(solved-de);
            ++counts.scalar_solves;
        }
    };
    if(o.order==Schedule::LE||o.order==Schedule::LEL) {
        material();energy();if(o.order==Schedule::LEL)material();
    } else {
        energy();material();if(o.order==Schedule::ELE)energy();
    }
    z.x+=base+(r.n*dl+r.g*de)/m;z.lambda+=dl;z.eta+=de;
    return finite(z.x)&&std::isfinite(z.lambda)&&std::isfinite(z.eta);
}
State exact(double v0,double t) {
    return {{1+.25*std::cos(10*t)+(v0/10)*std::sin(10*t),0},
            {-2.5*std::sin(10*t)+v0*std::cos(10*t),0}};
}
void validate() {
    Model m;m.mass=1.7;m.length=1.2;m.alpha=.035;m.acceleration={.2,-.3};
    const State old{{1.1,.7},{.15,-.12}};
    const Context c=context(m,old,.08,.4,energy(m,old)-.03);
    const Iterate initial{{1.06,.76},-.014,.003};
    for(Schedule s:{Schedule::LE,Schedule::EL,Schedule::LEL,Schedule::ELE})for(bool retain:{false,true}) {
        Iterate a=initial,b=initial;Counters ca,cb;
        Options o{s,Relax::EnergyOnly,1,false,retain,1};
        if(!update(c,a,o,ca)||!schedule_update(c,b,s,retain,cb)
           ||norm(a.x-b.x)>1e-13||std::abs(a.lambda-b.lambda)>1e-13||std::abs(a.eta-b.eta)>1e-13)
            throw std::runtime_error("energy-only omega=1 parity failed");
    }
}
void run(std::ostream& out,double v0,const Options& o,double h=1./60,int steps=600,bool annotate=false) {
    Model m;State state{{1.25,0},{v0,0}};
    const double E0=energy(m,state);
    double rmsx=0,rmsv=0,maxH=0,maxC=0,maxR=0,path=0,refpath=0;
    double minx=state.x.x,maxx=state.x.x;
    long frozen=0;int completed=0;std::string status="completed";
    Counters counts;
    for(int step=1;step<=steps;++step) {
        const Context c=context(m,state,h,0,E0);Iterate z{c.predictor};
        if(o.half_seed) {
            const V2 n=state.x/norm(state.x);
            const V2 force=n*(-(norm(state.x)-m.length)/m.alpha);
            z.x+=force*(.5*h*h/m.mass);
        }
        bool ok=true;
        for(int i=0;i<o.sweeps;++i) {
            ok=update(c,z,o,counts);
            const State attempt{z.x,(z.x-state.x)/h};
            const double E=energy(m,attempt);
            if(!ok||!finite(z.x)||!std::isfinite(E)) {status="undefined";ok=false;break;}
            if(norm(attempt.x)>1e6||norm(attempt.v)>1e6||std::abs(E-E0)>1e12*E0
                ||std::abs(z.lambda)>1e18||std::abs(z.eta)>1e18) {
                status="guarded";ok=false;break;
            }
        }
        if(!ok)break;
        const State next{z.x,(z.x-state.x)/h};const State ref=exact(v0,step*h),prevref=exact(v0,(step-1)*h);
        const Rows r=evaluate(c,z);if(!r.valid){status="undefined";break;}
        ++completed;
        const double dx=norm(next.x-state.x);
        path+=dx;refpath+=norm(ref.x-prevref.x);
        if(dx<1e-10)++frozen;
        rmsx+=dot(next.x-ref.x,next.x-ref.x);rmsv+=dot(next.v-ref.v,next.v-ref.v);
        maxH=std::max(maxH,std::abs(energy(m,next)-E0));maxC=std::max(maxC,std::abs(r.rC));
        maxR=std::max(maxR,norm(r.R)/(h*h));minx=std::min(minx,next.x.x);maxx=std::max(maxx,next.x.x);
        state=next;
    }
    if(completed){rmsx=std::sqrt(rmsx/completed);rmsv=std::sqrt(rmsv/completed);}
    if(status=="completed"&&counts.incompatible)status="completed_incompatible";
    out<<(v0==0?"rest":"radial")<<','<<name(o.order)<<','
       <<(o.relax==Relax::WholeSweep?"whole_sweep":"energy_only")<<','<<o.omega<<','
       <<(o.half_seed?"half_force_seed":"ordinary_seed")<<','<<o.retain<<','<<o.sweeps<<','
       <<completed<<','<<status<<','<<rmsx<<','<<rmsv<<','<<maxH/E0<<','<<maxC<<','<<maxR<<','
       <<path<<','<<(refpath>0?path/refpath:0)<<','<<minx<<','<<maxx<<','<<frozen<<','
       <<counts.singular<<','<<counts.incompatible;
    if(annotate)out<<','<<h<<','<<steps<<','<<steps*h;
    out<<'\n';
}
} // namespace
int main(int argc,char** argv) {
    try {
        validate();
        std::ofstream out(argc>1?argv[1]:"results/spring_repair/existing_controls.csv");
        if(!out)throw std::runtime_error("cannot open result CSV");
        out<<std::setprecision(17);
        const std::string header="case,order,relaxation,omega,initialization,retain_rx,sweeps,completed,status,rms_position,rms_velocity,max_relative_energy_error,"
              "max_material_residual,max_rx,path_length,path_ratio,min_x,max_x,frozen_steps,singular_updates,incompatible_updates";
        out<<header<<'\n';
        int runs=0;
        for(double v0:{0.,-.7})for(Schedule s:{Schedule::LE,Schedule::EL,Schedule::LEL,Schedule::ELE})
        for(Relax relax:{Relax::WholeSweep,Relax::EnergyOnly})for(double omega:{.1,.25,.5,1.}) {
            if(relax==Relax::EnergyOnly&&omega==1)continue;
            for(bool half:{false,true})for(bool retain:{false,true})for(int sweeps:{5,20}) {
                run(out,v0,{s,relax,omega,half,retain,sweeps});++runs;
            }
        }
        std::cout<<"Completed "<<runs<<" control runs; omega=1 parity checks passed\n";
        std::ofstream targeted(argc>2?argv[2]:"results/spring_repair/existing_controls_targeted.csv");
        if(!targeted)throw std::runtime_error("cannot open targeted result CSV");
        targeted<<std::setprecision(17)<<header<<",h,requested_steps,duration\n";
        int targeted_runs=0;
        for(double v0:{0.,-.7})for(double h:{1./30,1./60,1./120,1./240})
            for(int sweeps:{5,20,80})for(double duration:{10.,100.}) {
                run(targeted,v0,{Schedule::EL,Relax::EnergyOnly,.5,true,false,sweeps},
                    h,static_cast<int>(std::lround(duration/h)),true);
                ++targeted_runs;
            }
        for(double v0:{0.,-.7})for(double duration:{10.,100.}) {
            run(targeted,v0,{Schedule::EL,Relax::EnergyOnly,.5,true,true,20},
                1./60,static_cast<int>(std::lround(duration*60)),true);
            ++targeted_runs;
        }
        std::cout<<"Completed "<<targeted_runs<<" targeted consistency/long-run controls\n";
    } catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}
