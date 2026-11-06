#include "cases.hpp"
#include "stability.hpp"
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
using namespace agenda;
int main(int argc,char** argv){try{
    const std::filesystem::path dir=argc>1?argv[1]:"results";std::filesystem::create_directories(dir);std::ofstream out(dir/"stability.csv");out<<std::setprecision(17);
    out<<"case,method,hz,stiffness_factor,light_mass_ratio,perturbation,duration,status,completed,steps,max_energy_ratio,max_absolute_energy_error,attempts,substeps,rejections,deepest\n";
    int guarded_failures=0,unguarded_growth=0,guarded_runs=0;
    for(const auto& original:cases())if(original.name=="spring_release"||original.name=="spring_orbit"||original.name=="pendulum_rest"||original.name=="double_large")
    for(int hz:{240,960})for(double stiffness:{1.,100.,10000.})for(double ratio:{1.,100.})for(double perturbation:{0.,1e-8})for(int method=0;method<3;++method){
        auto test=original;for(int i=0;i<test.model.count;++i)if(test.model.alpha[i]>0)test.model.alpha[i]/=stiffness;test.model.mass[test.model.count-1]/=ratio;
        const double c=std::cos(perturbation),sn=std::sin(perturbation);for(int i=0;i<test.model.count;++i){const auto x=test.initial.x[i],v=test.initial.v[i];test.initial.x[i]={c*x.x-sn*x.y,sn*x.x+c*x.y};test.initial.v[i]={c*v.x-sn*v.y,sn*v.x+c*v.y};}
        const auto& m=test.model;const double E0=energy(m,test.initial),minimum=minimum_energy(m),scale=E0-minimum,h=1./hz;State state=test.initial;const int steps=hz;int completed=0;std::string status="completed";double max_ratio=1,max_error=0;StableReport d;
        for(int k=1;k<=steps;++k){State next;bool ok=true;try{if(method==0){const auto ctx=context(m,state,h,E0);next=endpoint(ctx,baseline(ctx).x);}else if(method==1){Report r;next=midpoint_step(context(m,state,h,E0),4,true,r,scale);}else ok=guarded_midpoint(m,state,h,E0,scale,next,d);}catch(const std::exception&){ok=false;}
            if(!ok||!valid(next.x,m.count)||!valid(next.v,m.count)||!std::isfinite(energy(m,next))){status=method==2?"retry_budget_exhausted":"nonfinite";break;}
            const double ratioE=(energy(m,next)-minimum)/scale;max_ratio=std::max(max_ratio,ratioE);max_error=std::max(max_error,std::abs(ratioE-1));
            if(ratioE>100){status="energy_growth_limit";break;}state=next;completed=k;
        }
        if(method==2){++guarded_runs;if(status!="completed")++guarded_failures;}if(method==1&&(status!="completed"||max_ratio>1.01))++unguarded_growth;
        out<<test.name<<','<<(method==0?"xpbd":method==1?"midpoint_4":"guarded_midpoint")<<','<<hz<<','<<stiffness<<','<<ratio<<','<<perturbation<<",1,"<<status<<','<<completed<<','<<steps<<','<<max_ratio<<','<<max_error<<','<<d.attempts<<','<<d.accepted<<','<<d.rejected<<','<<d.deepest<<'\n';
    }
    std::cout<<"Guarded runs "<<guarded_runs<<", retry failures "<<guarded_failures<<"; fixed-four-pass runs with >1% energy growth or failure "<<unguarded_growth<<'\n';
    std::ofstream equivalent(dir/"acceptance_equivalence.csv");equivalent<<std::setprecision(17)<<"case,hz,duration,steps,rejections,max_state_difference,max_rx,max_constitutive\n";
    for(const auto& test:cases()){State state=test.initial;const double target=energy(test.model,state),scale=target-minimum_energy(test.model);StableReport d;double difference=0,rx=0,material=0;const int steps=480*60;
        for(int k=0;k<steps;++k){Report r;const auto fixed=midpoint_step(context(test.model,state,1./480,target),4,true,r,scale);State accepted;
            if(!guarded_midpoint(test.model,state,1./480,target,scale,accepted,d))throw std::runtime_error("nominal guard retry failed");
            for(int i=0;i<test.model.count;++i){difference=std::max(difference,norm(fixed.x[i]-accepted.x[i]));difference=std::max(difference,norm(fixed.v[i]-accepted.v[i]));}rx=std::max(rx,r.rx);material=std::max(material,r.constitutive);state=accepted;}
        equivalent<<test.name<<",480,60,"<<steps<<','<<d.rejected<<','<<difference<<','<<rx<<','<<material<<'\n';
        if(d.rejected||difference)throw std::runtime_error("nominal four-pass trajectory changed under acceptance checks");
    }
    std::cout<<"All seven nominal 60-second trajectories exactly unchanged by guard checks\n";
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
