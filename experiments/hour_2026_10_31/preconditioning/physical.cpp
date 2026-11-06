#include "solver.hpp"
#include <fstream>
#include <iomanip>
#include <iostream>

struct PhysicalReference {std::vector<State> states;bool certified=false;double dx=INFINITY,dv=INFINITY,de=INFINITY;int substeps=0;};
PhysicalReference physical_reference(const Problem&p,double h,int steps){
    PhysicalReference out;double omega=std::sqrt(4*p.m.stiffness.maxCoeff()/p.m.mass.minCoeff());
    int substeps=std::max(8,int(std::ceil(h*omega/.08)));double vscale=std::sqrt(2*p.escale/p.mass.sum());
    auto integrate=[&](int count){std::vector<State> states{p.initial};for(int j=0;j<steps;++j){State s=reference_step(p.m,states.back(),h,count);if(!compose::finite(s))break;states.push_back(s);}return states;};
    auto coarse=integrate(substeps);
    for(int attempt=0;attempt<5;++attempt){substeps*=2;auto fine=integrate(substeps);out.dx=out.dv=out.de=0;out.substeps=substeps;out.states=fine;
        if(fine.size()!=size_t(steps+1)||coarse.size()!=fine.size())return out;
        for(int k=1;k<=steps;++k){out.dx=std::max(out.dx,p.msqrt.cwiseProduct(fine[k].x-coarse[k].x).norm()/std::sqrt(p.mass.sum()));out.dv=std::max(out.dv,p.msqrt.cwiseProduct(fine[k].v-coarse[k].v).norm()/std::sqrt(p.mass.sum())/vscale);out.de=std::max(out.de,std::abs(energy(p.m,fine[k])-energy(p.m,p.initial))/p.escale);}
        if(out.dx<1e-6&&out.dv<1e-5&&out.de<1e-6){out.certified=true;break;}coarse=std::move(fine);
    }return out;
}

int main(){try{
    Eigen::setNbThreads(1);std::ofstream out("physical.csv"),refs("physical_references.csv"),traces("branch_traces.csv");
    out<<std::setprecision(12)<<"N,regime,contrast,motion,h,predictor,requested,completed,reference_certified,reference_dx,reference_dv,rms_position,rms_velocity,max_mode,max_phase,max_energy,max_residual,cpu_ms\n";
    refs<<std::setprecision(12)<<"N,regime,contrast,motion,h,steps,substeps,certified,position_gap,velocity_gap_scaled,energy_gap,cpu_ms\n";
    traces<<std::setprecision(12)<<"N,regime,contrast,motion,h,predictor,time,residual,energy_error,position_error,velocity_error,mode_error,phase_error,tip_x,tip_y,reference_tip_x,reference_tip_y\n";
    struct Case {int N;std::string regime,motion;double contrast,h;};
    std::vector<Case> cases={{8,"both100","mixed",100,.04},{8,"both100","mixed",100,.02},{8,"both100","mixed",100,.01},{8,"both100","mixed",100,.005},{32,"both100","mixed",100,.04},{32,"both100","mixed",100,.01},{8,"stiff100","large",100,.1},{8,"both100","large",100,.1},{32,"both100","mixed",1000,.04}};
    for(auto c:cases){Problem p=problem(c.N,c.regime,c.motion,c.contrast);int steps=int(std::round(.4/c.h));double start=cpu_us();auto ref=physical_reference(p,c.h,steps);double refcpu=(cpu_us()-start)/1000;
        refs<<c.N<<','<<c.regime<<','<<c.contrast<<','<<c.motion<<','<<c.h<<','<<steps<<','<<ref.substeps<<','<<ref.certified<<','<<ref.dx<<','<<ref.dv<<','<<ref.de<<','<<refcpu<<'\n';refs.flush();
        for(std::string predictor_name:{"ballistic","acceleration","secant","continuation"}){State s=p.initial,prev=s;int done=0;double sumx=0,sumv=0,maxmode=0,maxphase=0,maxe=0,maxr=0,cpu=0;Metric metric(p,c.h,"direct_block");
            for(int k=0;k<steps;++k){start=cpu_us();Step r=predictor_name=="continuation"?continuation(p,s,c.h):solve(p,s,c.h,metric,0,predictor(p,s,c.h,predictor_name,k?&prev:nullptr),true);cpu+=cpu_us()-start;if(!r.converged||!r.finite)break;prev=s;s=r.state;++done;maxr=std::max(maxr,r.residual);double en=std::abs(energy(p.m,s)-energy(p.m,p.initial))/p.escale;maxe=std::max(maxe,en);
                const State&reference=ref.states.at(k+1);double dx=p.msqrt.cwiseProduct(s.x-reference.x).norm()/std::sqrt(p.mass.sum()),dv=p.msqrt.cwiseProduct(s.v-reference.v).norm()/std::sqrt(p.mass.sum());auto mode=mode_error(p,s,reference);sumx+=dx*dx;sumv+=dv*dv;maxmode=std::max(maxmode,mode.relative);maxphase=std::max(maxphase,mode.phase);
                traces<<c.N<<','<<c.regime<<','<<c.contrast<<','<<c.motion<<','<<c.h<<','<<predictor_name<<','<<(k+1)*c.h<<','<<r.residual<<','<<en<<','<<dx<<','<<dv<<','<<mode.relative<<','<<mode.phase<<','<<s.x[s.x.size()-2]<<','<<s.x[s.x.size()-1]<<','<<reference.x[reference.x.size()-2]<<','<<reference.x[reference.x.size()-1]<<'\n';
            }out<<c.N<<','<<c.regime<<','<<c.contrast<<','<<c.motion<<','<<c.h<<','<<predictor_name<<','<<steps<<','<<done<<','<<ref.certified<<','<<ref.dx<<','<<ref.dv<<','<<(done?std::sqrt(sumx/done):NAN)<<','<<(done?std::sqrt(sumv/done):NAN)<<','<<(done?maxmode:NAN)<<','<<(done?maxphase:NAN)<<','<<maxe<<','<<maxr<<','<<cpu/1000<<'\n';out.flush();traces.flush();
        }std::cout<<"physical N="<<c.N<<" "<<c.regime<<" contrast="<<c.contrast<<" h="<<c.h<<" certified="<<ref.certified<<" substeps="<<ref.substeps<<std::endl;
    }return 0;
}catch(const std::exception&e){std::cerr<<e.what()<<std::endl;return 1;}}
