#define main physical_benchmark_main
#include "physical.cpp"
#undef main
int main(){try{
    std::ofstream out("phase_refinement.csv"),rounded("phase_refinement_roundoff.csv");std::string header="h,steps,completed,omega_max,h_omega_max,reference_certified,reference_dx,reference_dv,position_error,velocity_error,mode_error,phase_error,energy_error,max_residual,max_tolerance\n";out<<std::setprecision(12)<<header;rounded<<std::setprecision(12)<<header;
    Problem p=problem(8,"both100","mixed",100);p.initial.x=p.equilibrium+.005*(p.initial.x-p.equilibrium);p.initial.v*=.005;p.escale=std::max(1e-12,energy(p.m,p.initial)-minimum_energy(p.m));
    for(double h:{.002,.001,.0005,.00025,.000125,.0000625,.00003125,.000015625}){int steps=int(std::round(.02/h));auto ref=physical_reference(p,h,steps);
      for(bool cancellation_scaled:{false,true}){State s=p.initial;Metric metric(p,h,"direct_block");int done=0;double maxr=0,maxe=0,maxtol=0;
        for(int k=0;k<steps;++k){double tol=1e-10;if(cancellation_scaled)tol=std::max(tol,64*std::numeric_limits<double>::epsilon()*(2*p.msqrt.cwiseProduct(s.x).norm()+h*p.msqrt.cwiseProduct(s.v).norm())/(h*std::sqrt(2*p.escale)));maxtol=std::max(maxtol,tol);Step r=solve(p,s,h,metric,0,predictor(p,s,h,"ballistic"),true,0,"transformed","fixed",tol);if(!r.converged||!r.finite)break;s=r.state;++done;maxr=std::max(maxr,r.residual);maxe=std::max(maxe,std::abs(energy(p.m,s)-energy(p.m,p.initial))/p.escale);}
        double dx=NAN,dv=NAN;ModeError mode;if(done==steps){dx=p.msqrt.cwiseProduct(s.x-ref.states.back().x).norm()/std::sqrt(p.mass.sum());dv=p.msqrt.cwiseProduct(s.v-ref.states.back().v).norm()/std::sqrt(p.mass.sum());mode=mode_error(p,s,ref.states.back());}
        auto&stream=cancellation_scaled?rounded:out;stream<<h<<','<<steps<<','<<done<<','<<p.frequencies.maxCoeff()<<','<<h*p.frequencies.maxCoeff()<<','<<ref.certified<<','<<ref.dx<<','<<ref.dv<<','<<dx<<','<<dv<<','<<(done==steps?mode.relative:NAN)<<','<<(done==steps?mode.phase:NAN)<<','<<maxe<<','<<maxr<<','<<maxtol<<'\n';stream.flush();
      }
    }return 0;
}catch(const std::exception&e){std::cerr<<e.what()<<std::endl;return 1;}}
