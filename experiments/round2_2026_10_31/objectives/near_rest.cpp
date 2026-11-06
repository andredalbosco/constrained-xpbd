#include "objectives.hpp"
#include <fstream>
#include <iomanip>
#include <iostream>
using namespace objectives2;
int main(){std::ofstream out("near_rest.csv");out<<std::setprecision(17)<<"variant,perturbation,hz,beta,steps,completed,unconverged,max_relative_energy,max_relative_state_error,max_balance,max_stationarity,status\n";
 for(double perturb:{1.,1e-2,1e-4,1e-6,0.})for(int hz:{80,320})for(auto&q:variants()){
  if(q.id!="original_midpoint_h2"&&q.id!="original_midpoint_h4"&&q.id!="quarter_b1"&&q.id!="midpotential_b1"&&q.id!="ls_dg_E"&&q.id!="projection_verlet_E")continue;
  auto c=make_case("scalar_release",1,1,perturb);State s=c.initial,ref=s;double E0=energy(c.model,s),scale=std::max(E0,1e-20),maxE=0,maxerr=0,bal=0,stat=0;int completed=0,unc=0;std::string status="complete";
  for(int k=0;k<hz;++k){auto r=step(c.model,s,1./hz,q,scale,32);if(!r.finite){status="solver_failure";break;}s=r.state;ref=reference_step(c.model,ref,1./hz,64);maxE=std::max(maxE,std::abs(energy(c.model,s)-E0)/scale);maxerr=std::max(maxerr,((s.x-ref.x).norm()+.1*(s.v-ref.v).norm())/std::max(.25*perturb,1e-10));bal=std::max(bal,r.balance);stat=std::max(stat,r.stationarity);unc+=!r.converged;++completed;}
  out<<q.id<<','<<perturb<<','<<hz<<','<<beta_for(q,c.model,1./hz)<<','<<hz<<','<<completed<<','<<unc<<','<<maxE<<','<<maxerr<<','<<bal<<','<<stat<<','<<status<<'\n';
 }
 std::cout<<"60 near-rest trajectories completed\n";
}
