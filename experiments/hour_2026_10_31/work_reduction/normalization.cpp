#include "reduction.hpp"
#include <fstream>
#include <filesystem>
#include <iomanip>
#include <iostream>
using namespace work_reduction;
namespace fs=std::filesystem;
// Exact symmetry reduction of N disconnected identical copies with total-E
// normalization: per-copy gamma becomes gamma/N. Compensating gamma by N
// before total normalization restores the original per-copy objective.
int main(int argc,char**argv){
 fs::path out=argc>1?argv[1]:"normalization_results";fs::create_directories(out);
 std::ofstream csv(out/"trajectories.csv");csv<<std::setprecision(16)<<"case,target,copies,normalization,effective_gamma,h,steps,nonconverged,energy_error,position_rms,reference_gap,evaluations\n";
 for(std::string name:{"radial_orbit","elastic_mixed","chain_mixed"})for(std::string target:{"frozen_power","implicit_midpoint_dg"})for(int copies:{1,2,8,32})for(std::string norm:copies==1?std::vector<std::string>{"local"}:std::vector<std::string>{"total","compensated"}){
  auto c=make_case(name);State s=c.initial,ref=s,ref2=s;double h=1./30,E0=energy(c.model,s),E=std::max(1e-6,E0-minimum_energy(c.model)),maxE=0,sumx=0,gap=0;int bad=0,eval=0;
  Config cfg=configuration(target);cfg.gamma=norm=="total"?1./copies:1;
  std::ofstream trace(out/(name+"_"+target+"_"+std::to_string(copies)+"_"+norm+".csv"));trace<<std::setprecision(16)<<"t";for(int k=0;k<dof(c.model);++k)trace<<",x"<<k;trace<<'\n';
  for(int i=0;i<150;++i){CrossProblem full(c.model,s,h,cfg,E,target,true);Reduced p(full);V q=p.midpoint_initial();auto a=joint_solve(p,128,&q);s=p.state(a.z);if(!finite(s))return 2;if(!a.converged)++bad;eval+=a.evaluations;maxE=std::max(maxE,std::abs(energy(c.model,s)-E0)/E);ref=reference_step(c.model,ref,h,128);ref2=reference_step(c.model,ref2,h,256);sumx+=(s.x-ref2.x).squaredNorm()/c.model.count;gap=std::max(gap,(ref.x-ref2.x).norm()+(ref.v-ref2.v).norm()/std::sqrt(E));trace<<(i+1)*h;for(int k=0;k<dof(c.model);++k)trace<<','<<s.x[k];trace<<'\n';}
  csv<<name<<','<<target<<','<<copies<<','<<norm<<','<<cfg.gamma<<','<<h<<",150,"<<bad<<','<<maxE<<','<<std::sqrt(sumx/150)<<','<<gap<<','<<eval<<'\n';csv.flush();std::cout<<name<<' '<<target<<' '<<copies<<' '<<norm<<'\n';
 }
}
