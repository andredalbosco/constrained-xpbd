#include "numerics.hpp"
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
using namespace hour_local;
int main(int argc,char**argv){std::filesystem::path out=argc>1?argv[1]:"results";std::filesystem::create_directories(out);std::ofstream csv(out/"rod_reaction_control.csv");csv<<std::setprecision(16)<<"reaction,h,duration,steps,inexact_steps,max_balance,max_stationarity,max_energy_error,rms_position,final_position_error,final_velocity_error\n";
 for(bool reaction:{false,true})for(int hz:{30,60,120,240}){auto c=make_case("rod_kick");const auto&m=c.model;State s=c.initial,ref=s;double h=1./hz,E0=energy(m,s),E=E0-minimum_energy(m),maxP=0,maxS=0,maxE=0,error=0;int bad=0;for(int i=0;i<2*hz;++i){Problem p(m,s,h,E);if(!reaction)p.aold=masses(m).cwiseInverse().cwiseProduct(force(m,s.x));Options o;auto r=solve(p,o);if(!r.converged)++bad;maxP=std::max(maxP,r.primal);maxS=std::max(maxS,r.dual_residual);s=p.state(r.z);ref=reference_step(m,ref,h,128);error+=(s.x-ref.x).squaredNorm();maxE=std::max(maxE,std::abs(energy(m,s)-E0)/E);}csv<<reaction<<','<<h<<",2,"<<2*hz<<','<<bad<<','<<maxP<<','<<maxS<<','<<maxE<<','<<std::sqrt(error/(2*hz))<<','<<(s.x-ref.x).norm()<<','<<(s.v-ref.v).norm()<<'\n';}
 std::cout<<"rod reaction control complete\n";
}
