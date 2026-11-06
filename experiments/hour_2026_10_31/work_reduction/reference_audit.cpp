#include "../local_solver/problem.hpp"
#include <fstream>
#include <filesystem>
#include <iomanip>
#include <iostream>
using namespace hour_local;
int main(int argc,char**argv){std::filesystem::path out=argc>1?argv[1]:"rollout_results";std::filesystem::create_directories(out);std::ofstream csv(out/"references.csv");csv<<std::setprecision(16)<<"N,stiffness,mass_ratio,h,position_gap,scaled_velocity_gap,energy_error,certified\n";for(int N:{8,32})for(double k:{1.,10.})for(double h:{1./30,1./60}){auto c=chain_case(N,k,10);State a=c.initial,b=a;double E0=energy(c.model,a),E=std::max(1e-6,E0-minimum_energy(c.model)),xgap=0,vgap=0,energyerr=0;for(int i=0;i<std::lround(1/h);++i){a=reference_step(c.model,a,h,128);b=reference_step(c.model,b,h,256);xgap=std::max(xgap,(a.x-b.x).norm());vgap=std::max(vgap,(a.v-b.v).norm()/std::sqrt(E));energyerr=std::max(energyerr,std::abs(energy(c.model,b)-E0)/E);}bool pass=xgap+vgap<1e-5&&energyerr<1e-7;csv<<N<<','<<k<<",10,"<<h<<','<<xgap<<','<<vgap<<','<<energyerr<<','<<pass<<'\n';std::cout<<N<<' '<<k<<' '<<h<<' '<<pass<<std::endl;}}
