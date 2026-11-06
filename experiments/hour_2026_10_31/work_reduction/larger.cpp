#include "analytic_reduction.hpp"
#include "reduction.hpp"
#include <fstream>
#include <filesystem>
#include <iostream>
#include <iomanip>
namespace fs=std::filesystem;
using namespace hour_local;
void validate(const fs::path&out){std::ofstream csv(out/"analytic_validation.csv");csv<<"case,objective,gradient,hessian,pass\n";for(std::string name:{"scalar_mixed","radial_orbit","elastic_mixed","chain_mixed"}){auto c=make_case(name);double E=std::max(1e-6,energy(c.model,c.initial)-minimum_energy(c.model));hour_local::Problem full(c.model,c.initial,1./60,E);work_reduction_fast::Problem fast(full);work_reduction::Reduced oracle(full);V q=fast.initial();for(int j=0;j<q.size();++j)q[j]+=.07*std::sin(j+1);auto a=fast.evaluate(q),b=oracle.evaluate(q);double o=std::abs(a.objective.a-b.objective.a),g=(a.objective.g-b.objective.g).norm(),H=(a.objective.H-b.objective.H).norm();bool pass=o<1e-10&&g<1e-10&&H<1e-10;csv<<name<<','<<o<<','<<g<<','<<H<<','<<pass<<'\n';if(!pass)throw std::runtime_error("analytic_reduction_oracle");}}
int main(int argc,char**argv){try{fs::path out=argc>1?argv[1]:"larger_results";fs::create_directories(out);validate(out);if(argc>2&&std::string(argv[2])=="validate")return 0;std::ofstream csv(out/"one_step.csv");csv<<std::setprecision(16)<<"N,stiffness,mass_ratio,h,form,budget,converged,primal,stationarity,full_stationarity,evaluations,ms,x_gap,v_gap,work_gap,kkt_size\n";
 for(int N:{2,8,32})for(double stiffness:{1.,10.})for(double mass:{1.,10.})for(double h:{1./60,1./30}){auto c=chain_case(N,stiffness,mass);double E=std::max(1e-6,energy(c.model,c.initial)-minimum_energy(c.model));hour_local::Problem f(c.model,c.initial,h,E);work_reduction_fast::Problem p(f);V q=p.initial(),z=p.expand(q);Solve ref=joint_solve(f,128,&z);State sref=f.state(ref.z);
 for(int budget:{16,64,128}){auto a=joint_solve(p,budget,&q);auto s=p.state(a.z);V expanded=p.expand(a.z);csv<<N<<','<<stiffness<<','<<mass<<','<<h<<",reduced,"<<budget<<','<<a.converged<<','<<a.primal<<','<<a.dual_residual<<','<<stationarity(f.evaluate(expanded),p.full_dual(a))<<','<<a.evaluations<<','<<a.milliseconds<<','<<(s.x-sref.x).norm()<<','<<(s.v-sref.v).norm()<<','<<(expanded.tail(f.p.wcount)-ref.z.tail(f.p.wcount)).norm()<<','<<p.p.n+1<<'\n';}
 csv<<N<<','<<stiffness<<','<<mass<<','<<h<<",full,128,"<<ref.converged<<','<<ref.primal<<','<<ref.dual_residual<<','<<ref.dual_residual<<','<<ref.evaluations<<','<<ref.milliseconds<<",0,0,0,"<<f.p.n+f.evaluate(ref.z).rows.size()<<'\n';csv.flush();std::cout<<N<<' '<<stiffness<<' '<<mass<<' '<<h<<std::endl;
 }return 0;}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}
