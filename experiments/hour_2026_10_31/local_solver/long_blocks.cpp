#include "numerics.hpp"
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <thread>
#include <mutex>
#include <atomic>
using namespace hour_local;
struct Job{int N,budget;std::string scenario,partition;double rho;};
int main(int argc,char**argv){try{std::filesystem::path out=argc>1?argv[1]:"results";bool stable=argc>2&&std::string(argv[2])=="stable";std::filesystem::create_directories(out);std::vector<Job>jobs;for(int N:{8,32})for(std::string scenario:{"nominal","mass10"})for(std::string partition:{"vertex_full","edge_full"})for(double rho:{10.,100.})jobs.push_back({N,N==8?4096:16384,scenario,partition,rho});for(int N:{2,8,32})for(std::string scenario:{"nominal","mass10"})for(double rho:{10.,100.})jobs.push_back({N,512,scenario,"global",rho});
 std::ofstream csv(out/(stable?"long_block_stable.csv":"long_block_comparison.csv"));csv<<std::setprecision(16)<<"N,scenario,partition,rho,budget,converged,primal,stationarity,objective_gap,scaled_state_difference,reduced_curvature,evaluations,visits,dual_sweeps,compute_ms\n";std::mutex mutex;std::atomic<int>next{0};std::vector<std::thread>threads;
 for(int t=0;t<4;++t)threads.emplace_back([&](){for(;;){int i=next.fetch_add(1);if(i>=int(jobs.size()))break;auto j=jobs[i];auto c=chain_case(j.N,1,j.scenario=="mass10"?10:1);double E=energy(c.model,c.initial)-minimum_energy(c.model);Problem p(c.model,c.initial,1./60,E);Options q;auto ref=solve(p,q);Options o;o.solver="AL";o.rho=j.rho;o.budget=j.budget;o.partition=j.partition;o.stable=stable;auto r=solve(p,o);double curvature=r.converged?reduced_curvature(p,r):0;std::ostringstream row;row<<std::setprecision(16)<<j.N<<','<<j.scenario<<','<<j.partition<<','<<j.rho<<','<<j.budget<<','<<r.converged<<','<<r.primal<<','<<r.dual_residual<<','<<r.objective-ref.objective<<','<<(r.z-ref.z).norm()<<','<<curvature<<','<<r.evaluations<<','<<r.visits<<','<<r.sweeps<<','<<r.milliseconds<<'\n';std::lock_guard<std::mutex>lock(mutex);csv<<row.str();csv.flush();std::cout<<j.N<<' '<<j.scenario<<' '<<j.partition<<' '<<j.rho<<" converged="<<r.converged<<" eval="<<r.evaluations<<'\n';}});for(auto&t:threads)t.join();return 0;}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}
