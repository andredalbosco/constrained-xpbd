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
struct Job{int N,budget;std::string scenario,partition,order;double rho;};
int main(int argc,char**argv){try{std::filesystem::path out=argc>1?argv[1]:"results";std::filesystem::create_directories(out);std::vector<Job>jobs;for(int N:{2,8,32})for(std::string scenario:{"nominal","stiff10","mass10"})for(std::string partition:{"vertex_full","edge_full"})for(double rho:{1.,10.,100.})for(int budget:{256,1024})jobs.push_back({N,budget,scenario,partition,"forward",rho});
 std::ofstream csv(out/"block_comparison.csv");csv<<std::setprecision(16)<<"N,scenario,partition,order,rho,budget,converged,primal,stationarity,objective_gap,scaled_state_difference,reduced_curvature,evaluations,visits,dual_sweeps,compute_ms,reference_evaluations,reference_curvature\n";std::mutex mutex;std::atomic<int>next{0};std::vector<std::thread>threads;
 for(int t=0;t<4;++t)threads.emplace_back([&](){for(;;){int i=next.fetch_add(1);if(i>=int(jobs.size()))break;auto j=jobs[i];auto c=chain_case(j.N,j.scenario=="stiff10"?10:1,j.scenario=="mass10"?10:1);double E=energy(c.model,c.initial)-minimum_energy(c.model);Problem p(c.model,c.initial,1./60,E);Options q;auto ref=solve(p,q);if(!ref.converged)throw std::runtime_error("reference_failed");Options o;o.solver="AL";o.rho=j.rho;o.budget=j.budget;o.partition=j.partition;o.order=j.order;auto r=solve(p,o);double rc=r.converged?reduced_curvature(p,r):0,refc=reduced_curvature(p,ref);std::ostringstream row;row<<std::setprecision(16)<<j.N<<','<<j.scenario<<','<<j.partition<<','<<j.order<<','<<j.rho<<','<<j.budget<<','<<r.converged<<','<<r.primal<<','<<r.dual_residual<<','<<r.objective-ref.objective<<','<<(r.z-ref.z).norm()<<','<<rc<<','<<r.evaluations<<','<<r.visits<<','<<r.sweeps<<','<<r.milliseconds<<','<<ref.evaluations<<','<<refc<<'\n';std::lock_guard<std::mutex>lock(mutex);csv<<row.str();csv.flush();if(i%12==0)std::cout<<i<<'/'<<jobs.size()<<'\n';}});for(auto&t:threads)t.join();std::cout<<jobs.size()<<" block comparisons complete\n";return 0;}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}
