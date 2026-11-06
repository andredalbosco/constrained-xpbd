#include "incremental.hpp"
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
using namespace hour_local;
int main(int argc,char**argv){try{std::filesystem::path out=argc>1?argv[1]:"results";std::filesystem::create_directories(out);std::ofstream csv(out/"incremental_validation.csv");csv<<std::setprecision(16)<<"case,probe,objective_difference_error,constraint_difference_error,passed\n";std::vector<Case>cases;for(std::string name:{"scalar_release","radial_orbit","rod_kick","elastic_mixed","chain_mixed"})cases.push_back(make_case(name));cases.push_back(chain_case(8,1,10));for(auto&c:cases){double E=energy(c.model,c.initial)-minimum_energy(c.model);Problem p(c.model,c.initial,1./60,E);V initial=p.initial();for(int probe=0;probe<20;++probe){V a=initial,b=initial;for(int i=0;i<a.size();++i){a[i]+=.1*std::sin(i+.3*probe);b[i]+=.1*std::cos(2*i+.2*probe);}auto ea=p.evaluate(a),eb=p.evaluate(b);auto delta=difference(p,a,b,ea);double f=std::abs(delta.objective-(eb.objective.a-ea.objective.a)),g=(delta.rows-(residuals(eb)-residuals(ea))).cwiseAbs().maxCoeff();bool pass=std::max(f,g)<1e-10;csv<<c.name<<','<<probe<<','<<f<<','<<g<<','<<pass<<'\n';if(!pass)throw std::runtime_error("incremental_difference_mismatch");}}std::cout<<"120 exact-difference comparisons pass\n";return 0;}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}
