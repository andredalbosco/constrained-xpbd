#include "objectives.hpp"
#include "solvers.hpp"
#include <fstream>
#include <filesystem>
#include <iomanip>
#include <iostream>
using namespace compose;
using namespace compose::independent_detail;
namespace fs=std::filesystem;
struct Combined {
 const objectives2::Objective&a;const objectives2::Objective&b;bool global;
 struct Metadata{int n;};Metadata p;
 Combined(const objectives2::Objective&aa,const objectives2::Objective&bb,bool g):a(aa),b(bb),global(g),p{aa.p.n+bb.p.n}{}
 Jet embed(const Jet&j,int offset)const{Jet r(j.a,p.n);r.g.segment(offset,j.g.size())=j.g;r.H.block(offset,offset,j.g.size(),j.g.size())=j.H;return r;}
 Evaluation evaluate(const V&z)const{auto x=a.evaluate(z.head(a.p.n)),y=b.evaluate(z.tail(b.p.n));Evaluation e;e.objective=embed(x.objective,0)+embed(y.objective,a.p.n);Jet total(0,p.n);for(const auto&row:x.rows){Jet r=embed(row.value,0);if(global&&row.kind==Global)total=total+r;else e.rows.push_back({r,row.kind,row.index});}for(const auto&row:y.rows){Jet r=embed(row.value,a.p.n);if(global&&row.kind==Global)total=total+r;else e.rows.push_back({r,row.kind,row.index+1});}if(global)e.rows.push_back({total,Global,0});return e;}
 V initial()const{V z(p.n);z.head(a.p.n)=a.initial();z.tail(b.p.n)=b.initial();return z;}
};
int main(int argc,char**argv){try{fs::path out=argc>1?argv[1]:"disconnected_results";fs::create_directories(out);std::ofstream csv(out/"trajectories.csv");csv<<std::setprecision(16)<<"first,second,objective,closure,h,steps,nonconverged,total_energy_error,first_energy_error,second_energy_error,position_rms,reference_gap\n";
 for(auto pair:std::vector<std::pair<std::string,std::string>>{{"scalar_mixed","radial_orbit"},{"radial_orbit","elastic_mixed"},{"scalar_mixed","chain_mixed"}})for(std::string objective:{"original","quarter","projection"})for(bool global:{false,true})for(double h:{1./30,1./120}){
  auto a=make_case(pair.first),b=make_case(pair.second);State sa=a.initial,sb=b.initial,ra=sa,rb=sb,ra2=sa,rb2=sb;double E0a=energy(a.model,sa),E0b=energy(b.model,sb),Ea=std::max(1e-6,E0a-minimum_energy(a.model)),Eb=std::max(1e-6,E0b-minimum_energy(b.model)),E=Ea+Eb,emax=0,amax=0,bmax=0,sumx=0,gap=0;int bad=0,done=0;
  objectives2::Variant variant;variant.objective=objective;variant.energy=true;variant.beta=1;std::ofstream trace(out/(pair.first+"_"+pair.second+"_"+objective+"_"+(global?"global":"component")+"_"+std::to_string(h)+".csv"));trace<<std::setprecision(16)<<"t,first_energy_change,second_energy_change,total_energy_change\n";
  for(int i=0;i<std::lround(5/h);++i){objectives2::Objective fa(a.model,sa,h,variant,E),fb(b.model,sb,h,variant,E);Combined p(fa,fb,global);auto s=local_work::joint_solve(p,128);if(!s.converged)++bad;sa=fa.state(fa.evaluate(s.z.head(fa.p.n)));sb=fb.state(fb.evaluate(s.z.tail(fb.p.n)));if(!finite(sa)||!finite(sb))break;ra=reference_step(a.model,ra,h,64);rb=reference_step(b.model,rb,h,64);ra2=reference_step(a.model,ra2,h,128);rb2=reference_step(b.model,rb2,h,128);gap=std::max(gap,(ra.x-ra2.x).norm()+(rb.x-rb2.x).norm()+(ra.v-ra2.v).norm()/std::sqrt(Ea)+(rb.v-rb2.v).norm()/std::sqrt(Eb));double da=energy(a.model,sa)-E0a,db=energy(b.model,sb)-E0b;emax=std::max(emax,std::abs(da+db)/E);amax=std::max(amax,std::abs(da)/Ea);bmax=std::max(bmax,std::abs(db)/Eb);sumx+=((sa.x-ra2.x).squaredNorm()+(sb.x-rb2.x).squaredNorm())/(a.model.count+b.model.count);++done;trace<<(i+1)*h<<','<<da<<','<<db<<','<<da+db<<'\n';}
  csv<<pair.first<<','<<pair.second<<','<<objective<<','<<(global?"global":"component")<<','<<h<<','<<done<<','<<bad<<','<<emax<<','<<amax<<','<<bmax<<','<<std::sqrt(sumx/std::max(1,done))<<','<<gap<<'\n';csv.flush();std::cout<<pair.first<<' '<<pair.second<<' '<<objective<<' '<<global<<' '<<h<<std::endl;
 }return 0;}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}
