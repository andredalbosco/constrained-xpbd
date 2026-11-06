#include "correction.hpp"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>

using namespace energy_ordering;
namespace cs=correction_schedules;
namespace {
struct LD2 {long double x=0,y=0;};
long double length(LD2 a){return std::hypotl(a.x,a.y);}
LD2 subtract(LD2 a,LD2 b){return {a.x-b.x,a.y-b.y};}
LD2 as_long(V2 a){return {a.x,a.y};}
void require(bool condition,const char* message){if(!condition)throw std::runtime_error(message);}
V2 naive(const Rows& r,double mass,cs::Proposal z){return (r.n*z.lambda+r.g*z.eta)/mass;}
V2 stable_material_closed(const Rows& r,double mass,double eta){
    const double normal=(-r.rC+r.at*r.p*eta)/(mass*r.A);
    const double transverse=r.q*eta/mass;
    return r.n*normal+tangent(r.n)*transverse;
}
Rows fixture(double mass,double at,double p,double q,double C,double H){
    Rows r;r.n={1,0};r.g={p,q};r.at=at;r.p=p;r.q=q;r.rC=r.C=C;r.rE=r.H=H;
    r.A=1/mass+at;r.B=p/mass;r.D=(p*p+q*q)/mass;
    r.schur=q*q/mass+(at/r.A)*p*p/mass;r.gradient_scale=1+norm(r.g);
    return r;
}
struct Reference {long double lambda=0,eta=0;LD2 dx;};
Reference exact_reference(const Rows& r,double mass,double blend,bool fallback){
    const long double m=mass,at=r.at,p=r.p,q=r.q,C=r.rC,H=r.rE,A=1/m+at;
    long double eta=0;
    if(!fallback){const long double schur=q*q/m+(at/A)*p*p/m;eta=(-H+(p/m)*C/A)/schur;}
    eta*=blend;
    const long double lambda=(-C-(p/m)*eta)/A;
    return {lambda,eta,{(-C+at*p*eta)/(m*A),q*eta/m}};
}

void probe_rows(const std::string& filename){
    std::ofstream out(filename);require(bool(out),"cannot open assembly CSV");out<<std::setprecision(17);
    out<<"mass,dimensionless_at,transverse_gradient,material_residual,energy_residual,blend,rank_fallback,lambda,eta,"
          "naive_dx,naive_dy,stable_dx,stable_dy,stored_pair_long_dx,stored_pair_long_dy,oracle_dx,oracle_dy,"
          "naive_error,stable_error,stored_pair_error,naive_material_error,stable_material_error,"
          "naive_energy_error,stable_energy_error,stable_top_absolute_error,stable_top_relative_error\n";
    long count=0,fallbacks=0;long double max_naive=0,max_stable=0,max_pair=0,max_naive_scaled=0,max_stable_scaled=0,max_top=0,max_top_relative=0;
    Context c;c.h=1./8;c.model.alpha=.01;c.target=1;
    for(double mass:{.125,1.,8.})for(int exponent:{0,8,16,24,32,36,38,40,44,48,50,52})
    for(double transverse:{0.,std::ldexp(1.,-12),std::ldexp(1.,-24)})for(double C:{0.,.125})for(double blend:{1.,.5,.001}){
        const double at=std::ldexp(1.,-exponent)/mass,p=2*mass,q=transverse*mass;
        const Rows r=fixture(mass,at,p,q,C,.3*mass);c.model.mass=mass;
        cs::Options o;o.schedule=cs::Schedule::Joint;cs::Diagnostics d;const auto raw=cs::propose(c,r,o,0,d);
        const double baseline=-r.rC/r.A;
        const cs::Proposal accepted{blend==1?raw.lambda:baseline+blend*(raw.lambda-baseline),blend*raw.eta};
        const V2 dx_naive=naive(r,mass,accepted),dx_stable=stable_material_closed(r,mass,accepted.eta);
        const LD2 stored_pair{(static_cast<long double>(accepted.lambda)+static_cast<long double>(p)*accepted.eta)/mass,
                              static_cast<long double>(q)*accepted.eta/mass};
        const Reference oracle=exact_reference(r,mass,blend,d.rank_fallbacks!=0);
        const long double naive_error=length(subtract(as_long(dx_naive),oracle.dx));
        const long double stable_error=length(subtract(as_long(dx_stable),oracle.dx));
        const long double pair_error=length(subtract(stored_pair,oracle.dx));
        const long double material_naive=std::abs(static_cast<long double>(C)+dx_naive.x+static_cast<long double>(at)*accepted.lambda);
        const long double material_stable=std::abs(static_cast<long double>(C)+dx_stable.x+static_cast<long double>(at)*accepted.lambda);
        const long double expected_energy=d.rank_fallbacks?(r.rE+r.B*baseline):(1-blend)*(r.rE+r.B*baseline);
        const long double energy_naive=std::abs(static_cast<long double>(r.rE)+static_cast<long double>(p)*dx_naive.x+static_cast<long double>(q)*dx_naive.y-expected_energy);
        const long double energy_stable=std::abs(static_cast<long double>(r.rE)+static_cast<long double>(p)*dx_stable.x+static_cast<long double>(q)*dx_stable.y-expected_energy);
        const long double top_error=mass*length(subtract(as_long(dx_stable),stored_pair));
        const long double top_scale=1+std::abs(static_cast<long double>(accepted.lambda))+length({static_cast<long double>(p)*accepted.eta,static_cast<long double>(q)*accepted.eta});
        max_naive=std::max(max_naive,naive_error);max_stable=std::max(max_stable,stable_error);max_pair=std::max(max_pair,pair_error);
        max_naive_scaled=std::max(max_naive_scaled,naive_error/(1+length(oracle.dx)));max_stable_scaled=std::max(max_stable_scaled,stable_error/(1+length(oracle.dx)));
        max_top=std::max(max_top,top_error);max_top_relative=std::max(max_top_relative,top_error/top_scale);
        require(stable_error<=128*std::numeric_limits<double>::epsilon()*(1+length(oracle.dx)),"stable assembly differs from long-double exact row solve");
        require(top_error<=32*std::numeric_limits<double>::epsilon()*top_scale,"stable assembly is not backward consistent with rounded multipliers");
        ++count;fallbacks+=d.rank_fallbacks;
        out<<mass<<','<<at*mass<<','<<q<<','<<C<<','<<r.rE<<','<<blend<<','<<d.rank_fallbacks<<','<<accepted.lambda<<','<<accepted.eta<<','
           <<dx_naive.x<<','<<dx_naive.y<<','<<dx_stable.x<<','<<dx_stable.y<<','<<stored_pair.x<<','<<stored_pair.y<<','<<oracle.dx.x<<','<<oracle.dx.y<<','
           <<naive_error<<','<<stable_error<<','<<pair_error<<','<<material_naive<<','<<material_stable<<','<<energy_naive<<','<<energy_stable<<','<<top_error<<','<<top_error/top_scale<<'\n';
    }
    std::cout<<std::setprecision(12)<<"Assembly rows="<<count<<", rank fallbacks="<<fallbacks
             <<"\nmax absolute errors: naive="<<max_naive<<", stable="<<max_stable<<", long-double stored pair="<<max_pair
             <<"\nmax scaled errors: naive="<<max_naive_scaled<<", stable="<<max_stable_scaled
             <<"\nrounded-pair top residual: max absolute="<<max_top<<", max relative="<<max_top_relative<<'\n';
}

void scope_checks(){
    // Material-ending candidates and any blend with the baseline retain the
    // row required by the stable formula. Energy-ending candidates need not.
    const Rows r=fixture(1,1,2,3,.25,.5);Context c;c.h=.125;
    for(auto schedule:{cs::Schedule::EL,cs::Schedule::LEL,cs::Schedule::Joint}){
        cs::Options o;o.schedule=schedule;cs::Diagnostics d;const auto raw=cs::propose(c,r,o,0,d);
        for(double blend:{0.,.125,.5,1.}){
            const cs::Proposal z{-r.rC/r.A+blend*(raw.lambda+r.rC/r.A),blend*raw.eta};
            require(norm(naive(r,1,z)-stable_material_closed(r,1,z.eta))<1e-15,"stable formula should match material-closing candidates and blends");
        }
    }
    for(auto schedule:{cs::Schedule::LE,cs::Schedule::ELE,cs::Schedule::Jacobi}){
        cs::Options o;o.schedule=schedule;cs::Diagnostics d;const auto z=cs::propose(c,r,o,0,d);
        require(norm(naive(r,1,z)-stable_material_closed(r,1,z.eta))>1e-3,"scope fixture must expose invalid use on non-material-closing candidates");
        const double residual=r.rC+r.A*z.lambda+r.B*z.eta;
        const V2 generalized=stable_material_closed(r,1,z.eta)+r.n*(residual/r.A);
        require(norm(naive(r,1,z)-generalized)<1e-15,"adding material residual yields the general identity");
    }
}

void release_control(){
    // Stable arithmetic cannot repair the exact Joint method's frozen release.
    Model m;const State initial{{1.25,0},{}};const double E0=energy(m,initial),h=1./480;
    State naive_state=initial,stable_state=initial;double max_difference=0,naive_path=0,stable_path=0;
    for(int step=0;step<4800;++step){
        for(bool use_stable:{false,true}){
            State& state=use_stable?stable_state:naive_state;Context c=context(m,state,h,0,E0);Iterate z{c.predictor};
            const Rows r=evaluate(c,z);cs::Options o;o.schedule=cs::Schedule::Joint;cs::Diagnostics d;const auto proposal=cs::propose(c,r,o,step,d);
            const V2 dx=use_stable?stable_material_closed(r,m.mass,proposal.eta):naive(r,m.mass,proposal);
            const V2 next=z.x+dx;const double distance=norm(next-state.x);
            (use_stable?stable_path:naive_path)+=distance;
            state={next,(next-state.x)/h};
        }
        max_difference=std::max(max_difference,norm(naive_state.x-stable_state.x));
    }
    require(max_difference<1e-12&&naive_path<1e-12&&stable_path<1e-12,"release control should remain frozen for both equivalent joint assemblies");
    std::cout<<"10s spring-release control: naive path="<<naive_path<<", stable path="<<stable_path<<", max position difference="<<max_difference<<'\n';
}
} // namespace
int main(int argc,char** argv){
    try{scope_checks();probe_rows(argc>1?argv[1]:"/tmp/xpbd-correction-schedules/assembly_probe.csv");release_control();}
    catch(const std::exception& e){std::cerr<<"FAIL: "<<e.what()<<'\n';return 1;}
}
