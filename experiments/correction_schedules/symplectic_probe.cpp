#include "correction.hpp"
#include <array>
#include <complex>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
using namespace energy_ordering;
enum class ProbeMethod { Predictor, Cold, Warm };
const char* name(ProbeMethod m){return m==ProbeMethod::Predictor?"old_force_predictor":m==ProbeMethod::Cold?"cold_xpbd":"warm_plus_xpbd";}
template<class T> std::array<T,2> scalar(std::array<T,2> z,double h,double mass,double stiffness,double gravity,ProbeMethod method) {
    const T oldq=z[0],oldp=z[1],oldlambda=-h*h*stiffness*oldq;
    T q=oldq+h*oldp/mass+h*h*gravity,lambda=0.;
    if(method!=ProbeMethod::Cold)q+=oldlambda/mass;
    if(method==ProbeMethod::Warm)lambda=oldlambda;
    if(method!=ProbeMethod::Predictor)q+=(-(q+lambda/(stiffness*h*h))/(1/mass+1/(stiffness*h*h)))/mass;
    return {q,mass*(q-oldq)/h};
}
template<class T> std::array<T,4> radial(std::array<T,4> z,double h,double mass,double stiffness,double gravity,ProbeMethod method) {
    using std::sqrt;
    const T oldr=sqrt(z[0]*z[0]+z[1]*z[1]),oldlambda=-h*h*stiffness*(oldr-1.);
    T x=z[0]+h*z[2]/mass,y=z[1]+h*z[3]/mass-h*h*gravity,lambda=0.;
    if(method!=ProbeMethod::Cold){x+=z[0]/oldr*oldlambda/mass;y+=z[1]/oldr*oldlambda/mass;}
    if(method==ProbeMethod::Warm)lambda=oldlambda;
    if(method!=ProbeMethod::Predictor) {
        const T r=sqrt(x*x+y*y),dl=-((r-1.)+lambda/(stiffness*h*h))/(1/mass+1/(stiffness*h*h));
        x+=x/r*dl/mass;y+=y/r*dl/mass;
    }
    return {x,y,mass*(x-z[0])/h,mass*(y-z[1])/h};
}
template<size_t N> using Matrix=std::array<std::array<double,N>,N>;
template<size_t N> double determinant(Matrix<N> a) {
    double det=1;
    for(size_t i=0;i<N;++i) {
        size_t p=i;for(size_t j=i+1;j<N;++j)if(std::abs(a[j][i])>std::abs(a[p][i]))p=j;
        if(p!=i){std::swap(a[p],a[i]);det=-det;}
        det*=a[i][i];
        if(a[i][i]==0)return 0;
        for(size_t j=i+1;j<N;++j){const double r=a[j][i]/a[i][i];for(size_t k=i+1;k<N;++k)a[j][k]-=r*a[i][k];}
    }
    return det;
}
template<size_t N> double defect(const Matrix<N>& a) {
    double largest=0;
    for(size_t i=0;i<N;++i)for(size_t j=0;j<N;++j) {
        double entry=0;
        for(size_t k=0;k<N/2;++k)entry+=a[k][i]*a[k+N/2][j]-a[k+N/2][i]*a[k][j];
        const double expected=j==i+N/2?1:i==j+N/2?-1:0;
        largest=std::max(largest,std::abs(entry-expected));
    }
    return largest;
}
template<size_t N,class F> Matrix<N> jacobian(const std::array<double,N>& state,F map) {
    Matrix<N> out{};const double eps=1e-30;
    for(size_t j=0;j<N;++j) {
        std::array<std::complex<double>,N> z;
        for(size_t i=0;i<N;++i)z[i]=state[i];
        z[j]+=std::complex<double>(0,eps);const auto next=map(z);
        for(size_t i=0;i<N;++i)out[i][j]=next[i].imag()/eps;
    }
    return out;
}
void require(bool ok,const char* message){if(!ok)throw std::runtime_error(message);}
void reconstruction_checks() {
    Model m;m.mass=2;m.acceleration={0,-9.81};State state{{.7,-1.1},{.4,.2}};
    const double h=1./480;const Context c=context(m,state,h,0,energy(m,state));
    Iterate previous{c.predictor};Counters count;
    require(material_update(c,previous,count),"old material update failed");
    const Rows r=evaluate(c,previous);
    require(norm(r.R)<1e-14,"cold single radial solution must satisfy top equation");
    State next{previous.x,(previous.x-state.x)/h};const Context cnext=context(m,next,h,0,energy(m,next));
    const V2 reused=cnext.predictor+r.n*(previous.lambda/m.mass);
    const V2 translated=previous.x+(cnext.predictor-c.predictor);
    require(norm(reused-translated)<1e-14,"previous solution correction did not reconstruct on the new prediction");
}
}
int main(int argc,char** argv) {
    try {
        const std::filesystem::path dir=argc>1?argv[1]:"symplectic_results";
        std::filesystem::create_directories(dir);
        std::ofstream out(dir/"jacobians.csv");out<<std::setprecision(17);
        out<<"case,method,h,mass,stiffness,gravity,determinant,symplectic_defect,reference_error,q0,q1,p0,p1\n";
        reconstruction_checks();size_t rows=0;double maxpredictor=0,maxparity=0;
        for(double h:{1./240,1./480,1./960})for(double mass:{.3,1.,4.})for(double gravity:{0.,9.81})
        for(auto method:{ProbeMethod::Predictor,ProbeMethod::Cold,ProbeMethod::Warm}) {
            const double k=100,a=h*h*k/mass;
            const auto j=jacobian<2>({.25,.7},[&](auto z){return scalar(z,h,mass,k,gravity,method);});
            Matrix<2> expected=method==ProbeMethod::Predictor?Matrix<2>{{{1-a,h/mass},{-h*k,1}}}:
                Matrix<2>{{{1/(1+a),h/mass/(1+a)},{-h*k/(1+a),1/(1+a)}}};
            double err=0;for(int i=0;i<2;++i)for(int l=0;l<2;++l)err=std::max(err,std::abs(j[i][l]-expected[i][l]));
            require(err<1e-11,"scalar map Jacobian disagrees with closed-form matrix");
            const double target=method==ProbeMethod::Predictor?1:1/(1+a);
            require(std::abs(determinant(j)-target)<1e-11,"scalar phase area factor mismatch");
            if(method==ProbeMethod::Predictor)require(defect(j)<1e-11,"scalar predictor is not symplectic");
            else require(defect(j)>1e-6,"XPBD should fail exact symplecticity on scalar spring");
            out<<"scalar,"<<name(method)<<','<<h<<','<<mass<<','<<k<<','<<gravity<<','<<determinant(j)<<','<<defect(j)<<','<<err<<",0.25,0,0.7,0\n";++rows;
            for(const auto& z:std::vector<std::array<double,4>>{{1.2,.3,.4,-.2},{.4,-.85,-.7,.5},{1.1,0,0,.8}}) {
                const auto jr=jacobian<4>(z,[&](auto s){return radial(s,h,mass,k,gravity,method);});
                const auto value=radial(z,h,mass,k,gravity,method);
                Model model;model.mass=mass;model.alpha=1/k;model.acceleration={0,-gravity};
                State old{{z[0],z[1]},{z[2]/mass,z[3]/mass}};const Context ctx=context(model,old,h,0,energy(model,old));
                Iterate actual{ctx.predictor};const double lold=-h*h*k*(norm(old.x)-1);
                if(method!=ProbeMethod::Cold)actual.x+=old.x/norm(old.x)*(lold/mass);
                if(method==ProbeMethod::Warm)actual.lambda=lold;
                if(method!=ProbeMethod::Predictor) {
                    correction_schedules::Options options;options.schedule=correction_schedules::Schedule::Baseline;
                    correction_schedules::Diagnostics d;
                    require(correction_schedules::update(ctx,actual,options,0,d),"actual production update failed");
                }
                const V2 p=(actual.x-old.x)*(mass/h);
                const double parity=std::max({std::abs(actual.x.x-value[0]),std::abs(actual.x.y-value[1]),std::abs(p.x-value[2]),std::abs(p.y-value[3])});
                maxparity=std::max(maxparity,parity);require(parity<1e-10,"complex-step map does not match actual XPBD update");
                if(method==ProbeMethod::Predictor){maxpredictor=std::max(maxpredictor,defect(jr));require(defect(jr)<1e-10,"nonlinear elastic predictor is not symplectic");}
                out<<"radial2d,"<<name(method)<<','<<h<<','<<mass<<','<<k<<','<<gravity<<','<<determinant(jr)<<','<<defect(jr)<<','<<parity<<','<<z[0]<<','<<z[1]<<','<<z[2]<<','<<z[3]<<'\n';++rows;
            }
        }
        std::cout<<"PASS "<<rows<<" canonical Jacobian checks; max nonlinear predictor defect="<<maxpredictor
            <<", max parity error versus production="<<maxparity<<"; warm reconstruction exact in the tested r_x=0 case\n";
        return 0;
    }catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}
