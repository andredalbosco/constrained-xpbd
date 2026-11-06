#pragma once
#include "../correction_schedules/correction.hpp"
#include <array>
#include <stdexcept>
#include <string>
#include <vector>

namespace traditional {
using energy_ordering::V2;
using energy_ordering::dot;
using energy_ordering::norm;
using energy_ordering::finite;
using energy_ordering::tangent;
using Vec=std::array<V2,2>;
using Matrix=std::array<std::array<double,3>,3>;
struct Model {
    int count=1;
    std::array<double,2> mass{1,1},length{1,1},alpha{.01,.01};
    V2 gravity{};
};
struct State {Vec x{},v{};};
struct Context {Model model;State old;Vec predictor{};double h=1./480,target=0;};
struct Iterate {Vec x{};std::array<double,2> lambda{};double eta=0;};
struct Rows {
    std::array<Vec,2> J{};Vec g{},objective{},R{};
    std::array<double,2> C{},rC{},at{};
    Matrix S{};double H=0,gradient_scale=1;bool valid=true;
};
enum class Schedule {Baseline,EL,LE,Jacobi,Joint,LEL,ELE};
enum class Policy {Raw,Direction,RepairCap,RawEnergyGate,RawEnergyGateCap};
struct Options {Schedule schedule=Schedule::Baseline;Policy policy=Policy::Raw;bool fresh=false;double fraction=.5,cap=.1;};
struct Diagnostics {
    long evaluations=0,scalar_solves=0,blocks=0,singular=0,rank_fallbacks=0,caps=0,skips=0,repair_scaled=0;
    double max_ratio=0,max_prediction_error=0,last_blend=1,last_extra=0,last_H0=0,last_predicted_H=0;
};
inline double inner(const Vec& a,const Vec& b,int n){double s=0;for(int i=0;i<n;++i)s+=dot(a[i],b[i]);return s;}
inline double mass_norm(const Model& m,const Vec& a){double s=0;for(int i=0;i<m.count;++i)s+=m.mass[i]*dot(a[i],a[i]);return std::sqrt(s);}
inline double weighted(const Model& m,const Vec& a,const Vec& b){double s=0;for(int i=0;i<m.count;++i)s+=dot(a[i],b[i])/m.mass[i];return s;}
inline Vec minv(const Model& m,Vec a){for(int i=0;i<m.count;++i)a[i]=a[i]/m.mass[i];return a;}
inline bool valid(const Vec& a,int n){for(int i=0;i<n;++i)if(!finite(a[i]))return false;return true;}
inline double cross(V2 a,V2 b){return a.x*b.y-a.y*b.x;}
inline double energy(const Model& m,const State& s) {
    double E=0;for(int i=0;i<m.count;++i){const double C=norm(s.x[i]-(i?s.x[i-1]:V2{}))-m.length[i];
        E+=.5*m.mass[i]*dot(s.v[i],s.v[i])-m.mass[i]*dot(m.gravity,s.x[i]);if(m.alpha[i]>0)E+=C*C/(2*m.alpha[i]);}
    return E;
}
inline double minimum_energy(const Model& m) {
    double V=0,down=0;for(int j=m.count-1;j>=0;--j){down+=m.mass[j];const double T=down*norm(m.gravity);V-=T*m.length[j]+.5*m.alpha[j]*T*T;}return V;
}
inline Vec potential_gradient(const Model& m,const Vec& x) {
    Vec u{};for(int j=0;j<m.count;++j)if(m.alpha[j]>0){const V2 edge=x[j]-(j?x[j-1]:V2{});const double r=norm(edge);
        const V2 f=edge*((r-m.length[j])/(m.alpha[j]*r));u[j]+=f;if(j)u[j-1]+=f*(-1);}
    return u;
}
inline Vec acceleration(const Model& m,const Vec& x){Vec a=potential_gradient(m,x);for(int i=0;i<m.count;++i)a[i]=m.gravity-a[i]/m.mass[i];return a;}
inline Context context(const Model& m,State old,double h,double target){Context c{m,old,{},h,target};for(int i=0;i<m.count;++i)c.predictor[i]=old.x[i]+old.v[i]*h+m.gravity*(h*h);return c;}
inline Rows evaluate(const Context& c,const Iterate& z) {
    Rows r;const auto& m=c.model;const int n=m.count;
    for(int j=0;j<n;++j) {
        const V2 edge=z.x[j]-(j?z.x[j-1]:V2{});const double length=norm(edge);
        if(!(length>1e-14*m.length[j])){r.valid=false;return r;}
        r.J[j][j]=edge/length;if(j)r.J[j][j-1]=r.J[j][j]*(-1);
        r.C[j]=length-m.length[j];r.at[j]=m.alpha[j]/(c.h*c.h);r.rC[j]=r.C[j]+r.at[j]*z.lambda[j];
    }
    const Vec u=potential_gradient(m,z.x);
    State endpoint;endpoint.x=z.x;
    for(int i=0;i<n;++i) {
        endpoint.v[i]=(z.x[i]-c.old.x[i])/c.h;
        r.g[i]=(z.x[i]-(c.old.x[i]+m.gravity*(c.h*c.h)))*(m.mass[i]/(c.h*c.h))+u[i];
        r.objective[i]=(z.x[i]-c.predictor[i])*(m.mass[i]/(c.h*c.h))+u[i];
        r.R[i]=(z.x[i]-c.predictor[i])*m.mass[i]-r.g[i]*z.eta;
        for(int j=0;j<n;++j)r.R[i]+=r.J[j][i]*(-z.lambda[j]);
        r.gradient_scale+=norm((z.x[i]-c.old.x[i])*(m.mass[i]/(c.h*c.h)))+norm(u[i])+m.mass[i]*norm(m.gravity);
    }
    r.H=energy(m,endpoint)-c.target;
    for(int j=0;j<n;++j){for(int k=0;k<n;++k)r.S[j][k]=weighted(m,r.J[j],r.J[k])+(j==k?r.at[j]:0);
        r.S[j][n]=r.S[n][j]=weighted(m,r.J[j],r.g);}
    r.S[n][n]=weighted(m,r.g,r.g);
    r.valid=valid(r.g,n)&&valid(r.R,n)&&std::isfinite(r.H)&&std::isfinite(r.S[n][n]);return r;
}
inline Vec displacement(const Model& m,const Rows& r,const std::array<double,3>& d) {
    Vec x{};for(int i=0;i<m.count;++i){x[i]=r.g[i]*d[m.count];for(int j=0;j<m.count;++j)x[i]+=r.J[j][i]*d[j];x[i]=x[i]/m.mass[i];}return x;
}
inline std::array<double,3> material_pass(const Model& m,const Rows& r,std::array<double,3> d,Diagnostics* stats=nullptr) {
    for(int i=0;i<m.count;++i){double rhs=-r.rC[i];for(int j=0;j<=m.count;++j)if(j!=i)rhs-=r.S[i][j]*d[j];d[i]=rhs/r.S[i][i];if(stats)++stats->scalar_solves;}return d;
}
inline bool normalized_solve(const Matrix& s,const std::array<double,3>& rhs,int count,std::array<double,3>& answer) {
    Matrix L{};std::array<double,3> scale{},z{};
    for(int i=0;i<count;++i){if(!(s[i][i]>0))return false;scale[i]=std::sqrt(s[i][i]);}
    for(int i=0;i<count;++i)for(int j=0;j<=i;++j){double v=s[i][j]/(scale[i]*scale[j]);for(int k=0;k<j;++k)v-=L[i][k]*L[j][k];
        if(i==j){if(!(v>1e-12))return false;L[i][j]=std::sqrt(v);}else L[i][j]=v/L[j][j];}
    for(int i=0;i<count;++i){double v=rhs[i]/scale[i];for(int j=0;j<i;++j)v-=L[i][j]*z[j];z[i]=v/L[i][i];}
    for(int i=count-1;i>=0;--i){double v=z[i];for(int j=i+1;j<count;++j)v-=L[j][i]*answer[j];answer[i]=v/L[i][i];}
    for(int i=0;i<count;++i)answer[i]/=scale[i];return true;
}
inline double motion_scale(const Context& c){return c.h*mass_norm(c.model,c.old.v)+c.h*c.h*mass_norm(c.model,acceleration(c.model,c.old.x));}
inline bool fresh_material(const Context& c,Iterate& z,Diagnostics& d) {
    for(int j=0;j<c.model.count;++j){const Rows r=evaluate(c,z);++d.evaluations;if(!r.valid)return false;
        const double dl=-r.rC[j]/r.S[j][j];for(int i=0;i<c.model.count;++i)z.x[i]+=r.J[j][i]*(dl/c.model.mass[i]);z.lambda[j]+=dl;++d.scalar_solves;}
    return valid(z.x,c.model.count);
}
inline bool update(const Context& c,Iterate& z,const Options& o,Diagnostics& stats) {
    if(o.fresh)return fresh_material(c,z,stats);
    const Rows r=evaluate(c,z);++stats.evaluations;if(!r.valid)return false;
    const auto& m=c.model;const int n=m.count;
    const auto base=material_pass(m,r,{});auto d=std::array<double,3>{};
    const Vec dx0=displacement(m,r,base);
    const double gradient_tolerance=64*std::numeric_limits<double>::epsilon()*r.gradient_scale;
    const bool zero=std::sqrt(inner(r.g,r.g,n))<=gradient_tolerance;
    auto L=[&](){d=material_pass(m,r,d,&stats);};
    auto E=[&](){
        if(zero){d[n]=0;++stats.singular;return;}
        double rhs=-r.H;for(int j=0;j<n;++j)rhs-=r.S[n][j]*d[j];d[n]=rhs/r.S[n][n];++stats.scalar_solves;
        if(o.policy==Policy::RawEnergyGate||o.policy==Policy::RawEnergyGateCap){
            Vec raw=minv(m,r.g);for(int i=0;i<n;++i)raw[i]=raw[i]*d[n];
            if(inner(r.objective,raw,n)>0){d[n]=0;++stats.skips;}
            else if(o.policy==Policy::RawEnergyGateCap&&mass_norm(m,raw)>o.cap*motion_scale(c)){
                d[n]*=o.cap*motion_scale(c)/mass_norm(m,raw);++stats.caps;}
        }
    };
    if(o.schedule==Schedule::Baseline)L();
    else if(o.schedule==Schedule::Joint){
        std::array<double,3> rhs{};for(int j=0;j<n;++j)rhs[j]=-r.rC[j];rhs[n]=-r.H;++stats.blocks;
        if(zero||!normalized_solve(r.S,rhs,n+1,d)){d=base;++stats.rank_fallbacks;if(zero)++stats.singular;}
    } else if(o.schedule==Schedule::Jacobi){for(int j=0;j<n;++j){d[j]=-r.rC[j]/r.S[j][j];++stats.scalar_solves;}
        if(zero)++stats.singular;else{d[n]=-r.H/r.S[n][n];++stats.scalar_solves;}}
    else if(o.schedule==Schedule::EL||o.schedule==Schedule::ELE){E();L();if(o.schedule==Schedule::ELE)E();}
    else{L();E();if(o.schedule==Schedule::LEL)L();}
    Vec candidate=displacement(m,r,d),extra{};for(int i=0;i<n;++i)extra[i]=candidate[i]-dx0[i];
    double blend=1;const double H0=r.H+inner(r.g,dx0,n),change=inner(r.g,extra,n);
    if(o.policy==Policy::Direction&&inner(r.objective,extra,n)>0){blend=0;++stats.skips;}
    if(o.policy==Policy::RepairCap){
        const double tol=128*std::numeric_limits<double>::epsilon()*(1+std::abs(c.target)+std::abs(r.H));
        if(std::abs(H0)<=tol||H0*change>=0||change==0){blend=0;++stats.skips;}
        else if(std::abs(change)>o.fraction*std::abs(H0)){blend=o.fraction*std::abs(H0)/std::abs(change);++stats.repair_scaled;}
        const double size=mass_norm(m,extra),scale=motion_scale(c);
        if(blend*size>o.cap*scale){blend=o.cap*scale/size;++stats.caps;}
        if(scale>0)stats.max_ratio=std::max(stats.max_ratio,blend*size/scale);
    }
    for(int j=0;j<=n;++j)d[j]=base[j]+blend*(d[j]-base[j]);
    const Vec dx=displacement(m,r,d);for(int i=0;i<n;++i)z.x[i]+=dx[i];for(int j=0;j<n;++j)z.lambda[j]+=d[j];z.eta+=d[n];
    if(!valid(z.x,n)||!std::isfinite(z.eta))return false;
    stats.last_H0=H0;stats.last_blend=blend;stats.last_extra=blend*mass_norm(m,extra);stats.last_predicted_H=r.H+inner(r.g,dx,n);
    State endpoint;endpoint.x=z.x;for(int i=0;i<n;++i)endpoint.v[i]=(z.x[i]-c.old.x[i])/c.h;
    const double actual=energy(m,endpoint)-c.target;stats.max_prediction_error=std::max(stats.max_prediction_error,std::abs(actual-stats.last_predicted_H));
    return std::isfinite(actual);
}
} // namespace traditional
