#include "cases.hpp"
#include "stability.hpp"
#include <iostream>
#include <random>
using namespace agenda;
void require(bool v,const char* message){if(!v)throw std::runtime_error(message);}
void near(double a,double b,double tol,const char* message){if(!std::isfinite(a)||!std::isfinite(b)||std::abs(a-b)>tol*(1+std::abs(b)))throw std::runtime_error(std::string(message)+": "+std::to_string(a)+" vs "+std::to_string(b));}
Context fixture(){Model m;m.count=2;m.mass={2,3};m.alpha={.03,.02};m.gravity={0,-2};State s;s.x={V2{.75,-1},V2{1.3,-1.7}};s.v={V2{.3,.5},V2{-.2,.1}};return context(m,s,.125,energy(m,s));}
void derivatives(){const auto c=fixture();Vec x=c.predictor;x[0].x+=.02;x[1].y-=.03;
    for(auto estimator:{Estimate::Lagged,Estimate::Column,Estimate::Row,Estimate::Secant}){
        const auto values=work_values(c,x,estimator);const auto l=layout(c.model,c.old);
        for(int k=0;k<4;++k){Vec plus=x,minus=x;double& p=k%2?plus[k/2].y:plus[k/2].x;double& q=k%2?minus[k/2].y:minus[k/2].x;const double eps=1e-6;p+=eps;q-=eps;
            const auto a=work_values(c,plus,estimator),b=work_values(c,minus,estimator);
            for(int e=0;e<l.ne;++e)near(values.w[e].d[k],(a.w[e].v-b.w[e].v)/(2*eps),2e-8,"work derivative finite difference");
            for(int i=0;i<2;++i)near(values.dk[i].d[k],(a.dk[i].v-b.dk[i].v)/(2*eps),2e-8,"kinetic derivative finite difference");
            for(int j=0;j<l.np;++j)near(values.du[j].d[k],(a.du[j].v-b.du[j].v)/(2*eps),2e-8,"potential derivative finite difference");}
        if(estimator==Estimate::Column||estimator==Estimate::Secant)for(int j=0;j<l.np;++j){AD defect=values.du[j];for(int e=0;e<l.ne;++e)if(l.edge[e].potential==j)defect=defect+values.w[e];near(defect.v,0,1e-13,"column work identity");for(double d:defect.d)near(d,0,1e-12,"column gradient identity");}
        if(estimator==Estimate::Row)for(int i=0;i<2;++i){AD defect=values.dk[i];for(int e=0;e<l.ne;++e)if(l.edge[e].body==i)defect=defect-values.w[e];near(defect.v,0,1e-13,"row work identity");for(double d:defect.d)near(d,0,1e-12,"row gradient identity");}
    }
}
Array gaussian(Dense s,Array b,int n){std::array<std::array<long double,N+1>,N> a{};for(int i=0;i<n;++i){for(int j=0;j<n;++j)a[i][j]=s[i][j];a[i][n]=b[i];}
    for(int k=0;k<n;++k){int p=k;for(int i=k+1;i<n;++i)if(std::abs(a[i][k])>std::abs(a[p][k]))p=i;require(std::abs(a[p][k])>1e-18,"oracle singular");std::swap(a[k],a[p]);for(int i=k+1;i<n;++i){const long double f=a[i][k]/a[k][k];for(int j=k;j<=n;++j)a[i][j]-=f*a[k][j];}}
    Array x{};for(int i=n-1;i>=0;--i){long double v=a[i][n];for(int j=i+1;j<n;++j)v-=a[i][j]*x[j];x[i]=double(v/a[i][i]);}return x;
}
void schur_and_scaling(){const auto c=fixture();Array t{},z{};t[0]=.02;z[0]=.03;Vec x=c.predictor;x[0].x+=.01;WorkOptions o;const auto a=assemble_work(c,x,t,z,{},o,1);Array rhs{},d{},scaled{},scale{};
    for(int i=0;i<a.nr;++i){rhs[i]=-a.residual[i];for(int k=0;k<a.nv;++k)rhs[i]+=a.G[i][k]*a.inverse[k]*a.primal[k];scale[i]=std::sqrt(a.S[i][i]);}
    const auto exact=gaussian(a.S,rhs,a.nr);
    for(int pass=0;pass<20000;++pass)for(int i=0;i<a.nr;++i){double b=rhs[i],bs=rhs[i]/scale[i];for(int j=0;j<a.nr;++j)if(i!=j){b-=a.S[i][j]*d[j];bs-=a.S[i][j]/(scale[i]*scale[j])*scaled[j];}d[i]=b/a.S[i][i];scaled[i]=bs;}
    for(int i=0;i<a.nr;++i){near(d[i],exact[i],1e-9,"GS against independent pivoted oracle");near(scaled[i]/scale[i],d[i],1e-10,"diagonal scaling preserves GS iterate");}
    Array dp{};for(int k=0;k<a.nv;++k){double v=-a.primal[k];for(int i=0;i<a.nr;++i)v+=a.G[i][k]*exact[i];dp[k]=v*a.inverse[k];}
    for(int i=0;i<a.nr;++i){double residual=a.residual[i];for(int k=0;k<a.nv;++k)residual+=a.G[i][k]*dp[k];if(i<c.model.count)residual+=c.model.alpha[i]/(c.h*c.h)*exact[i];near(residual,0,1e-10,"full linearized constraint after Schur reconstruction");}
}
void dual_gradients(){const auto c=fixture();for(bool EL:{true,false})for(bool keep:{true,false}){Report a,b;const auto x=energy_step(c,EL,keep,0,.1,a),y=energy_step(c,EL,keep,1,.1,b);for(int i=0;i<2;++i)near(norm(x.x[i]-y.x[i]),0,2e-13,"residual-complete dualization preserves step");}}
void baseline_and_midpoint(){Model m;State s;s.x[0]={1.25,0};s.v[0]={.3,0};const double h=1./480,k=100,E0=energy(m,s);const auto c=context(m,s,h,E0);const auto b=endpoint(c,baseline(c).x);
    const double v=(s.v[0].x-h*k*(s.x[0].x-1))/(1+h*h*k);near(b.v[0].x,v,1e-12,"XPBD exact backward Euler spring");near(b.x[0].x,s.x[0].x+h*v,1e-14,"XPBD position oracle");
    Report r;const auto next=midpoint_step(c,1,true,r);const double q=((1-h*h*k/4)*.25+h*.3)/(1+h*h*k/4);near(next.x[0].x,1+q,1e-14,"midpoint analytic position");near(next.v[0].x,.3-h*k*(q+.25)/2,2e-12,"midpoint analytic velocity");
    for(int n=0;n<3000;++n){Report rr;s=midpoint_step(context(m,s,h,E0),4,true,rr);}near(energy(m,s),E0,1e-9,"midpoint spring physical energy over 3000 steps");
}
void reservoirs_and_cap(){auto c=fixture();const double E0=energy(c.model,c.old),scale=E0-minimum_energy(c.model);Reservoir reservoir;
    for(int n=0;n<100;++n){Report r;WorkOptions o;o.repay_time=.25;const auto next=work_step(c,o,scale,reservoir,r);near(energy(c.model,next)+total(reservoir),E0,2e-11,"physical plus reservoir ledger");require(r.max_extra_ratio<=.1+1e-12,"work correction cap");c=context(c.model,next,c.h,E0);}
}
void midpoint_multibody_identities(){
    for(const auto& test:cases()){const auto& m=test.model;State s=test.initial;const double E0=energy(m,s),scale=E0-minimum_energy(m),h=1./480;
        for(int k=0;k<50;++k){Report r;const auto next=midpoint_step(context(m,s,h,E0),8,true,r,scale);near(r.identity_error,0,3e-12,"kinetic defect equals retained primal residual identity");near(energy(m,next),energy(m,s),5e-12,"multibody midpoint one-step physical energy");
            double change=0,external=0;for(int i=0;i<m.count;++i){change+=m.mass[i]*(cross(next.x[i],next.v[i])-cross(s.x[i],s.v[i]));external+=h*m.mass[i]*cross((next.x[i]+s.x[i])*.5,m.gravity);}
            near(change,external,5e-12,"discrete midpoint angular momentum/torque identity");s=next;}
    }
}
void acceptance_checks(){
    for(auto test:cases()){const auto& m=test.model;const double E0=energy(m,test.initial),scale=E0-minimum_energy(m);State s=test.initial;
        for(int k=0;k<100;++k){State accepted;StableReport guard;Report r;const auto expected=midpoint_step(context(m,s,1./480,E0),4,true,r,scale);require(guarded_midpoint(m,s,1./480,E0,scale,accepted,guard),"guarded baseline step completes");require(guard.rejected==0,"standard small steps need no adaptive rejection");for(int i=0;i<m.count;++i)near(norm(accepted.x[i]-expected.x[i]),0,0,"guard preserves accepted midpoint trajectory");s=accepted;}}
    auto test=cases().back();test.model.alpha={1e-6,1e-6};test.model.mass[1]=.01;const double E0=energy(test.model,test.initial),scale=E0-minimum_energy(test.model);State accepted;StableReport d;
    require(guarded_midpoint(test.model,test.initial,1./240,E0,scale,accepted,d),"stiff guarded step completes");require(d.rejected>0,"stiff step exercises actual adaptive rejection");require(std::abs(energy(test.model,accepted)-E0)<=1e-6*scale,"accepted stiff state respects global physical energy bound");
    State untouched;untouched.x[0]={123,456};StableReport failed;require(!guarded_midpoint(test.model,test.initial,1./240,E0,scale,untouched,failed,12),"retry budget is a visible failure");near(untouched.x[0].x,123,0,"failed solve does not commit an invalid state");
}
int main(){try{derivatives();schur_and_scaling();dual_gradients();baseline_and_midpoint();reservoirs_and_cap();midpoint_multibody_identities();acceptance_checks();std::cout<<"PASS: work derivatives/marginals; Schur oracle/scaling; dualization; baseline/midpoint oracles; reservoir ledger/cap; physical identities; adaptive acceptance/rejection and no invalid commit\n";}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
