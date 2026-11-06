#include "solver.hpp"
#include "../research_agenda/cases.hpp"
#include <iostream>
using namespace implicit;
void require(bool value,const char* message){if(!value)throw std::runtime_error(message);}
void near(double a,double b,double tol,const char* message){if(!std::isfinite(a)||!std::isfinite(b)||std::abs(a-b)>tol*(1+std::abs(b)))throw std::runtime_error(std::string(message)+": "+std::to_string(a)+" vs "+std::to_string(b));}
Array gaussian(Dense s,Array b,int n){std::array<std::array<long double,N+1>,N> a{};for(int i=0;i<n;++i){for(int j=0;j<n;++j)a[i][j]=s[i][j];a[i][n]=b[i];}
    for(int k=0;k<n;++k){int p=k;for(int i=k+1;i<n;++i)if(std::abs(a[i][k])>std::abs(a[p][k]))p=i;require(std::abs(a[p][k])>1e-24L,"singular oracle");std::swap(a[k],a[p]);for(int i=k+1;i<n;++i){const long double f=a[i][k]/a[k][k];for(int j=k;j<=n;++j)a[i][j]-=f*a[k][j];}}
    Array x{};for(int i=n-1;i>=0;--i){long double v=a[i][n];for(int j=i+1;j<n;++j)v-=a[i][j]*x[j];x[i]=double(v/a[i][i]);}return x;
}
void spring_oracle(){Model m;State old;old.x[0]={1.25,0};const double h=1./480,scale=energy(m,old);const auto ctx=context(m,old,h,scale);Config cfg;cfg.EL=false;StepReport r;
    const auto p=energy_proposal(ctx,cfg,scale,r);const double A=1+m.alpha[0]/(h*h),g=25,beta=scale/(h*h*cfg.strength),dl=-.25/A,de=-g*dl/(g*g+beta);
    near(p.x[0].x,1.25+dl+g*de,1e-15,"hand-derived compliant LE scalar solve");
    cfg.EL=true;StepReport rr;const auto el=energy_proposal(ctx,cfg,scale,rr);near(el.x[0].x,1.25+dl,1e-15,"from-rest EL soft row leaves first material step intact");
}
void angular_derivative(){const auto original=cases().back();const auto ctx=context(original.model,original.initial,.125,energy(original.model,original.initial));Vec x=ctx.predictor;const auto row=angular_row(ctx,x);
    for(int i=0;i<2;++i)for(int axis=0;axis<2;++axis){Vec plus=x,minus=x;double& p=axis?plus[i].y:plus[i].x;double& m=axis?minus[i].y:minus[i].x;p+=1e-6;m-=1e-6;const double fd=(angular_row(ctx,plus).value-angular_row(ctx,minus).value)/(2e-6);near(fd,axis?row.gradient[i].y:row.gradient[i].x,2e-9,"angular balance linear-row derivative");}
}
void primal_dual_gauss_newton(){auto test=cases().back();test.model.mass={2,3};const auto ctx=context(test.model,test.initial,.125,energy(test.model,test.initial));const auto r=evaluate(ctx,Iterate{ctx.predictor});const int nx=4;const double scale=energy(test.model,test.initial)-minimum_energy(test.model);
    for(double softness:{-1.,0.,.1}){const bool constrained=softness>=0;Config cfg;cfg.inner=20000;cfg.angular=constrained;cfg.angular_softness=std::max(0.,softness);const double beta=scale/(ctx.h*ctx.h*cfg.strength);Dense A{};Array b{};
        for(int i=0;i<2;++i)for(int axis=0;axis<2;++axis){const int row=2*i+axis;A[row][row]=test.model.mass[i];const double g=axis?r.g[i].y:r.g[i].x;b[row]=-g*r.H/beta;
            for(int k=0;k<nx;++k){const double gk=k%2?r.g[k/2].y:r.g[k/2].x;A[row][k]+=g*gk/beta;}
            for(int j=0;j<2;++j){const double J=axis?r.J[j][i].y:r.J[j][i].x;b[row]-=J*r.C[j]/r.at[j];for(int k=0;k<nx;++k){const double Jk=k%2?r.J[j][k/2].y:r.J[j][k/2].x;A[row][k]+=J*Jk/r.at[j];}}
        }
        int n=nx;if(constrained){const auto l=angular_row(ctx,ctx.predictor);if(softness==0){for(int k=0;k<nx;++k)A[k][nx]=A[nx][k]=k%2?l.gradient[k/2].y:l.gradient[k/2].x;b[nx]=-l.value;n++;}
            else{const double betaL=softness*weighted(test.model,l.gradient,l.gradient);for(int i=0;i<nx;++i){const double gi=i%2?l.gradient[i/2].y:l.gradient[i/2].x;b[i]-=gi*l.value/betaL;for(int j=0;j<nx;++j){const double gj=j%2?l.gradient[j/2].y:l.gradient[j/2].x;A[i][j]+=gi*gj/betaL;}}}}
        const auto exact=gaussian(A,b,n);StepReport stats;const auto p=energy_proposal(ctx,cfg,scale,stats);
        for(int k=0;k<nx;++k)near((k%2?p.x[k/2].y-ctx.predictor[k/2].y:p.x[k/2].x-ctx.predictor[k/2].x),exact[k],2e-11,"dual GS equals independent primal positive Gauss-Newton system");
    }
}
void work_primal_dual_oracle(){auto test=cases().back();const auto ctx=context(test.model,test.initial,.125,energy(test.model,test.initial));const double scale=ctx.target-minimum_energy(test.model);Config cfg;cfg.family=Family::Work;cfg.inner=20000;cfg.repay=.05;
    Reservoir reservoir;reservoir.k[0]=.02;reservoir.u[1]=-.03;WorkOptions o;o.strength=cfg.work_strength;o.repay_time=cfg.repay;const auto a=assemble_work(ctx,ctx.predictor,{}, {},reservoir,o,scale);const auto material=evaluate(ctx,Iterate{ctx.predictor});Dense A{};Array b{};const double beta=scale/(ctx.h*ctx.h*cfg.strength);
    for(int k=0;k<a.nv;++k)A[k][k]=1/a.inverse[k];
    for(int i=0;i<a.nr;++i){const double compliance=i<test.model.count?material.at[i]:beta;for(int k=0;k<a.nv;++k){b[k]-=a.G[i][k]*a.residual[i]/compliance;for(int j=0;j<a.nv;++j)A[k][j]+=a.G[i][k]*a.G[i][j]/compliance;}}
    const auto expected=gaussian(A,b,a.nv);StepReport report;const auto p=work_proposal(ctx,cfg,scale,reservoir,report);for(int k=0;k<4;++k)near(k%2?p.x[k/2].y-ctx.predictor[k/2].y:p.x[k/2].x-ctx.predictor[k/2].x,expected[k],2e-10,"soft work dual GS equals independent positive primal GN system");
    for(int e=0;e<a.graph.ne;++e)near(p.t[e],expected[4+e],2e-10,"soft work primal t reconstruction");
}
void nominal_bounds_and_ledger(){
    for(const auto& test:cases())for(const auto& cfg:configurations()){State old=test.initial;Reservoir reservoir;const double target=energy(test.model,old),scale=target-minimum_energy(test.model);
        for(int step=0;step<20;++step){StepReport report;State next;const auto ctx=context(test.model,old,1./480,target);require(advance(ctx,cfg,scale,reservoir,next,report),"nominal implicit step advances safely");require(energy(test.model,next)<=target+1.01e-10*scale,"global accepted-state energy ceiling");require(report.cap_ratio<=.100001,"actual correction cap");
            if(cfg.family==Family::Work||cfg.family==Family::Baseline)near(energy(test.model,next)+total(reservoir),target,2e-11,"work-defined reservoir ledger after accepted correction");
            for(int i=0;i<test.model.count;++i){near(norm(next.v[i]-(next.x[i]-old.x[i])*480),0,1e-12,"endpoint velocity reconstruction retained");if(test.model.alpha[i]==0)near(norm(next.x[i])-test.model.length[i],0,1e-12,"rod projection remains feasible");}
            old=next;
        }
    }
    const auto test=cases()[0];const double E0=energy(test.model,test.initial),scale=E0-minimum_energy(test.model);auto ctx=context(test.model,test.initial,1./480,E0-scale*10);State unchanged;unchanged.x[0]={123,456};Reservoir reservoir;StepReport r;Config baseline_cfg;baseline_cfg.family=Family::Baseline;
    require(!advance(ctx,baseline_cfg,scale,reservoir,unchanged,r),"infeasible energy ceiling fails visibly");near(unchanged.x[0].x,123,0,"no invalid state commit on failure");near(total(reservoir),0,0,"reservoir not mutated on failed step");
}
int main(){try{spring_oracle();angular_derivative();primal_dual_gauss_newton();work_primal_dual_oracle();nominal_bounds_and_ledger();std::cout<<"PASS: scalar compliance; angular derivative; energy/work primal-dual GN oracles; accepted bound/cap; reservoirs; endpoint kinematics; safe failure\n";}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
