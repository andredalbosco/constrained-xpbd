#include "traditional_cases/solver.hpp"
#include <complex>
#include <iostream>
#include <stdexcept>
#include <string>

namespace t=traditional;
using t::V2;
using t::Vec;
namespace {
void require(bool condition,const std::string& message){if(!condition)throw std::runtime_error(message);}
void near(double actual,double expected,double tolerance,const std::string& message){
    if(!std::isfinite(actual)||!std::isfinite(expected)||std::abs(actual-expected)>tolerance)
        throw std::runtime_error(message+"; actual="+std::to_string(actual)+", expected="+std::to_string(expected));
}
void near(V2 actual,V2 expected,double tolerance,const std::string& message){near(t::norm(actual-expected),0,tolerance,message);}
void near(const Vec& actual,const Vec& expected,int count,double tolerance,const std::string& message){for(int i=0;i<count;++i)near(actual[i],expected[i],tolerance,message+" particle "+std::to_string(i));}
double independent_energy(const t::Model& m,const t::State& s){
    double result=0;
    for(int i=0;i<m.count;++i){
        const V2 edge=s.x[i]-(i?s.x[i-1]:V2{});const double C=std::hypot(edge.x,edge.y)-m.length[i];
        result+=.5*m.mass[i]*(s.v[i].x*s.v[i].x+s.v[i].y*s.v[i].y)-m.mass[i]*(m.gravity.x*s.x[i].x+m.gravity.y*s.x[i].y);
        if(m.alpha[i]>0)result+=C*C/(2*m.alpha[i]);
    }
    return result;
}
Vec independent_elastic_gradient(const t::Model& m,const Vec& x){
    Vec gradient{};
    for(int j=0;j<m.count;++j){
        if(m.alpha[j]==0)continue;
        const V2 edge=x[j]-(j?x[j-1]:V2{});const double radius=std::hypot(edge.x,edge.y);
        const V2 force=edge*((radius-m.length[j])/(m.alpha[j]*radius));
        gradient[j]+=force;if(j)gradient[j-1]+=force*(-1);
    }
    return gradient;
}
struct IndependentRows {std::array<Vec,2> J{};Vec g{},objective{};std::array<double,2> C{},rC{},at{};t::Matrix S{};double H=0;};
IndependentRows independent_rows(const t::Context& c,const t::Iterate& z){
    IndependentRows r;const auto& m=c.model;const double h2=c.h*c.h;
    for(int j=0;j<m.count;++j){const V2 edge=z.x[j]-(j?z.x[j-1]:V2{});const double radius=std::hypot(edge.x,edge.y);
        r.J[j][j]=edge/radius;if(j)r.J[j][j-1]=r.J[j][j]*(-1);
        r.C[j]=radius-m.length[j];r.at[j]=m.alpha[j]/h2;r.rC[j]=r.C[j]+r.at[j]*z.lambda[j];}
    const Vec u=independent_elastic_gradient(m,z.x);t::State endpoint;endpoint.x=z.x;
    for(int i=0;i<m.count;++i){
        r.g[i]=(z.x[i]-c.old.x[i])*(m.mass[i]/h2)+u[i]-m.gravity*m.mass[i];
        r.objective[i]=(z.x[i]-c.predictor[i])*(m.mass[i]/h2)+u[i];
        endpoint.v[i]=(z.x[i]-c.old.x[i])/c.h;
    }
    r.H=independent_energy(m,endpoint)-c.target;
    for(int j=0;j<m.count;++j){
        for(int k=0;k<m.count;++k){for(int i=0;i<m.count;++i)r.S[j][k]+=t::dot(r.J[j][i],r.J[k][i])/m.mass[i];if(j==k)r.S[j][k]+=r.at[j];}
        for(int i=0;i<m.count;++i)r.S[j][m.count]+=t::dot(r.J[j][i],r.g[i])/m.mass[i];r.S[m.count][j]=r.S[j][m.count];
    }
    for(int i=0;i<m.count;++i)r.S[m.count][m.count]+=t::dot(r.g[i],r.g[i])/m.mass[i];
    return r;
}
t::Context coupled_fixture(){
    t::Model m;m.count=2;m.mass={2,4};m.length={1,.75};m.alpha={1./32,1./16};m.gravity={0,-2};
    t::State old;old.x={V2{.6875,1.0625},V2{1.515625,-.015625}};old.v={V2{.5,-.25},V2{-.125,.375}};
    return t::context(m,old,1./8,7.0625);
}
t::Iterate apply(const t::Context& c,t::Iterate z,const t::Options& o,t::Diagnostics& d){require(t::update(c,z,o,d),"traditional update failed");return z;}
t::State step(const t::Model& m,const t::State& old,double h,const t::Options& o,t::Diagnostics& diagnostics){
    const auto c=t::context(m,old,h,t::energy(m,old));const t::Iterate start{c.predictor};
    require(start.lambda[0]==0&&start.lambda[1]==0&&start.eta==0,"each baseline substep cold-starts all multipliers");
    const auto z=apply(c,start,o,diagnostics);t::State next;next.x=z.x;for(int i=0;i<m.count;++i)next.v[i]=(next.x[i]-old.x[i])/h;return next;
}

void hand_derived_rows(){
    const auto c=coupled_fixture();const t::Iterate z{c.predictor};const auto r=t::evaluate(c,z);
    near(z.x[0],{.75,1},0,"fixture predictor first mass");near(z.x[1],{1.5,0},0,"fixture predictor second mass");
    near(r.J[0][0],{.6,.8},0,"anchor spring Jacobian");near(r.J[0][1],{},0,"anchor spring sparse zero");
    near(r.J[1][0],{-.6,.8},0,"second spring first incidence");near(r.J[1][1],{.6,-.8},0,"second spring second incidence");
    near(r.C[0],.25,0,"first extension");near(r.C[1],.5,0,"second extension");
    near(r.g[0],{8,8.8},2e-15,"first endpoint energy gradient");near(r.g[1],{.8,5.6},2e-15,"second endpoint energy gradient");
    near(r.objective[0],{0,12.8},2e-15,"first implicit objective gradient");near(r.objective[1],{4.8,-6.4},2e-15,"second implicit objective gradient");
    const t::Matrix expected{{{{2.5,.14,5.92}},{{.14,4.75,.12}},{{5.92,.12,78.72}}}};
    for(int i=0;i<3;++i)for(int j=0;j<3;++j)near(r.S[i][j],expected[i][j],3e-14,"hand-derived Schur entry");
    near(r.H,.5,0,"fixture energy residual");
    const auto acceleration=t::acceleration(c.model,z.x);
    near(acceleration[0],{0,-8.4},2e-15,"first double-spring acceleration");near(acceleration[1],{-1.2,-.4},2e-15,"second double-spring acceleration");
    const V2 f1=acceleration[0]*2,f2=acceleration[1]*4;
    near(f1+f2,{-4.8,-18.4},5e-15,"internal second-spring force cancels from total");
    near(t::cross(z.x[0],f1)+t::cross(z.x[1],f2),-15,5e-15,"central springs exert zero total torque about anchor");
}

void finite_difference_derivatives(){
    const auto c=coupled_fixture();const t::Iterate z{c.predictor};const auto r=t::evaluate(c,z);constexpr double epsilon=1e-6;
    auto potential=[&](const Vec& x){return independent_energy(c.model,{x,{}});};
    auto objective=[&](const Vec& x){
        double value=0;for(int i=0;i<c.model.count;++i)value+=.5*c.model.mass[i]*t::dot(x[i]-c.predictor[i],x[i]-c.predictor[i])/(c.h*c.h);
        for(int j=0;j<c.model.count;++j){const double C=t::norm(x[j]-(j?x[j-1]:V2{}))-c.model.length[j];value+=C*C/(2*c.model.alpha[j]);}return value;
    };
    auto endpoint_energy=[&](const Vec& x){t::State s;s.x=x;for(int i=0;i<c.model.count;++i)s.v[i]=(x[i]-c.old.x[i])/c.h;return independent_energy(c.model,s);};
    for(int particle=0;particle<2;++particle)for(int axis=0;axis<2;++axis){
        Vec plus=z.x,minus=z.x;double& p=axis?plus[particle].y:plus[particle].x;double& q=axis?minus[particle].y:minus[particle].x;p+=epsilon;q-=epsilon;
        const double g=(endpoint_energy(plus)-endpoint_energy(minus))/(2*epsilon);
        near(g,axis?r.g[particle].y:r.g[particle].x,2e-8,"endpoint energy gradient finite difference");
        const double F=(objective(plus)-objective(minus))/(2*epsilon);
        near(F,axis?r.objective[particle].y:r.objective[particle].x,2e-8,"implicit objective gradient finite difference");
        const auto force=t::acceleration(c.model,z.x);const double gradV=(potential(plus)-potential(minus))/(2*epsilon);
        near(gradV,-c.model.mass[particle]*(axis?force[particle].y:force[particle].x),2e-8,"force equals negative physical potential gradient");
        for(int j=0;j<2;++j){const double cp=t::norm(plus[j]-(j?plus[j-1]:V2{}))-c.model.length[j];const double cm=t::norm(minus[j]-(j?minus[j-1]:V2{}))-c.model.length[j];
            near((cp-cm)/(2*epsilon),axis?r.J[j][particle].y:r.J[j][particle].x,2e-10,"constraint Jacobian incidence finite difference");}
    }
}

void static_equilibrium(){
    for(double m1:{.5,1.,2.})for(double m2:{.75,1.,4.}){
        t::Model m;m.count=2;m.mass={m1,m2};m.alpha={.01,.025};m.length={1,.75};m.gravity={0,-9.81};
        const double r1=m.length[0]+m.alpha[0]*(m1+m2)*9.81,r2=m.length[1]+m.alpha[1]*m2*9.81;
        t::State equilibrium;equilibrium.x={V2{0,-r1},V2{0,-r1-r2}};
        const auto a=t::acceleration(m,equilibrium.x);near(a[0],{},1e-13,"static first mass physical force balance");near(a[1],{},1e-13,"static second mass physical force balance");
        const double minimum=-9.81*((m1+m2)*m.length[0]+m2*m.length[1])-.5*9.81*9.81*(m.alpha[0]*(m1+m2)*(m1+m2)+m.alpha[1]*m2*m2);
        near(t::energy(m,equilibrium),minimum,3e-14*(1+std::abs(minimum)),"analytical hanging minimum energy");near(t::minimum_energy(m),minimum,2e-14*(1+std::abs(minimum)),"model minimum energy formula");
    }
    // Document finite-sweep equilibrium error instead of falsely demanding a
    // coupled cold-start GS sweep solve its complete material block exactly.
    t::Model m;m.count=2;m.gravity={0,-9.81};t::State old;old.x={V2{0,-1.1962},V2{0,-2.2943}};
    double previous_error=1;
    for(int hz:{240,480,960}){
        const double h=1./hz;const auto c=t::context(m,old,h,t::energy(m,old));const auto r=independent_rows(c,{c.predictor});
        const double l1=-r.rC[0]/r.S[0][0],l2=(-r.rC[1]-r.S[1][0]*l1)/r.S[1][1];
        Vec expected=c.predictor;for(int i=0;i<2;++i)expected[i]+=(r.J[0][i]*l1+r.J[1][i]*l2)/m.mass[i];
        t::Options o;o.fresh=true;t::Diagnostics d;const auto actual=step(m,old,h,o,d);
        near(actual.x,expected,2,5e-16,"static one-sweep error matches exact scalar GS oracle");
        const double error=t::norm(actual.x[0]-old.x[0]);require(error>0&&error<previous_error/8,"static finite-sweep position error decreases faster than cubic");previous_error=error;
    }
}

void fresh_frozen_baseline(){
    const auto c=coupled_fixture();const t::Iterate start{c.predictor};
    const double l1=-.1,length2=std::sqrt(1.53),l2=-(length2-.75)/4.75;
    const V2 n1{.6,.8},n2{.78/length2,-.96/length2};
    Vec fresh_expected=start.x;fresh_expected[0]+=n1*(l1/2);fresh_expected[0]+=n2*(-l2/2);fresh_expected[1]+=n2*(l2/4);
    t::Options o;o.fresh=true;t::Diagnostics d;const auto fresh=apply(c,start,o,d);
    near(fresh.x,fresh_expected,2,3e-16,"ordinary XPBD refreshes the second spring after the first correction");
    near(fresh.lambda[0],l1,0,"first fresh lambda");near(fresh.lambda[1],l2,3e-17,"second fresh lambda");require(d.evaluations==2&&d.scalar_solves==2,"fresh baseline evaluates each constraint once");
    const double frozen_l2=-.486/4.75;Vec frozen_expected=start.x;
    frozen_expected[0]+=(n1*l1-V2{.6,-.8}*frozen_l2)/2;frozen_expected[1]+=V2{.6,-.8}*(frozen_l2/4);
    o.fresh=false;d={};const auto frozen=apply(c,start,o,d);
    near(frozen.x,frozen_expected,2,3e-16,"frozen baseline uses one linearized coupled GS pass");near(frozen.lambda[1],frozen_l2,3e-17,"frozen second lambda includes first-row Schur response");
    require(d.evaluations==1&&d.scalar_solves==2,"frozen baseline has one evaluation and two scalar solves");
    require(t::norm(frozen.x[0]-fresh.x[0])>1e-4,"noncollinear fixture distinguishes baseline variants");
}

std::array<double,3> gaussian_oracle(const t::Matrix& input,const std::array<double,3>& rhs,int count){
    std::array<std::array<long double,4>,3> a{};for(int i=0;i<count;++i){for(int j=0;j<count;++j)a[i][j]=input[i][j];a[i][count]=rhs[i];}
    for(int k=0;k<count;++k){int pivot=k;for(int i=k+1;i<count;++i)if(std::abs(a[i][k])>std::abs(a[pivot][k]))pivot=i;require(std::abs(a[pivot][k])>1e-24L,"Gaussian oracle singular fixture");std::swap(a[k],a[pivot]);
        for(int i=k+1;i<count;++i){const long double factor=a[i][k]/a[k][k];for(int j=k;j<=count;++j)a[i][j]-=factor*a[k][j];}}
    std::array<long double,3> x{};for(int i=count-1;i>=0;--i){long double v=a[i][count];for(int j=i+1;j<count;++j)v-=a[i][j]*x[j];x[i]=v/a[i][i];}
    return {static_cast<double>(x[0]),static_cast<double>(x[1]),static_cast<double>(x[2])};
}
void normalized_solver_audit(){
    const auto c=coupled_fixture();const auto r=independent_rows(c,{c.predictor});
    for(int count:{2,3})for(double a:{.001,1.,1000.})for(double b:{.01,1.,100.}){
        const std::array<double,3> diagonal{a,b,2.};t::Matrix S{};std::array<double,3> rhs{},answer{};
        for(int i=0;i<count;++i){rhs[i]=(-.2+.3*i)*diagonal[i];for(int j=0;j<count;++j)S[i][j]=r.S[i][j]*diagonal[i]*diagonal[j];}
        require(t::normalized_solve(S,rhs,count,answer),"normalized SPD solve rejected valid fixture");const auto expected=gaussian_oracle(S,rhs,count);
        for(int i=0;i<count;++i)near(answer[i],expected[i],3e-12*(1+std::abs(expected[i])),"normalized solve matches independent pivoted Gaussian elimination");
    }
    t::Matrix rank_one{{{{1,2,0}},{{2,4,0}},{{0,0,1}}}};std::array<double,3> answer{};
    require(!t::normalized_solve(rank_one,{1,2,0},2,answer),"dependent normalized block explicitly rejected");
    rank_one[1][1]=3;require(!t::normalized_solve(rank_one,{1,2,0},2,answer),"indefinite normalized block explicitly rejected");
    rank_one[1][1]=0;require(!t::normalized_solve(rank_one,{1,2,0},2,answer),"zero diagonal block explicitly rejected");
}

std::array<double,3> schedule_oracle(const IndependentRows& r,int n,t::Schedule schedule){
    std::array<double,3> d{};auto L=[&](){for(int i=0;i<n;++i){long double rhs=-r.rC[i];for(int j=0;j<=n;++j)if(j!=i)rhs-=static_cast<long double>(r.S[i][j])*d[j];d[i]=static_cast<double>(rhs/r.S[i][i]);}};
    auto E=[&](){long double rhs=-r.H;for(int j=0;j<n;++j)rhs-=static_cast<long double>(r.S[n][j])*d[j];d[n]=static_cast<double>(rhs/r.S[n][n]);};
    switch(schedule){
    case t::Schedule::Baseline:L();break;
    case t::Schedule::EL:E();L();break;
    case t::Schedule::LE:L();E();break;
    case t::Schedule::LEL:L();E();L();break;
    case t::Schedule::ELE:E();L();E();break;
    case t::Schedule::Jacobi:for(int j=0;j<n;++j)d[j]=-r.rC[j]/r.S[j][j];d[n]=-r.H/r.S[n][n];break;
    case t::Schedule::Joint:{std::array<double,3> rhs{};for(int j=0;j<n;++j)rhs[j]=-r.rC[j];rhs[n]=-r.H;d=gaussian_oracle(r.S,rhs,n+1);break;}
    }return d;
}
void frozen_schedule_oracles(){
    const auto c=coupled_fixture();
    for(bool shifted:{false,true}){
        t::Iterate start{c.predictor};if(shifted){start.x[0]+=V2{.01,-.02};start.x[1]+=V2{-.03,.015};start.lambda={-.01,.02};start.eta=.003;}
        const auto r=independent_rows(c,start);
        for(auto schedule:{t::Schedule::Baseline,t::Schedule::EL,t::Schedule::LE,t::Schedule::Jacobi,t::Schedule::Joint,t::Schedule::LEL,t::Schedule::ELE}){
            t::Options o;o.schedule=schedule;t::Diagnostics stats;const auto actual=apply(c,start,o,stats);const auto d=schedule_oracle(r,2,schedule);
            Vec expected=start.x;for(int i=0;i<2;++i)expected[i]+=(r.J[0][i]*d[0]+r.J[1][i]*d[1]+r.g[i]*d[2])/c.model.mass[i];
            near(actual.x,expected,2,4e-15,"frozen schedule matches independent three-row oracle with dropped r_x");
            for(int j=0;j<2;++j)near(actual.lambda[j],start.lambda[j]+d[j],2e-15,"schedule material multiplier");near(actual.eta,start.eta+d[2],2e-15,"schedule energy multiplier");
            require(stats.evaluations==1,"frozen schedule evaluates one block");
            if(schedule==t::Schedule::EL||schedule==t::Schedule::LEL||schedule==t::Schedule::Joint){
                double last=r.rC[1];for(int j=0;j<3;++j)last+=r.S[1][j]*d[j];near(last,0,2e-14,"material-ending sweep closes last material row");}
            if(schedule==t::Schedule::Joint){for(int i=0;i<3;++i){double residual=i==2?r.H:r.rC[i];for(int j=0;j<3;++j)residual+=r.S[i][j]*d[j];near(residual,0,3e-14,"joint closes all three frozen rows");}}
        }
    }
}

void spring_baseline_closed_form(){
    t::Model m;const t::State initial{{V2{1.25,0},V2{}},{V2{0,0},V2{}}};constexpr double omega=10;
    for(int hz:{240,480,960})for(bool fresh:{false,true}){
        const double h=1./hz;t::State s=initial;t::Options o;o.fresh=fresh;t::Diagnostics stats;
        for(int step_index=0;step_index<hz;++step_index)s=step(m,s,h,o,stats);
        const std::complex<double> expected=.25*std::pow(std::complex<double>(1,-h*omega)/(1+h*h*omega*omega),hz);
        near(s.x[0],{1+expected.real(),0},3e-12,"single spring baseline exact backward-Euler position");near(s.v[0],{omega*expected.imag(),0},3e-11,"single spring baseline exact backward-Euler velocity");
        near(t::energy(m,s)/t::energy(m,initial),std::pow(1+h*h*omega*omega,-hz),2e-11,"single spring baseline expected numerical damping");
        require(stats.evaluations==hz&&stats.scalar_solves==hz,"single-spring baseline has one scalar solve per cold substep");
    }
}

t::State rk4(const t::Model& m,const t::State& state,double h){
    auto derivative=[&](const t::State& s){return t::State{s.v,t::acceleration(m,s.x)};};
    auto add=[&](const t::State& s,const t::State& k,double factor){t::State out=s;for(int i=0;i<m.count;++i){out.x[i]+=k.x[i]*factor;out.v[i]+=k.v[i]*factor;}return out;};
    const auto a=derivative(state),b=derivative(add(state,a,h/2)),c=derivative(add(state,b,h/2)),d=derivative(add(state,c,h));
    t::State out=state;for(int i=0;i<m.count;++i){out.x[i]+=(a.x[i]+b.x[i]*2+c.x[i]*2+d.x[i])*(h/6);out.v[i]+=(a.v[i]+b.v[i]*2+c.v[i]*2+d.v[i])*(h/6);}return out;
}
void double_spring_normal_mode(){
    t::Model m;m.count=2;m.gravity={0,-9.81};const double phi=(1+std::sqrt(5.))/2,omega=10/phi;
    t::State initial;initial.x={V2{0,-1.1962+.01},V2{0,-2.2943+phi*.01}};
    auto exact=[&](double time){t::State s;const double q=.01*std::cos(omega*time),v=-.01*omega*std::sin(omega*time);s.x={V2{0,-1.1962+q},V2{0,-2.2943+phi*q}};s.v={V2{0,v},V2{0,phi*v}};return s;};
    t::State reference=initial;const double E0=t::energy(m,initial);double reference_energy_error=0;
    for(int i=0;i<10000;++i){reference=rk4(m,reference,1e-4);reference_energy_error=std::max(reference_energy_error,std::abs(t::energy(m,reference)-E0));}
    const auto expected=exact(1);near(reference.x,expected.x,2,2e-12,"RK4 double-spring reference matches analytic normal mode");near(reference.v,expected.v,2,2e-11,"RK4 reference velocity matches analytic normal mode");
    require(reference_energy_error<2e-12,"independent double-spring reference conserves physical energy");
    double previous_rms=1;
    for(int hz:{240,480,960}){
        t::State fresh=initial,frozen=initial;t::Options a;a.fresh=true;t::Options b;t::Diagnostics da,db;double sum=0;
        for(int i=1;i<=hz;++i){fresh=step(m,fresh,1./hz,a,da);frozen=step(m,frozen,1./hz,b,db);const auto truth=exact(double(i)/hz);
            for(int j=0;j<2;++j)sum+=t::dot(fresh.x[j]-truth.x[j],fresh.x[j]-truth.x[j]);
            near(fresh.x,frozen.x,2,1e-12,"collinear fresh and frozen baseline trajectories agree");}
        const double rms=std::sqrt(sum/hz);require(rms<previous_rms*.65,"double-spring baseline converges as cold substeps shrink");previous_rms=rms;
        require(da.evaluations==2*hz&&db.evaluations==hz&&da.scalar_solves==2*hz&&db.scalar_solves==2*hz,"double baseline evaluation/solve cost matches one GS sweep");
    }
}

void pendulum_baseline(){
    t::Model m;m.alpha[0]=0;m.gravity={0,-9.81};const double theta0=.65;
    t::State initial;initial.x[0]={std::sin(theta0),-std::cos(theta0)};
    double theta=theta0,velocity=0;constexpr int reference_steps=20000;const double dt=1./reference_steps;
    auto angular_acceleration=[](double angle){return -9.81*std::sin(angle);};
    for(int i=0;i<reference_steps;++i){const double a=velocity,b=angular_acceleration(theta),c=velocity+dt*b/2,d=angular_acceleration(theta+dt*a/2),e=velocity+dt*d/2,f=angular_acceleration(theta+dt*c/2),g=velocity+dt*f,k=angular_acceleration(theta+dt*e);theta+=dt*(a+2*c+2*e+g)/6;velocity+=dt*(b+2*d+2*f+k)/6;}
    const V2 reference{std::sin(theta),-std::cos(theta)};double previous_error=1;
    for(int hz:{240,480,960}){
        t::State s=initial;t::Options o;o.fresh=true;t::Diagnostics diagnostics;
        for(int i=0;i<hz;++i){s=step(m,s,1./hz,o,diagnostics);near(t::norm(s.x[0]),1,5e-16,"ordinary hard-rod XPBD preserves length each substep");}
        const double error=t::norm(s.x[0]-reference);require(error<previous_error*.7,"baseline pendulum trajectory converges with smaller cold substeps");previous_error=error;
        require(diagnostics.evaluations==hz&&diagnostics.scalar_solves==hz,"pendulum baseline one scalar projection per substep");
    }
}

double independent_mass_norm(const t::Model& m,const Vec& v){double squared=0;for(int i=0;i<m.count;++i)squared+=m.mass[i]*t::dot(v[i],v[i]);return std::sqrt(squared);}
double independent_motion_scale(const t::Context& c){
    Vec acceleration=independent_elastic_gradient(c.model,c.old.x);
    for(int i=0;i<c.model.count;++i)acceleration[i]=c.model.gravity-acceleration[i]/c.model.mass[i];
    return c.h*independent_mass_norm(c.model,c.old.v)+c.h*c.h*independent_mass_norm(c.model,acceleration);
}
Vec independent_displacement(const t::Model& m,const IndependentRows& r,const std::array<double,3>& d){
    Vec dx{};for(int i=0;i<m.count;++i){dx[i]=r.g[i]*d[m.count];for(int j=0;j<m.count;++j)dx[i]+=r.J[j][i]*d[j];dx[i]=dx[i]/m.mass[i];}return dx;
}

void repair_cap_guards(){
    auto coupled=coupled_fixture();t::Model spring;spring.alpha[0]=1./64;
    t::State released;released.x[0]={1.25,0};const auto release=t::context(spring,released,1./8,2);
    t::State moving;moving.x[0]={1.125,0};moving.v[0]={-1,0};const auto second=t::context(spring,moving,1./8,2);
    long capped=0,skipped=0,accepted=0;
    for(const auto& original:{coupled,release,second})for(double target_offset:{-1.,0.,1.}){
        auto c=original;c.target+=target_offset;const int n=c.model.count;const t::Iterate start{c.predictor};const auto r=independent_rows(c,start);
        const auto base=schedule_oracle(r,n,t::Schedule::Baseline);const Vec baseline_dx=independent_displacement(c.model,r,base);
        double H0=r.H;for(int i=0;i<n;++i)H0+=t::dot(r.g[i],baseline_dx[i]);
        const double motion=independent_motion_scale(c);near(t::motion_scale(c),motion,2e-14*(1+motion),"cap uses old physical mass-weighted velocity and acceleration");
        for(auto schedule:{t::Schedule::EL,t::Schedule::LE,t::Schedule::Jacobi,t::Schedule::Joint}){
            const auto raw=schedule_oracle(r,n,schedule);const Vec candidate=independent_displacement(c.model,r,raw);Vec extra{};double change=0;
            for(int i=0;i<n;++i){extra[i]=candidate[i]-baseline_dx[i];change+=t::dot(r.g[i],extra[i]);}
            const double size=independent_mass_norm(c.model,extra);
            const double tolerance=128*std::numeric_limits<double>::epsilon()*(1+std::abs(c.target)+std::abs(r.H));
            double blend=1;
            if(std::abs(H0)<=tolerance||H0*change>=0||change==0)blend=0;
            else if(std::abs(change)>.5*std::abs(H0))blend=.5*std::abs(H0)/std::abs(change);
            if(blend*size>.1*motion)blend=.1*motion/size;
            std::array<double,3> expected_increment{};for(int j=0;j<=n;++j)expected_increment[j]=base[j]+blend*(raw[j]-base[j]);
            const Vec expected_dx=independent_displacement(c.model,r,expected_increment);
            t::Options o;o.schedule=schedule;o.policy=t::Policy::RepairCap;t::Diagnostics d;const auto actual=apply(c,start,o,d);
            Vec actual_extra{},actual_dx{};
            for(int i=0;i<n;++i){actual_dx[i]=actual.x[i]-start.x[i];actual_extra[i]=actual_dx[i]-baseline_dx[i];}
            near(actual_dx,expected_dx,n,5e-15,"repair-cap accepted position reconstructs blended multiplier pair");
            for(int j=0;j<n;++j)near(actual.lambda[j],expected_increment[j],3e-15,"repair-cap blends each material multiplier toward baseline");
            near(actual.eta,expected_increment[n],3e-15,"repair-cap scales energy multiplier consistently");
            near(d.last_blend,blend,3e-14,"independent repair-cap blend");near(d.last_H0,H0,3e-14,"repair reference is energy after frozen baseline");
            require(independent_mass_norm(c.model,actual_extra)<=.1*motion+1e-14,"repair-cap actual extra respects mass-norm displacement limit");
            const double predicted_change=blend*change;
            require(H0*predicted_change<=1e-13,"accepted repair has non-worsening predicted sign");
            require(std::abs(predicted_change)<=.5*std::abs(H0)+1e-13,"accepted repair corrects at most half baseline linear residual");
            near(d.last_predicted_H,H0+predicted_change,3e-13,"reported linear energy prediction includes accepted correction");
            require(d.evaluations==1,"repair-cap keeps one frozen linearization");
            capped+=d.caps;skipped+=d.skips;accepted+=blend>0;
        }
    }
    require(capped>0&&skipped>0&&accepted>0,"repair fixtures exercise caps, rejections, and accepted corrections");
}

void raw_energy_direction_and_cap(){
    auto coupled=coupled_fixture();t::Model spring;spring.alpha[0]=1./64;
    t::State old;old.x[0]={1.25,0};const auto released=t::context(spring,old,1./8,2);
    t::State moving;moving.x[0]={1.125,0};moving.v[0]={-1,0};const auto second=t::context(spring,moving,1./8,2);
    long rejected=0,capped=0;
    for(const auto& original:{coupled,released,second})for(double target_offset:{-1.,0.,1.})for(bool cap:{false,true}){
        auto c=original;c.target+=target_offset;const int n=c.model.count;const t::Iterate start{c.predictor};const auto r=independent_rows(c,start);
        std::array<double,3> d{};d[n]=-r.H/r.S[n][n];Vec raw{};double work=0;
        for(int i=0;i<n;++i){raw[i]=r.g[i]*(d[n]/c.model.mass[i]);work+=t::dot(r.objective[i],raw[i]);}
        const double size=independent_mass_norm(c.model,raw),motion=independent_motion_scale(c);
        const bool oppose=work>0;const bool clip=!oppose&&cap&&size>.1*motion;
        if(oppose)d[n]=0;else if(clip)d[n]*=.1*motion/size;
        // The material GS pass must see the accepted/clipped energy multiplier.
        for(int j=0;j<n;++j){double rhs=-r.rC[j];for(int k=0;k<=n;++k)if(k!=j)rhs-=r.S[j][k]*d[k];d[j]=rhs/r.S[j][j];}
        const Vec expected_dx=independent_displacement(c.model,r,d);
        t::Options o;o.schedule=t::Schedule::EL;o.policy=cap?t::Policy::RawEnergyGateCap:t::Policy::RawEnergyGate;t::Diagnostics stats;
        const auto actual=apply(c,start,o,stats);Vec actual_dx{};for(int i=0;i<n;++i)actual_dx[i]=actual.x[i]-start.x[i];
        near(actual_dx,expected_dx,n,5e-15,"raw energy gate/cap precedes dependent material GS pass");
        for(int j=0;j<n;++j)near(actual.lambda[j],d[j],3e-15,"raw energy controls induce matching material multiplier");near(actual.eta,d[n],3e-15,"accepted raw energy multiplier");
        require(stats.skips==int(oppose)&&stats.caps==int(clip),"raw gate and cap counters distinguish rejection from scaling");
        Vec accepted_raw{};double accepted_work=0;for(int i=0;i<n;++i){accepted_raw[i]=r.g[i]*(actual.eta/c.model.mass[i]);accepted_work+=t::dot(r.objective[i],accepted_raw[i]);}
        require(accepted_work<=1e-13,"raw-energy gate rejects objective-opposed increments");
        if(cap)require(independent_mass_norm(c.model,accepted_raw)<=.1*motion+1e-14,"optional raw cap limits energy displacement before material response");
        rejected+=stats.skips;capped+=stats.caps;
    }
    require(rejected>0&&capped>0,"raw energy fixtures exercise uphill rejection and partial clipping");
}
} // namespace
int main(){
    try{
        hand_derived_rows();finite_difference_derivatives();static_equilibrium();fresh_frozen_baseline();normalized_solver_audit();frozen_schedule_oracles();
        spring_baseline_closed_form();double_spring_normal_mode();pendulum_baseline();
        repair_cap_guards();raw_energy_direction_and_cap();
        std::cout<<"PASS: independent double-spring forces/Jacobians/Schur block; static equilibrium; fresh/frozen XPBD; normalized solve; frozen schedules; spring/pendulum/double-spring baseline convergence; repair/raw-energy guard properties\n";
    }catch(const std::exception& e){std::cerr<<"FAIL: "<<e.what()<<'\n';return 1;}
}
