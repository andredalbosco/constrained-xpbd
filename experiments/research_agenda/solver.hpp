#pragma once
// Small, isolated experiments. All solves below are linearized scalar GS visits.
#include "../traditional_cases/solver.hpp"
#include <algorithm>
#include <limits>

namespace agenda {
using namespace traditional;
constexpr int N=12;
using Array=std::array<double,N>;
using Dense=std::array<Array,N>;
struct AD {
    double v=0;std::array<double,4> d{};
    AD()=default;AD(double value):v(value){}
    static AD variable(double value,int i){AD a(value);a.d[i]=1;return a;}
};
inline AD operator+(AD a,const AD& b){a.v+=b.v;for(int k=0;k<4;++k)a.d[k]+=b.d[k];return a;}
inline AD operator-(AD a,const AD& b){a.v-=b.v;for(int k=0;k<4;++k)a.d[k]-=b.d[k];return a;}
inline AD operator-(AD a){return AD(0)-a;}
inline AD operator*(const AD& a,const AD& b){AD c(a.v*b.v);for(int k=0;k<4;++k)c.d[k]=a.d[k]*b.v+a.v*b.d[k];return c;}
inline AD operator/(const AD& a,const AD& b){AD c(a.v/b.v);for(int k=0;k<4;++k)c.d[k]=(a.d[k]-c.v*b.d[k])/b.v;return c;}
inline AD root(const AD& a){AD b(std::sqrt(a.v));for(int k=0;k<4;++k)b.d[k]=a.d[k]/(2*b.v);return b;}
struct D2 {AD x{},y{};};
inline D2 operator+(D2 a,D2 b){return {a.x+b.x,a.y+b.y};}
inline D2 operator-(D2 a,D2 b){return {a.x-b.x,a.y-b.y};}
inline D2 operator*(D2 a,AD b){return {a.x*b,a.y*b};}
inline AD dp(D2 a,D2 b){return a.x*b.x+a.y*b.y;}
inline AD length(D2 a){return root(dp(a,a));}
inline D2 constant(V2 a){return {AD(a.x),AD(a.y)};}

enum class Estimate {Lagged,Column,Row,Secant};
struct Edge {int body=0,potential=0;V2 force{};};
struct WorkLayout {int np=0,ne=0;std::array<int,4> kind{},body{};std::array<Edge,6> edge{};};
inline WorkLayout layout(const Model& m,const State& old){
    WorkLayout l;l.np=m.count;
    for(int j=0;j<m.count;++j){
        l.kind[j]=0;l.body[j]=j;const V2 e=old.x[j]-(j?old.x[j-1]:V2{});const double r=norm(e);const V2 normal=e/r;
        double scale=m.alpha[j]>0?-(r-m.length[j])/m.alpha[j]:-m.mass[j]*(dot(old.v[j],old.v[j])/r+dot(m.gravity,normal));
        l.edge[l.ne++]={j,j,normal*scale};if(j)l.edge[l.ne++]={j-1,j,normal*(-scale)};
    }
    if(norm(m.gravity)>0)for(int i=0;i<m.count;++i){const int p=l.np++;l.kind[p]=1;l.body[p]=i;l.edge[l.ne++]={i,p,m.gravity*m.mass[i]};}
    return l;
}
struct WorkValues {std::array<AD,2> dk{};std::array<AD,4> du{};std::array<AD,6> w{};};
inline WorkValues work_values(const Context& c,const Vec& x,Estimate estimate){
    const auto& m=c.model;const auto l=layout(m,c.old);WorkValues z;std::array<D2,2> p{},displacement{};
    for(int i=0;i<m.count;++i){p[i]={AD::variable(x[i].x,2*i),AD::variable(x[i].y,2*i+1)};displacement[i]=p[i]-constant(c.old.x[i]);
        z.dk[i]=dp(displacement[i],displacement[i])*(m.mass[i]/(2*c.h*c.h))-.5*m.mass[i]*dot(c.old.v[i],c.old.v[i]);}
    for(int j=0;j<l.np;++j){const int i=l.body[j];if(l.kind[j])z.du[j]=-dp(displacement[i],constant(m.gravity))*m.mass[i];
        else if(m.alpha[i]>0){const AD C=length(p[i]-(i?p[i-1]:D2{}))-m.length[i];const double oldC=norm(c.old.x[i]-(i?c.old.x[i-1]:V2{}))-m.length[i];z.du[j]=(C*C-oldC*oldC)/(2*m.alpha[i]);}}
    for(int e=0;e<l.ne;++e){const auto& q=l.edge[e];D2 force=constant(q.force);
        if(estimate==Estimate::Secant&&!l.kind[q.potential]){const int j=q.potential;const D2 en=p[j]-(j?p[j-1]:D2{});const V2 eo=c.old.x[j]-(j?c.old.x[j-1]:V2{});const AD rn=length(en);const double ro=norm(eo);
            AD scale=m.alpha[j]>0?-(rn+ro-2*m.length[j])/(2*m.alpha[j]):dot(q.force,q.body==j?eo/ro:eo*(-1/ro));
            force=(en+constant(eo))*(scale/(rn+ro));if(q.body!=j)force=force*AD(-1);
        }
        z.w[e]=dp(force,displacement[q.body]);
    }
    if(estimate==Estimate::Column)for(int j=0;j<l.np;++j){AD sum;int count=0;for(int e=0;e<l.ne;++e)if(l.edge[e].potential==j){sum=sum+z.w[e];++count;}
        const AD repair=(-z.du[j]-sum)/count;for(int e=0;e<l.ne;++e)if(l.edge[e].potential==j)z.w[e]=z.w[e]+repair;}
    if(estimate==Estimate::Row)for(int i=0;i<m.count;++i){AD sum;int count=0;for(int e=0;e<l.ne;++e)if(l.edge[e].body==i){sum=sum+z.w[e];++count;}
        const AD repair=(z.dk[i]-sum)/count;for(int e=0;e<l.ne;++e)if(l.edge[e].body==i)z.w[e]=z.w[e]+repair;}
    return z;
}
struct Reservoir {std::array<double,2> k{};std::array<double,4> u{};};
inline double total(const Reservoir& r){double v=0;for(double x:r.k)v+=x;for(double x:r.u)v+=x;return v;}
struct Report {
    long visits=0,evaluations=0,singular=0,caps=0;double rx=0,balance=0,work_defect=0,ledger_error=0,kinematic=0;
    double max_extra_ratio=0,constitutive=0,identity_error=0;std::array<double,6> w{};std::array<double,2> lambda{};double nu=0;
};
inline Iterate baseline(const Context& c){Iterate z{c.predictor};Diagnostics d;Options o;o.fresh=true;if(!update(c,z,o,d))throw std::runtime_error("baseline failed");return z;}
inline State endpoint(const Context& c,const Vec& x){State s;s.x=x;for(int i=0;i<c.model.count;++i)s.v[i]=(x[i]-c.old.x[i])/c.h;return s;}
inline double proposal_blend(const Context& c,const Vec& proposal,const Vec& base,double cap,Report& r,bool repair){
    Vec extra{};for(int i=0;i<c.model.count;++i)extra[i]=proposal[i]-base[i];double blend=1;
    if(repair){const auto initial=evaluate(c,Iterate{c.predictor});Vec dx{};for(int i=0;i<c.model.count;++i)dx[i]=base[i]-c.predictor[i];
        const double H0=initial.H+inner(initial.g,dx,c.model.count),change=inner(initial.g,extra,c.model.count);
        const double tol=128*std::numeric_limits<double>::epsilon()*(1+std::abs(c.target)+std::abs(initial.H));
        if(std::abs(H0)<=tol||H0*change>=0||change==0)blend=0;else if(std::abs(change)>.5*std::abs(H0))blend=.5*std::abs(H0)/std::abs(change);}
    const double size=mass_norm(c.model,extra),scale=motion_scale(c);
    if(cap>0&&size*blend>cap*scale){blend=cap*scale/size;++r.caps;}
    if(scale>0)r.max_extra_ratio=std::max(r.max_extra_ratio,blend*size/scale);return blend;
}

// Rebuild every visited row. keep_rx changes the primal residual treatment.
// dual=1 is residual-complete and identical to the analytic gradient; dual=2
// deliberately omits that residual correction, as an approximation control.
inline State energy_step(const Context& c,bool EL,bool keep_rx,int dual,double cap,Report& stats){
    Iterate z{c.predictor};const auto& m=c.model;
    auto visit=[&](int j){auto r=evaluate(c,z);++stats.evaluations;if(!r.valid)throw std::runtime_error("invalid endpoint row");
        if(dual){Vec g{};for(int i=0;i<m.count;++i){g[i]=(z.x[i]-c.old.x[i])*(m.mass[i]/(c.h*c.h))-m.gravity*m.mass[i];
                for(int k=0;k<m.count;++k)if(m.alpha[k]>0){g[i]+=r.J[k][i]*(-z.lambda[k]/(c.h*c.h));if(dual==1)g[i]+=r.J[k][i]*(r.rC[k]/m.alpha[k]);}}
            r.g=g;for(int i=0;i<m.count;++i){r.R[i]=(z.x[i]-c.predictor[i])*m.mass[i]-r.g[i]*z.eta;for(int k=0;k<m.count;++k)r.R[i]+=r.J[k][i]*(-z.lambda[k]);}}
        Vec R=keep_rx?r.R:Vec{},g=j==m.count?r.g:r.J[j];
        const double denom=weighted(m,g,g)+(j==m.count?0:r.at[j]);
        if(j==m.count){double tolerance=64*std::numeric_limits<double>::epsilon()*r.gradient_scale;for(int i=0;i<m.count;++i)tolerance+=4*std::numeric_limits<double>::epsilon()*m.mass[i]/(c.h*c.h)*(norm(z.x[i])+norm(c.old.x[i]+m.gravity*(c.h*c.h)));
            if(std::sqrt(inner(g,g,m.count))<=tolerance){++stats.singular;return;}}
        const double delta=(-(j==m.count?r.H:r.rC[j])+weighted(m,g,R))/denom;
        for(int i=0;i<m.count;++i)z.x[i]+=(g[i]*delta-R[i])/m.mass[i];if(j==m.count)z.eta+=delta;else z.lambda[j]+=delta;++stats.visits;
    };
    if(EL)visit(m.count);for(int j=0;j<m.count;++j)visit(j);if(!EL)visit(m.count);
    if(cap>0){const auto b=baseline(c);const double f=proposal_blend(c,z.x,b.x,cap,stats,true);for(int i=0;i<m.count;++i){z.x[i]=b.x[i]+(z.x[i]-b.x[i])*f;z.lambda[i]=b.lambda[i]+(z.lambda[i]-b.lambda[i])*f;}z.eta*=f;}
    stats.lambda=z.lambda;stats.nu=z.eta;const auto r=evaluate(c,z);stats.rx=mass_norm(m,minv(m,r.R));return endpoint(c,z.x);
}

struct WorkOptions {Estimate estimate=Estimate::Lagged;double strength=1,cap=.1,repay_time=0;int sweeps=1;bool kinetic_first=false;};
struct WorkAssembly {int nv=0,nr=0;Array residual{},inverse{},primal{};Dense G{},S{};WorkValues values{};WorkLayout graph{};};
inline WorkAssembly assemble_work(const Context& c,const Vec& x,const Array& t,const Array& z,const Reservoir& reservoirs,const WorkOptions& o,double energy_scale){
    WorkAssembly a;const auto& m=c.model;const int nx=2*m.count;a.graph=layout(m,c.old);a.nv=nx+a.graph.ne;a.nr=m.count+a.graph.np+m.count;
    a.values=work_values(c,x,o.estimate);const auto material=evaluate(c,Iterate{x});if(!material.valid)throw std::runtime_error("collapsed work edge");
    for(int i=0;i<m.count;++i){a.inverse[2*i]=a.inverse[2*i+1]=1/m.mass[i];a.primal[2*i]=m.mass[i]*(x[i].x-c.predictor[i].x);a.primal[2*i+1]=m.mass[i]*(x[i].y-c.predictor[i].y);}
    const double penalty=o.strength*c.h*c.h/energy_scale;
    for(int e=0;e<a.graph.ne;++e){a.inverse[nx+e]=1/penalty;a.primal[nx+e]=penalty*t[e];}
    for(int j=0;j<m.count;++j){a.residual[j]=material.C[j]+material.at[j]*z[j];for(int i=0;i<m.count;++i){a.G[j][2*i]=material.J[j][i].x;a.G[j][2*i+1]=material.J[j][i].y;}}
    const double repay=o.repay_time>0?-std::expm1(-c.h/o.repay_time):0;
    for(int p=0;p<a.graph.np;++p){AD b=a.values.du[p];const int row=m.count+p;for(int e=0;e<a.graph.ne;++e)if(a.graph.edge[e].potential==p){b=b+a.values.w[e]+t[e];a.G[row][nx+e]=1;}
        a.residual[row]=b.v-repay*reservoirs.u[p];for(int k=0;k<nx;++k)a.G[row][k]=b.d[k];}
    for(int i=0;i<m.count;++i){AD b=a.values.dk[i];const int row=m.count+a.graph.np+i;for(int e=0;e<a.graph.ne;++e)if(a.graph.edge[e].body==i){b=b-a.values.w[e]-t[e];a.G[row][nx+e]=-1;}
        a.residual[row]=b.v-repay*reservoirs.k[i];for(int k=0;k<nx;++k)a.G[row][k]=b.d[k];}
    for(int k=0;k<a.nv;++k)for(int j=0;j<a.nr;++j)a.primal[k]-=a.G[j][k]*z[j];
    for(int i=0;i<a.nr;++i)for(int j=0;j<a.nr;++j){for(int k=0;k<a.nv;++k)a.S[i][j]+=a.G[i][k]*a.inverse[k]*a.G[j][k];if(i==j&&i<m.count)a.S[i][j]+=material.at[i];}
    return a;
}
inline State work_step(const Context& c,const WorkOptions& o,double scale,Reservoir& reservoir,Report& report){
    const auto& m=c.model;Array t{},z{},delta{},rhs{};Vec x=c.predictor;const auto a=assemble_work(c,x,t,z,reservoir,o,scale);++report.evaluations;
    for(int i=0;i<a.nr;++i){rhs[i]=-a.residual[i];for(int k=0;k<a.nv;++k)rhs[i]+=a.G[i][k]*a.inverse[k]*a.primal[k];}
    std::array<int,N> order{};int ns=0;for(int j=0;j<m.count;++j)order[ns++]=j;
    auto U=[&](){for(int p=0;p<a.graph.np;++p)order[ns++]=m.count+p;};auto K=[&](){for(int i=0;i<m.count;++i)order[ns++]=m.count+a.graph.np+i;};
    if(o.kinetic_first){K();U();}else{U();K();}
    for(int pass=0;pass<o.sweeps;++pass)for(int q=0;q<a.nr;++q){const int i=order[q];double b=rhs[i];for(int j=0;j<a.nr;++j)if(j!=i)b-=a.S[i][j]*delta[j];
        if(!(a.S[i][i]>0)){++report.singular;continue;}delta[i]=b/a.S[i][i];++report.visits;}
    Array dp{};for(int k=0;k<a.nv;++k){double v=-a.primal[k];for(int j=0;j<a.nr;++j)v+=a.G[j][k]*delta[j];dp[k]=v*a.inverse[k];}
    for(int i=0;i<m.count;++i)x[i]+={dp[2*i],dp[2*i+1]};for(int e=0;e<a.graph.ne;++e)t[e]=dp[2*m.count+e];z=delta;
    if(o.cap>0){const auto b=baseline(c);const double blend=proposal_blend(c,x,b.x,o.cap,report,false);for(int i=0;i<m.count;++i){x[i]=b.x[i]+(x[i]-b.x[i])*blend;z[i]=b.lambda[i]+(z[i]-b.lambda[i])*blend;}
        for(int j=m.count;j<a.nr;++j)z[j]*=blend;for(int e=0;e<a.graph.ne;++e)t[e]*=blend;}
    const auto values=work_values(c,x,o.estimate);const double ledger_before=energy(m,c.old)+total(reservoir);
    for(int i=0;i<m.count;++i){double defect=values.dk[i].v;for(int e=0;e<a.graph.ne;++e)if(a.graph.edge[e].body==i)defect-=values.w[e].v+t[e];reservoir.k[i]-=defect;report.work_defect=std::max(report.work_defect,std::abs(defect)/scale);}
    for(int p=0;p<a.graph.np;++p){double defect=values.du[p].v;for(int e=0;e<a.graph.ne;++e)if(a.graph.edge[e].potential==p)defect+=values.w[e].v+t[e];reservoir.u[p]-=defect;report.work_defect=std::max(report.work_defect,std::abs(defect)/scale);}
    for(int e=0;e<a.graph.ne;++e)report.w[e]=values.w[e].v+t[e];
    const auto next=endpoint(c,x);report.ledger_error=std::abs(energy(m,next)+total(reservoir)-ledger_before)/scale;
    const auto final=assemble_work(c,x,t,z,Reservoir{},o,scale);++report.evaluations;
    for(int j=m.count;j<final.nr;++j)report.balance=std::max(report.balance,std::abs(final.residual[j])/scale);
    for(int i=0;i<m.count;++i)report.lambda[i]=z[i];return next;
}

// Midpoint kinematics + discrete spring gradient. The material Jacobian is the
// endpoint derivative, whereas the force direction is a discrete gradient.
// We retain r_x and drop derivatives of the force direction (curvature).
inline State midpoint_step(const Context& c,int iterations,bool keep_rx,Report& report,double energy_scale=1){
    const auto& m=c.model;Vec q{},predictor{};for(int i=0;i<m.count;++i)q[i]=predictor[i]=c.old.x[i]+c.old.v[i]*(c.h/2)+m.gravity*(c.h*c.h/4);
    std::array<double,2> lambda{};
    for(int it=0;it<iterations;++it){std::array<Vec,2> right{},left{};std::array<double,2> residual{},at{},dl{};Vec R{},dq{};
        for(int i=0;i<m.count;++i)R[i]=(q[i]-predictor[i])*m.mass[i];
        for(int j=0;j<m.count;++j){const V2 eo=c.old.x[j]-(j?c.old.x[j-1]:V2{}),em=q[j]-(j?q[j-1]:V2{}),en=em*2-eo;const double ro=norm(eo),rn=norm(en);
            if(!(rn>1e-14))throw std::runtime_error("collapsed midpoint spring");V2 force,grad;
            if(m.alpha[j]>0){force=em*(2/(ro+rn));grad=en/rn;at[j]=4*m.alpha[j]/(c.h*c.h);residual[j]=.5*(ro+rn)-m.length[j]+at[j]*lambda[j];}
            else{force=em;grad=en;residual[j]=.25*(ro*ro+rn*rn-2*m.length[j]*m.length[j]);}
            right[j][j]=force;left[j][j]=grad;if(j){right[j][j-1]=force*(-1);left[j][j-1]=grad*(-1);}
            for(int i=0;i<m.count;++i)R[i]+=right[j][i]*(-lambda[j]);
        }
        if(keep_rx)for(int i=0;i<m.count;++i)dq[i]=R[i]*(-1/m.mass[i]);
        for(int j=0;j<m.count;++j){const double den=weighted(m,left[j],right[j])+at[j];if(std::abs(den)<1e-20){++report.singular;continue;}
            dl[j]=-(residual[j]+inner(left[j],dq,m.count))/den;for(int i=0;i<m.count;++i)dq[i]+=right[j][i]*(dl[j]/m.mass[i]);++report.visits;}
        for(int j=0;j<m.count;++j)lambda[j]+=dl[j];for(int i=0;i<m.count;++i)q[i]+=dq[i];++report.evaluations;
    }
    State s;for(int i=0;i<m.count;++i){s.x[i]=q[i]*2-c.old.x[i];s.v[i]=(s.x[i]-c.old.x[i])*(2/c.h)-c.old.v[i];report.kinematic=std::max(report.kinematic,norm(s.x[i]-c.old.x[i]-(s.v[i]+c.old.v[i])*(c.h/2)));}
    Vec R{},impulse{};for(int i=0;i<m.count;++i){R[i]=(q[i]-predictor[i])*m.mass[i];impulse[i]=m.gravity*(c.h*m.mass[i]);}
    for(int j=0;j<m.count;++j){const V2 eo=c.old.x[j]-(j?c.old.x[j-1]:V2{}),en=s.x[j]-(j?s.x[j-1]:V2{});const double ro=norm(eo),rn=norm(en);
        const V2 force=m.alpha[j]>0?(en+eo)/(rn+ro):(en+eo)*.5;
        double residual=m.alpha[j]>0?.5*(ro+rn)-m.length[j]+4*m.alpha[j]/(c.h*c.h)*lambda[j]:.25*(ro*ro+rn*rn-2*m.length[j]*m.length[j]);
        report.constitutive=std::max(report.constitutive,std::abs(residual));
        R[j]+=force*(-lambda[j]);impulse[j]+=force*(4*lambda[j]/c.h);if(j){R[j-1]+=force*lambda[j];impulse[j-1]+=force*(-4*lambda[j]/c.h);}
        double potential_change=0;if(m.alpha[j]>0)potential_change=((rn-m.length[j])*(rn-m.length[j])-(ro-m.length[j])*(ro-m.length[j]))/(2*m.alpha[j]);
        const double spring_work=dot(force,en-eo)*4*lambda[j]/(c.h*c.h);report.work_defect=std::max(report.work_defect,std::abs(potential_change+spring_work)/energy_scale);
    }
    for(int i=0;i<m.count;++i){const V2 average=(s.v[i]+c.old.v[i])*.5;const double defect=.5*m.mass[i]*(dot(s.v[i],s.v[i])-dot(c.old.v[i],c.old.v[i]))-dot(average,impulse[i]);
        report.work_defect=std::max(report.work_defect,std::abs(defect)/energy_scale);report.identity_error=std::max(report.identity_error,std::abs(defect-dot(average,R[i])*4/c.h)/energy_scale);}
    report.rx=mass_norm(m,minv(m,R));report.lambda=lambda;return s;
}
} // namespace agenda
