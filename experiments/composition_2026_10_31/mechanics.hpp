#pragma once
#include "common.hpp"
#include <functional>
#include <numeric>
#include <limits>

namespace compose {
namespace mechanical_detail {
using Residual = std::function<V(const V&)>;
using Jacobian = std::function<Mat(const V&)>;
inline Mat jacobian(const Residual& f,const V& z,StepResult& out) {
    V r=f(z); ++out.evaluations; Mat j(r.size(),z.size());
    for(int k=0;k<z.size();++k){double e=2e-6*std::max(1.0,std::abs(z[k]));V a=z,b=z;a[k]+=e;b[k]-=e;j.col(k)=(f(a)-f(b))/(2*e);out.evaluations+=2;}
    return j;
}
inline std::vector<int> range(int a,int b){std::vector<int> r;for(int k=a;k<b;++k)r.push_back(k);return r;}
// Joint Newton, or nonlinear minimum-norm row-block residual projections.
// The latter are equation-block sweeps, not coordinate descent on an objective.
inline V root(V z,const Residual& f,const Config& c,StepResult& out,
              std::vector<std::vector<int>> blocks={},const Jacobian& derivative={}) {
    const int limit=c.iterations<=0?80:c.iterations; Mat frozen;
    auto getj=[&](const V& p){if(!c.relinearize&&frozen.size())return frozen;Mat j=derivative?derivative(p):jacobian(f,p,out);if(!c.relinearize)frozen=j;return j;};
    bool joint=c.schedule=="joint"||blocks.empty();
    if(c.reverse)std::reverse(blocks.begin(),blocks.end());
    const bool symmetric=c.schedule.find("symmetric")!=std::string::npos;
    if(symmetric){auto back=blocks;std::reverse(back.begin(),back.end());blocks.insert(blocks.end(),back.begin(),back.end());}
    V initial=f(z);++out.evaluations;const double tol=2e-10*(1+initial.norm());
    for(int it=0;it<limit;++it){V r=f(z);++out.evaluations;if(!r.allFinite()||!z.allFinite()){out.valid=false;out.reason="nonfinite_newton";break;}if(r.norm()<tol)break;
        ++out.iterations;
        if(joint){Mat j=getj(z);V dz=j.completeOrthogonalDecomposition().solve(-r);++out.visits;if(!dz.allFinite()){out.valid=false;out.reason="singular_newton";break;}
            double alpha=1;V trial=z+dz;double best=std::numeric_limits<double>::infinity();for(int ls=0;ls<14;++ls){V rr=f(trial);++out.evaluations;best=rr.norm();if(std::isfinite(best)&&(best<r.norm()||alpha<1e-4))break;alpha*=.5;trial=z+alpha*dz;}z=trial;
        }else for(const auto& rows:blocks){r=f(z);++out.evaluations;Mat j=getj(z),a(rows.size(),z.size());V b(rows.size());for(int k=0;k<int(rows.size());++k){a.row(k)=j.row(rows[k]);b[k]=r[rows[k]];}
            V dz=a.completeOrthogonalDecomposition().solve(-b);++out.visits;if(!dz.allFinite()){out.valid=false;out.reason="singular_block";break;}
            double alpha=1;V trial=z+dz;for(int ls=0;ls<12;++ls){V rr=f(trial);++out.evaluations;double score=0;for(int k:rows)score+=rr[k]*rr[k];if(std::isfinite(score)&&(score<b.squaredNorm()||alpha<1e-3))break;alpha*=.5;trial=z+alpha*dz;}z=trial;
        }
        if(z.norm()>1e12){out.valid=false;out.reason="divergent_newton";break;}
    }
    if(c.iterations<=0){V finalr=f(z);++out.evaluations;if(!finalr.allFinite()||finalr.norm()>10*tol){out.valid=false;if(out.reason.empty())out.reason="converged_iteration_limit";}}
    return z;
}
inline Mat gradient(const Model&m,const V&x){Mat g=Mat::Zero(m.count,dof(m));for(int j=0;j<m.count;++j){V n=edge_direction(m,edge(m,x,j));g.block(j,j*m.dim,1,m.dim)=n.transpose();if(j)g.block(j,(j-1)*m.dim,1,m.dim)=-n.transpose();}return g;}
inline V constraints(const Model&m,const V&x){V c(m.count);for(int j=0;j<m.count;++j)c[j]=edge_length(m,edge(m,x,j))-m.length[j];return c;}
inline Mat constraint_hessian(const Model&m,const V&x,const V&q){Mat H=Mat::Zero(dof(m),dof(m));if(m.dim==1)return H;for(int j=0;j<m.count;++j){V e=edge(m,x,j);double r=e.norm();V n=e/r;Mat a=q[j]*(Mat::Identity(m.dim,m.dim)-n*n.transpose())/r;H.block(j*m.dim,j*m.dim,m.dim,m.dim)+=a;if(j){H.block((j-1)*m.dim,(j-1)*m.dim,m.dim,m.dim)+=a;H.block(j*m.dim,(j-1)*m.dim,m.dim,m.dim)-=a;H.block((j-1)*m.dim,j*m.dim,m.dim,m.dim)-=a;}}return H;}
inline V gravity_force(const Model&m){V f(dof(m));for(int i=0;i<m.count;++i)f.segment(i*m.dim,m.dim)=m.mass[i]*m.gravity;return f;}
inline double escale(double s){return std::max(1e-10,std::abs(s));}
inline void finish(const Model&m,const State&old,double h,double scale,StepResult&out){if(!finite(out.state)){out.valid=false;out.reason="nonfinite_state";}out.balance=std::abs(energy(m,out.state)-energy(m,old))/escale(scale);out.kinematic=(out.state.v-(out.state.x-old.x)/h).norm();}

inline StepResult xpbd(const Model&m,const State&old,double h,const Config&c,double scale){
    StepResult out;const int n=dof(m),e=m.count;const V mass=masses(m),inv=mass.cwiseInverse();const double mref=m.mass.mean();
    V pred=old.x+h*old.v+h*h*inv.cwiseProduct(gravity_force(m));V z=V::Zero(n+e);z.head(n)=pred;
    auto physical_rx=[&](const V&a)->V{return (a.head(n)-pred)/h-mref*inv.cwiseProduct(gradient(m,a.head(n)).transpose()*a.tail(e));};
    Residual f=[&](const V&a){V r(n+e);r.head(n)=c.keep_rx?physical_rx(a):V::Zero(n);r.tail(e)=constraints(m,a.head(n))/h;if(!m.rigid)r.tail(e)+=(mref/(h*h))*a.tail(e).cwiseQuotient(m.stiffness);return r;};
    Jacobian j=[&](const V&a){Mat G=gradient(m,a.head(n)),J=Mat::Zero(n+e,n+e);J.topLeftCorner(n,n)=Mat::Identity(n,n)/h;if((c.curvature=="full"||c.curvature=="exact"))J.topLeftCorner(n,n)-=mref*inv.asDiagonal()*constraint_hessian(m,a.head(n),a.tail(e));J.topRightCorner(n,e)=-mref*inv.asDiagonal()*G.transpose();J.bottomLeftCorner(e,n)=G/h;if(!m.rigid)J.bottomRightCorner(e,e)=(mref/(h*h))*m.stiffness.cwiseInverse().asDiagonal();return J;};
    // The local XPBD control is the standard accumulated-multiplier edge sweep.
    if(c.schedule!="joint"){
        V x=pred,lambda=V::Zero(e);int lim=c.iterations<=0?100:c.iterations;std::vector<int> order=range(0,e);if(c.reverse)std::reverse(order.begin(),order.end());if(c.schedule.find("symmetric")!=std::string::npos){auto rev=order;std::reverse(rev.begin(),rev.end());order.insert(order.end(),rev.begin(),rev.end());}
        for(int it=0;it<lim;++it){double maxupdate=0;for(int k:order){Mat G=gradient(m,x);double C=constraints(m,x)[k],alpha=m.rigid?0:1/(h*h*m.stiffness[k]);V g=G.row(k).transpose();double dl=(-C-alpha*lambda[k])/(g.dot(inv.cwiseProduct(g))+alpha);x+=dl*inv.cwiseProduct(g);lambda[k]+=dl;maxupdate=std::max(maxupdate,std::abs(dl));++out.visits;++out.evaluations;}++out.iterations;if(maxupdate<1e-12)break;}z.head(n)=x;z.tail(e)=lambda/(h*mref);
    }else z=root(z,f,c,out,{},j);
    out.state={z.head(n),(z.head(n)-old.x)/h};out.rx=physical_rx(z).norm();finish(m,old,h,scale,out);return out;
}

inline V edge_force(const Model&m,const State&old,const V&x,int j,const std::string&law){
    V a=edge(m,old.x,j),b=edge(m,x,j),na=edge_direction(m,a),nb=edge_direction(m,b);double ra=edge_length(m,a),rb=edge_length(m,b);double qa=-m.stiffness[j]*(ra-m.length[j]),qb=-m.stiffness[j]*(rb-m.length[j]);
    if(law=="start"||law=="startpoint")return qa*na;
    if(law=="end"||law=="endpoint")return qb*nb;
    if(law=="force_average"||law=="average_force")return .5*(qa*na+qb*nb);
    if(law=="separate_average")return .25*(qa+qb)*(na+nb);
    V dg=m.dim==1?V::Ones(1):V((a+b)/(ra+rb));
    if(law=="dg_start")return qa*dg;
    if(law=="dg_end")return qb*dg;
    if(law=="dg_average")return .5*(qa+qb)*dg;
    throw std::runtime_error("unsupported_mechanical_work:"+law);
}
inline V rod_direction(const Model&m,const State&old,const V&x,int j,const std::string&law){V a=edge(m,old.x,j),b=edge(m,x,j);if(law=="start"||law=="startpoint")return edge_direction(m,a);if(law=="end"||law=="endpoint")return edge_direction(m,b);if(law=="force_average"||law=="average_force"||law=="separate_average")return .5*(edge_direction(m,a)+edge_direction(m,b));if(law=="dg_start"||law=="dg_end"||law=="dg_average")return (a+b)/(a.norm()+b.norm());throw std::runtime_error("unsupported_rod_force_direction:"+law);}
inline StepResult force_momentum(const Model&m,const State&old,double h,const Config&c,double scale){
    StepResult out;int n=dof(m),nr=m.rigid?m.count:0;double mref=m.mass.mean();V mass=masses(m),inv=mass.cwiseInverse();
    V z=V::Zero(2*n+nr);z.head(n)=old.x+h*old.v;z.segment(n,n)=old.v;
    auto forces=[&](const V&a,std::vector<V>*parts=nullptr){V f=gravity_force(m);for(int j=0;j<m.count;++j){V fj=m.rigid?(mref/h)*a[2*n+j]*rod_direction(m,old,a.head(n),j,c.work):edge_force(m,old,a.head(n),j,c.work);f.segment(j*m.dim,m.dim)+=fj;if(j)f.segment((j-1)*m.dim,m.dim)-=fj;if(parts)parts->push_back(fj);}return f;};
    Residual f=[&](const V&a){V r(2*n+nr);r.head(n)=(a.head(n)-old.x)/h-.5*(a.segment(n,n)+old.v);r.segment(n,n)=a.segment(n,n)-old.v-h*inv.cwiseProduct(forces(a));if(nr)r.tail(nr)=constraints(m,a.head(n))/h;return r;};
    std::vector<std::vector<int>> blocks;
    if(c.schedule=="separate"||c.schedule=="separate_symmetric"){blocks={range(0,n),range(n,2*n)};if(nr)blocks.push_back(range(2*n,2*n+nr));}
    else for(int i=0;i<m.count;++i){auto b=range(i*m.dim,(i+1)*m.dim);auto vv=range(n+i*m.dim,n+(i+1)*m.dim);b.insert(b.end(),vv.begin(),vv.end());if(nr)b.push_back(2*n+i);blocks.push_back(b);}
    z=root(z,f,c,out,blocks);out.state={z.head(n),z.segment(n,n)};V r=f(z);++out.evaluations;out.rx=r.norm();
    std::vector<V> parts;V total=forces(z,&parts),mean=.5*(out.state.v+old.v);out.impulse=(mass.cwiseProduct(out.state.v-old.v)-h*total).norm();
    double springwork=0;for(int j=0;j<m.count;++j)springwork+=h*parts[j].dot(edge(m,mean,j));
    double gravwork=h*gravity_force(m).dot(mean);out.work_defect=std::abs(potential(m,out.state.x)-potential(m,old.x)+springwork+gravwork)/escale(scale);
    finish(m,old,h,scale,out);return out;
}

inline StepResult verlet_projection(const Model&m,const State&old,double h,const Config&c,double scale,bool tangent){
    StepResult out;V inv=masses(m).cwiseInverse(),vhalf=old.v+.5*h*inv.cwiseProduct(force(m,old.x));V trial=old.x+h*vhalf,x=trial;
    if(m.rigid){Mat G0=gradient(m,old.x);auto position=[&](const V&q)->V{return trial+inv.cwiseProduct(G0.transpose()*q);};Residual r=[&](const V&q)->V{return constraints(m,position(q))/h;};Jacobian j=[&](const V&q)->Mat{return gradient(m,position(q))*inv.asDiagonal()*G0.transpose()/h;};std::vector<std::vector<int>> blocks;for(int k=0;k<m.count;++k)blocks.push_back({k});V q=root(V::Zero(m.count),r,c,out,blocks,j);x=position(q);vhalf=(x-old.x)/h;out.rx=constraints(m,x).norm();}
    V v=vhalf+.5*h*inv.cwiseProduct(force(m,x));if(m.rigid&&tangent){Mat G=gradient(m,x),A=G*inv.asDiagonal()*G.transpose();V q=A.completeOrthogonalDecomposition().solve(-G*v);v+=inv.cwiseProduct(G.transpose()*q);++out.visits;out.rx=std::max(out.rx,(G*v).norm());}
    out.state={x,v};finish(m,old,h,scale,out);return out;
}

// FEPR-style energy projection only: fixed tentative target, no claim to reproduce
// Dinev et al.'s full angular/linear momentum constraints or auxiliary variables.
inline StepResult fepr_projection(const Model&m,const State&old,double h,const Config&c,double scale){
    Config base=c;base.family="xpbd";base.schedule="joint";base.keep_rx=true;base.relinearize=true;base.iterations=16;
    StepResult tentative=xpbd(m,old,h,base,scale),out;int n=dof(m),nc=1+(m.rigid&&c.tangency?2*m.count:0);V mass=masses(m);double mref=m.mass.mean(),es=escale(scale),vs=std::sqrt(es/mref);V target(2*n);target<<tentative.state.x/h,tentative.state.v;
    auto decode=[&](const V&y){return State{h*y.head(n),y.segment(n,n)};};
    Residual g=[&](const V&y){State s=decode(y);V a(nc);a[0]=(energy(m,s)-energy(m,old))/es;if(nc>1){a.segment(1,m.count)=constraints(m,s.x)/(h*vs);a.tail(m.count)=gradient(m,s.x)*s.v/vs;}return a;};
    auto dg=[&](const V&y){State s=decode(y);Mat J=Mat::Zero(nc,2*n);J.block(0,0,1,n)=(-h*force(m,s.x)/es).transpose();J.block(0,n,1,n)=mass.cwiseProduct(s.v).transpose()/es;
        if(nc>1){Mat G=gradient(m,s.x);J.block(1,0,m.count,n)=G/vs;J.block(1+m.count,n,m.count,n)=G/vs;for(int k=0;k<m.count;++k){V q=V::Zero(m.count);q[k]=1;J.block(1+m.count+k,0,1,n)=(h*constraint_hessian(m,s.x,q)*s.v/vs).transpose();}}return J;};
    V D(2*n);D<<mass/mref,mass/(mref*c.beta);V z=V::Zero(2*n+nc);z.head(2*n)=target;
    Residual physical=[&](const V&a){V y=a.head(2*n),r(2*n+nc);Mat G=dg(y);r.head(2*n)=D.cwiseProduct(y-target)+G.transpose()*a.tail(nc);r.tail(nc)=g(y);++out.evaluations;return r;};
    Residual f=[&](const V&a){V r=physical(a);if(!c.keep_rx)r.head(2*n).setZero();return r;};
    Jacobian approx=[&](const V&a){Mat G=dg(a.head(2*n)),J=Mat::Zero(2*n+nc,2*n+nc);J.topLeftCorner(2*n,2*n)=D.asDiagonal();J.topRightCorner(2*n,nc)=G.transpose();J.bottomLeftCorner(nc,2*n)=G;return J;};
    std::vector<std::vector<int>> blocks;blocks.push_back(range(0,2*n));for(int i=0;i<nc;++i)blocks.push_back({2*n+i});
    Jacobian exact=[&](const V&a){return jacobian(physical,a,out);};
    z=root(z,f,c,out,blocks,(c.curvature=="full"||c.curvature=="exact")?exact:approx);out.state=decode(z.head(2*n));out.rx=physical(z).head(2*n).norm();finish(m,old,h,scale,out);out.iterations+=tentative.iterations;out.visits+=tentative.visits;out.evaluations+=tentative.evaluations;return out;
}
// Per-impulse work allocations differ between common-final and sequential
// midpoint velocities, although both total correctly when momentum closes.
struct ImpulseAudit {double common_defect=0,sequential_defect=0,allocation_gap=0,closure=0;};
inline ImpulseAudit audit_impulses(const Model&m,const V&before,const V&after,const std::vector<V>&impulses){
    ImpulseAudit a;V mass=masses(m),inv=mass.cwiseInverse(),sum=V::Zero(before.size()),running=before,mean=.5*(before+after);double common=0,sequential=0;
    for(const V&I:impulses){V next=running+inv.cwiseProduct(I);double wc=I.dot(mean),ws=I.dot(.5*(running+next));common+=wc;sequential+=ws;a.allocation_gap=std::max(a.allocation_gap,std::abs(wc-ws));running=next;sum+=I;}
    double dk=.5*(mass.array()*(after.array().square()-before.array().square())).sum();a.common_defect=std::abs(dk-common);a.sequential_defect=std::abs(dk-sequential);a.closure=(mass.cwiseProduct(after-before)-sum).norm();return a;
}
} // namespace mechanical_detail

inline StepResult mechanics_step(const Model&m,const State&old,double h,const Config&cfg,double scale){
    try{
        if(!(h>0)||!(cfg.beta>0))throw std::runtime_error("nonpositive_step_or_beta");
        if(cfg.family=="xpbd")return mechanical_detail::xpbd(m,old,h,cfg,scale);
        if(cfg.family=="force_momentum"||cfg.family=="midpoint_dg")return mechanical_detail::force_momentum(m,old,h,cfg,scale);
        if(cfg.family=="shake")return mechanical_detail::verlet_projection(m,old,h,cfg,scale,false);
        if(cfg.family=="rattle")return mechanical_detail::verlet_projection(m,old,h,cfg,scale,true);
        if(cfg.family=="fepr_projection"||cfg.family=="fepr")return mechanical_detail::fepr_projection(m,old,h,cfg,scale);
        StepResult out;out.state=old;out.valid=false;out.reason="unknown_mechanics_family";return out;
    }catch(const std::exception&e){StepResult out;out.state=old;out.valid=false;out.reason=e.what();return out;}
}
} // namespace compose
