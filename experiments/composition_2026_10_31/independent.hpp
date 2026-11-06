#pragma once
#include "common.hpp"
#include <Eigen/SVD>
#include <numeric>
#include <limits>

// Independent-velocity constrained optimization. The AD implementation keeps
// exact objective/constraint derivatives, including the chosen work estimate.
// No force-based kinetic estimate is treated as an identity: momentum closure
// is reported independently below.
namespace compose { namespace independent_detail {
struct Jet {
    double a=0; V g; Mat H;
    Jet()=default;
    Jet(double value,int n):a(value),g(V::Zero(n)),H(Mat::Zero(n,n)){}
    static Jet variable(double value,int n,int i){Jet r(value,n);r.g[i]=1;return r;}
};
inline Jet operator+(const Jet&a,const Jet&b){Jet r;r.a=a.a+b.a;r.g=a.g+b.g;r.H=a.H+b.H;return r;}
inline Jet operator-(const Jet&a,const Jet&b){Jet r;r.a=a.a-b.a;r.g=a.g-b.g;r.H=a.H-b.H;return r;}
inline Jet operator-(const Jet&a){Jet r;r.a=-a.a;r.g=-a.g;r.H=-a.H;return r;}
inline Jet operator*(const Jet&a,const Jet&b){Jet r;r.a=a.a*b.a;r.g=b.a*a.g+a.a*b.g;r.H=b.a*a.H+a.a*b.H+a.g*b.g.transpose()+b.g*a.g.transpose();return r;}
inline Jet operator*(const Jet&a,double b){Jet r;r.a=a.a*b;r.g=a.g*b;r.H=a.H*b;return r;}
inline Jet operator*(double a,const Jet&b){return b*a;}
inline Jet operator/(const Jet&a,double b){return a*(1/b);}
inline Jet operator+(const Jet&a,double b){Jet r=a;r.a+=b;return r;}
inline Jet operator+(double a,const Jet&b){return b+a;}
inline Jet operator-(const Jet&a,double b){return a+(-b);}
inline Jet operator-(double a,const Jet&b){return -b+a;}
inline Jet inv(const Jet&a){double q=1/a.a;Jet r;r.a=q;r.g=-q*q*a.g;r.H=2*q*q*q*(a.g*a.g.transpose())-q*q*a.H;return r;}
inline Jet operator/(const Jet&a,const Jet&b){return a*inv(b);}
inline Jet square_root(const Jet&a){if(a.a<=1e-24)throw std::runtime_error("collapsed_edge");double q=std::sqrt(a.a);Jet r;r.a=q;r.g=a.g/(2*q);r.H=a.H/(2*q)-(a.g*a.g.transpose())/(4*q*q*q);return r;}
using JV=std::vector<Jet>;
inline Jet dot(const JV&a,const JV&b){Jet r(0,a[0].g.size());for(int i=0;i<int(a.size());++i)r=r+a[i]*b[i];return r;}
inline JV plus(const JV&a,const JV&b){JV r=a;for(int i=0;i<int(a.size());++i)r[i]=r[i]+b[i];return r;}
inline JV minus(const JV&a,const JV&b){JV r=a;for(int i=0;i<int(a.size());++i)r[i]=r[i]-b[i];return r;}
inline JV times(const JV&a,const Jet&b){JV r=a;for(auto&x:r)x=x*b;return r;}
inline JV times(const JV&a,double b){JV r=a;for(auto&x:r)x=x*b;return r;}
inline JV segment(const JV&a,int start,int size){return JV(a.begin()+start,a.begin()+start+size);}
inline JV jedge(const Model&m,const JV&x,int j){JV r=segment(x,j*m.dim,m.dim);if(j)r=minus(r,segment(x,(j-1)*m.dim,m.dim));return r;}
inline Jet length(const Model&m,const JV&e){return m.dim==1?e[0]:square_root(dot(e,e));}
inline JV direction(const Model&m,const JV&e){if(m.dim==1)return JV{Jet(1,e[0].g.size())};return times(e,inv(length(m,e)));}
enum Kind {Kinetic,Storage,GravityStorage,Position,Tangency,Global};
struct Row {Jet value;Kind kind;int index;};
struct Evaluation {Jet objective;std::vector<Row> rows;JV x,v,work,force;double estimate_defect=0;};
struct Problem {
    const Model&m;const State&old;double h,E,vs,xs;const Config&cfg;int d,n,wcount;bool hard,penalized,global;
    V old_parts;std::vector<int> work_child,work_parent,work_gravity;
    Problem(const Model&model,const State&s,double dt,const Config&c,double scale):m(model),old(s),h(dt),E(std::max(scale,1e-10)),cfg(c){
        d=dof(m);vs=std::sqrt(E/m.mass.mean());xs=h*vs;hard=cfg.kinematics=="hard_endpoint";penalized=cfg.treatment=="penalized";global=cfg.treatment=="global"||cfg.treatment=="free";
        wcount=0;work_child.resize(m.count);work_parent.assign(m.count,-1);work_gravity.assign(m.count,-1);
        for(int j=0;j<m.count;++j){work_child[j]=wcount++;if(j)work_parent[j]=wcount++;}
        if(has_gravity(m))for(int i=0;i<m.count;++i)work_gravity[i]=wcount++;
        n=d*(hard?1:2)+(penalized?wcount:0);old_parts=potential_parts(m,old.x);
    }
    JV constants(const V&x)const {JV a;for(int i=0;i<x.size();++i)a.emplace_back(x[i],n);return a;}
    // Instantaneous acceleration-level rod reaction. This is a force estimate,
    // not the impulse produced by the constrained optimization.
    JV coefficients(const JV&x,const JV&v)const {
        JV q(m.count,Jet(0,n));
        if(!m.rigid){for(int j=0;j<m.count;++j)q[j]=-m.stiffness[j]*(length(m,jedge(m,x,j))-m.length[j]);return q;}
        std::vector<JV> B(m.count,JV(d,Jet(0,n)));
        JV b(m.count,Jet(0,n));
        for(int j=0;j<m.count;++j){JV e=jedge(m,x,j),u=jedge(m,v,j),normal=direction(m,e);Jet r=length(m,e),speed=dot(normal,u);
            for(int k=0;k<m.dim;++k){B[j][j*m.dim+k]=normal[k];if(j)B[j][(j-1)*m.dim+k]=-normal[k];}
            b[j]=-(dot(u,u)-speed*speed)/r;
            for(int i=0;i<m.count;++i)for(int k=0;k<m.dim;++k)b[j]=b[j]-B[j][i*m.dim+k]*m.gravity[k];
        }
        std::vector<JV>A(m.count,JV(m.count,Jet(0,n)));
        for(int j=0;j<m.count;++j)for(int k=0;k<m.count;++k){for(int i=0;i<d;++i)A[j][k]=A[j][k]+B[j][i]*B[k][i]/m.mass[i/m.dim];if(j==k)A[j][k]=A[j][k]+1e-12;}
        for(int k=0;k<m.count;++k){for(int j=k+1;j<m.count;++j){Jet a=A[j][k]/A[k][k];for(int l=k;l<m.count;++l)A[j][l]=A[j][l]-a*A[k][l];b[j]=b[j]-a*b[k];}}
        for(int j=m.count-1;j>=0;--j){Jet a=b[j];for(int k=j+1;k<m.count;++k)a=a-A[j][k]*q[k];q[j]=a/A[j][j];}return q;
    }
    Evaluation evaluate(const V&z)const {
        Evaluation out;out.objective=Jet(0,n);out.x.resize(d);out.v.resize(d);JV oldx=constants(old.x),oldv=constants(old.v);
        for(int i=0;i<d;++i){out.x[i]=old.x[i]+xs*Jet::variable(z[i],n,i);out.v[i]=hard?vs*Jet::variable(z[i],n,i):vs*Jet::variable(z[d+i],n,d+i);}
        JV up(m.count+(has_gravity(m)?m.count:0),Jet(0,n)),kin(m.count,Jet(0,n));
        for(int i=0;i<d;++i){Jet a=(out.x[i]-old.x[i])/h-old.v[i];out.objective=out.objective+.5*m.mass[i/m.dim]/E*a*a;
            // Both targets have unit coefficient on the unknown endpoint v,
            // so equal beta means equal velocity compliance in the two modes.
            if(!hard){Jet r=cfg.kinematics=="soft_midpoint"?out.v[i]+old.v[i]-2*(out.x[i]-old.x[i])/h:out.v[i]-(out.x[i]-old.x[i])/h;out.objective=out.objective+.5*m.mass[i/m.dim]/(std::max(cfg.beta,1e-12)*E)*r*r;}
            kin[i/m.dim]=kin[i/m.dim]+.5*m.mass[i/m.dim]*(out.v[i]*out.v[i]-old.v[i]*old.v[i]);
        }
        for(int j=0;j<m.count;++j){if(!m.rigid){Jet a=length(m,jedge(m,out.x,j))-m.length[j];up[j]=.5*m.stiffness[j]*a*a;}out.objective=out.objective+up[j]/E;}
        if(has_gravity(m))for(int i=0;i<m.count;++i){for(int k=0;k<m.dim;++k)up[m.count+i]=up[m.count+i]-m.mass[i]*m.gravity[k]*out.x[i*m.dim+k];out.objective=out.objective+up[m.count+i]/E;}
        JV q0=coefficients(oldx,oldv),q1=coefficients(out.x,out.v);out.force=JV(d,Jet(0,n));out.work=JV(wcount,Jet(0,n));
        auto component_work=[&](const JV&f0,const JV&f1,const JV&fwork,int i)->Jet{
            JV dx=minus(segment(out.x,i*m.dim,m.dim),segment(oldx,i*m.dim,m.dim));
            if(cfg.work=="frozen_power")return h*dot(f0,segment(oldv,i*m.dim,m.dim));
            if(cfg.work=="kinetic_start"||cfg.work=="kinetic_end")return .5*h*dot(cfg.work=="kinetic_start"?f0:f1,plus(segment(out.v,i*m.dim,m.dim),segment(oldv,i*m.dim,m.dim)));
            return dot(fwork,dx);
        };
        for(int j=0;j<m.count;++j){JV e0=jedge(m,oldx,j),e1=jedge(m,out.x,j),n0=direction(m,e0),n1=direction(m,e1),f0=times(n0,q0[j]),f1=times(n1,q1[j]),fw=f1;
            if(cfg.work=="start")fw=f0;
            else if(cfg.work=="force_average")fw=times(plus(f0,f1),.5);
            else if(cfg.work=="separate_average")fw=times(plus(n0,n1),.25*(q0[j]+q1[j]));
            else if(cfg.work=="dg_start"||cfg.work=="dg_end"||cfg.work=="dg_average"){
                JV dg=m.dim==1?JV{Jet(1,n)}:times(plus(e0,e1),inv(length(m,e0)+length(m,e1)));
                Jet q=cfg.work=="dg_start"?q0[j]:(cfg.work=="dg_end"?q1[j]:.5*(q0[j]+q1[j]));fw=times(dg,q);
            }
            out.work[work_child[j]]=component_work(f0,f1,fw,j);
            if(j)out.work[work_parent[j]]=component_work(times(f0,-1),times(f1,-1),times(fw,-1),j-1);
            for(int k=0;k<m.dim;++k){out.force[j*m.dim+k]=out.force[j*m.dim+k]+f1[k];if(j)out.force[(j-1)*m.dim+k]=out.force[(j-1)*m.dim+k]-f1[k];}
        }
        if(has_gravity(m))for(int i=0;i<m.count;++i){JV fg(m.dim,Jet(0,n));for(int k=0;k<m.dim;++k){fg[k]=Jet(m.mass[i]*m.gravity[k],n);out.force[i*m.dim+k]=out.force[i*m.dim+k]+fg[k];}out.work[work_gravity[i]]=component_work(fg,fg,fg,i);}
        for(int j=0;j<m.count;++j){Jet e=up[j]-old_parts[j]+out.work[work_child[j]];if(j)e=e+out.work[work_parent[j]];out.estimate_defect=std::max(out.estimate_defect,std::abs(e.a)/E);}
        if(penalized)for(int j=0;j<wcount;++j){Jet t=Jet::variable(z[d*(hard?1:2)+j],n,d*(hard?1:2)+j);out.work[j]=out.work[j]+E*t;out.objective=out.objective+.5*cfg.gamma*t*t;}
        auto row=[&](const Jet&a,Kind kind,int i){out.rows.push_back({a,kind,i});};
        if(global){Jet e(0,n);for(auto&k:kin)e=e+k;for(int j=0;j<int(up.size());++j)e=e+up[j]-old_parts[j];row(e/E,Global,0);}
        else {
            for(int i=0;i<m.count;++i){Jet e=kin[i]-out.work[work_child[i]];if(i+1<m.count)e=e-out.work[work_parent[i+1]];if(has_gravity(m))e=e-out.work[work_gravity[i]];row(e/E,Kinetic,i);}
            for(int j=0;j<m.count;++j){Jet e=up[j]-old_parts[j]+out.work[work_child[j]];if(j)e=e+out.work[work_parent[j]];
                // For compliant quadratic springs DG-average gives an identity.
                // For rigid rods it is proportional to the positional residual,
                // not an identity off the constraint manifold, so retain it.
                if(!(cfg.treatment=="prescribed"&&!m.rigid&&cfg.work=="dg_average"))row(e/E,Storage,j);
            }
            if(has_gravity(m))for(int i=0;i<m.count;++i){bool spatial=cfg.work!="frozen_power"&&cfg.work!="kinetic_start"&&cfg.work!="kinetic_end";
                if(!(cfg.treatment=="prescribed"&&spatial))row((up[m.count+i]-old_parts[m.count+i]+out.work[work_gravity[i]])/E,GravityStorage,i);
            }
        }
        if(m.rigid)for(int j=0;j<m.count;++j){JV e=jedge(m,out.x,j);row((length(m,e)-m.length[j])/xs,Position,j);if(cfg.tangency)row(dot(direction(m,e),jedge(m,out.v,j))/vs,Tangency,j);}
        return out;
    }
    V initial()const {
        V x=old.x+h*old.v,mm=masses(m);
        if(m.rigid){x+=h*h*mm.cwiseInverse().cwiseProduct(force(m,old.x));for(int j=0;j<m.count;++j){V e=edge(m,x,j);V anchor=j?V(x.segment((j-1)*m.dim,m.dim)):V::Zero(m.dim);x.segment(j*m.dim,m.dim)=anchor+m.length[j]*edge_direction(m,e);}}
        else for(int it=0;it<6;++it){V g=mm.cwiseProduct((x-old.x-h*old.v)/(h*h))-force(m,x);Mat H=(mm/(h*h)).asDiagonal();
            for(int j=0;j<m.count;++j){V e=edge(m,x,j),normal=edge_direction(m,e);double r=edge_length(m,e),c=r-m.length[j];Mat B=m.stiffness[j]*(normal*normal.transpose());if(m.dim>1)B+=m.stiffness[j]*c/r*(Mat::Identity(m.dim,m.dim)-normal*normal.transpose());
                H.block(j*m.dim,j*m.dim,m.dim,m.dim)+=B;if(j){H.block((j-1)*m.dim,(j-1)*m.dim,m.dim,m.dim)+=B;H.block(j*m.dim,(j-1)*m.dim,m.dim,m.dim)-=B;H.block((j-1)*m.dim,j*m.dim,m.dim,m.dim)-=B;}}
            V dx=H.ldlt().solve(-g);if(!dx.allFinite())break;x+=dx;if(dx.norm()<1e-12)break;
        }
        V z=V::Zero(n);z.head(d)=(x-old.x)/xs;if(!hard){V v=(x-old.x)/h;if(cfg.kinematics=="soft_midpoint")v=2*v-old.v;z.segment(d,d)=v/vs;}return z;
    }
};
inline std::vector<std::vector<int>> groups(const Evaluation&e,const Problem&p){
    std::vector<std::vector<int>> out;std::vector<int> all(e.rows.size());std::iota(all.begin(),all.end(),0);
    if(p.global||p.cfg.schedule=="joint")out.push_back(all);
    else if(p.cfg.schedule=="constraint_block"){
        for(int j=0;j<p.m.count;++j){std::vector<int>b;for(int r=0;r<int(e.rows.size());++r){const auto&a=e.rows[r];bool yes=(a.kind==Kinetic&&(a.index==j||(j&&a.index==j-1)))||((a.kind==Storage||a.kind==Position||a.kind==Tangency||a.kind==GravityStorage)&&a.index==j);if(yes)b.push_back(r);}if(!b.empty())out.push_back(b);}
    }else if(p.cfg.schedule=="separate"||p.cfg.schedule=="symmetric"){
        std::vector<int>a,b;for(int r=0;r<int(e.rows.size());++r)(e.rows[r].kind==Kinetic?b:a).push_back(r);if(!a.empty())out.push_back(a);if(!b.empty())out.push_back(b);if(p.cfg.schedule=="symmetric"&&!a.empty())out.push_back(a);
    }else if(p.cfg.schedule=="scalar_cyclic"){for(int r:all)out.push_back({r});}
    else throw std::runtime_error("unknown_schedule");
    if(p.cfg.reverse)std::reverse(out.begin(),out.end());return out;
}
inline V residuals(const Evaluation&e){V g(e.rows.size());for(int i=0;i<g.size();++i)g[i]=e.rows[i].value.a;return g;}
inline Mat jacobian(const Evaluation&e){Mat J(e.rows.size(),e.objective.g.size());for(int i=0;i<J.rows();++i)J.row(i)=e.rows[i].value.g.transpose();return J;}
inline double stationarity(const Evaluation&e,const V&lambda){return (e.objective.g+jacobian(e).transpose()*lambda).norm();}
} // independent_detail

inline StepResult independent_step(const Model&m,const State&old,double h,const Config&cfg,double energy_scale){
    using namespace independent_detail;StepResult result;result.state=old;
    try{
        Problem p(m,old,h,cfg,energy_scale);V z=p.initial();Evaluation current=p.evaluate(z);++result.evaluations;Evaluation frozen=current;V lambda=V::Zero(current.rows.size());auto blocks=groups(current,p);
        for(int it=0;it<cfg.iterations;++it){result.iterations=it+1;
            for(const auto&block:blocks){++result.visits;const Evaluation&der=cfg.relinearize?current:frozen;Mat Jall=jacobian(der);V g=residuals(current);Mat H=der.objective.H;
                if(cfg.curvature=="exact")for(int r=0;r<lambda.size();++r)H+=lambda[r]*der.rows[r].value.H;
                H=.5*(H+H.transpose()).eval();H.diagonal().array()+=1e-10;int q=block.size();Mat A=Mat::Zero(p.n+q,p.n+q);A.topLeftCorner(p.n,p.n)=H;V rhs=V::Zero(p.n+q);
                if(cfg.keep_rx)rhs.head(p.n)=-(current.objective.g+Jall.transpose()*lambda);
                for(int k=0;k<q;++k){A.block(0,p.n+k,p.n,1)=Jall.row(block[k]).transpose();A.block(p.n+k,0,1,p.n)=Jall.row(block[k]);rhs[p.n+k]=-g[block[k]];}
                Eigen::JacobiSVD<Mat> svd(A,Eigen::ComputeThinU|Eigen::ComputeThinV);svd.setThreshold(1e-11);V delta=svd.solve(rhs);if(!delta.allFinite())throw std::runtime_error("nonfinite_kkt");
                // Backtrack the active KKT residual. Shared multipliers are
                // retained across overlapping edge blocks.
                double before=0;for(int r:block)before+=g[r]*g[r];if(cfg.keep_rx)before+=std::pow(stationarity(current,lambda),2);
                double step=std::min(1.,10./std::max(10.,delta.head(p.n).norm()));bool accepted=false;V bestz=z,bestl=lambda;Evaluation best=current;double bestscore=before;
                for(int line=0;line<9;++line){V trialz=z+step*delta.head(p.n),triall=lambda;for(int k=0;k<q;++k)triall[block[k]]+=step*delta[p.n+k];
                    try{Evaluation trial=p.evaluate(trialz);++result.evaluations;double score=0;for(int r:block)score+=trial.rows[r].value.a*trial.rows[r].value.a;if(cfg.keep_rx)score+=std::pow(stationarity(trial,triall),2);
                        if(std::isfinite(score)&&score<bestscore){bestscore=score;bestz=trialz;bestl=triall;best=std::move(trial);accepted=true;}if(score<=before*(1-1e-4*step)||score<1e-22)break;
                    }catch(const std::exception&){}step*=.5;
                }
                if(accepted){z=std::move(bestz);lambda=std::move(bestl);current=std::move(best);}
            }
            if(residuals(current).lpNorm<Eigen::Infinity>()<1e-9&&(!cfg.keep_rx||stationarity(current,lambda)<1e-8))break;
        }
        result.state.x.resize(p.d);result.state.v.resize(p.d);V impulse(p.d),kr(p.d);
        for(int i=0;i<p.d;++i){result.state.x[i]=current.x[i].a;result.state.v[i]=current.v[i].a;double u=(result.state.x[i]-old.x[i])/h;kr[i]=cfg.kinematics=="soft_midpoint"?result.state.v[i]+old.v[i]-2*u:result.state.v[i]-u;impulse[i]=m.mass[i/m.dim]*(result.state.v[i]-old.v[i])-h*current.force[i].a;}
        result.kinematic=kr.norm()/p.vs;result.impulse=impulse.norm()/(m.mass.mean()*p.vs);result.work_defect=current.estimate_defect;result.rx=stationarity(current,lambda);
        for(const auto&r:current.rows)if(r.kind!=Position&&r.kind!=Tangency)result.balance=std::max(result.balance,std::abs(r.value.a));
        if(!finite(result.state)||result.state.x.norm()>1e6||result.state.v.norm()>1e7){result.valid=false;result.reason="nonfinite_or_unbounded";}
        else if(result.balance>1e-6)result.reason="balance_residual";
        if(cfg.treatment=="free"&&result.reason.empty())result.reason="free_work_global_energy_alias";
    }catch(const std::exception&e){result.valid=false;result.reason=e.what();}
    return result;
}
} // compose
