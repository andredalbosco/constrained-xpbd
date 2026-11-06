#pragma once
#include "model.hpp"
#include <Eigen/SVD>
#include <ctime>
#include <functional>
#include <numeric>

using M2=Eigen::Matrix2d;
using M4=Eigen::Matrix4d;
using V2=Eigen::Vector2d;
inline double cpu_us(){return 1e6*std::clock()/CLOCKS_PER_SEC;}

struct Stats {
    int residuals=0,jacobians=0,products=0,model_products=0,outer=0,rejects=0,
        fallbacks=0,direct_solves=0,stages=0,metric_updates=0;
    double setup_us=0,linear_us=0;
    Stats& operator+=(const Stats&o){
        residuals+=o.residuals;jacobians+=o.jacobians;products+=o.products;
        model_products+=o.model_products;outer+=o.outer;rejects+=o.rejects;
        fallbacks+=o.fallbacks;direct_solves+=o.direct_solves;stages+=o.stages;metric_updates+=o.metric_updates;
        setup_us+=o.setup_us;linear_us+=o.linear_us;return *this;
    }
};

// A generic nonsymmetric 2x2 block-tridiagonal matrix. The spring Jacobian is
// block tridiagonal even though an individual block need not be symmetric.
struct Blocks {
    std::vector<M2> d,l,u;
    explicit Blocks(int N=0):d(N,M2::Zero()),l(std::max(0,N-1),M2::Zero()),u(l){}
    V apply(const V&v,bool transpose=false) const {
        int N=d.size();V r=V::Zero(2*N);
        for(int i=0;i<N;++i){
            if(transpose)r.segment<2>(2*i)=d[i].transpose()*v.segment<2>(2*i);
            else r.segment<2>(2*i)=d[i]*v.segment<2>(2*i);
            if(i){if(transpose)r.segment<2>(2*i)+=u[i-1].transpose()*v.segment<2>(2*i-2);
                  else r.segment<2>(2*i)+=l[i-1]*v.segment<2>(2*i-2);}
            if(i+1<N){if(transpose)r.segment<2>(2*i)+=l[i].transpose()*v.segment<2>(2*i+2);
                      else r.segment<2>(2*i)+=u[i]*v.segment<2>(2*i+2);}
        }return r;
    }
    Mat dense()const {int N=d.size();Mat A=Mat::Zero(2*N,2*N);for(int i=0;i<N;++i){A.block<2,2>(2*i,2*i)=d[i];if(i+1<N){A.block<2,2>(2*i,2*i+2)=u[i];A.block<2,2>(2*i+2,2*i)=l[i];}}return A;}
    V solve(const V&b,bool&ok)const {
        int N=d.size();std::vector<M2> c(N);std::vector<V2> z(N);ok=true;
        for(int i=0;i<N;++i){M2 a=d[i];V2 r=b.segment<2>(2*i);if(i){a-=l[i-1]*c[i-1];r-=l[i-1]*z[i-1];}
            if(!a.allFinite()||std::abs(a.determinant())<1e-14*std::max(1e-30,a.squaredNorm())){ok=false;return V::Zero(b.size());}
            M2 inv=a.inverse();z[i]=inv*r;if(i+1<N)c[i]=inv*u[i];
        }V x(b.size());for(int i=N-1;i>=0;--i){x.segment<2>(2*i)=z[i];if(i+1<N)x.segment<2>(2*i)-=c[i]*x.segment<2>(2*i+2);}ok=x.allFinite();return x;
    }
};

struct Eval {V r;Blocks J;bool valid=true;};
inline Eval evaluate(const Problem&p,const State&old,const V&x,double h,bool jac,Stats&stats){
    ++stats.residuals;if(jac)++stats.jacobians;int N=p.m.count;Eval e;e.J=Blocks(jac?N:0);V F=V::Zero(2*N);
    for(int j=0;j<N;++j){F.segment<2>(2*j)=p.m.mass[j]*p.m.gravity;if(jac)e.J.d[j]=p.m.mass[j]*M2::Identity();}
    for(int j=0;j<N;++j){V2 a=old.x.segment<2>(2*j),b=x.segment<2>(2*j);if(j){a-=old.x.segment<2>(2*j-2);b-=x.segment<2>(2*j-2);}double ra=a.norm(),rb=b.norm(),den=ra+rb;
        if(ra<1e-12||rb<1e-12){e.valid=false;e.r=V::Constant(2*N,NAN);return e;}
        double t=1-2*p.m.length[j]/den;V2 fj=-.5*p.m.stiffness[j]*t*(a+b);F.segment<2>(2*j)+=fj;if(j)F.segment<2>(2*j-2)-=fj;
        if(jac){M2 aJ=.25*h*h*p.m.stiffness[j]*(t*M2::Identity()+(2*p.m.length[j]/(den*den))*(a+b)*(b/rb).transpose());e.J.d[j]+=aJ;if(j){e.J.d[j-1]+=aJ;e.J.l[j-1]-=aJ;e.J.u[j-1]-=aJ;}}
    }e.r=p.mass.cwiseProduct(x-old.x-h*old.v)-.5*h*h*F;e.valid=e.r.allFinite();return e;
}

// Every metric is used as L J Q dy = -L R, dx = Q dy. It changes only the
// linear solve, never the force, physical mass, energy, or acceptance norm.
struct Metric {
    std::string name;const Problem*p;double h;V diagonal,graph_d,graph_l,coverage;
    std::vector<M2> vertex,chol_d,chol_l;std::vector<M4> patches;double initial_setup_us=0;
    int refresh_period=1,update_calls=0;bool cache_per_step=false;
    Metric(const Problem&pr,double step,std::string id):name(id),p(&pr),h(step){
        double start=cpu_us();int N=pr.m.count;V D=pr.m.mass,off=V::Zero(std::max(0,N-1));
        for(int i=0;i<N;++i){double k=.25*h*h*pr.m.stiffness[i];D[i]+=k;if(i){D[i-1]+=k;off[i-1]-=k;}}
        diagonal=V::Ones(2*N);
        if(name=="diagonal")for(int i=0;i<N;++i)diagonal.segment<2>(2*i).setConstant(1/std::sqrt(D[i]));
        else if(name!="none"&&name!="graph"&&name!="direct_block")diagonal=pr.misqrt;
        if(name=="graph"){graph_d=V::Zero(N);graph_l=off;for(int i=0;i<N;++i){double a=D[i];if(i){graph_l[i-1]/=graph_d[i-1];a-=graph_l[i-1]*graph_l[i-1];}graph_d[i]=std::sqrt(a);}}
        coverage=V::Ones(N);if((name=="edge_rows"||name=="edge_right")&&N>1){coverage.setConstant(2);coverage[0]=coverage[N-1]=1;coverage=coverage.cwiseSqrt().cwiseInverse();}
        initial_setup_us=cpu_us()-start;
    }
    template<int D> static Eigen::Matrix<double,D,D> balance(const Eigen::Matrix<double,D,D>&A,bool inverse=false){
        using MD=Eigen::Matrix<double,D,D>;Eigen::JacobiSVD<MD> svd(A,Eigen::ComputeFullU|Eigen::ComputeFullV);
        auto sig=svd.singularValues();double floor=.05*std::max(1.,sig[D/2]);
        if(inverse){auto w=sig.cwiseMax(floor).cwiseInverse();return svd.matrixV()*w.asDiagonal()*svd.matrixU().transpose();}
        auto w=(sig.array().square()+floor*floor).sqrt().inverse().matrix().eval();return svd.matrixU()*w.asDiagonal()*svd.matrixU().transpose();
    }
    void update(const Blocks&J,Stats&s){
        if(update_calls++%refresh_period!=0)return;
        if(name=="curvature_graph"){
            ++s.metric_updates;
            double start=cpu_us();int N=J.d.size();std::vector<M2> edges(N),diags(N),lower(std::max(0,N-1));
            for(int i=1;i<N;++i)edges[i]=-J.u[i-1];edges[0]=J.d[0]-p->m.mass[0]*M2::Identity();if(N>1)edges[0]-=edges[1];
            for(int i=0;i<N;++i){M2 sym=.5*(edges[i]+edges[i].transpose());Eigen::SelfAdjointEigenSolver<M2> eig(sym);edges[i]=eig.eigenvectors()*eig.eigenvalues().cwiseMax(0).asDiagonal()*eig.eigenvectors().transpose();diags[i]=p->m.mass[i]*M2::Identity();}
            for(int i=0;i<N;++i){diags[i]+=edges[i];if(i){diags[i-1]+=edges[i];lower[i-1]=-edges[i];}}
            chol_d.resize(N);chol_l.resize(std::max(0,N-1));
            for(int i=0;i<N;++i){M2 S=diags[i];if(i){chol_l[i-1]=chol_d[i-1].triangularView<Eigen::Lower>().solve(lower[i-1].transpose()).transpose();S-=chol_l[i-1]*chol_l[i-1].transpose();}Eigen::LLT<M2> factor(S);if(factor.info()!=Eigen::Success)throw std::runtime_error("positive graph factor failed");chol_d[i]=factor.matrixL();}
            s.setup_us+=cpu_us()-start;return;
        }
        if(name!="vertex_rows"&&name!="vertex_jacobi"&&name!="edge_rows"&&name!="vertex_right"&&name!="edge_right")return;
        ++s.metric_updates;
        double start=cpu_us();int N=J.d.size();
        if((name=="edge_rows"||name=="edge_right")&&N>1){patches.resize(N-1);for(int i=0;i<N-1;++i){M4 A;double a=p->m.mass[i],b=p->m.mass[i+1];A.block<2,2>(0,0)=J.d[i]/a;A.block<2,2>(2,2)=J.d[i+1]/b;A.block<2,2>(0,2)=J.u[i]/std::sqrt(a*b);A.block<2,2>(2,0)=J.l[i]/std::sqrt(a*b);patches[i]=balance<4>(A);}}
        else{vertex.resize(N);for(int i=0;i<N;++i)vertex[i]=balance<2>(M2(J.d[i]/p->m.mass[i]),name=="vertex_jacobi"||name=="vertex_right");}
        s.setup_us+=cpu_us()-start;
    }
    V local(V r)const{
        if((name=="edge_rows"||name=="edge_right")&&p->m.count>1){V q=V::Zero(r.size());for(int i=0;i<p->m.count;++i)r.segment<2>(2*i)*=coverage[i];for(int i=0;i<p->m.count-1;++i)q.segment<4>(2*i)+=patches[i]*r.segment<4>(2*i);for(int i=0;i<p->m.count;++i)q.segment<2>(2*i)*=coverage[i];return q;}
        for(int i=0;i<p->m.count;++i)r.segment<2>(2*i)=vertex[i]*V2(r.segment<2>(2*i));return r;
    }
    V right(const V&v)const{
        if(name=="curvature_graph"){V r=v;for(int i=p->m.count-1;i>=0;--i){if(i+1<p->m.count)r.segment<2>(2*i)-=chol_l[i].transpose()*r.segment<2>(2*i+2);r.segment<2>(2*i)=chol_d[i].transpose().triangularView<Eigen::Upper>().solve(V2(r.segment<2>(2*i)));}return r;}
        if(name=="vertex_right"||name=="edge_right")return diagonal.cwiseProduct(local(v));
        if(name!="graph")return diagonal.cwiseProduct(v);
        V r=v;for(int i=p->m.count-1;i>=0;--i){if(i+1<p->m.count)r.segment<2>(2*i)-=graph_l[i]*r.segment<2>(2*i+2);r.segment<2>(2*i)/=graph_d[i];}return r;
    }
    V left(const V&v)const{
        if(name=="curvature_graph"){V r=v;for(int i=0;i<p->m.count;++i){if(i)r.segment<2>(2*i)-=chol_l[i-1]*r.segment<2>(2*i-2);r.segment<2>(2*i)=chol_d[i].triangularView<Eigen::Lower>().solve(V2(r.segment<2>(2*i)));}return r;}
        if(name=="graph"){V r=v;for(int i=0;i<p->m.count;++i){if(i)r.segment<2>(2*i)-=graph_l[i-1]*r.segment<2>(2*i-2);r.segment<2>(2*i)/=graph_d[i];}return r;}
        V r=diagonal.cwiseProduct(v);
        if(name=="vertex_rows"||name=="vertex_jacobi"||name=="edge_rows")return local(r);return r;
    }
};

inline V gmres(const Blocks&J,const V&r,Metric&metric,int budget,double relative,Stats&s,const std::string&policy="transformed"){
    V b=-metric.left(r);int n=b.size(),kmax=std::min(n,budget);double beta=b.norm();if(beta<1e-30||kmax==0)return V::Zero(n);
    Mat Q=Mat::Zero(n,kmax+1),H=Mat::Zero(kmax+1,kmax),physical;V physical_r;
    if(policy!="transformed"){physical=Mat::Zero(n,kmax);physical_r=metric.p->misqrt.cwiseProduct(r);}
    V cs=V::Zero(kmax),sn=V::Zero(kmax),g=V::Zero(kmax+1);Q.col(0)=b/beta;g[0]=beta;int used=0;
    for(int k=0;k<kmax;++k){V raw=J.apply(metric.right(Q.col(k))),w=metric.left(raw);++s.products;if(policy!="transformed")physical.col(k)=metric.p->misqrt.cwiseProduct(raw);
        for(int j=0;j<=k;++j){H(j,k)=Q.col(j).dot(w);w-=H(j,k)*Q.col(j);}for(int j=0;j<=k;++j){double t=Q.col(j).dot(w);H(j,k)+=t;w-=t*Q.col(j);}
        H(k+1,k)=w.norm();bool breakdown=H(k+1,k)<1e-14;if(!breakdown)Q.col(k+1)=w/H(k+1,k);
        for(int j=0;j<k;++j){double a=cs[j]*H(j,k)+sn[j]*H(j+1,k);H(j+1,k)=-sn[j]*H(j,k)+cs[j]*H(j+1,k);H(j,k)=a;}
        double den=std::hypot(H(k,k),H(k+1,k));cs[k]=den>0?H(k,k)/den:1;sn[k]=den>0?H(k+1,k)/den:0;H(k,k)=cs[k]*H(k,k)+sn[k]*H(k+1,k);H(k+1,k)=0;g[k+1]=-sn[k]*g[k];g[k]=cs[k]*g[k];used=k+1;
        bool stop=std::abs(g[k+1])<=relative*beta;
        if(policy!="transformed"){V y;if(policy=="physical_min")y=physical.leftCols(used).colPivHouseholderQr().solve(-physical_r);else y=H.topLeftCorner(used,used).triangularView<Eigen::Upper>().solve(g.head(used));stop=(physical_r+physical.leftCols(used)*y).norm()<=relative*physical_r.norm();}
        if(stop||breakdown)break;
    }if(!used)return V::Zero(n);V y;if(policy=="physical_min")y=physical.leftCols(used).colPivHouseholderQr().solve(-physical_r);else y=H.topLeftCorner(used,used).triangularView<Eigen::Upper>().solve(g.head(used));return metric.right(Q.leftCols(used)*y);
}

struct Step {State state;Stats stats;bool finite=true,converged=false;std::string failure;double residual=0;};
inline Step solve(const Problem&p,const State&old,double h,Metric&metric,int budget,const V&guess,bool direct=false,int chunk_override=0,const std::string&linear_policy="transformed",const std::string&forcing="fixed",double tolerance_override=0){
    if(metric.cache_per_step)metric.update_calls=0;
    Step out;V x=guess;const double scale=std::max(1e-12,h*std::sqrt(2*p.escale));
    int limit=budget==0?80:12,maxproducts=budget==0?50000:budget,maxeval=budget==0?900:80;double tol=tolerance_override>0?tolerance_override:(budget==0?1e-10:1e-8);
    auto norm=[&](const V&r){return p.misqrt.cwiseProduct(r).norm()/scale;};double previous_residual=0,eta=.5;
    for(int k=0;k<limit&&out.stats.residuals<maxeval;++k){
        Eval e=evaluate(p,old,x,h,true,out.stats);if(!e.valid){out.finite=false;out.failure="invalid_edge";break;}out.residual=norm(e.r);if(out.residual<tol){out.converged=true;break;}if(!direct&&out.stats.products>=maxproducts)break;++out.stats.outer;
        V dx;double start=cpu_us();
        if(direct){bool ok;dx=e.J.solve(-e.r,ok);++out.stats.direct_solves;out.stats.setup_us+=cpu_us()-start;if(!ok){out.failure="singular_block_factor";break;}}
        else{metric.update(e.J,out.stats);int chunk=chunk_override?chunk_override:(budget==0?std::min(160,int(x.size())):std::max(4,budget/4));double relative=budget==0?1e-12:1e-3;
            if(forcing=="loose")relative=.1;
            if(forcing=="ew2"){if(previous_residual>0){double safeguard=.9*eta*eta;eta=.9*std::pow(out.residual/previous_residual,2);if(safeguard>.1)eta=std::max(eta,safeguard);eta=std::max(1e-4,std::min(.5,eta));}relative=eta;}
            start=cpu_us();dx=gmres(e.J,e.r,metric,std::min(chunk,maxproducts-out.stats.products),relative,out.stats,linear_policy);out.stats.linear_us+=cpu_us()-start;}
        if(!dx.allFinite()){out.finite=false;out.failure="nonfinite_update";break;}
        V weighted=p.misqrt.cwiseProduct(e.r),Jdx=e.J.apply(dx);++out.stats.model_products;double phi=.5*weighted.squaredNorm(),slope=weighted.dot(p.misqrt.cwiseProduct(Jdx));
        if(slope>=0||dx.norm()<1e-18){V g=e.J.apply(p.misqrt.cwiseProduct(weighted),true);V Jg=e.J.apply(g);out.stats.model_products+=2;double den=p.misqrt.cwiseProduct(Jg).squaredNorm();if(den<1e-40){out.failure="stationary_nonroot";break;}dx=-g*(g.squaredNorm()/den);Jdx=e.J.apply(dx);++out.stats.model_products;slope=weighted.dot(p.misqrt.cwiseProduct(Jdx));++out.stats.fallbacks;}
        bool accepted=false;double alpha=1;V candidate=x;
        for(int ls=0;ls<16&&out.stats.residuals<maxeval;++ls){candidate=x+alpha*dx;Eval next=evaluate(p,old,candidate,h,false,out.stats);double actual=next.valid?phi-.5*p.misqrt.cwiseProduct(next.r).squaredNorm():-INFINITY;if(next.valid&&actual>=-1e-4*alpha*slope){accepted=true;break;}alpha*=.5;++out.stats.rejects;}
        if(!accepted){out.failure="line_search_stall";break;}previous_residual=out.residual;x=candidate;
    }
    Eval final=evaluate(p,old,x,h,false,out.stats);out.residual=final.valid?norm(final.r):INFINITY;out.converged=out.residual<tol;out.state={x,2*(x-old.x)/h-old.v};out.finite=out.finite&&compose::finite(out.state)&&final.valid;if(!out.converged&&out.failure.empty())out.failure="budget_or_iteration_limit";return out;
}

inline V predictor(const Problem&p,const State&old,double h,const std::string&name,const State*previous=nullptr){
    V x=old.x+h*old.v;if(name=="ballistic")return x;V correction;
    if(name=="secant"&&previous)correction=.5*h*(old.v-previous->v);
    else correction=.5*h*h*p.mass.cwiseInverse().cwiseProduct(force(p.m,old.x));
    // Clip only the guess to an energy-based displacement scale. No force or
    // equation is clipped, and Newton remains free to leave this guess region.
    double norm=p.msqrt.cwiseProduct(correction).norm(),cap=2*h*std::sqrt(2*p.escale);
    if(norm>cap)correction*=cap/norm;return x+correction;
}

inline Step continuation(const Problem&p,const State&old,double h){
    Step total;total.state=old;double tau=0,increment=.125;V lastx=old.x;V tangent=old.v;
    for(int stage=0;stage<100&&tau<1-1e-12;++stage){double next=std::min(1.,tau+increment),dt=h*next;Metric metric(p,dt,"direct_block");V guess=lastx+(next-tau)*h*tangent;
        Step trial=solve(p,old,dt,metric,0,guess,true);total.stats+=trial.stats;++total.stats.stages;
        if(!trial.finite||!trial.converged){increment*=.5;++total.stats.rejects;if(increment<1./4096){total.failure="continuation_minimum_increment";break;}continue;}
        tau=next;lastx=trial.state.x;total.state=trial.state;total.residual=trial.residual;
        if(tau>=1-1e-12){total.converged=true;break;}
        Eval e=evaluate(p,old,lastx,dt,true,total.stats);V F=2*(p.mass.cwiseProduct(lastx-old.x-dt*old.v)-e.r)/(dt*dt);bool ok;double start=cpu_us();tangent=e.J.solve(p.mass.cwiseProduct(old.v)+dt*F,ok);++total.stats.direct_solves;total.stats.setup_us+=cpu_us()-start;if(!ok)tangent=trial.state.v;
        increment=std::min(.25,increment*1.5);
    }
    if(!total.converged){total.failure=total.failure.empty()?"continuation_stage_limit":total.failure;total.finite=false;}
    return total;
}
