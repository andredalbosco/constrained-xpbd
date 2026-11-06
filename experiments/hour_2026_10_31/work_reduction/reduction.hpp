#pragma once
#include "cross_problem.hpp"
#include <Eigen/Cholesky>
namespace work_reduction {
using namespace local_work;
using namespace compose::independent_detail;
struct Reduced {
    const CrossProblem&full;
    struct Metadata {int n,d;};Metadata p;
    Mat B,P,Lplus;V target;
    std::vector<int> energy_rows,geometry_rows;
    std::string extension;
    Reduced(const CrossProblem&f,std::string ext="balanced"):full(f),p{2*f.p.d,f.p.d},extension(ext){
        if(!full.is_frozen)throw std::runtime_error("frozen_targets_only");
        V z=full.initial();auto e=full.evaluate(z);for(int i=0;i<int(e.rows.size());++i){auto kind=e.rows[i].kind;if(kind==Kinetic||kind==Storage||kind==GravityStorage)energy_rows.push_back(i);else geometry_rows.push_back(i);}
        B.resize(energy_rows.size(),full.p.wcount);for(int i=0;i<B.rows();++i)B.row(i)=-e.rows[energy_rows[i]].value.g.tail(full.p.wcount).transpose();
        if(B.fullPivLu().rank()!=B.rows()-1)throw std::runtime_error("one_energy_row_requires_connected_graph");
        P=B.completeOrthogonalDecomposition().pseudoInverse();Lplus=P.transpose()*P;
        if(extension=="grounded_first"||extension=="grounded_last"){
            int omit=extension=="grounded_first"?0:B.rows()-1;Mat K(B.rows()-1,B.cols()),select=Mat::Zero(B.rows()-1,B.rows());int j=0;for(int i=0;i<B.rows();++i)if(i!=omit){K.row(j)=B.row(i);select(j,i)=1;++j;}P=K.completeOrthogonalDecomposition().pseudoInverse()*select;
        }
        target=full.frozen/full.p.E;
    }
    Jet restrict_jet(const Jet&a)const {Jet b(a.a,p.n);b.g=a.g.head(p.n);b.H=a.H.topLeftCorner(p.n,p.n);return b;}
    struct Map {V z;Mat derivative;std::vector<Mat> work_hessians;JV changes;};
    Map map(const V&q)const {
        V z=V::Zero(full.p.n);z.head(p.n)=q;auto e=full.evaluate(z);Map a;a.z=z;a.derivative=Mat::Zero(full.p.n,p.n);a.derivative.topRows(p.n).setIdentity();V delta(B.rows());Mat grad(B.rows(),p.n);
        for(int i=0;i<B.rows();++i){Jet change=restrict_jet(e.rows[energy_rows[i]].value);a.changes.push_back(change);delta[i]=change.a;grad.row(i)=change.g.transpose();}
        a.z.tail(full.p.wcount)=target+P*(delta-B*target);a.derivative.bottomRows(full.p.wcount)=P*grad;
        for(int j=0;j<full.p.wcount;++j){Mat H=Mat::Zero(p.n,p.n);for(int i=0;i<B.rows();++i)H+=P(j,i)*a.changes[i].H;a.work_hessians.push_back(H);}return a;
    }
    Evaluation evaluate(const V&q)const {
        Map a=map(q);auto e=full.evaluate(a.z);Evaluation r;r.objective=Jet(e.objective.a,p.n);r.objective.g=a.derivative.transpose()*e.objective.g;r.objective.H=a.derivative.transpose()*e.objective.H*a.derivative;
        for(int j=0;j<full.p.wcount;++j)r.objective.H+=e.objective.g[p.n+j]*a.work_hessians[j];
        Jet total(0,p.n);for(const auto&u:a.changes)total=total+u;r.rows.push_back({total,Global,0});
        for(int i:geometry_rows)r.rows.push_back({restrict_jet(e.rows[i].value),e.rows[i].kind,e.rows[i].index});
        for(auto x:e.x)r.x.push_back(restrict_jet(x));for(auto v:e.v)r.v.push_back(restrict_jet(v));
        for(int j=0;j<full.p.wcount;++j){Jet w(a.z[p.n+j]*full.p.E,p.n);w.g=full.p.E*a.derivative.row(p.n+j).transpose();w.H=full.p.E*a.work_hessians[j];r.work.push_back(w);}return r;
    }
    V initial()const{return full.initial().head(p.n);}
    State state(const V&q)const {return full.state(map(q).z);}
    V full_dual(const Solve&s)const {
        auto a=map(s.z);V residual=a.z.tail(full.p.wcount)-target;
        V y=P.transpose()*(full.cfg.gamma*residual);y.array()+=s.dual[0];
        V dual=V::Zero(energy_rows.size()+geometry_rows.size());for(int i=0;i<int(energy_rows.size());++i)dual[energy_rows[i]]=y[i];for(int i=0;i<int(geometry_rows.size());++i)dual[geometry_rows[i]]=s.dual[1+i];return dual;
    }
    V midpoint_initial()const {
        const auto&m=full.p.m;const auto&s=full.p.old;double h=full.p.h;int d=full.p.d;V M=masses(m);Mat H=M.asDiagonal();
        for(int j=0;j<m.count;++j){V a=V::Zero(d),n=edge_direction(m,edge(m,s.x,j));a.segment(j*m.dim,m.dim)=n;if(j)a.segment((j-1)*m.dim,m.dim)=-n;H+=.25*h*h*m.stiffness[j]*a*a.transpose();}
        V dx=H.ldlt().solve(h*M.cwiseProduct(s.v)+.5*h*h*force(m,s.x));V q(p.n);q.head(d)=dx/full.p.xs;q.tail(d)=(2*dx/h-s.v)/full.p.vs;return q;
    }
};
inline Case chain(int count,double mass_ratio=1,double stiffness=100){Case c=make_case("chain_mixed");c.name="chain"+std::to_string(count)+(mass_ratio>1?"_contrast":"");auto&m=c.model;m.count=count;m.mass=V::Ones(count);m.length=V::Ones(count);m.stiffness=V::Constant(count,stiffness);for(int i=0;i<count;++i)m.mass[i]=i%2?mass_ratio:1;m.mass/=m.mass.mean();c.initial={V::Zero(2*count),V::Zero(2*count)};V prev=V::Zero(2);double downstream=m.mass.sum();for(int i=0;i<count;++i){double angle=.15*std::sin(.7*(i+1)),r=1+downstream*9.81/m.stiffness[i]+.02*(i%2?1:-1);V n(2);n<<std::sin(angle),-std::cos(angle);prev+=r*n;c.initial.x.segment(2*i,2)=prev;c.initial.v[2*i]=.2*std::sin(i+1);c.initial.v[2*i+1]=.1*std::cos(i+1);downstream-=m.mass[i];}return c;}
}
