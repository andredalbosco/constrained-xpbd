#pragma once
#include "../local_solver/problem.hpp"
namespace work_reduction_fast {
using namespace hour_local;
// Exact elimination for frozen targets and a uniform quadratic work metric.
// The energy incidence graph must be connected. This implementation stores P
// densely; a tree admits linear-time forward/adjoint flow applications instead.
struct Problem {
    const hour_local::Problem& full;
    struct Metadata {int n,d;}; Metadata p;
    Mat B,P; V target;
    std::vector<int> energy_rows,geometry_rows;
    Problem(const hour_local::Problem& f):full(f),p{2*f.p.d,f.p.d}{
        auto e=f.evaluate(f.initial());
        for(int i=0;i<int(e.rows.size());++i){auto k=e.rows[i].kind;if(k==Kinetic||k==Storage||k==GravityStorage)energy_rows.push_back(i);else geometry_rows.push_back(i);}
        B.resize(energy_rows.size(),f.p.wcount);
        for(int i=0;i<B.rows();++i)B.row(i)=-e.rows[energy_rows[i]].value.g.tail(f.p.wcount).transpose();
        if(B.fullPivLu().rank()!=B.rows()-1)throw std::runtime_error("one_energy_row_requires_connected_graph");
        P=B.completeOrthogonalDecomposition().pseudoInverse();target=f.frozen/f.p.E;
    }
    Jet restricted(const Jet&a)const{Jet b(a.a,p.n);b.g=a.g.head(p.n);b.H=a.H.topLeftCorner(p.n,p.n);return b;}
    struct Data {Evaluation zero; V delta,work; Mat G;};
    Data data(const V&q)const{V z=V::Zero(full.p.n);z.head(p.n)=q;Data a;a.zero=full.evaluate(z);a.delta.resize(B.rows());a.G.resize(B.rows(),p.n);for(int i=0;i<B.rows();++i){const auto&row=a.zero.rows[energy_rows[i]].value;a.delta[i]=row.a;a.G.row(i)=row.g.head(p.n).transpose();}a.work=target+P*(a.delta-B*target);return a;}
    Evaluation evaluate(const V&q)const{
        auto a=data(q);Evaluation e; e.objective=restricted(a.zero.objective);V w=a.work-target,dual=full.cfg.gamma*P.transpose()*w;Mat J=P*a.G;
        e.objective.a+=.5*full.cfg.gamma*(w.squaredNorm()-target.squaredNorm());
        e.objective.g+=a.G.transpose()*dual;
        e.objective.H+=full.cfg.gamma*J.transpose()*J;
        Jet total(0,p.n);for(int i=0;i<B.rows();++i){Jet row=restricted(a.zero.rows[energy_rows[i]].value);total=total+row;e.objective.H+=dual[i]*row.H;}
        e.rows.push_back({total,independent_detail::Global,0});
        for(int i:geometry_rows)e.rows.push_back({restricted(a.zero.rows[i].value),a.zero.rows[i].kind,a.zero.rows[i].index});
        // x/v/work are value-only diagnostics, as in the analytic full oracle.
        e.x=a.zero.x;e.v=a.zero.v;for(int j=0;j<a.work.size();++j){Jet wj;wj.a=full.p.E*a.work[j];e.work.push_back(wj);}return e;
    }
    V expand(const V&q)const{auto a=data(q);V z(full.p.n);z.head(p.n)=q;z.tail(full.p.wcount)=a.work;return z;}
    V initial()const{return full.initial().head(p.n);}
    State state(const V&q)const{return {full.p.old.x+full.p.xs*q.head(p.d),full.p.vs*q.tail(p.d)};}
    V full_dual(const Solve&s)const{auto a=data(s.z);V y=full.cfg.gamma*P.transpose()*(a.work-target);y.array()+=s.dual[0];V out=V::Zero(energy_rows.size()+geometry_rows.size());for(int i=0;i<int(energy_rows.size());++i)out[energy_rows[i]]=y[i];for(int i=0;i<int(geometry_rows.size());++i)out[geometry_rows[i]]=s.dual[i+1];return out;}
};
}
