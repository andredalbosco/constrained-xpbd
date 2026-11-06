#pragma once
#include "solvers.hpp"
namespace local_work {
// A separate objective variant; all work coordinates, penalties and balances
// are inherited unchanged. Only the inertial position target and U weight vary.
struct CrossProblem : ExplicitProblem {
    bool quarter;
    V aold;
    CrossProblem(const Model&m,const State&s,double h,const Config&cfg,double E,const std::string&target,bool q)
      : ExplicitProblem(m,s,h,cfg,E,target),quarter(q),aold(masses(m).cwiseInverse().cwiseProduct(force(m,s.x))) {}
    Evaluation evaluate(const V&z)const {
        auto e=ExplicitProblem::evaluate(z);
        if(!quarter)return e;
        using namespace independent_detail;
        for(int i=0;i<p.d;++i){
            Jet r=(e.x[i]-p.old.x[i])/p.h-p.old.v[i];
            Jet rq=r-.25*p.h*aold[i];
            e.objective=e.objective+.5*p.m.mass[i/p.m.dim]/p.E*(rq*rq-r*r);
        }
        Jet U(0,p.n);
        for(int j=0;j<p.m.count;++j){
            if(!p.m.rigid){Jet c=length(p.m,jedge(p.m,e.x,j))-p.m.length[j];U=U+.5*p.m.stiffness[j]*c*c;}
            for(int k=0;k<p.m.dim;++k)U=U-p.m.mass[j]*p.m.gravity[k]*e.x[j*p.m.dim+k];
        }
        e.objective=e.objective-.75/p.E*U;
        return e;
    }
};

}
