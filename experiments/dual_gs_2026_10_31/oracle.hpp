#pragma once
// Independent dense oracle for the material-eliminated hour objective.
// This header intentionally does not include or call the multiplier GS engine.
#ifndef DUAL_GS_HOUR_PROBLEM_HEADER
#define DUAL_GS_HOUR_PROBLEM_HEADER "/home/blacknoir/Projects/constrained-xpbd/experiments/hour_2026_10_31/local_solver/problem.hpp"
#endif
#include DUAL_GS_HOUR_PROBLEM_HEADER
#include <chrono>
#include <limits>
#include <string>

namespace dual_gs {
using compose::V;
using compose::Mat;
using compose::Model;
using compose::State;

struct OracleOptions {
    int evaluation_budget=768;
    int max_iterations=96;
    double primal_tolerance=1e-10;
    double stationarity_tolerance=1e-9;
    double rank_tolerance=1e-10;
    double curvature_tolerance=1e-8;
    bool lm_fallback=true;
};
struct OracleResult {
    State state;
    V z, y_unlifted, y_lifted;
    double objective=std::numeric_limits<double>::quiet_NaN();
    double primal=std::numeric_limits<double>::infinity();
    double stationarity=std::numeric_limits<double>::infinity();
    double residual=std::numeric_limits<double>::infinity();
    double curvature=std::numeric_limits<double>::quiet_NaN();
    double smallest_jacobian_singular_value=std::numeric_limits<double>::quiet_NaN();
    double milliseconds=0;
    int evaluations=0, iterations=0, fallbacks=0;
    int jacobian_rank=0, constraint_count=0, tangent_dimension=0;
    bool converged=false, regular=false, curvature_certified=false;
    bool local_minimum_certified=false;
    std::string reason="not_started";
};

namespace oracle_detail {
using Evaluation=compose::independent_detail::Evaluation;
inline V constraints(const Evaluation&e){return compose::independent_detail::residuals(e);}
inline Mat row_jacobian(const Evaluation&e){return compose::independent_detail::jacobian(e);}
inline V stationarity_vector(const Evaluation&e,const V&y){
    return e.objective.g-row_jacobian(e).transpose()*y;
}
inline V residual(const Evaluation&e,const V&y){
    V r(e.objective.g.size()+y.size());
    r.head(e.objective.g.size())=stationarity_vector(e,y);
    r.tail(y.size())=-constraints(e);
    return r;
}
inline Mat lagrangian_hessian(const Evaluation&e,const V&y){
    Mat H=e.objective.H;
    for(int i=0;i<y.size();++i)H-=y[i]*e.rows[i].value.H;
    return .5*(H+H.transpose()).eval();
}
inline Mat newton_matrix(const Evaluation&e,const V&y){
    const int n=e.objective.g.size(),q=y.size();
    Mat J=row_jacobian(e),A=Mat::Zero(n+q,n+q);
    A.topLeftCorner(n,n)=lagrangian_hessian(e,y);
    A.topRightCorner(n,q)=-J.transpose();
    A.bottomLeftCorner(q,n)=-J;
    return A;
}
inline bool finite(const Evaluation&e){
    if(!std::isfinite(e.objective.a)||!e.objective.g.allFinite()||!e.objective.H.allFinite())return false;
    for(const auto&r:e.rows)if(!std::isfinite(r.value.a)||!r.value.g.allFinite()||!r.value.H.allFinite())return false;
    return true;
}
inline void measure(const Evaluation&e,OracleResult&r,const OracleOptions&o){
    r.objective=e.objective.a;
    r.primal=constraints(e).lpNorm<Eigen::Infinity>();
    r.stationarity=stationarity_vector(e,r.y_unlifted).norm();
    r.residual=std::hypot(constraints(e).norm(),r.stationarity);
    r.converged=finite(e)&&r.z.allFinite()&&r.y_unlifted.allFinite()
        &&r.primal<=o.primal_tolerance&&r.stationarity<=o.stationarity_tolerance;
}
// Stable public order: material/length, kinetic, storage, gravity, tangency.
inline V lift_dual(const hour_local::Problem&p,const Evaluation&e,const V&z,const V&y){
    const auto&m=p.p.m;
    const int N=m.count, gravity=compose::has_gravity(m)?N:0;
    V lifted=V::Zero(3*N+gravity+(m.rigid?N:0));
    if(!m.rigid){
        const State state=p.state(z);
        const double mref=m.mass.mean();
        for(int j=0;j<N;++j){
            if(!(m.stiffness[j]>0))throw std::runtime_error("oracle_lift_requires_positive_stiffness");
            const double C=compose::edge_length(m,compose::edge(m,state.x,j))-m.length[j];
            const double alpha=4*mref/(p.p.h*p.p.h*m.stiffness[j]);
            lifted[j]=-(C/p.p.xs)/alpha;
        }
    }
    for(int i=0;i<int(e.rows.size());++i){
        const auto&row=e.rows[i];int target=-1;
        if(row.kind==compose::independent_detail::Position)target=row.index;
        else if(row.kind==compose::independent_detail::Kinetic)target=N+row.index;
        else if(row.kind==compose::independent_detail::Storage)target=2*N+row.index;
        else if(row.kind==compose::independent_detail::GravityStorage)target=3*N+row.index;
        else if(row.kind==compose::independent_detail::Tangency)target=3*N+gravity+row.index;
        else throw std::runtime_error("oracle_unrecognized_constraint_row");
        lifted[target]=y[i];
    }
    return lifted;
}
inline void certify_curvature(const Evaluation&e,OracleResult&r,const OracleOptions&o){
    const Mat J=row_jacobian(e);
    Eigen::JacobiSVD<Mat> svd(J,Eigen::ComputeFullV);
    svd.setThreshold(o.rank_tolerance);
    r.constraint_count=J.rows();r.jacobian_rank=svd.rank();
    r.tangent_dimension=J.cols()-r.jacobian_rank;
    r.regular=r.jacobian_rank==r.constraint_count;
    if(svd.singularValues().size())r.smallest_jacobian_singular_value=svd.singularValues().tail(1)[0];
    if(r.tangent_dimension){
        const Mat Z=svd.matrixV().rightCols(r.tangent_dimension);
        Mat R=Z.transpose()*lagrangian_hessian(e,r.y_unlifted)*Z;
        R=.5*(R+R.transpose()).eval();
        Eigen::SelfAdjointEigenSolver<Mat> eig(R);
        if(eig.info()==Eigen::Success)r.curvature=eig.eigenvalues()[0];
        r.curvature_certified=r.converged&&r.regular&&std::isfinite(r.curvature)
            &&r.curvature>o.curvature_tolerance;
    }
    r.local_minimum_certified=r.converged&&r.regular
        &&(r.tangent_dimension==0||r.curvature_certified);
}
} // namespace oracle_detail

inline OracleResult solve_oracle(const Model&m,const State&old,double h,double E,
                                 const OracleOptions&o={}){
    const auto start=std::chrono::steady_clock::now();
    OracleResult out;
    try{
        if(!(h>0)||!(E>0))throw std::runtime_error("oracle_requires_positive_h_and_E");
        hour_local::Problem p(m,old,h,E);
        out.z=p.initial("midpoint");
        auto e=p.evaluate(out.z);++out.evaluations;
        out.y_unlifted=V::Zero(e.rows.size());
        out.reason="evaluation_budget";
        while(out.evaluations<o.evaluation_budget&&out.iterations<o.max_iterations){
            oracle_detail::measure(e,out,o);
            if(out.converged){out.reason="converged";break;}
            if(!oracle_detail::finite(e)){out.reason="nonfinite_evaluation";break;}
            const V r=oracle_detail::residual(e,out.y_unlifted);
            const Mat A=oracle_detail::newton_matrix(e,out.y_unlifted);
            Eigen::JacobiSVD<Mat> svd(A,Eigen::ComputeThinU|Eigen::ComputeThinV);
            svd.setThreshold(1e-12);
            V delta=svd.solve(-r);
            const double before=r.squaredNorm();
            bool accepted=false;
            auto attempt=[&](const V&direction,double step,int lines){
                if(!direction.allFinite())return false;
                for(int line=0;line<lines&&out.evaluations<o.evaluation_budget;++line){
                    V z=out.z+step*direction.head(p.p.n);
                    V y=out.y_unlifted+step*direction.tail(out.y_unlifted.size());
                    try{
                        auto trial=p.evaluate(z);++out.evaluations;
                        if(oracle_detail::finite(trial)){
                            const double score=oracle_detail::residual(trial,y).squaredNorm();
                            if(score<before*(1-1e-4*step)||score<1e-24){
                                out.z=std::move(z);out.y_unlifted=std::move(y);
                                e=std::move(trial);return true;
                            }
                        }
                    }catch(const std::exception&){++out.evaluations;}
                    step*=.5;
                }
                return false;
            };
            accepted=attempt(delta,std::min(1.,20./std::max(20.,delta.head(p.p.n).norm())),20);
            if(!accepted&&o.lm_fallback){
                const Mat normal=A.transpose()*A;
                const double scale=std::max(1.,normal.diagonal().cwiseAbs().maxCoeff());
                const V right=-A.transpose()*r;
                for(int k=0;k<10&&!accepted&&out.evaluations<o.evaluation_budget;++k){
                    Mat regularized=normal;
                    regularized.diagonal().array()+=scale*std::pow(10.,k-10.);
                    delta=regularized.ldlt().solve(right);++out.fallbacks;
                    accepted=attempt(delta,1.,8);
                }
            }
            ++out.iterations;
            if(!accepted){out.reason="kkt_line_search_stalled";break;}
        }
        oracle_detail::measure(e,out,o);
        if(out.converged)out.reason="converged";
        else if(out.iterations>=o.max_iterations)out.reason="iteration_limit";
        out.state=p.state(out.z);
        out.y_lifted=oracle_detail::lift_dual(p,e,out.z,out.y_unlifted);
        oracle_detail::certify_curvature(e,out,o);
    }catch(const std::exception&e){out.reason=e.what();out.converged=false;}
    out.milliseconds=std::chrono::duration<double,std::milli>(std::chrono::steady_clock::now()-start).count();
    return out;
}
inline OracleResult solve_oracle(const compose::Case&c,double h,double E,const OracleOptions&o={}){
    return solve_oracle(c.model,c.initial,h,E,o);
}
} // namespace dual_gs
