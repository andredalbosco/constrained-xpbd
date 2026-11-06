#pragma once
#include "analytic_reduction.hpp"
#include <filesystem>
#include <fstream>
#include <iostream>
#include <iomanip>
namespace fs=std::filesystem;
using namespace hour_local;
struct ScaledReduced : work_reduction_fast::Problem {
 double scale;
 ScaledReduced(const hour_local::Problem&f,double s):work_reduction_fast::Problem(f),scale(s){}
 Evaluation evaluate(const V&q)const{auto e=work_reduction_fast::Problem::evaluate(q);e.rows[0].value=e.rows[0].value/scale;return e;}
 V full_dual_scaled(Solve s)const{s.dual[0]/=scale;return full_dual(s);}
 Mat gauss_newton(const V&q,const Evaluation&e)const{auto a=data(q);V y=full.cfg.gamma*P.transpose()*(a.work-target);Mat H=e.objective.H;for(int i=0;i<B.rows();++i)H-=y[i]*a.zero.rows[energy_rows[i]].value.H.topLeftCorner(p.n,p.n);return H;}
};
template<class P,class Curvature> Solve solve(const P&p,int budget,Curvature curvature,const V*dual_seed=nullptr,const V*primal_seed=nullptr){
 auto start=std::chrono::steady_clock::now();Solve out;out.z=primal_seed?*primal_seed:p.initial();auto e=p.evaluate(out.z);++out.evaluations;out.dual=dual_seed?*dual_seed:V::Zero(e.rows.size());
 while(out.evaluations<budget){finish(e,out);if(out.converged)break;Mat J=jacobian(e),H=curvature(out.z,e);for(int i=0;i<out.dual.size();++i)H+=out.dual[i]*e.rows[i].value.H;int n=p.p.n,q=out.dual.size();Mat A=Mat::Zero(n+q,n+q);A.topLeftCorner(n,n)=.5*(H+H.transpose()).eval();A.topRightCorner(n,q)=J.transpose();A.bottomLeftCorner(q,n)=J;V rhs(n+q);rhs.head(n)=-(e.objective.g+J.transpose()*out.dual);rhs.tail(q)=-residuals(e);Eigen::JacobiSVD<Mat> svd(A,Eigen::ComputeThinU|Eigen::ComputeThinV);svd.setThreshold(1e-12);V delta=svd.solve(rhs);++out.visits;if(!delta.allFinite()){out.reason="nonfinite_kkt";break;}double before=rhs.squaredNorm(),step=std::min(1.,20./std::max(20.,delta.head(n).norm()));bool accepted=false;
 for(int line=0;line<18&&out.evaluations<budget;++line){V z=out.z+step*delta.head(n),dual=out.dual+step*delta.tail(q);try{auto trial=p.evaluate(z);++out.evaluations;double score=residuals(trial).squaredNorm()+std::pow(stationarity(trial,dual),2);if(score<before*(1-1e-4*step)||score<1e-20){out.z=z;out.dual=dual;e=std::move(trial);accepted=true;break;}}catch(const std::exception&){++out.evaluations;}step*=.5;}if(!accepted){out.reason="kkt_line_search_stalled";break;}
 }finish(e,out);out.milliseconds=std::chrono::duration<double,std::milli>(std::chrono::steady_clock::now()-start).count();return out;
}
struct LiftedFull : hour_local::Problem {
 V z0;
 LiftedFull(const Model&m,const State&s,double h,double E):hour_local::Problem(m,s,h,E){work_reduction_fast::Problem p(*this);z0=p.expand(p.initial());}
 V initial()const{return z0;}
};
V newton_step(const Evaluation&e,const V&dual){int n=e.objective.g.size(),q=e.rows.size();Mat H=e.objective.H,J=jacobian(e);for(int i=0;i<q;++i)H+=dual[i]*e.rows[i].value.H;Mat K=Mat::Zero(n+q,n+q);K.topLeftCorner(n,n)=H;K.topRightCorner(n,q)=J.transpose();K.bottomLeftCorner(q,n)=J;V rhs(n+q);rhs.head(n)=-(e.objective.g+J.transpose()*dual);rhs.tail(q)=-residuals(e);return K.fullPivLu().solve(rhs);}
