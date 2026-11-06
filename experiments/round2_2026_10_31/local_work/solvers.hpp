#pragma once
#include "predictions.hpp"
#include <chrono>
#include <numeric>

namespace local_work {
using independent_detail::Evaluation;using independent_detail::Jet;
using independent_detail::residuals;using independent_detail::jacobian;using independent_detail::stationarity;
struct Solve {V z,dual;double objective=0,primal=0,dual_residual=0,rho=0,milliseconds=0;int evaluations=0,visits=0,sweeps=0,regularizations=0;bool converged=false;std::string reason="budget";};
inline void finish(const Evaluation&e,Solve&s){s.objective=e.objective.a;s.primal=residuals(e).lpNorm<Eigen::Infinity>();s.dual_residual=stationarity(e,s.dual);s.converged=s.primal<1e-8&&s.dual_residual<1e-7;if(s.converged)s.reason="converged";}
template<class ProblemType> inline double reduced_curvature(const ProblemType&p,const Solve&s){
    auto e=p.evaluate(s.z);Mat J=jacobian(e);Eigen::JacobiSVD<Mat> svd(J,Eigen::ComputeFullV);svd.setThreshold(1e-10);int nullity=p.p.n-svd.rank();if(!nullity)return 0;Mat Z=svd.matrixV().rightCols(nullity),H=e.objective.H;for(int i=0;i<s.dual.size();++i)H+=s.dual[i]*e.rows[i].value.H;Mat reduced=Z.transpose()*H*Z;reduced=.5*(reduced+reduced.transpose()).eval();return Eigen::SelfAdjointEigenSolver<Mat>(reduced).eigenvalues()[0];
}
template<class ProblemType> inline Solve joint_solve(const ProblemType&p,int budget,const V*initial=nullptr){
    auto start=std::chrono::steady_clock::now();Solve out;out.z=initial?*initial:p.initial();auto e=p.evaluate(out.z);++out.evaluations;out.dual=V::Zero(e.rows.size());
    while(out.evaluations<budget){finish(e,out);if(out.converged)break;Mat J=jacobian(e),H=e.objective.H;for(int i=0;i<out.dual.size();++i)H+=out.dual[i]*e.rows[i].value.H;int n=p.p.n,q=out.dual.size();Mat A=Mat::Zero(n+q,n+q);A.topLeftCorner(n,n)=.5*(H+H.transpose()).eval();A.topRightCorner(n,q)=J.transpose();A.bottomLeftCorner(q,n)=J;V rhs(n+q);rhs.head(n)=-(e.objective.g+J.transpose()*out.dual);rhs.tail(q)=-residuals(e);
        Eigen::JacobiSVD<Mat> svd(A,Eigen::ComputeThinU|Eigen::ComputeThinV);svd.setThreshold(1e-12);V delta=svd.solve(rhs);++out.visits;++out.sweeps;if(!delta.allFinite()){out.reason="nonfinite_kkt";break;}
        double before=rhs.squaredNorm(),step=std::min(1.,20./std::max(20.,delta.head(n).norm()));bool accepted=false;
        for(int line=0;line<18&&out.evaluations<budget;++line){V z=out.z+step*delta.head(n),dual=out.dual+step*delta.tail(q);try{auto trial=p.evaluate(z);++out.evaluations;double score=residuals(trial).squaredNorm()+std::pow(stationarity(trial,dual),2);if(score<before*(1-1e-4*step)||score<1e-20){out.z=z;out.dual=dual;e=std::move(trial);accepted=true;break;}}catch(const std::exception&){++out.evaluations;}step*=.5;}
        if(!accepted){out.reason="kkt_line_search_stalled";break;}
    }finish(e,out);out.milliseconds=std::chrono::duration<double,std::milli>(std::chrono::steady_clock::now()-start).count();return out;
}
template<class ProblemType> inline std::vector<std::vector<int>> primal_blocks(const ProblemType&p,std::string order){
    std::vector<std::vector<int>> blocks;int d=p.p.d,dim=p.p.m.dim,offset=2*d;
    for(int j=0;j<p.p.m.count;++j){std::vector<int>b;for(int i:{j,j-1})if(i>=0)for(int k=0;k<dim;++k)b.push_back(i*dim+k);b.push_back(offset+p.p.work_child[j]);if(j)b.push_back(offset+p.p.work_parent[j]);blocks.push_back(b);}
    for(int i=0;i<p.p.m.count;++i){std::vector<int>b;for(int k=0;k<dim;++k)b.push_back(d+i*dim+k);b.push_back(offset+p.p.work_child[i]);if(i+1<p.p.m.count)b.push_back(offset+p.p.work_parent[i+1]);if(has_gravity(p.p.m))b.push_back(offset+p.p.work_gravity[i]);blocks.push_back(b);}
    if(order=="reverse")std::reverse(blocks.begin(),blocks.end());
    if(order=="symmetric"){auto extra=blocks;std::reverse(extra.begin(),extra.end());blocks.insert(blocks.end(),extra.begin(),extra.end());}return blocks;
}
inline Jet augmented(const Evaluation&e,const V&dual,double rho){Jet a=e.objective;for(int i=0;i<dual.size();++i)a=a+dual[i]*e.rows[i].value+.5*rho*e.rows[i].value*e.rows[i].value;return a;}
template<class ProblemType> inline Solve local_al(const ProblemType&p,int budget,double initial_rho,bool adaptive,const std::string&order){
    auto start=std::chrono::steady_clock::now();Solve out;out.z=p.initial();auto e=p.evaluate(out.z);++out.evaluations;out.dual=V::Zero(e.rows.size());out.rho=initial_rho;auto blocks=primal_blocks(p,order);double previous_primal=residuals(e).norm();
    while(out.evaluations<budget){finish(e,out);if(out.converged)break;
        // Three primal sweeps at fixed shared multipliers before one dual update.
        for(int inner=0;inner<3&&out.evaluations<budget;++inner)for(const auto&block:blocks){if(out.evaluations>=budget)break;Jet L=augmented(e,out.dual,out.rho);int q=block.size();V g(q);Mat H(q,q);for(int i=0;i<q;++i){g[i]=L.g[block[i]];for(int j=0;j<q;++j)H(i,j)=L.H(block[i],block[j]);}++out.visits;
            H=.5*(H+H.transpose()).eval();Eigen::SelfAdjointEigenSolver<Mat> eigen(H);double floor=1e-6*std::max(1.,H.diagonal().cwiseAbs().maxCoeff()),shift=std::max(0.,floor-eigen.eigenvalues()[0]);if(shift>0){H.diagonal().array()+=shift;++out.regularizations;}V direction=H.ldlt().solve(-g);if(!direction.allFinite()){out.reason="nonfinite_block";continue;}double step=std::min(1.,5./std::max(5.,direction.norm())),slope=g.dot(direction);
            if(g.norm()<1e-12){continue;}
            for(int line=0;line<16&&out.evaluations<budget;++line){V z=out.z;for(int i=0;i<q;++i)z[block[i]]+=step*direction[i];try{auto trial=p.evaluate(z);++out.evaluations;double value=augmented(trial,out.dual,out.rho).a;if(std::isfinite(value)&&value<=L.a+1e-4*step*slope){out.z=z;e=std::move(trial);break;}}catch(const std::exception&){++out.evaluations;}step*=.5;}
        }
        ++out.sweeps;V g=residuals(e);out.dual+=out.rho*g;
        // A progress-based penalty schedule, not an ADMM dual-residual rule.
        if(adaptive&&out.sweeps%5==0){double primal=g.norm();if(primal>.7*previous_primal)out.rho=std::min(1e5,2*out.rho);previous_primal=primal;}
        if(!out.z.allFinite()||!out.dual.allFinite()){out.reason="nonfinite";break;}
        if(out.visits>100000){out.reason="stationary_primal_blocks";break;}
    }
    finish(e,out);out.milliseconds=std::chrono::duration<double,std::milli>(std::chrono::steady_clock::now()-start).count();return out;
}
}
