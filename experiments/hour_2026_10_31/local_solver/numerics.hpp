#pragma once
#include "incremental.hpp"
#include <chrono>
namespace hour_local {
struct Options {std::string id="kkt_midpoint",solver="kkt",initial="midpoint",order="forward";int budget=128;double rho=10;bool warm=false;std::string partition="split";bool stable=false;};
struct Result : local_work::Solve {int fallback=0;};
inline V kkt_residual(const Evaluation&e,const V&dual){V r(e.objective.g.size()+dual.size());r.head(e.objective.g.size())=e.objective.g+jacobian(e).transpose()*dual;r.tail(dual.size())=residuals(e);return r;}
inline Mat kkt_matrix(const Evaluation&e,const V&dual){int n=e.objective.g.size(),q=dual.size();Mat J=jacobian(e),H=e.objective.H;for(int i=0;i<q;++i)H+=dual[i]*e.rows[i].value.H;Mat K=Mat::Zero(n+q,n+q);K.topLeftCorner(n,n)=.5*(H+H.transpose()).eval();K.topRightCorner(n,q)=J.transpose();K.bottomLeftCorner(q,n)=J;return K;}
inline Result joint(const Problem&p,const Options&o,const V*seed=nullptr,const V*dual_seed=nullptr){
    auto start=std::chrono::steady_clock::now();Result s;s.z=seed?*seed:p.initial(o.initial);auto e=p.evaluate(s.z);s.evaluations=1;s.dual=dual_seed&&dual_seed->size()==int(e.rows.size())?*dual_seed:V::Zero(e.rows.size());int n=p.p.n,q=s.dual.size();
    while(s.evaluations<o.budget){finish(e,s);if(s.converged)break;V r=kkt_residual(e,s.dual);Mat K=kkt_matrix(e,s.dual);V delta;double before=r.squaredNorm();bool accepted=false;
        if(o.solver=="merit"){
            Mat H=K.topLeftCorner(n,n),J=jacobian(e);Eigen::SelfAdjointEigenSolver<Mat>eig(H);double shift=std::max(0.,1e-4-eig.eigenvalues().minCoeff());H.diagonal().array()+=shift;if(shift)++s.regularizations;K.topLeftCorner(n,n)=H;V rhs=-r;rhs.head(n)=-e.objective.g;delta=K.fullPivLu().solve(rhs);V lambda=delta.tail(q);delta.tail(q)=lambda-s.dual;double penalty=std::max(10.,1.5*lambda.lpNorm<Eigen::Infinity>()),merit=e.objective.a+penalty*residuals(e).lpNorm<1>();double slope=e.objective.g.dot(delta.head(n))-penalty*residuals(e).lpNorm<1>();double step=std::min(1.,10./std::max(10.,delta.head(n).norm()));
            for(int line=0;line<20&&s.evaluations<o.budget;++line){V z=s.z+step*delta.head(n);auto trial=p.evaluate(z);++s.evaluations;if(trial.objective.a+penalty*residuals(trial).lpNorm<1>()<=merit+1e-4*step*std::min(-1e-16,slope)){s.z=z;s.dual+=step*delta.tail(q);e=std::move(trial);accepted=true;break;}step*=.5;}
        }else{
            delta=K.fullPivLu().solve(-r);double step=std::min(1.,20./std::max(20.,delta.head(n).norm()));
            for(int line=0;line<18&&s.evaluations<o.budget;++line){V z=s.z+step*delta.head(n),dual=s.dual+step*delta.tail(q);auto trial=p.evaluate(z);++s.evaluations;double score=kkt_residual(trial,dual).squaredNorm();if(score<before*(1-1e-4*step)||score<1e-20){s.z=z;s.dual=dual;e=std::move(trial);accepted=true;break;}step*=.5;}
            if(!accepted&&o.solver=="lm"){
                // Safeguard the residual equations, not the constrained objective.
                // A residual stationary point is still reported as unconverged.
                Mat normal=K.transpose()*K;double scale=std::max(1.,normal.diagonal().cwiseAbs().maxCoeff());V right=-K.transpose()*r;
                for(int attempt=0;attempt<9&&s.evaluations<o.budget;++attempt){Mat A=normal;A.diagonal().array()+=scale*std::pow(10.,attempt-8.);delta=A.ldlt().solve(right);V z=s.z+delta.head(n),dual=s.dual+delta.tail(q);auto trial=p.evaluate(z);++s.evaluations;++s.fallback;if(kkt_residual(trial,dual).squaredNorm()<before*(1-1e-8)){s.z=z;s.dual=dual;e=std::move(trial);accepted=true;break;}}
            }
        }
        ++s.visits;++s.sweeps;if(!accepted){s.reason="line_search_stalled";break;}if(!s.z.allFinite()||!s.dual.allFinite()){s.reason="nonfinite";break;}
    }
    finish(e,s);s.milliseconds=std::chrono::duration<double,std::milli>(std::chrono::steady_clock::now()-start).count();return s;
}
inline double aug_value(const Evaluation&e,const V&y,double rho){double L=e.objective.a;for(int i=0;i<y.size();++i){double g=e.rows[i].value.a;L+=y[i]*g+.5*rho*g*g;}return L;}
inline std::vector<std::vector<int>> blocks_for(const Problem&p,const Options&o){
    if(o.partition=="global"){std::vector<int>b(p.p.n);std::iota(b.begin(),b.end(),0);return {b};}
    if(o.partition=="split")return primal_blocks(p,o.order);
    std::vector<std::vector<int>> blocks;int dim=p.p.m.dim,d=p.p.d,offset=2*d;
    for(int j=0;j<p.p.m.count;++j){std::vector<int>b;auto vertex=[&](int i){for(int k=0;k<dim;++k){b.push_back(i*dim+k);b.push_back(d+i*dim+k);}b.push_back(offset+p.p.work_child[i]);if(i+1<p.p.m.count)b.push_back(offset+p.p.work_parent[i+1]);if(p.p.work_gravity[i]>=0)b.push_back(offset+p.p.work_gravity[i]);};vertex(j);if(o.partition=="edge_full"&&j)vertex(j-1);std::sort(b.begin(),b.end());b.erase(std::unique(b.begin(),b.end()),b.end());blocks.push_back(b);}
    if(o.order=="reverse")std::reverse(blocks.begin(),blocks.end());if(o.order=="symmetric"){auto b=blocks;std::reverse(b.begin(),b.end());blocks.insert(blocks.end(),b.begin(),b.end());}return blocks;
}
inline Result local(const Problem&p,const Options&o,const V*seed=nullptr,const V*dual_seed=nullptr){
    auto start=std::chrono::steady_clock::now();Result s;s.z=seed?*seed:p.initial(o.initial);auto e=p.evaluate(s.z);s.evaluations=1;s.dual=dual_seed&&dual_seed->size()==int(e.rows.size())?*dual_seed:V::Zero(e.rows.size());s.rho=o.rho;auto blocks=blocks_for(p,o);
    while(s.evaluations<o.budget){finish(e,s);if(s.converged)break;
        for(int inner=0;inner<3&&s.evaluations<o.budget;++inner)for(const auto&block:blocks){if(s.evaluations>=o.budget)break;int b=block.size();V gradient(b);Mat H(b,b);for(int i=0;i<b;++i){gradient[i]=e.objective.g[block[i]];for(int j=0;j<b;++j)H(i,j)=e.objective.H(block[i],block[j]);}
            for(int r=0;r<int(e.rows.size());++r){const auto&row=e.rows[r].value;double factor=s.dual[r]+s.rho*row.a;for(int i=0;i<b;++i){gradient[i]+=factor*row.g[block[i]];for(int j=0;j<b;++j)H(i,j)+=factor*row.H(block[i],block[j])+s.rho*row.g[block[i]]*row.g[block[j]];}}
            ++s.visits;if(gradient.norm()<1e-12)continue;H=.5*(H+H.transpose()).eval();Eigen::SelfAdjointEigenSolver<Mat>eig(H);double floor=1e-6*std::max(1.,H.diagonal().cwiseAbs().maxCoeff()),shift=std::max(0.,floor-eig.eigenvalues().minCoeff());if(shift){H.diagonal().array()+=shift;++s.regularizations;}V direction=H.ldlt().solve(-gradient);double step=std::min(1.,5./std::max(5.,direction.norm())),before=aug_value(e,s.dual,s.rho),slope=gradient.dot(direction);
            for(int line=0;line<16&&s.evaluations<o.budget;++line){V z=s.z;for(int i=0;i<b;++i)z[block[i]]+=step*direction[i];auto trial=p.evaluate(z);++s.evaluations;bool accepted=o.stable?augmented_difference(p,s.z,z,s.dual,s.rho,e)<=1e-4*step*slope:aug_value(trial,s.dual,s.rho)<=before+1e-4*step*slope;if(accepted){s.z=z;e=std::move(trial);break;}step*=.5;}
        }
        s.dual+=s.rho*residuals(e);++s.sweeps;if(s.visits>100000){s.reason="stationary_primal_blocks";break;}
    }
    finish(e,s);s.milliseconds=std::chrono::duration<double,std::milli>(std::chrono::steady_clock::now()-start).count();return s;
}
inline Result solve(const Problem&p,const Options&o,const V*dual_seed=nullptr){
    if(o.solver=="AL")return local(p,o,nullptr,dual_seed);
    if(o.solver!="continuation")return joint(p,o,nullptr,dual_seed);
    Result out;V seed=p.initial(o.initial),dual;int remaining=o.budget;
    for(double fraction:{0.,.5,1.}){Problem stage(p.p.m,p.p.old,p.p.h,p.p.E,fraction);Options q=o;q.solver="kkt";q.budget=std::max(1,remaining-(fraction==1?0:8));auto current=joint(stage,q,&seed,dual.size()?&dual:nullptr);remaining-=current.evaluations;out.evaluations+=current.evaluations;out.visits+=current.visits;out.sweeps+=current.sweeps;out.milliseconds+=current.milliseconds;seed=current.z;dual=current.dual;if(fraction==1){int evals=out.evaluations,visits=out.visits,sweeps=out.sweeps;double ms=out.milliseconds;out=current;out.evaluations=evals;out.visits=visits;out.sweeps=sweeps;out.milliseconds=ms;}}
    return out;
}
inline std::vector<Options> small_options(){return {
    {"kkt_BE","kkt","backward_euler","forward",128,10,false},
    {"kkt_midpoint","kkt","midpoint","forward",128,10,false},
    {"kkt_midpoint_warm","kkt","midpoint","forward",128,10,true},
    {"continuation_BE","continuation","backward_euler","forward",128,10,false},
    {"residual_LM_BE","lm","backward_euler","forward",128,10,false},
    {"merit_BE","merit","backward_euler","forward",128,10,false},
    {"AL_BE_512","AL","backward_euler","forward",512,10,false},
    {"AL_midpoint_512","AL","midpoint","forward",512,10,false},
    {"AL_midpoint_warm_512","AL","midpoint","forward",512,10,true}
};}
} // namespace hour_local
