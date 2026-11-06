#pragma once
#include "cross_problem.hpp"

namespace hour_local {
using namespace local_work;
using independent_detail::Kind;
using independent_detail::Row;
using independent_detail::Kinetic;
using independent_detail::Storage;
using independent_detail::GravityStorage;
using independent_detail::Position;
using independent_detail::Tangency;

// The old AD evaluator remains an independent derivative oracle. The numerical
// experiments use the same explicit-W objective assembled analytically below.
struct Problem : CrossProblem {
    double fraction=1;
    Problem(const Model&m,const State&s,double h,double E,double blend=1)
      : CrossProblem(m,s,h,config(m),E,"implicit_midpoint_dg",true),fraction(blend){
        aold=components(m,s,false).acceleration;
    }
    static Config config(const Model&m){auto c=configuration("implicit_midpoint_dg","soft_midpoint");c.tangency=m.rigid;return c;}
    Evaluation evaluate(const V&z)const {
        const auto&m=p.m;const int d=p.d,n=p.n,dim=m.dim,offset=2*d;const double E=p.E,xs=p.xs,vs=p.vs;
        V x=p.old.x+xs*z.head(d),v=vs*z.segment(d,d);
        Evaluation e;e.objective=Jet(0,n);e.x.resize(d);e.v.resize(d);e.work.resize(p.wcount);
        // Only values are needed in these diagnostic fields. Objective and row
        // derivatives below carry the complete exact derivatives.
        for(int i=0;i<d;++i){e.x[i].a=x[i];e.v[i].a=v[i];}
        auto addrow=[&](double value,Kind kind,int index)->Jet&{e.rows.push_back({Jet(value,n),kind,index});return e.rows.back().value;};
        for(int i=0;i<d;++i){double mass=m.mass[i/dim],r=(x[i]-p.old.x[i])/p.h-p.old.v[i]-.25*fraction*p.h*aold[i],rv=v[i]+p.old.v[i]-2*(x[i]-p.old.x[i])/p.h,c=mass/E;
            e.objective.a+=.5*c*(r*r+rv*rv/cfg.beta);e.objective.g[i]+=c*vs*(r-2*rv/cfg.beta);e.objective.g[d+i]+=c*vs*rv/cfg.beta;
            e.objective.H(i,i)+=c*vs*vs*(1+4/cfg.beta);e.objective.H(d+i,d+i)+=c*vs*vs/cfg.beta;e.objective.H(i,d+i)-=2*c*vs*vs/cfg.beta;e.objective.H(d+i,i)-=2*c*vs*vs/cfg.beta;
        }
        for(int j=0;j<p.wcount;++j){double r=z[offset+j]-frozen[j]/E;e.work[j].a=E*z[offset+j];e.objective.a+=.5*cfg.gamma*r*r;e.objective.g[offset+j]=cfg.gamma*r;e.objective.H(offset+j,offset+j)=cfg.gamma;}
        for(int i=0;i<m.count;++i){double K=.5*m.mass[i]*(v.segment(i*dim,dim).squaredNorm()-p.old.v.segment(i*dim,dim).squaredNorm())/E;Jet&r=addrow(K,Kinetic,i);
            for(int k=0;k<dim;++k){int vi=d+i*dim+k;r.g[vi]=m.mass[i]*v[i*dim+k]*vs/E;r.H(vi,vi)=m.mass[i]*vs*vs/E;}
            for(int j:{p.work_child[i],i+1<m.count?p.work_parent[i+1]:-1,p.work_gravity[i]})if(j>=0){r.a-=z[offset+j];r.g[offset+j]=-1;}
        }
        const double weight=1-.75*fraction;
        for(int j=0;j<m.count;++j){V delta=edge(m,x,j),normal=edge_direction(m,delta);double len=edge_length(m,delta),C=len-m.length[j],U=m.rigid?0:.5*m.stiffness[j]*C*C;V gradient=V::Zero(dim);Mat H=Mat::Zero(dim,dim);
            if(!m.rigid){gradient=m.stiffness[j]*C*normal;H=m.stiffness[j]*(normal*normal.transpose());if(dim>1)H+=m.stiffness[j]*C/len*(Mat::Identity(dim,dim)-normal*normal.transpose());}
            Jet&r=addrow((U-p.old_parts[j])/E,Storage,j);r.a+=z[offset+p.work_child[j]];r.g[offset+p.work_child[j]]=1;if(j){r.a+=z[offset+p.work_parent[j]];r.g[offset+p.work_parent[j]]=1;}
            e.objective.a+=weight*U/E;
            for(int a=0;a<(j?2:1);++a){int ia=(j-a)*dim;double sa=a?-1:1;r.g.segment(ia,dim)+=sa*xs/E*gradient;e.objective.g.segment(ia,dim)+=weight*sa*xs/E*gradient;
                for(int b=0;b<(j?2:1);++b){int ib=(j-b)*dim;double sb=b?-1:1;Mat local=sa*sb*xs*xs/E*H;r.H.block(ia,ib,dim,dim)+=local;e.objective.H.block(ia,ib,dim,dim)+=weight*local;}}
        }
        if(has_gravity(m))for(int i=0;i<m.count;++i){V gradient=-m.mass[i]*m.gravity;double U=gradient.dot(x.segment(i*dim,dim));e.objective.a+=weight*U/E;e.objective.g.segment(i*dim,dim)+=weight*xs/E*gradient;Jet&r=addrow((U-p.old_parts[m.count+i])/E+z[offset+p.work_gravity[i]],GravityStorage,i);r.g.segment(i*dim,dim)=xs/E*gradient;r.g[offset+p.work_gravity[i]]=1;}
        if(m.rigid){if(m.count!=1)throw std::runtime_error("only_single_rigid_pendulum_supported");V normal=x/x.norm();double len=x.norm();Mat P=Mat::Identity(dim,dim)-normal*normal.transpose();Jet&r=addrow((len-m.length[0])/xs,Position,0);r.g.head(dim)=normal;r.H.topLeftCorner(dim,dim)=xs/len*P;
            if(cfg.tangency){double nv=normal.dot(v);Jet&t=addrow(nv/vs,Tangency,0);t.g.head(dim)=xs/(vs*len)*P*v;t.g.segment(d,dim)=normal;Mat cross=xs/len*P;t.H.block(0,d,dim,dim)=cross;t.H.block(d,0,dim,dim)=cross;Mat B=(-nv*Mat::Identity(dim,dim)-v*normal.transpose()-normal*v.transpose()+3*nv*normal*normal.transpose());t.H.topLeftCorner(dim,dim)=xs*xs/(vs*len*len)*B;}
        }
        return e;
    }
    V initial(const std::string&kind="midpoint")const {
        if(kind=="backward_euler")return ExplicitProblem::initial();
        State s=stiffness_aware_state(p.m,p.old,p.h,true);V z=V::Zero(p.n);z.head(p.d)=(s.x-p.old.x)/p.xs;z.segment(p.d,p.d)=s.v/p.vs;z.tail(p.wcount)=frozen/p.E;return z;
    }
    State state(const V&z)const{return {p.old.x+p.xs*z.head(p.d),p.vs*z.segment(p.d,p.d)};}
    Evaluation oracle(const V&z)const {
        if(fraction!=1)throw std::runtime_error("oracle_only_final_objective");return CrossProblem::evaluate(z);
    }
};
inline Case chain_case(int N,double stiffness=1,double mass_ratio=1){
    Case c;c.name="chain"+std::to_string(N);auto&m=c.model;m.name=c.name;m.count=N;m.dim=2;m.mass=V::Ones(N);for(int i=0;i<N;++i)if(i%2)m.mass[i]=mass_ratio;m.length=V::Constant(N,1./N);m.stiffness=V::Constant(N,100.*N*stiffness);m.gravity=V::Zero(2);m.gravity[1]=-9.81;c.initial={V::Zero(2*N),V::Zero(2*N)};
    V x=V::Zero(2),v=V::Zero(2);double downstream=m.mass.sum();
    for(int i=0;i<N;++i){double angle=.35*std::sin(.7*(i+1)),r=m.length[i]+downstream*9.81/m.stiffness[i]+.03/N;V normal(2),tangent(2);normal<<std::sin(angle),-std::cos(angle);tangent<<std::cos(angle),std::sin(angle);x+=r*normal;v+=(.12/N)*tangent;c.initial.x.segment(2*i,2)=x;c.initial.v.segment(2*i,2)=v;downstream-=m.mass[i];}return c;
}
} // namespace hour_local
