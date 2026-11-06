#pragma once
#include "independent.hpp"
#include <string>

namespace local_work {
using namespace compose;
struct Incidence {int body,potential;V f,fdot;};
struct Components {std::vector<Incidence> entries;V acceleration;};
inline Components components(const Model&m,const State&s,bool derivative=true){
    Components out;V total=V::Zero(dof(m));
    for(int i=0;i<m.count;++i)total.segment(i*m.dim,m.dim)=m.mass[i]*m.gravity;
    std::vector<V> normals(m.count),edges(m.count),relative(m.count);std::vector<double> q(m.count),r(m.count),speed(m.count);
    for(int j=0;j<m.count;++j){edges[j]=edge(m,s.x,j);relative[j]=edge(m,s.v,j);r[j]=edge_length(m,edges[j]);normals[j]=edge_direction(m,edges[j]);speed[j]=normals[j].dot(relative[j]);
        if(m.rigid){if(m.count!=1)throw std::runtime_error("rigid_prediction_requires_single_pendulum");q[j]=-m.mass[j]*((relative[j].squaredNorm()-speed[j]*speed[j])/r[j]+normals[j].dot(m.gravity));}
        else q[j]=-m.stiffness[j]*(r[j]-m.length[j]);
        V f=q[j]*normals[j];total.segment(j*m.dim,m.dim)+=f;if(j)total.segment((j-1)*m.dim,m.dim)-=f;
    }
    out.acceleration=masses(m).cwiseInverse().cwiseProduct(total);
    for(int j=0;j<m.count;++j){V ndot=m.dim==1?V::Zero(1):V((relative[j]-normals[j]*speed[j])/r[j]);double qdot=0;
        if(derivative){if(!m.rigid)qdot=-m.stiffness[j]*speed[j];else{V ar=edge(m,out.acceleration,j);double sdot=ndot.dot(relative[j])+normals[j].dot(ar),tangent2=relative[j].squaredNorm()-speed[j]*speed[j];qdot=-m.mass[j]*((2*relative[j].dot(ar)-2*speed[j]*sdot)/r[j]-tangent2*speed[j]/(r[j]*r[j])+ndot.dot(m.gravity));}}
        V f=q[j]*normals[j],fd=derivative?V(qdot*normals[j]+q[j]*ndot):V::Zero(m.dim);out.entries.push_back({j,j,f,fd});if(j)out.entries.push_back({j-1,j,-f,-fd});
    }
    if(has_gravity(m))for(int i=0;i<m.count;++i)out.entries.push_back({i,m.count+i,m.mass[i]*m.gravity,V::Zero(m.dim)});
    return out;
}
inline V powers(const Model&m,const State&s){auto c=components(m,s,false);V w(c.entries.size());for(int j=0;j<w.size();++j)w[j]=c.entries[j].f.dot(s.v.segment(c.entries[j].body*m.dim,m.dim));return w;}
inline State predicted_state(const Model&m,const State&old,double h,bool project){
    auto c=components(m,old);State s{old.x+h*old.v+.5*h*h*c.acceleration,old.v+h*c.acceleration};
    if(project&&m.rigid){double L=m.length[0],angle=std::atan2(old.x[0],-old.x[1]),omega=(old.x[0]*old.v[1]-old.x[1]*old.v[0])/(L*L),alpha=(m.gravity[0]*std::cos(angle)+m.gravity[1]*std::sin(angle))/L;
        double a=angle+h*omega+.5*h*h*alpha,w=omega+h*alpha;s.x<<L*std::sin(a),-L*std::cos(a);s.v<<L*w*std::cos(a),L*w*std::sin(a);
    }return s;
}
inline State stiffness_aware_state(const Model&m,const State&old,double h,bool midpoint){
    if(m.rigid)return predicted_state(m,old,h,true);
    V mass=masses(m);Mat B=Mat::Zero(dof(m),dof(m));
    for(int j=0;j<m.count;++j){V b=V::Zero(dof(m)),n=edge_direction(m,edge(m,old.x,j));b.segment(j*m.dim,m.dim)=n;if(j)b.segment((j-1)*m.dim,m.dim)=-n;B+=m.stiffness[j]*b*b.transpose();}
    const double curvature=midpoint?.25:1,acceleration=midpoint?.5:1;
    Mat A=mass.asDiagonal();A+=curvature*h*h*B;
    V rhs=h*mass.cwiseProduct(old.v)+acceleration*h*h*force(m,old.x),dx=A.ldlt().solve(rhs);
    if(!dx.allFinite())throw std::runtime_error("start_predictor_linear_solve_failed");
    return {old.x+dx,midpoint?V(2*dx/h-old.v):V(dx/h)};
}
inline V work_prediction(const Model&m,const State&old,double h,const std::string&law){
    auto c=components(m,old);V w=V::Zero(c.entries.size());
    if(law=="frozen_power"||law=="acceleration_chord"||law=="taylor2"){
        for(int j=0;j<w.size();++j){const auto&a=c.entries[j];V v=old.v.segment(a.body*m.dim,m.dim),acc=c.acceleration.segment(a.body*m.dim,m.dim);w[j]=h*a.f.dot(v);if(law!="frozen_power")w[j]+=.5*h*h*a.f.dot(acc);if(law=="taylor2")w[j]+=.5*h*h*a.fdot.dot(v);}return w;
    }
    if(law=="predicted_midpoint"){State mid=predicted_state(m,old,h/2,false);auto f=components(m,mid,false);for(int j=0;j<w.size();++j)w[j]=h*f.entries[j].f.dot(mid.v.segment(f.entries[j].body*m.dim,m.dim));return w;}
    if(law=="predicted_dg"||law=="projected_dg"||law=="implicit_dg"||law=="implicit_midpoint_dg"){
        State end=(law=="implicit_dg"||law=="implicit_midpoint_dg")?stiffness_aware_state(m,old,h,law=="implicit_midpoint_dg"):predicted_state(m,old,h,law=="projected_dg");auto f=components(m,end,false);int incidence=0;
        for(int j=0;j<m.count;++j){V eo=edge(m,old.x,j),en=edge(m,end.x,j),no=edge_direction(m,eo),nn=edge_direction(m,en),dg=m.dim==1?V::Ones(1):V((en+eo)/(en.norm()+eo.norm()));double qo=c.entries[incidence].f.dot(no),qn=f.entries[incidence].f.dot(nn);V fw=.5*(qo+qn)*dg;
            w[incidence++]=fw.dot(end.x.segment(j*m.dim,m.dim)-old.x.segment(j*m.dim,m.dim));if(j)w[incidence++]=-fw.dot(end.x.segment((j-1)*m.dim,m.dim)-old.x.segment((j-1)*m.dim,m.dim));}
        if(has_gravity(m))for(int i=0;i<m.count;++i)w[incidence++]=m.mass[i]*m.gravity.dot(end.x.segment(i*m.dim,m.dim)-old.x.segment(i*m.dim,m.dim));return w;
    }
    throw std::runtime_error("unknown_frozen_prediction_"+law);
}
struct IntegratedWork {State end;V work;};
inline IntegratedWork integrate_work(const Model&m,const State&old,double h,int subdivisions){
    IntegratedWork r{old,V::Zero(components(m,old,false).entries.size())};double dt=h/subdivisions;
    for(int i=0;i<subdivisions;++i){State mid=reference_step(m,r.end,dt/2,1),end=reference_step(m,mid,dt/2,1);r.work+=dt/6*(powers(m,r.end)+4*powers(m,mid)+powers(m,end));r.end=std::move(end);}return r;
}
inline const std::vector<std::string>&prediction_names(){static std::vector<std::string> names={"frozen_power","acceleration_chord","taylor2","predicted_midpoint","predicted_dg","projected_dg","implicit_dg","implicit_midpoint_dg"};return names;}

// Convert the old estimate-correction coordinates to true explicit W/E
// coordinates. Constraint balances then contain W only, regardless of target.
struct ExplicitProblem {
    Config cfg;independent_detail::Problem p;std::string target;V frozen;bool is_frozen;
    ExplicitProblem(const Model&m,const State&s,double h,const Config&config,double scale,std::string law):cfg(config),p(m,s,h,cfg,scale),target(std::move(law)),is_frozen(target.rfind("live_",0)!=0){
        if(!p.penalized||p.hard)throw std::runtime_error("explicit_work_requires_penalized_soft_kinematics");
        if(is_frozen)frozen=work_prediction(m,s,h,target);
    }
    independent_detail::Evaluation evaluate(const V&z)const {
        using namespace independent_detail;auto e=p.evaluate(z);int offset=2*p.d;JV estimate(p.wcount),actual(p.wcount),chosen(p.wcount);
        for(int j=0;j<p.wcount;++j){Jet a=Jet::variable(z[offset+j],p.n,offset+j);actual[j]=p.E*a;estimate[j]=e.work[j]-actual[j];chosen[j]=is_frozen?Jet(frozen[j],p.n):estimate[j];e.objective=e.objective-.5*cfg.gamma*a*a+.5*cfg.gamma*(a-chosen[j]/p.E)*(a-chosen[j]/p.E);e.work[j]=actual[j];}
        for(auto&row:e.rows){int i=row.index;if(row.kind==Kinetic){row.value=row.value+estimate[p.work_child[i]]/p.E;if(i+1<p.m.count)row.value=row.value+estimate[p.work_parent[i+1]]/p.E;if(has_gravity(p.m))row.value=row.value+estimate[p.work_gravity[i]]/p.E;}
            else if(row.kind==Storage){row.value=row.value-estimate[p.work_child[i]]/p.E;if(i)row.value=row.value-estimate[p.work_parent[i]]/p.E;}
            else if(row.kind==GravityStorage)row.value=row.value-estimate[p.work_gravity[i]]/p.E;
        }return e;
    }
    V initial()const {V z=p.initial();auto e=p.evaluate(z);for(int j=0;j<p.wcount;++j)z[2*p.d+j]=(is_frozen?frozen[j]:e.work[j].a)/p.E;return z;}
    State state(const V&z)const {auto e=evaluate(z);State s{V(p.d),V(p.d)};for(int i=0;i<p.d;++i){s.x[i]=e.x[i].a;s.v[i]=e.v[i].a;}return s;}
};
inline Config configuration(std::string target,std::string kin="soft_midpoint"){Config c;c.work=target=="live_start"?"start":target=="live_dg"?"dg_average":"frozen_power";c.treatment="penalized";c.kinematics=kin;c.curvature="exact";c.beta=1;c.gamma=1;c.tangency=false;return c;}
}
