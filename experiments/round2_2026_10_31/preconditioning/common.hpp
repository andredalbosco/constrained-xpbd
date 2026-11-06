#pragma once
#include <Eigen/Dense>
#include <algorithm>
#include <cmath>
#include <string>
#include <vector>
#include <stdexcept>

namespace compose {
using V=Eigen::VectorXd;
using Mat=Eigen::MatrixXd;
struct Model {std::string name; int count=1,dim=2; bool rigid=false; V mass,length,stiffness,gravity;};
struct State {V x,v;};
struct Config {
    std::string id,family,work="dg_average",treatment="global",kinematics="soft_endpoint",schedule="joint",curvature="drop";
    double beta=1,gamma=1; int iterations=8;
    bool tangency=false,keep_rx=true,relinearize=true,reverse=false;
};
struct StepResult {
    State state; bool valid=true; std::string reason;
    int iterations=0,visits=0,evaluations=0;
    double balance=0,rx=0,kinematic=0,impulse=0,work_defect=0;
};
inline int dof(const Model&m){return m.count*m.dim;}
inline V masses(const Model&m){V a(dof(m));for(int i=0;i<m.count;++i)a.segment(i*m.dim,m.dim).setConstant(m.mass[i]);return a;}
inline V edge(const Model&m,const V&x,int j){V e=x.segment(j*m.dim,m.dim);if(j)e-=x.segment((j-1)*m.dim,m.dim);return e;}
// The scalar environment is an affine spring, C=x-l, not an absolute-distance spring.
inline double edge_length(const Model&m,const V&e){return m.dim==1?e[0]:e.norm();}
inline V edge_direction(const Model&m,const V&e){if(m.dim==1)return V::Ones(1);double n=e.norm();if(n<1e-13)throw std::runtime_error("collapsed_edge");return e/n;}
inline bool has_gravity(const Model&m){return m.gravity.squaredNorm()>0;}
inline V potential_parts(const Model&m,const V&x){
    V u=V::Zero(m.count+(has_gravity(m)?m.count:0));
    for(int j=0;j<m.count;++j)if(!m.rigid){double c=edge_length(m,edge(m,x,j))-m.length[j];u[j]=.5*m.stiffness[j]*c*c;}
    if(has_gravity(m))for(int i=0;i<m.count;++i)u[m.count+i]=-m.mass[i]*m.gravity.dot(x.segment(i*m.dim,m.dim));
    return u;
}
inline double potential(const Model&m,const V&x){return potential_parts(m,x).sum();}
inline V force(const Model&m,const V&x){V f=V::Zero(dof(m));for(int i=0;i<m.count;++i)f.segment(i*m.dim,m.dim)=m.mass[i]*m.gravity;
    if(!m.rigid)for(int j=0;j<m.count;++j){V e=edge(m,x,j);V fj=-m.stiffness[j]*(edge_length(m,e)-m.length[j])*edge_direction(m,e);f.segment(j*m.dim,m.dim)+=fj;if(j)f.segment((j-1)*m.dim,m.dim)-=fj;}return f;}
inline double energy(const Model&m,const State&s){return .5*(masses(m).array()*s.v.array().square()).sum()+potential(m,s.x);}
inline double minimum_energy(const Model&m){double e=0,weight=0;for(int j=m.count-1;j>=0;--j){weight+=m.mass[j]*m.gravity.norm();e-=weight*m.length[j];if(!m.rigid)e-=.5*weight*weight/m.stiffness[j];}return e;}
inline bool finite(const State&s){return s.x.allFinite()&&s.v.allFinite();}
inline double angular(const Model&m,const State&s){if(m.dim!=2)return 0;double a=0;for(int i=0;i<m.count;++i)a+=m.mass[i]*(s.x[2*i]*s.v[2*i+1]-s.x[2*i+1]*s.v[2*i]);return a;}
inline double torque(const Model&m,const State&s){if(m.dim!=2)return 0;double a=0;for(int i=0;i<m.count;++i)a+=m.mass[i]*(s.x[2*i]*m.gravity[1]-s.x[2*i+1]*m.gravity[0]);return a;}
inline State reference_step(const Model&m,State s,double h,int substeps=32){
    const double dt=h/substeps;
    if(m.rigid){if(m.count!=1||m.dim!=2)throw std::runtime_error("unsupported_rigid_reference");
        double q=std::atan2(s.x[0],-s.x[1]),w=(s.x[0]*s.v[1]-s.x[1]*s.v[0])/(m.length[0]*m.length[0]);
        auto a=[&](double t){return (m.gravity[0]*std::cos(t)+m.gravity[1]*std::sin(t))/m.length[0];};
        for(int k=0;k<substeps;++k){double k1q=w,k1w=a(q),k2q=w+.5*dt*k1w,k2w=a(q+.5*dt*k1q),k3q=w+.5*dt*k2w,k3w=a(q+.5*dt*k2q),k4q=w+dt*k3w,k4w=a(q+dt*k3q);q+=dt*(k1q+2*k2q+2*k3q+k4q)/6;w+=dt*(k1w+2*k2w+2*k3w+k4w)/6;}
        s.x<<m.length[0]*std::sin(q),-m.length[0]*std::cos(q);s.v<<m.length[0]*w*std::cos(q),m.length[0]*w*std::sin(q);return s;
    }
    const V inv=masses(m).cwiseInverse();auto acceleration=[&](const V&x)->V{return inv.cwiseProduct(force(m,x));};
    for(int k=0;k<substeps;++k){V k1x=s.v,k1v=acceleration(s.x),k2x=s.v+.5*dt*k1v,k2v=acceleration(s.x+.5*dt*k1x),k3x=s.v+.5*dt*k2v,k3v=acceleration(s.x+.5*dt*k2x),k4x=s.v+dt*k3v,k4v=acceleration(s.x+dt*k3x);s.x+=dt*(k1x+2*k2x+2*k3x+k4x)/6;s.v+=dt*(k1v+2*k2v+2*k3v+k4v)/6;}
    return s;
}
struct Case {std::string name;Model model;State initial;};
inline Case make_case(const std::string&name,double stiffness_factor=1,double mass_ratio=1,double perturbation=1){
    Case c;c.name=name;Model&m=c.model;m.name=name;m.dim=name.find("scalar")==0?1:2;m.count=name.find("chain")==0?2:1;m.rigid=name.find("rod")==0;
    m.mass=V::Ones(m.count);m.mass[m.count-1]=mass_ratio;m.length=V::Ones(m.count);m.stiffness=V::Constant(m.count,100*stiffness_factor);m.gravity=V::Zero(m.dim);
    if(m.rigid||name.find("elastic")==0||m.count==2)m.gravity[m.dim-1]=-9.81;
    c.initial={V::Zero(dof(m)),V::Zero(dof(m))};auto&s=c.initial;
    if(m.dim==1){s.x[0]=1+.25*perturbation;if(name=="scalar_mixed")s.v[0]=.8*perturbation;return c;}
    if(name.find("radial")==0){s.x[0]=1+.1*perturbation;s.v[1]=.8*perturbation;if(name=="radial_release"){s.x[0]=1+.25*perturbation;s.v.setZero();}return c;}
    if(m.rigid){double a=.65*perturbation;s.x<<std::sin(a),-std::cos(a);if(name=="rod_kick")s.v<<1.2*perturbation*std::cos(a),1.2*perturbation*std::sin(a);return c;}
    double downstream=m.mass.sum();V prev=V::Zero(2),prevv=V::Zero(2);
    for(int j=0;j<m.count;++j){double radius=1+downstream*9.81/m.stiffness[j];double angle=(j?-.3:.45)*perturbation;double stretch=(j?.1:.05)*perturbation;if(name=="elastic_release"){angle=0;stretch=.2*perturbation;}if(name=="chain_large"){angle=(j?-.45:.65)*perturbation;stretch=(j?.2:.1)*perturbation;}
        V n(2),t(2);n<<std::sin(angle),-std::cos(angle);t<<std::cos(angle),std::sin(angle);V x=prev+(radius+stretch)*n;V v=prevv+perturbation*((j?.3:-.2)*n+(j?.8:.5)*t);if(name=="elastic_release"||name=="chain_large")v.setZero();s.x.segment(2*j,2)=x;s.v.segment(2*j,2)=v;prev=x;prevv=v;downstream-=m.mass[j];}
    return c;
}
inline std::vector<std::string> case_names(){return {"scalar_release","scalar_mixed","radial_orbit","radial_release","rod_release","rod_kick","elastic_mixed","elastic_release","chain_mixed","chain_large"};}
} // namespace compose
