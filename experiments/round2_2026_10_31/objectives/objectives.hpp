#pragma once
#include "independent.hpp"
#include "mechanics.hpp"
namespace objectives2 {
using namespace compose;
using namespace compose::independent_detail;
struct Variant {std::string id,objective="original",kin="soft_midpoint",schedule="constant",target="verlet";double beta=1;bool energy=true,angular=false;};
inline std::vector<Variant> variants(){
 return {
 {"original_endpoint_b1","original","soft_endpoint","constant","verlet",1,true,false},
 {"original_midpoint_b1","original","soft_midpoint","constant","verlet",1,true,false},
 {"original_midpoint_b01","original","soft_midpoint","constant","verlet",.1,true,false},
 {"original_midpoint_b001","original","soft_midpoint","constant","verlet",.01,true,false},
 {"original_midpoint_h2","original","soft_midpoint","h2","verlet",1,true,false},
 {"original_midpoint_h2_c01","original","soft_midpoint","h2","verlet",.1,true,false},
 {"original_midpoint_h4","original","soft_midpoint","h4","verlet",1,true,false},
 {"quarter_b1","quarter","soft_midpoint","constant","verlet",1,true,false},
 {"quarter_h2","quarter","soft_midpoint","h2","verlet",1,true,false},
 {"half_half_bad","half_half","soft_midpoint","constant","verlet",1,true,false},
 {"midpotential_b1","midpotential","soft_midpoint","constant","verlet",1,true,false},
 {"midpotential_h2","midpotential","soft_midpoint","h2","verlet",1,true,false},
 {"projection_xpbd_E","projection","soft_midpoint","constant","xpbd",1,true,false},
 {"projection_verlet_E","projection","soft_midpoint","constant","verlet",1,true,false},
 {"projection_verlet_EL","projection","soft_midpoint","constant","verlet",1,true,true},
 {"projection_midpoint_E","projection","soft_midpoint","constant","midpoint",1,true,false},
 {"projection_midpoint_EL","projection","soft_midpoint","constant","midpoint",1,true,true},
 {"ls_midpoint","ls_midpoint","soft_midpoint","constant","verlet",1,false,false},
 {"ls_midpoint_E","ls_midpoint","soft_midpoint","constant","verlet",1,true,false},
 {"ls_dg","ls_dg","soft_midpoint","constant","verlet",1,false,false},
 {"ls_dg_E","ls_dg","soft_midpoint","constant","verlet",1,true,false},
 {"ls_dg_EL","ls_dg","soft_midpoint","constant","verlet",1,true,true}
 };
}
inline double omega_bound(const Model&m){if(m.rigid)return std::sqrt(std::max(m.gravity.norm()/m.length.minCoeff(),1.));return std::sqrt(4*m.stiffness.maxCoeff()/m.mass.minCoeff());}
inline double beta_for(const Variant&q,const Model&m,double h){double s=h*omega_bound(m);return std::max(1e-12,q.beta*(q.schedule=="h2"?s*s:(q.schedule=="h4"?s*s*s*s:1.)));}
inline State target_state(const Variant&q,const Model&m,const State&old,double h,double E){
 Config c;c.iterations=32;c.curvature="exact";c.schedule="joint";
 if(q.target=="xpbd"){c.family="xpbd";return mechanics_step(m,old,h,c,E).state;}
 if(q.target=="midpoint"){c.family="force_momentum";c.work="dg_average";c.kinematics="hard_midpoint";c.tangency=false;return mechanics_step(m,old,h,c,E).state;}
 if(m.rigid){c.family="rattle";c.tangency=true;return mechanics_step(m,old,h,c,E).state;}
 V invm=masses(m).cwiseInverse();V a=invm.cwiseProduct(force(m,old.x));State s;s.x=old.x+h*old.v+.5*h*h*a;s.v=old.v+.5*h*(a+invm.cwiseProduct(force(m,s.x)));return s;
}
struct Objective {
 const Model&m;const State&old;double h,E;Variant variant;Config cfg;Problem p;State target;V aold;bool reactions;
 Objective(const Model&mm,const State&s,double dt,const Variant&q,double scale):m(mm),old(s),h(dt),E(scale),variant(q),cfg(configuration(q,mm,dt)),p(mm,s,dt,cfg,scale),target(target_state(q,mm,s,dt,scale)){
  reactions=m.rigid&&q.objective.find("ls_")==0;if(reactions)p.n+=m.count;aold=masses(m).cwiseInverse().cwiseProduct(force(m,old.x));
  if(m.rigid){JV q0=p.coefficients(p.constants(old.x),p.constants(old.v));for(int j=0;j<m.count;++j){V f=q0[j].a*edge_direction(m,edge(m,old.x,j));aold.segment(j*m.dim,m.dim)+=f/m.mass[j];if(j)aold.segment((j-1)*m.dim,m.dim)-=f/m.mass[j-1];}}
 }
 static Config configuration(const Variant&q,const Model&m,double h){Config c;c.family="independent";c.kinematics=q.kin;c.treatment="global";c.work="dg_average";c.beta=beta_for(q,m,h);c.curvature="exact";c.iterations=32;c.tangency=m.rigid;return c;}
 Jet potential_jet(const JV&x)const {Jet u(0,p.n);for(int j=0;j<m.count;++j){if(!m.rigid){Jet c=length(m,jedge(m,x,j))-m.length[j];u=u+.5*m.stiffness[j]*c*c;}for(int k=0;k<m.dim;++k)u=u-m.mass[j]*m.gravity[k]*x[j*m.dim+k];}return u;}
 JV mechanical_force(const JV&x,const V&z)const {
  JV f(p.d,Jet(0,p.n)),oldx=p.constants(old.x);for(int i=0;i<p.d;++i)f[i]=Jet(m.mass[i/m.dim]*m.gravity[i%m.dim],p.n);
  bool dg=variant.objective=="ls_dg";
  for(int j=0;j<m.count;++j){JV e0=jedge(m,oldx,j),e1=jedge(m,x,j),em=times(plus(e0,e1),.5),g=direction(m,em);Jet strength(0,p.n);
   if(m.rigid){strength=Jet::variable(z[2*p.d+j],p.n,2*p.d+j)*(m.mass.mean()*p.vs/h);g=times(plus(e0,e1),inv(length(m,e0)+length(m,e1)));}
   else if(dg){strength=-.5*m.stiffness[j]*(length(m,e0)+length(m,e1)-2*m.length[j]);g=m.dim==1?JV{Jet(1,p.n)}:times(plus(e0,e1),inv(length(m,e0)+length(m,e1)));}
   else strength=-m.stiffness[j]*(length(m,em)-m.length[j]);
   JV fj=times(g,strength);for(int k=0;k<m.dim;++k){f[j*m.dim+k]=f[j*m.dim+k]+fj[k];if(j)f[(j-1)*m.dim+k]=f[(j-1)*m.dim+k]-fj[k];}
  }return f;
 }
 Evaluation evaluate(const V&z)const {
  Evaluation e=p.evaluate(z);JV oldx=p.constants(old.x);
  if(variant.objective!="original"){
   e.objective=Jet(0,p.n);
   if(variant.objective=="projection"){
    for(int i=0;i<p.d;++i){Jet a=(e.x[i]-target.x[i])/h,b=e.v[i]-target.v[i];e.objective=e.objective+.5*m.mass[i/m.dim]/E*(a*a+b*b);}
   }else if(variant.objective.find("ls_")==0){
    JV f=mechanical_force(e.x,z);e.force=f;
    for(int i=0;i<p.d;++i){Jet rk=(e.x[i]-old.x[i])/h-.5*(e.v[i]+old.v[i]),rp=e.v[i]-old.v[i]-h/m.mass[i/m.dim]*f[i];e.objective=e.objective+.5*m.mass[i/m.dim]/E*(rk*rk+rp*rp);}
   }else{
    double a=variant.objective=="quarter"?.25:(variant.objective=="half_half"?.5:0),b=variant.objective=="quarter"?.25:(variant.objective=="half_half"?.5:1);
    for(int i=0;i<p.d;++i){Jet r=(e.x[i]-old.x[i])/h-old.v[i]-a*h*aold[i];Jet rv=e.v[i]+old.v[i]-2*(e.x[i]-old.x[i])/h;e.objective=e.objective+.5*m.mass[i/m.dim]/E*(r*r+rv*rv/cfg.beta);}
    e.objective=e.objective+b/E*potential_jet(variant.objective=="midpotential"?times(plus(oldx,e.x),.5):e.x);
   }
  }
  if(!variant.energy)e.rows.erase(std::remove_if(e.rows.begin(),e.rows.end(),[](const Row&r){return r.kind==Global;}),e.rows.end());
  if(variant.angular&&!has_gravity(m)&&m.dim==2){Jet L(0,p.n);for(int i=0;i<m.count;++i)L=L+m.mass[i]*(e.x[2*i]*e.v[2*i+1]-e.x[2*i+1]*e.v[2*i]);e.rows.push_back({(L-angular(m,old))/(m.mass.sum()*m.length.mean()*p.vs),Global,1});}
  return e;
 }
 V initial(bool original=false)const {V z=p.initial();if(!original){z.head(p.d)=(target.x-old.x)/p.xs;z.segment(p.d,p.d)=target.v/p.vs;}if(reactions){JV q=p.coefficients(p.constants(old.x),p.constants(old.v));for(int j=0;j<m.count;++j)z[2*p.d+j]=q[j].a*h/(m.mass.mean()*p.vs);}return z;}
 State state(const Evaluation&e)const {State s{V(p.d),V(p.d)};for(int i=0;i<p.d;++i){s.x[i]=e.x[i].a;s.v[i]=e.v[i].a;}return s;}
 V encode(const State&s)const {V z=initial();z.head(p.d)=(s.x-old.x)/p.xs;z.segment(p.d,p.d)=s.v/p.vs;return z;}
};
struct Solve {State state;bool finite=true,converged=false;double objective=0,balance=0,stationarity=0;int iterations=0,evaluations=0;std::string reason;};
inline Solve solve(const Objective&p,V z,int budget=24){Solve out;out.state=p.old;try{auto e=p.evaluate(z);++out.evaluations;V lambda=V::Zero(e.rows.size());
 for(int it=0;it<budget;++it){out.iterations=it+1;Mat J=jacobian(e),H=e.objective.H;V g=residuals(e);for(int r=0;r<lambda.size();++r)H+=lambda[r]*e.rows[r].value.H;H=.5*(H+H.transpose()).eval();H.diagonal().array()+=1e-11;
  V station=e.objective.g+J.transpose()*lambda;double score=station.squaredNorm()+g.squaredNorm();if(score<1e-20)break;int n=p.p.n,q=g.size();Mat K=Mat::Zero(n+q,n+q);K.topLeftCorner(n,n)=H;if(q){K.topRightCorner(n,q)=J.transpose();K.bottomLeftCorner(q,n)=J;}V rhs(n+q);rhs.head(n)=-station;if(q)rhs.tail(q)=-g;Eigen::JacobiSVD<Mat>svd(K,Eigen::ComputeThinU|Eigen::ComputeThinV);svd.setThreshold(1e-13);V dz=svd.solve(rhs);if(!dz.allFinite())throw std::runtime_error("nonfinite_kkt");double step=std::min(1.,10./std::max(10.,dz.head(n).norm()));bool accepted=false;
  for(int ls=0;ls<12;++ls){V nz=z+step*dz.head(n),nl=lambda+step*dz.tail(q);try{auto ne=p.evaluate(nz);++out.evaluations;double ns=(ne.objective.g+jacobian(ne).transpose()*nl).squaredNorm()+residuals(ne).squaredNorm();if(ns<score||ns<1e-22){z=nz;lambda=nl;e=std::move(ne);accepted=true;break;}}catch(const std::exception&){}step*=.5;}
  if(!accepted)break;
 }
 out.state=p.state(e);out.objective=e.objective.a;out.balance=e.rows.empty()?0:residuals(e).lpNorm<Eigen::Infinity>();out.stationarity=stationarity(e,lambda);out.converged=out.balance<1e-8&&out.stationarity<1e-7;out.finite=compose::finite(out.state);if(!out.converged)out.reason="residual";
 }catch(const std::exception&e){out.finite=false;out.reason=e.what();}return out;}
inline Solve step(const Model&m,const State&s,double h,const Variant&q,double E,int budget=24){Objective p(m,s,h,q,E);return solve(p,p.initial(q.objective=="original"),budget);}
} // namespace objectives2
