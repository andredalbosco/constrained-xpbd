#pragma once
#include "problem.hpp"
namespace hour_local {
struct Difference {double objective=0;V rows;};
// Exact algebraic differences, avoiding subtraction of the large constant
// gravitational potential or two nearly equal squared norms.
inline Difference difference(const Problem&p,const V&from,const V&to,const Evaluation&e){
    const auto&m=p.p.m;int d=p.p.d,dim=m.dim,offset=2*d;double E=p.p.E,xs=p.p.xs,vs=p.p.vs;V dz=to-from,dx=xs*dz.head(d),dv=vs*dz.segment(d,d),x=p.p.old.x+xs*from.head(d),v=vs*from.segment(d,d);Difference out;out.rows=V::Zero(e.rows.size());long double F=0;
    for(int i=0;i<d;++i){double r=vs*from[i]-p.p.old.v[i]-.25*p.fraction*p.p.h*p.aold[i],rv=vs*from[d+i]+p.p.old.v[i]-2*vs*from[i],dr=vs*dz[i],drv=vs*(dz[d+i]-2*dz[i]);F+=.5L*m.mass[i/dim]/E*(static_cast<long double>(dr)*(2*r+dr)+static_cast<long double>(drv)*(2*rv+drv)/p.cfg.beta);}
    V DU=V::Zero(m.count+(has_gravity(m)?m.count:0)),DL=V::Zero(m.count);
    for(int j=0;j<m.count;++j){V ec=edge(m,x,j),de=edge(m,dx,j);double r=edge_length(m,ec),rn=edge_length(m,V(ec+de)),dr=dim==1?de[0]:de.dot(2*ec+de)/(r+rn);DL[j]=dr;if(!m.rigid)DU[j]=.5*m.stiffness[j]*dr*(2*(r-m.length[j])+dr);}
    if(has_gravity(m))for(int i=0;i<m.count;++i)DU[m.count+i]=-m.mass[i]*m.gravity.dot(dx.segment(i*dim,dim));
    F+=(1-.75*p.fraction)/E*DU.sum();for(int j=0;j<p.p.wcount;++j){double a=from[offset+j]-p.frozen[j]/E,b=dz[offset+j];F+=.5L*p.cfg.gamma*b*(2*a+b);}out.objective=double(F);
    for(int r=0;r<int(e.rows.size());++r){const auto&row=e.rows[r];int i=row.index;double value=0;if(row.kind==Kinetic){value=.5*m.mass[i]/E*dv.segment(i*dim,dim).dot(2*v.segment(i*dim,dim)+dv.segment(i*dim,dim));for(int j:{p.p.work_child[i],i+1<m.count?p.p.work_parent[i+1]:-1,p.p.work_gravity[i]})if(j>=0)value-=dz[offset+j];}
      else if(row.kind==Storage){value=DU[i]/E+dz[offset+p.p.work_child[i]];if(i)value+=dz[offset+p.p.work_parent[i]];}
      else if(row.kind==GravityStorage)value=DU[m.count+i]/E+dz[offset+p.p.work_gravity[i]];
      else if(row.kind==Position)value=DL[i]/xs;
      else if(row.kind==Tangency){V n=x/x.norm(),dn=(dx-n*DL[i])/(x.norm()+DL[i]);value=(dn.dot(v)+n.dot(dv)+dn.dot(dv))/vs;}
      out.rows[r]=value;
    }return out;
}
inline double augmented_difference(const Problem&p,const V&from,const V&to,const V&dual,double rho,const Evaluation&e){auto delta=difference(p,from,to,e);long double value=delta.objective;for(int i=0;i<dual.size();++i)value+=(static_cast<long double>(dual[i])+rho*e.rows[i].value.a)*delta.rows[i]+.5L*rho*delta.rows[i]*delta.rows[i];return double(value);}
}
