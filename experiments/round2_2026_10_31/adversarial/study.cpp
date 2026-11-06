#include "mechanics.hpp"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <chrono>
#include <filesystem>
#include <sstream>
#include <omp.h>
using namespace compose;
struct Result {State s; bool ok=true;std::string reason;int iter=0,eval=0,leaves=1,reject=0;double residual=0;};
struct Scale {double l=1,v=1,e=1;};
Scale scale_for(const Case& c){Scale s;s.l=c.model.length.mean();s.e=std::max(1e-4,energy(c.model,c.initial)-minimum_energy(c.model));s.v=std::sqrt(s.e/c.model.mass.mean());return s;}
V encode(const State&s,const Scale&q){V y(s.x.size()*2);y<<s.x/q.l,s.v/q.v;return y;}
State decode(const V&y,const Scale&q){int n=y.size()/2;return {q.l*y.head(n),q.v*y.tail(n)};}
V ode(const Model&m,const V&y,const Scale&q){State s=decode(y,q);V a=masses(m).cwiseInverse().cwiseProduct(force(m,s.x));if(m.rigid){double r2=s.x.squaredNorm();if(r2<1e-20)throw std::runtime_error("zero_radius_stage");a-=s.x*((s.v.squaredNorm()+s.x.dot(a))/r2);}V r(y.size());r<<s.v/q.l,a/q.v;return r;}
V newton(V z,const mechanical_detail::Residual&f,int limit,Result&o){
    for(int it=0;it<limit;++it){V r=f(z);++o.eval;o.residual=r.norm();if(!r.allFinite()){o.ok=false;o.reason="nonfinite_residual";return z;}if(r.norm()<2e-11)return z;
        ++o.iter;Mat j(r.size(),z.size());for(int k=0;k<z.size();++k){double eps=1e-6*std::max(1.,std::abs(z[k]));V a=z,b=z;a[k]+=eps;b[k]-=eps;j.col(k)=(f(a)-f(b))/(2*eps);o.eval+=2;}
        V dz=j.completeOrthogonalDecomposition().solve(-r);if(!dz.allFinite()){o.ok=false;o.reason="nonfinite_increment";return z;}bool accepted=false;double alpha=1;for(int k=0;k<20;++k){V p=z+alpha*dz,rr=f(p);++o.eval;if(rr.allFinite()&&rr.norm()<r.norm()){z=p;accepted=true;break;}alpha*=.5;}if(!accepted){o.ok=false;o.reason="no_descent";return z;}
    }o.residual=f(z).norm();++o.eval;if(o.residual>2e-9){o.ok=false;o.reason="iteration_limit";}return z;
}
V invariants(const Model&m,const V&y,const Scale&q,double target,bool E){State s=decode(y,q);V r(E?3:2);r[0]=.5*(s.x.squaredNorm()/(q.l*q.l)-1);r[1]=s.x.dot(s.v)/(q.l*q.v);if(E)r[2]=(energy(m,s)-target)/q.e;return r;}
Mat inv_gradient(const Model&m,const V&y,const Scale&q,bool E){State s=decode(y,q);Mat a=Mat::Zero(E?3:2,4);a.block(0,0,1,2)=(s.x/q.l).transpose();a.block(1,0,1,2)=(s.v/q.v).transpose();a.block(1,2,1,2)=(s.x/q.l).transpose();if(E){a.block(2,0,1,2)=(-q.l*force(m,s.x)/q.e).transpose();a.block(2,2,1,2)=(q.v*masses(m).cwiseProduct(s.v)/q.e).transpose();}return a;}
Result projected_midpoint(const Model&m,const State&old,double h,const Scale&q,const std::string&method){
    Result o;V y0=encode(old,q);const bool E=method.find("energy")!=std::string::npos,sym=method.find("sym_")==0;const int k=E?3:2;const double target=energy(m,old);Mat G0=inv_gradient(m,y0,q,E);
    V yguess=y0+h*ode(m,y0,q);if(sym){V z=V::Zero(4+k);z.head(4)=yguess;mechanical_detail::Residual f=[&](const V&z)->V{V y=z.head(4),mu=z.tail(k),a=y0+G0.transpose()*mu,b=y-inv_gradient(m,y,q,E).transpose()*mu;V r(4+k);r.head(4)=b-a-h*ode(m,.5*(a+b),q);r.tail(k)=invariants(m,y,q,target,E);return r;};z=newton(z,f,64,o);o.s=decode(z.head(4),q);}
    else {mechanical_detail::Residual f=[&](const V&y)->V{return y-y0-h*ode(m,.5*(y+y0),q);};V target_y=newton(yguess,f,64,o);if(!o.ok){o.s=decode(target_y,q);return o;}if(method=="ode_midpoint"){o.s=decode(target_y,q);return o;}V z=V::Zero(4+k);z.head(4)=target_y;mechanical_detail::Residual p=[&](const V&z)->V{V y=z.head(4);V r(4+k);r.head(4)=y-target_y-inv_gradient(m,y,q,E).transpose()*z.tail(k);r.tail(k)=invariants(m,y,q,target,E);return r;};z=newton(z,p,64,o);o.s=decode(z.head(4),q);}
    return o;
}
Result raw_step(const Model&m,const State&old,double h,const Scale&q,std::string method){
    Result o;o.s=old;try{
    if(method=="angle_dg"){
        const double l=m.length[0],theta=std::atan2(old.x[0],-old.x[1]),w=angular(m,old)/(m.mass[0]*l*l);
        auto acceleration=[&](double d){double z=.5*d,sinc=std::abs(z)<1e-5?1-z*z/6+z*z*z*z/120:std::sin(z)/z;return (m.gravity[0]*std::cos(theta+z)+m.gravity[1]*std::sin(theta+z))*sinc/l;};
        mechanical_detail::Residual f=[&](const V&d)->V{V r(1);r[0]=d[0]-h*w-.5*h*h*acceleration(d[0]);return r;};
        V d(1);d[0]=h*w+.5*h*h*acceleration(0);d=newton(d,f,64,o);double t=theta+d[0],wn=w+h*acceleration(d[0]);o.s.x.resize(2);o.s.v.resize(2);o.s.x<<l*std::sin(t),-l*std::cos(t);o.s.v<<l*wn*std::cos(t),l*wn*std::sin(t);return o;
    }
    if(method.find("sym_")==0||method.find("post_")==0||method=="ode_midpoint")return projected_midpoint(m,old,h,q,method);
    if(method=="implicit_midpoint"){V y0=encode(old,q);mechanical_detail::Residual f=[&](const V&y)->V{return y-y0-h*ode(m,.5*(y+y0),q);};V y=newton(y0+h*ode(m,y0,q),f,64,o);o.s=decode(y,q);return o;}
    Config c;c.iterations=method=="dg4"?4:method=="dg16"?16:64;c.family="force_momentum";c.work="dg_average";c.schedule="joint";StepResult r;
    if(method.find("rattle")==0||method=="shake"){r=mechanical_detail::verlet_projection(m,old,h,c,q.e,method!="shake");}
    else r=mechanical_detail::force_momentum(m,old,h,c,q.e);
    o.s=r.state;o.ok=r.valid;o.reason=r.reason;o.iter=r.iterations;o.eval=r.evaluations;o.residual=r.rx/q.v;
    if(method.find("rattle")==0||method=="shake"){
        o.residual=mechanical_detail::constraints(m,o.s.x).norm()/q.l;
        if(method!="shake")o.residual=std::max(o.residual,(mechanical_detail::gradient(m,o.s.x)*o.s.v).norm()/q.v);
    }
    if(method=="dg_tangent"||method=="dg_tangent_energy"||method=="rattle_energy"){
       Mat G=mechanical_detail::gradient(m,o.s.x);V inv=masses(m).cwiseInverse();V correction=(G*inv.asDiagonal()*G.transpose()).completeOrthogonalDecomposition().solve(-G*o.s.v);o.s.v+=inv.cwiseProduct(G.transpose()*correction);
       if(method!="dg_tangent"){double needed=energy(m,old)-potential(m,o.s.x),kin=.5*(masses(m).array()*o.s.v.array().square()).sum();if(needed < -2e-11*q.e){o.ok=false;o.reason="negative_target_kinetic";return o;}needed=std::max(0.,needed);if(kin<1e-24&&needed>2e-11*q.e){o.ok=false;o.reason="zero_direction_for_energy";return o;}if(kin>1e-24)o.s.v*=std::sqrt(needed/kin);else o.s.v.setZero();}
    }
    if(!finite(o.s)){o.ok=false;o.reason="nonfinite_state";}
    }catch(const std::exception&e){o.ok=false;o.reason=e.what();}return o;
}
Result advance(const Model&m,const State&s,double h,const Scale&q,std::string method,int depth=0){
    if(method.find("accuracy")!=std::string::npos){
        const std::string base=method.find("rattle")==0?"rattle":"dg16";
        Result a=raw_step(m,s,h,q,base),b=raw_step(m,s,.5*h,q,base),c=b.ok?raw_step(m,b.s,.5*h,q,base):b;
        bool pass=a.ok&&b.ok&&c.ok&&std::max({a.residual,b.residual,c.residual})<2e-8;
        if(pass)pass=((a.s.x-c.s.x).norm()/q.l+(a.s.v-c.s.v).norm()/q.v)<1e-4*std::abs(h);
        int it=a.iter+b.iter+(b.ok?c.iter:0),ev=a.eval+b.eval+(b.ok?c.eval:0);
        if(pass){c.iter=it;c.eval=ev;c.leaves=2;return c;}
        if(depth>=14){c.ok=false;c.reason="accuracy_depth_limit";c.iter=it;c.eval=ev;return c;}
        Result l=advance(m,s,.5*h,q,method,depth+1);if(!l.ok){l.iter+=it;l.eval+=ev;++l.reject;return l;}
        Result r=advance(m,l.s,.5*h,q,method,depth+1);r.iter+=it+l.iter;r.eval+=ev+l.eval;r.leaves+=l.leaves;r.reject+=l.reject+1;return r;
    }
    bool adapt=method.find("adaptive")!=std::string::npos;std::string base=method;if(adapt)base=method.find("rattle")==0?"rattle":method.find("sym_")==0?"sym_energy":"dg16";
    Result a=raw_step(m,s,h,q,base);bool pass=a.ok&&a.residual<2e-8;if(m.rigid)pass=pass&&mechanical_detail::constraints(m,a.s.x).norm()/q.l<2e-8;
    if(!adapt)return a;
    if(pass)return a;if(depth>=10){a.ok=false;a.reason="adaptive_depth_limit:"+a.reason;return a;}
    Result b=advance(m,s,.5*h,q,method,depth+1);if(!b.ok){b.iter+=a.iter;b.eval+=a.eval;b.reject+=a.reject+1;return b;}Result c=advance(m,b.s,.5*h,q,method,depth+1);c.iter+=a.iter+b.iter;c.eval+=a.eval+b.eval;c.leaves+=b.leaves;c.reject+=a.reject+b.reject+1;return c;
}
Case rod(const std::string&name,double theta,double omega,double g=9.81){Case c=make_case("rod_release");c.name=name;c.model.gravity[1]=-g;c.initial.x<<std::sin(theta),-std::cos(theta);c.initial.v<<omega*std::cos(theta),omega*std::sin(theta);return c;}
std::vector<Case> cases(){std::vector<Case> a;
 for(auto n:case_names())a.push_back(make_case(n));
 a.push_back(rod("rod_separatrix_below",0,2*std::sqrt(9.81)*.999));a.push_back(rod("rod_separatrix_above",0,2*std::sqrt(9.81)*1.001));a.push_back(rod("rod_fast_spin",0,20));a.push_back(rod("rod_near_upright",3.13,0));a.push_back(rod("rod_free_spin",.2,8,0));
 Case c=make_case("radial_orbit");c.name="radial_compressed";c.initial.x<<.02,0;c.initial.v<<0,2;a.push_back(c);c.name="radial_fast";c.initial.x<<1.1,0;c.initial.v<<-2,20;a.push_back(c);c.name="radial_stiff";c.model.stiffness[0]=1e5;c.initial.x<<1.1,0;c.initial.v<<0,3;a.push_back(c);
 c=make_case("elastic_mixed",100);c.name="elastic_stiff";a.push_back(c);c=make_case("chain_large",10,100,2);c.name="chain_heavy";a.push_back(c);return a;}
struct Ref {std::vector<State>a,b;double gap=0,energy_gap=0;bool ok=true;};
Ref reference(const Case&c,double h,int steps){Ref r;State a=c.initial,b=a;Scale q=scale_for(c);r.a.push_back(a);r.b.push_back(b);double omega=std::sqrt(c.model.stiffness.maxCoeff()/c.model.mass.minCoeff());int sub=std::max(64,int(std::ceil(std::abs(h)*omega*40)));for(int i=0;i<steps;++i){try{a=reference_step(c.model,a,h,sub);b=reference_step(c.model,b,h,2*sub);}catch(...){r.ok=false;break;}if(!finite(a)||!finite(b)){r.ok=false;r.gap=std::numeric_limits<double>::infinity();break;}r.a.push_back(a);r.b.push_back(b);r.gap=std::max(r.gap,(a.x-b.x).norm()/q.l+(a.v-b.v).norm()/q.v);r.energy_gap=std::max(r.energy_gap,std::abs(energy(c.model,b)-energy(c.model,c.initial))/q.e);}r.ok=r.ok&&r.gap<1e-5&&r.energy_gap<1e-6;return r;}
std::vector<std::string> methods(bool rigid){std::vector<std::string> ms={"dg4","dg16","dg64","dg_adaptive","dg_accuracy"};if(!rigid)ms.push_back("implicit_midpoint");else {for(std::string m:{"shake","rattle","rattle_adaptive","rattle_accuracy","dg_tangent","dg_tangent_energy","rattle_energy","sym_geom","sym_energy","post_geom","post_energy","ode_midpoint","sym_energy_adaptive"})ms.push_back(m);}return ms;}
std::string run(const Case&c,const Ref&ref,double h,double duration,const std::string&method,const std::string&folder,bool trace){
 Scale q=scale_for(c);State s=c.initial;double E0=energy(c.model,s),L0=angular(c.model,s),maxE=0,maxg=0,maxt=0,maxkin=0,maximp=0,maxrx=0,maxrev=0,maxL=0,sumx=0,sumv=0;int steps=std::lround(duration/h),done=0,iter=0,eval=0,reject=0,leaves=0,unconv=0,matched=0;std::string status="complete",reason;
 std::ofstream tr;if(trace){tr.open(folder+"/trace_"+c.name+"_"+method+"_"+std::to_string(h)+".csv");tr<<std::setprecision(16)<<"t";for(int k=0;k<dof(c.model);++k)tr<<",x"<<k;for(int k=0;k<dof(c.model);++k)tr<<",v"<<k;tr<<",energy,reference_x0,reference_x1\n";}
 auto tic=std::chrono::steady_clock::now();
 for(int i=0;i<steps;++i){State old=s;Result o=advance(c.model,old,h,q,method);iter+=o.iter;eval+=o.eval;reject+=o.reject;leaves+=o.leaves;if(!o.ok||!finite(o.s)){status="failed";reason=o.reason;break;}s=o.s;++done;if(o.residual>2e-8)++unconv;maxrx=std::max(maxrx,o.residual);maxE=std::max(maxE,std::abs(energy(c.model,s)-E0)/q.e);maxL=std::max(maxL,std::abs(angular(c.model,s)-L0)/(c.model.mass.sum()*q.l*q.v));
 if(c.model.rigid){maxg=std::max(maxg,mechanical_detail::constraints(c.model,s.x).norm()/q.l);maxt=std::max(maxt,(mechanical_detail::gradient(c.model,s.x)*s.v).norm()/q.v);}
 maxkin=std::max(maxkin,((s.x-old.x)/h-.5*(s.v+old.v)).norm()/q.v);
 V f=mechanical_detail::gravity_force(c.model),mass=masses(c.model);if(c.model.rigid){V dir=mechanical_detail::rod_direction(c.model,old,s.x,0,"dg_average"),delta=mass.cwiseProduct(s.v-old.v)-h*f;double p=dir.squaredNorm()>1e-24?dir.dot(delta)/dir.squaredNorm():0;maximp=std::max(maximp,(delta-p*dir).norm()/(c.model.mass.mean()*q.v));}else{for(int j=0;j<c.model.count;++j){V fj=mechanical_detail::edge_force(c.model,old,s.x,j,"dg_average");f.segment(j*c.model.dim,c.model.dim)+=fj;if(j)f.segment((j-1)*c.model.dim,c.model.dim)-=fj;}maximp=std::max(maximp,(mass.cwiseProduct(s.v-old.v)-h*f).norm()/(c.model.mass.mean()*q.v));}
 if(i<32){Result back=advance(c.model,s,-h,q,method);double rev=back.ok?((back.s.x-old.x).norm()/q.l+(back.s.v-old.v).norm()/q.v):1e9;maxrev=std::max(maxrev,rev);}
 if(i+1<int(ref.b.size())){++matched;sumx+=(s.x-ref.b[i+1].x).squaredNorm()/c.model.count;sumv+=(s.v-ref.b[i+1].v).squaredNorm()/c.model.count;}
 if(trace){tr<<(i+1)*h;for(double z:s.x)tr<<","<<z;for(double z:s.v)tr<<","<<z;tr<<","<<energy(c.model,s)<<","<<ref.b[std::min(i+1,int(ref.b.size())-1)].x[0]<<","<<(dof(c.model)>1?ref.b[std::min(i+1,int(ref.b.size())-1)].x[1]:0)<<"\n";}
 if(s.x.norm()>1e6||maxE>1e8){status="runaway";reason="state_or_energy_limit";break;}
 }
 double ms=std::chrono::duration<double,std::milli>(std::chrono::steady_clock::now()-tic).count();std::ostringstream os;os<<std::setprecision(16)<<c.name<<","<<method<<","<<h<<","<<duration<<","<<done<<","<<steps<<","<<status<<","<<reason<<","<<ref.ok<<","<<ref.gap<<","<<ref.energy_gap<<","<<(matched?std::sqrt(sumx/matched):std::numeric_limits<double>::quiet_NaN())<<","<<(matched?std::sqrt(sumv/matched):std::numeric_limits<double>::quiet_NaN())<<","<<maxE<<","<<maxg<<","<<maxt<<","<<maxkin<<","<<maximp<<","<<maxrx<<","<<maxrev<<","<<maxL<<","<<unconv<<","<<iter<<","<<eval<<","<<reject<<","<<leaves<<","<<ms;return os.str();}
void validate(){
 Case c=make_case("scalar_release");Scale q=scale_for(c);auto a=raw_step(c.model,c.initial,.01,q,"dg64");auto b=raw_step(c.model,c.initial,.01,q,"implicit_midpoint");if(!a.ok||!b.ok||(a.s.x-b.s.x).norm()>1e-9||(a.s.v-b.s.v).norm()>1e-8)throw std::runtime_error("quadratic_midpoint_equivalence");
 for(auto meth:{"sym_geom","sym_energy","post_geom","post_energy"}){c=make_case("rod_kick");q=scale_for(c);a=raw_step(c.model,c.initial,.01,q,meth);if(!a.ok||invariants(c.model,encode(a.s,q),q,energy(c.model,c.initial),std::string(meth).find("energy")!=std::string::npos).norm()>1e-8)throw std::runtime_error("projection_invariants");if(std::string(meth).find("sym_")==0){b=raw_step(c.model,a.s,-.01,q,meth);if(!b.ok||(encode(b.s,q)-encode(c.initial,q)).norm()>1e-8)throw std::runtime_error("symmetric_reversal");}}
 c=rod("no_rattle_root",0,20,0);q=scale_for(c);a=raw_step(c.model,c.initial,.1,q,"rattle");if(a.residual<1e-4)throw std::runtime_error("missed_no_real_rattle_root");b=advance(c.model,c.initial,.1,q,"rattle_adaptive");if(!b.ok||b.reject==0||b.residual>2e-8)throw std::runtime_error("adaptive_rattle_recovery");
 for(auto rc:{make_case("rod_release"),rod("rotation",0,20)}){q=scale_for(rc);a=raw_step(rc.model,rc.initial,.25,q,"angle_dg");b=raw_step(rc.model,a.s,-.25,q,"angle_dg");if(!a.ok||!b.ok||std::abs(energy(rc.model,a.s)-energy(rc.model,rc.initial))/q.e>1e-8||(encode(b.s,q)-encode(rc.initial,q)).norm()>1e-8)throw std::runtime_error("angle_dg_energy_reversal");}
 // Analytic invariant derivatives, including off-manifold values.
 c=make_case("rod_kick");q=scale_for(c);V y=encode(c.initial,q);y[0]+=.1;y[2]-=.2;Mat G=inv_gradient(c.model,y,q,true),fd=G;for(int k=0;k<4;++k){V u=y,v=y;u[k]+=1e-6;v[k]-=1e-6;fd.col(k)=(invariants(c.model,u,q,0,true)-invariants(c.model,v,q,0,true))/(2e-6);}if((G-fd).norm()>1e-8)throw std::runtime_error("invariant_derivatives");
 std::cout<<"validation passed: quadratic equivalence, four projections, symmetric reversal, missing RATTLE root, adaptive recovery, analytic invariant gradients\n";
}
int main(int argc,char**argv){try{validate();if(argc>1&&std::string(argv[1])=="validate")return 0;std::string folder=argc>1?argv[1]:"results",mode=argc>2?argv[2]:"short";bool angles=mode.find("angles")!=std::string::npos;std::filesystem::create_directories(folder);auto cs=cases();std::vector<double> hs={1./120,1./30,1./10,1./4};double duration=10;if(mode.find("long")!=std::string::npos){duration=60;hs={1./120,1./30};}if(angles||duration==60){cs.erase(std::remove_if(cs.begin(),cs.end(),[](const Case&c){return !c.model.rigid;}),cs.end());}
 struct Batch{Case c;double h;};std::vector<Batch>batches;for(const auto&c:cs)for(double h:hs)batches.push_back({c,h});std::vector<std::string> rows(batches.size());
 #pragma omp parallel for schedule(dynamic)
 for(int j=0;j<int(batches.size());++j){auto&ba=batches[j];Ref ref=reference(ba.c,ba.h,std::lround(duration/ba.h));std::ostringstream out;for(const auto&method:angles?std::vector<std::string>{"angle_dg"}:methods(ba.c.model.rigid)){out<<run(ba.c,ref,ba.h,duration,method,folder,ba.h==1./30)<<"\n";}rows[j]=out.str();
 #pragma omp critical
 std::cout<<"completed "<<j+1<<"/"<<batches.size()<<" "<<ba.c.name<<" h="<<ba.h<<" ref="<<ref.ok<<std::endl;
 }
 std::ofstream out(folder+"/summary.csv");out<<"case,method,h,duration,steps_done,steps_requested,status,reason,reference_certified,reference_gap,reference_energy_gap,position_rms,velocity_rms,max_energy_relative,max_length_relative,max_tangency_scaled,max_midpoint_defect,max_dg_impulse_defect,max_solver_residual,max_local_reversal,max_angular_change,unconverged_steps,newton_iterations,residual_evaluations,rejected_attempts,accepted_substeps,wall_ms\n";for(auto&r:rows)out<<r;return 0;}catch(const std::exception&e){std::cerr<<e.what()<<"\n";return 1;}}
