#include "cross_problem.hpp"
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
using namespace local_work;
V encode(const CrossProblem&p,const State&s){V z=p.initial();z.head(p.p.d)=(s.x-p.p.old.x)/p.p.xs;z.segment(p.p.d,p.p.d)=s.v/p.p.vs;return z;}
int main(int argc,char**argv){std::filesystem::path dir=argc>1?argv[1]:"results";std::filesystem::create_directories(dir);std::ofstream csv(dir/"scalar_stalls.csv");csv<<std::setprecision(16)<<"target,step,old_x,old_v,initialization,budget,reason,converged,evaluations,primal,stationarity,objective,new_x,new_v,jacobian_min_singular,reduced_curvature\n";
 for(std::string target:{"frozen_power","implicit_midpoint_dg"}){auto c=make_case("scalar_release");const auto&m=c.model;State s=c.initial;double E=energy(m,s)-minimum_energy(m),h=1./60;auto cfg=configuration(target,"soft_midpoint");for(int step=0;step<60;++step){CrossProblem p(m,s,h,cfg,E,target,true);auto solution=joint_solve(p,128);if(!solution.converged){auto report=[&](const Solve&r,std::string initial,int budget){auto e=p.evaluate(r.z);State a=p.state(r.z);Eigen::JacobiSVD<Mat>svd(jacobian(e));csv<<target<<','<<step<<','<<s.x[0]<<','<<s.v[0]<<','<<initial<<','<<budget<<','<<r.reason<<','<<r.converged<<','<<r.evaluations<<','<<r.primal<<','<<r.dual_residual<<','<<r.objective<<','<<a.x[0]<<','<<a.v[0]<<','<<svd.singularValues().minCoeff()<<','<<reduced_curvature(p,r)<<'\n';};report(solution,"backward_euler",128);report(joint_solve(p,4096),"backward_euler",4096);
 State guess=stiffness_aware_state(m,s,h,true);guess.v=2*(guess.x-s.x)/h-s.v;V z=encode(p,guess);report(joint_solve(p,128,&z),"linear_midpoint",128);report(local_al(p,2048,10,false,"forward"),"local_AL",2048);
 // A feasible scalar circle parameterization gives an independent global scan.
 double A=std::sqrt(2*energy(m,s)/m.stiffness[0]),omega=std::sqrt(m.stiffness[0]/m.mass[0]),best=1e100,angle=0;
 auto at=[&](double t){State a=s;a.x[0]=m.length[0]+A*std::cos(t);a.v[0]=-omega*A*std::sin(t);V q=encode(p,a);q[2]=.5*m.mass[0]*(a.v.squaredNorm()-s.v.squaredNorm())/E;return q;};
 for(int i=0;i<8192;++i){double t=2*3.141592653589793*i/8192;double f=p.evaluate(at(t)).objective.a;if(f<best){best=f;angle=t;}}
 double lo=angle-2*3.141592653589793/8192,hi=angle+2*3.141592653589793/8192;for(int k=0;k<60;++k){double a=lo+(hi-lo)/3,b=hi-(hi-lo)/3;if(p.evaluate(at(a)).objective.a<p.evaluate(at(b)).objective.a)hi=b;else lo=a;}V q=at(.5*(lo+hi));report(joint_solve(p,128,&q),"feasible_circle_scan",128);
 }s=p.state(solution.z);}}
 std::cout<<"scalar stall diagnosis complete\n";
}
