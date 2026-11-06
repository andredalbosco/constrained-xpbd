#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

// Work-based central force at a midpoint stage. This is a new discrete force
// equation, not a solver of the old endpoint-energy-constrained objective.
// Each visit is one linearized scalar multiplier solve. The position Jacobian
// drops multiplier curvature; the current position residual is retained.
namespace {
constexpr double mass=1,stiffness=100,rest_length=1,omega=10;
constexpr double pi=3.14159265358979323846;
struct V { double x=0,y=0; };
V operator+(V a,V b){return {a.x+b.x,a.y+b.y};}
V operator-(V a,V b){return {a.x-b.x,a.y-b.y};}
V operator*(V a,double s){return {a.x*s,a.y*s};}
V operator/(V a,double s){return a*(1/s);}
double dot(V a,V b){return a.x*b.x+a.y*b.y;}
double cross(V a,V b){return a.x*b.y-a.y*b.x;}
double norm(V a){return std::hypot(a.x,a.y);}
bool finite(V a){return std::isfinite(a.x)&&std::isfinite(a.y);}
struct State { V x,v; };
struct Case { std::string name;State initial;bool tangential; };
const std::array<Case,3> cases{{
    {"spring_rest",{{1.25,0},{0,0}},false},
    {"spring_radial",{{1.25,0},{-.7,0}},false},
    {"spring_tangent",{{1.25,0},{0,1}},true}
}};
double potential(V x){const double c=norm(x)-rest_length;return .5*stiffness*c*c;}
double energy(State s){return .5*mass*dot(s.v,s.v)+potential(s.x);}
double angular(State s){return mass*cross(s.x,s.v);}
V force(V x){const double r=norm(x);if(!(r>0))throw std::runtime_error("undefined spring force direction");return x*(-stiffness*(1-rest_length/r));}
double wrapped(double a){return std::remainder(a,2*pi);}
double phase(State s,bool tangential){
    return tangential?std::atan2(s.x.y,s.x.x):std::atan2(-dot(s.v,s.x/norm(s.x))/omega,norm(s.x)-rest_length);
}
State reference_step(State s,double h){
    const V a1=force(s.x)/mass,q1=s.v;
    const V q2=s.v+a1*(h/2),a2=force(s.x+q1*(h/2))/mass;
    const V q3=s.v+a2*(h/2),a3=force(s.x+q2*(h/2))/mass;
    const V q4=s.v+a3*h,a4=force(s.x+q3*h)/mass;
    return {s.x+(q1+q2*2+q3*2+q4)*(h/6),s.v+(a1+a2*2+a3*2+a4)*(h/6)};
}
std::vector<State> reference(State s,double h,int steps,int subdivisions){
    std::vector<State> states{s};
    for(int i=0;i<steps;++i){for(int j=0;j<subdivisions;++j)s=reference_step(s,h/subdivisions);states.push_back(s);}
    return states;
}
struct Rows {V n,g,R,force;double C=0,rc=0,denominator=0;const char* failure="";};
bool evaluate(State old,V predictor,V z,double lambda,double h,Rows& out){
    const V xnew=z*2-old.x;
    const double rz=norm(z),rnew=norm(xnew),rold=norm(old.x),rbar=(rnew+rold)/2;
    if(!finite(z)){out.failure="nonfinite_stage";return false;}
    if(!(rbar>0)){out.failure="undefined_mean_radius";return false;}
    if(!(rz>0)){out.failure="undefined_stage_radius";return false;}
    if(!(rnew>0)){out.failure="undefined_endpoint_radius";return false;}
    if(!(rold>0)){out.failure="undefined_initial_radius";return false;}
    const double at=4/(stiffness*h*h),factor=1-rest_length/rbar;
    out.n=z/rz;
    out.C=rz*factor;
    out.g=out.n*factor+(xnew/rnew)*(rz*rest_length/(rbar*rbar));
    out.R=(z-predictor)*mass-out.n*lambda;
    out.rc=out.C+at*lambda;
    out.force=z*(-stiffness*factor);
    out.denominator=dot(out.g,out.n)/mass+at;
    if(!finite(out.g)||!finite(out.R)||!finite(out.force)||!std::isfinite(out.rc)||!std::isfinite(out.denominator)){
        out.failure="nonfinite_row";return false;
    }
    if(!(out.denominator>0)){out.failure="nonpositive_denominator";return false;}
    out.failure="";return true;
}
struct Step {
    State next;
    bool ok=true;
    std::string status="completed";
    double rc=0,rx=0,kinematic=0,force_defect=0,kinetic_work=0,potential_work=0;
    double balance_identity=0,angular_identity=0;
};
Step step(State old,double h,int sweeps){
    const V predictor=old.x+old.v*(h/2);
    V z=predictor;double lambda=0;Rows r;
    Step out;
    for(int i=0;i<sweeps;++i){
        if(!evaluate(old,predictor,z,lambda,h,r)){out.ok=false;out.status=r.failure;return out;}
        const V base=r.R*(-1/mass);
        const double dl=(-r.rc-dot(r.g,base))/r.denominator;
        z=z+base+r.n*(dl/mass);lambda+=dl;
        if(!finite(z)||!std::isfinite(lambda)||norm(z)>1e6||std::abs(lambda)>1e18){out.ok=false;out.status="guarded";return out;}
    }
    if(!evaluate(old,predictor,z,lambda,h,r)){out.ok=false;out.status=r.failure;return out;}
    out.next={z*2-old.x,(z-old.x)*(4/h)-old.v};
    if(!finite(out.next.v)||norm(out.next.v)>1e6){out.ok=false;out.status="guarded";return out;}
    const V dx=out.next.x-old.x,dv=out.next.v-old.v;
    const double W=dot(r.force,dx);
    const double dK=.5*mass*(dot(out.next.v,out.next.v)-dot(old.v,old.v));
    const double dU=potential(out.next.x)-potential(old.x);
    out.rc=std::abs(r.rc);out.rx=norm(r.R)*4/(h*h);
    out.kinematic=norm(dx-(out.next.v+old.v)*(h/2));
    out.force_defect=norm(dv*mass-r.force*h);
    out.kinetic_work=std::abs(dK-W);out.potential_work=std::abs(dU+W);
    // Exact algebraic defect decomposition; separates force-equation and row errors.
    out.balance_identity=std::abs(dK+dU-stiffness*r.rc*dot(r.n,dx)-4*dot(r.R,dx)/(h*h));
    out.angular_identity=std::abs(angular(out.next)-angular(old)-4*cross(z,r.R)/h);
    return out;
}
void require(bool ok,const char* message){if(!ok)throw std::runtime_error(message);}
void checks(std::ostream& out){
    const State old{{1.18,.17},{-.3,.4}};const V z{1.14,.19};const double h=.04,lambda=-.003;
    const V predictor=old.x+old.v*(h/2);Rows row;
    require(evaluate(old,predictor,z,lambda,h,row),"gradient fixture undefined");
    double gradient_error=0,work_error=0;
    for(V direction:{V{.6,.8},V{-.8,.6}}){
        constexpr double eps=1e-6;Rows plus,minus;
        require(evaluate(old,predictor,z+direction*eps,lambda,h,plus)&&evaluate(old,predictor,z-direction*eps,lambda,h,minus),"finite difference fixture undefined");
        gradient_error=std::max(gradient_error,std::abs((plus.C-minus.C)/(2*eps)-dot(row.g,direction)));
    }
    for(const State initial:{cases[0].initial,cases[1].initial,cases[2].initial,old}){
        for(V endpoint:{V{.8,.3},V{1.15,-.2},V{1.25,.01}}){
            const V middle=(initial.x+endpoint)/2;Rows q;
            require(evaluate(initial,initial.x,middle,0,h,q),"work fixture undefined");
            work_error=std::max(work_error,std::abs(potential(endpoint)-potential(initial.x)+dot(q.force,endpoint-initial.x)));
        }
    }
    const double dt=1./60,ratio=omega*dt/2,q0=.25;
    const Step one=step(cases[0].initial,dt,1);
    const double expected_q=q0*(1-ratio*ratio)/(1+ratio*ratio),expected_v=-q0*omega*2*ratio/(1+ratio*ratio);
    const double harmonic_error=std::max(std::abs(one.next.x.x-1-expected_q),std::abs(one.next.v.x-expected_v));
    out<<"work_scalar_gradient_finite_difference,"<<gradient_error<<",2e-9\n"
       <<"potential_work_identity,"<<work_error<<",1e-12\n"
       <<"harmonic_midpoint_one_step,"<<harmonic_error<<",1e-12\n";
    require(gradient_error<2e-9,"work scalar gradient incorrect");
    require(work_error<1e-12,"central averaged-force work identity failed");
    require(one.ok&&harmonic_error<1e-12,"harmonic midpoint one-step formula failed");
    Rows invalid;
    require(!evaluate({{1,0},{}},{1,0},{0,0},0,h,invalid)&&std::string(invalid.failure)=="undefined_stage_radius",
            "zero stage radius must be rejected explicitly");
    require(!evaluate({{0,0},{}},{0,0},{0,0},0,h,invalid)&&std::string(invalid.failure)=="undefined_mean_radius",
            "zero mean radius must be rejected explicitly");
    require(!evaluate({{1,0},{}},{1,0},{.5,0},0,h,invalid)&&std::string(invalid.failure)=="undefined_endpoint_radius",
            "zero endpoint radius must be rejected explicitly");
    require(!evaluate({{1,0},{}},{1,0},{.1,0},0,1,invalid)&&std::string(invalid.failure)=="nonpositive_denominator",
            "nonpositive scalar denominator must be rejected explicitly");
    out<<"explicit_radius_and_denominator_guards,0,0\n";
}
} // namespace

int main(int argc,char** argv){
try{
    const std::filesystem::path dir=argc>1?argv[1]:"results/spring_repair/work_midpoint";
    std::filesystem::create_directories(dir);
    std::ofstream summary(dir/"summary.csv"),traces(dir/"traces.csv"),validation(dir/"checks.csv");
    require(summary&&traces&&validation,"cannot open work-midpoint outputs");
    summary<<std::setprecision(17);traces<<std::setprecision(17);validation<<std::setprecision(17);
    validation<<"check,value,tolerance\n";checks(validation);
    summary<<"case,method,h,sweeps,duration,completed_steps,requested_steps,status,E0,L0,final_energy,max_relative_energy_error,max_angular_error,"
             "rms_position_error,max_position_error,final_phase_error,max_phase_error,path_length,max_material_residual,max_stationarity,"
             "max_midpoint_kinematic_residual,max_midpoint_force_residual,max_kinetic_work_defect,max_potential_work_defect,max_balance_identity_error,max_angular_identity_error,"
             "reference_refinement_error,reference_energy_error\n";
    traces<<"case,method,h,sweeps,duration,time,x,y,vx,vy,energy,angular_momentum,phase_error,ref_x,ref_y,ref_vx,ref_vy,ref_energy\n";
    int runs=0,failures=0;double worst_energy_converged=0,worst_angular=0,worst_work_identity=0,worst_balance_identity=0,worst_angular_identity=0;
    for(const Case& c:cases)for(double h:{1./60,1./240})for(double duration:{10.,100.}){
        if(duration==100&&(!c.tangential||h!=1./60))continue;
        const int count=static_cast<int>(std::lround(duration/h));
        const auto ref=reference(c.initial,h,count,64),fine=reference(c.initial,h,count,128);
        const double E0=energy(c.initial),L0=angular(c.initial);
        double refinement=0,ref_energy=0;
        for(int i=0;i<=count;++i){refinement=std::max(refinement,norm(ref[i].x-fine[i].x));ref_energy=std::max(ref_energy,std::abs(energy(ref[i])-E0));}
        require(refinement<1e-8&&ref_energy<1e-8,"work-midpoint reference not converged");
        for(int sweeps:{1,5,20}){
            State s=c.initial;double maxE=0,maxL=0,sumX=0,maxX=0,path=0,maxC=0,maxR=0,maxKin=0,maxForce=0,maxKW=0,maxUW=0,maxBI=0,maxAI=0;
            double old_phase=phase(s,c.tangential),old_ref_phase=old_phase,phase_difference=0,max_phase=0;
            int completed=0;std::string status="completed";
            for(int i=0;i<=count;++i){
                if(i){
                    const Step update=step(s,h,sweeps);
                    if(!update.ok){status=update.status;break;}
                    path+=norm(update.next.x-s.x);s=update.next;completed=i;
                    const double p=phase(s,c.tangential),pr=phase(ref[i],c.tangential);
                    phase_difference+=wrapped(p-old_phase)-wrapped(pr-old_ref_phase);old_phase=p;old_ref_phase=pr;
                    maxC=std::max(maxC,update.rc);maxR=std::max(maxR,update.rx);maxKin=std::max(maxKin,update.kinematic);
                    maxForce=std::max(maxForce,update.force_defect);maxKW=std::max(maxKW,update.kinetic_work);
                    maxUW=std::max(maxUW,update.potential_work);maxBI=std::max(maxBI,update.balance_identity);maxAI=std::max(maxAI,update.angular_identity);
                }
                const double error=norm(s.x-ref[i].x);sumX+=error*error;maxX=std::max(maxX,error);
                maxE=std::max(maxE,std::abs(energy(s)-E0));maxL=std::max(maxL,std::abs(angular(s)-L0));max_phase=std::max(max_phase,std::abs(phase_difference));
                traces<<c.name<<",work_midpoint,"<<h<<','<<sweeps<<','<<duration<<','<<i*h<<','<<s.x.x<<','<<s.x.y<<','<<s.v.x<<','<<s.v.y<<','
                      <<energy(s)<<','<<angular(s)<<','<<phase_difference<<','<<ref[i].x.x<<','<<ref[i].x.y<<','<<ref[i].v.x<<','<<ref[i].v.y<<','<<energy(ref[i])<<'\n';
            }
            summary<<c.name<<",work_midpoint,"<<h<<','<<sweeps<<','<<duration<<','<<completed<<','<<count<<','<<status<<','<<E0<<','<<L0<<','<<energy(s)<<','
                   <<maxE/E0<<','<<maxL<<','<<std::sqrt(sumX/(completed+1))<<','<<maxX<<','<<phase_difference<<','<<max_phase<<','<<path<<','
                   <<maxC<<','<<maxR<<','<<maxKin<<','<<maxForce<<','<<maxKW<<','<<maxUW<<','<<maxBI<<','<<maxAI<<','<<refinement<<','<<ref_energy<<'\n';
            ++runs;if(status!="completed")++failures;
            if(sweeps>=5)worst_energy_converged=std::max(worst_energy_converged,maxE/E0);
            worst_angular=std::max(worst_angular,maxL);worst_work_identity=std::max(worst_work_identity,maxUW);worst_balance_identity=std::max(worst_balance_identity,maxBI);worst_angular_identity=std::max(worst_angular_identity,maxAI);
        }
    }
    validation<<"converged_energy_relative,"<<worst_energy_converged<<",1e-9\n"
              <<"angular_absolute,"<<worst_angular<<",1e-9\n"
              <<"potential_work_absolute,"<<worst_work_identity<<",1e-11\n"
              <<"energy_defect_decomposition,"<<worst_balance_identity<<",1e-10\n"
              <<"angular_defect_decomposition,"<<worst_angular_identity<<",1e-10\n";
    require(failures==0,"one or more work-midpoint runs failed");
    require(worst_energy_converged<1e-9,"converged work-midpoint energy error");
    require(worst_angular<1e-9,"work-midpoint angular error");
    require(worst_work_identity<1e-11&&worst_balance_identity<1e-10,"work balance identity error");
    require(worst_angular_identity<1e-10,"angular balance identity error");
    std::cout<<"runs="<<runs<<" failures="<<failures<<" converged_energy_relative="<<worst_energy_converged
             <<" angular_absolute="<<worst_angular<<" potential_work="<<worst_work_identity<<" balance_identity="<<worst_balance_identity<<'\n';
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}
