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

// Independent single-particle, anchored, radial spring probe.
// No energy multiplier, damping, clipping, or velocity rescaling.
namespace {
constexpr double mass=1.0, stiffness=100.0, rest_length=1.0;
constexpr double omega=10.0, pi=3.14159265358979323846;
struct V {double x=0,y=0;};
V operator+(V a,V b){return {a.x+b.x,a.y+b.y};}
V operator-(V a,V b){return {a.x-b.x,a.y-b.y};}
V operator*(V a,double s){return {a.x*s,a.y*s};}
V operator/(V a,double s){return a*(1/s);}
double dot(V a,V b){return a.x*b.x+a.y*b.y;}
double cross(V a,V b){return a.x*b.y-a.y*b.x;}
double norm(V a){return std::hypot(a.x,a.y);}
struct State {V x,v;};
struct Case {std::string name; State initial; bool tangential;};
const std::array<Case,3> cases{{
    {"spring_rest",{{1.25,0},{0,0}},false},
    {"spring_radial",{{1.25,0},{-.7,0}},false},
    {"spring_tangent",{{1.25,0},{0,1}},true},
}};
double potential(V x){const double c=norm(x)-rest_length;return .5*stiffness*c*c;}
double energy(State s){return .5*mass*dot(s.v,s.v)+potential(s.x);}
double angular(State s){return mass*cross(s.x,s.v);}
V force(V x){const double r=norm(x);if(!(r>0))throw std::runtime_error("spring direction undefined");return x*(-stiffness*(r-rest_length)/r);}
double wrapped(double x){return std::remainder(x,2*pi);}
double phase(State s,bool tangent){
    if(tangent)return std::atan2(s.x.y,s.x.x);
    return std::atan2(-dot(s.v,s.x/norm(s.x))/omega,norm(s.x)-rest_length);
}
State reference_step(State s,double h){
    const V a1=force(s.x)/mass, q1=s.v;
    const V a2=force(s.x+q1*(h/2))/mass, q2=s.v+a1*(h/2);
    const V a3=force(s.x+q2*(h/2))/mass, q3=s.v+a2*(h/2);
    const V a4=force(s.x+q3*h)/mass, q4=s.v+a3*h;
    return {s.x+(q1+q2*2+q3*2+q4)*(h/6),s.v+(a1+a2*2+a3*2+a4)*(h/6)};
}
std::vector<State> reference(State initial,double h,int steps,int subdivisions){
    std::vector<State> result{initial};
    for(int i=0;i<steps;++i){
        for(int j=0;j<subdivisions;++j)initial=reference_step(initial,h/subdivisions);
        result.push_back(initial);
    }
    return result;
}
struct Step {
    State next;
    double constraint_residual=0,stationarity=0,lambda=0,chord_energy=0;
    double kinematic_residual=0,force_residual=0;
};
Step step(State old,double h,int sweeps,bool midpoint){
    const double tau=midpoint?h/2:h;
    const V predictor=old.x+old.v*tau;
    V z=predictor;
    double lambda=0;
    const double alpha_tilde=1/(stiffness*tau*tau);
    for(int i=0;i<sweeps;++i){
        const double r=norm(z);
        if(!(r>0))throw std::runtime_error("undefined constraint normal");
        const V n=z/r;
        const double dlambda=-(r-rest_length+alpha_tilde*lambda)/(1/mass+alpha_tilde);
        lambda+=dlambda;
        z=z+n*(dlambda/mass);
    }
    Step out;
    out.next=midpoint?State{z*2-old.x,(z-old.x)*(4/h)-old.v}:State{z,(z-old.x)/h};
    out.constraint_residual=norm(z)-rest_length+alpha_tilde*lambda;
    const V normal=z/norm(z);
    out.stationarity=norm((z-predictor)*mass-normal*lambda)/(tau*tau);
    out.lambda=lambda;
    out.chord_energy=energy({out.next.x,(out.next.x-old.x)/h});
    if(midpoint){
        out.kinematic_residual=norm(out.next.x-old.x-(out.next.v+old.v)*(h/2));
        out.force_residual=norm((out.next.v-old.v)*mass-force(z)*h);
    }
    return out;
}
void require(bool ok,const char* msg){if(!ok)throw std::runtime_error(msg);}
} // namespace

int main(int argc,char** argv){
try{
    const std::filesystem::path out=argc>1?argv[1]:"results/spring_repair/midpoint";
    std::filesystem::create_directories(out);
    std::ofstream summary(out/"summary.csv"),traces(out/"traces.csv"),checks(out/"checks.csv");
    require(summary&&traces&&checks,"cannot open outputs");
    summary<<std::setprecision(17);traces<<std::setprecision(17);checks<<std::setprecision(17);
    summary<<"case,method,h,sweeps,steps,E0,L0,final_energy,max_relative_energy_error,max_angular_error,"
              "rms_position_error,max_position_error,final_phase_error,max_phase_error,phase_kind,"
              "path_length,max_material_residual,max_stationarity,max_midpoint_kinematic_residual,"
              "max_midpoint_force_residual,max_chord_energy_error,reference_refinement_error,reference_energy_error\n";
    traces<<"case,method,h,sweeps,time,x,y,vx,vy,energy,angular_momentum,phase_error,ref_x,ref_y,ref_vx,ref_vy,ref_energy,chord_energy\n";
    checks<<"check,value,tolerance\n";
    double maximum_harmonic_energy=0,maximum_angular=0,maximum_single_update_residual=0;
    double maximum_reference=0,maximum_reference_energy=0,maximum_sweep_difference=0;
    double maximum_phase_formula_error=0;
    int runs=0;
    for(const Case& c:cases)for(double h:{1./60,1./240}){
        const int count=static_cast<int>(std::lround(10/h));
        const auto ref=reference(c.initial,h,count,64);
        const auto finer=reference(c.initial,h,count,128);
        const double E0=energy(c.initial),L0=angular(c.initial);
        double refinement=0,ref_energy_error=0;
        for(int i=0;i<=count;++i){
            refinement=std::max(refinement,norm(ref[i].x-finer[i].x));
            ref_energy_error=std::max(ref_energy_error,std::abs(energy(ref[i])-E0));
        }
        maximum_reference=std::max(maximum_reference,refinement);
        maximum_reference_energy=std::max(maximum_reference_energy,ref_energy_error);
        for(bool midpoint:{false,true}){
            std::vector<State> one_sweep;
            for(int sweeps:{1,5,20}){
                State s=c.initial;
                double max_dE=0,max_dL=0,sum_position_error=0,max_position_error=0;
                double old_phase=phase(s,c.tangential),old_ref_phase=old_phase,phase_difference=0,max_phase_error=0;
                double path=0,max_rC=0,max_rx=0,max_kinematic=0,max_force=0,max_chord_dE=0;
                const std::string method=midpoint?"midpoint_xpbd":"standard_xpbd";
                for(int i=0;i<=count;++i){
                    double chord_energy=E0;
                    if(i){
                        const Step z=step(s,h,sweeps,midpoint);
                        path+=norm(z.next.x-s.x);
                        s=z.next;chord_energy=z.chord_energy;
                        const double current_phase=phase(s,c.tangential),current_ref_phase=phase(ref[i],c.tangential);
                        phase_difference+=wrapped(current_phase-old_phase)-wrapped(current_ref_phase-old_ref_phase);
                        old_phase=current_phase;old_ref_phase=current_ref_phase;
                        max_rC=std::max(max_rC,std::abs(z.constraint_residual));
                        max_rx=std::max(max_rx,z.stationarity);
                        max_kinematic=std::max(max_kinematic,z.kinematic_residual);
                        max_force=std::max(max_force,z.force_residual);
                        max_chord_dE=std::max(max_chord_dE,std::abs(chord_energy-E0));
                    }
                    const double position_error=norm(s.x-ref[i].x);
                    max_position_error=std::max(max_position_error,position_error);
                    sum_position_error+=position_error*position_error;
                    max_dE=std::max(max_dE,std::abs(energy(s)-E0));
                    max_dL=std::max(max_dL,std::abs(angular(s)-L0));
                    max_phase_error=std::max(max_phase_error,std::abs(phase_difference));
                    traces<<c.name<<','<<method<<','<<h<<','<<sweeps<<','<<i*h<<','<<s.x.x<<','<<s.x.y<<','
                        <<s.v.x<<','<<s.v.y<<','<<energy(s)<<','<<angular(s)<<','<<phase_difference<<','
                        <<ref[i].x.x<<','<<ref[i].x.y<<','<<ref[i].v.x<<','<<ref[i].v.y<<','<<energy(ref[i])<<','<<chord_energy<<'\n';
                    if(sweeps==1)one_sweep.push_back(s);
                    else maximum_sweep_difference=std::max(maximum_sweep_difference,norm(s.x-one_sweep[i].x));
                }
                if(midpoint){
                    maximum_angular=std::max(maximum_angular,max_dL);
                    if(!c.tangential)maximum_harmonic_energy=std::max(maximum_harmonic_energy,max_dE/E0);
                    if(sweeps==1)maximum_single_update_residual=std::max(maximum_single_update_residual,max_rx);
                }
                if(!c.tangential){
                    const double discrete_phase=midpoint?2*std::atan(omega*h/2):std::atan(omega*h);
                    maximum_phase_formula_error=std::max(maximum_phase_formula_error,
                        std::abs(phase_difference-(count*discrete_phase-omega*count*h)));
                }
                summary<<c.name<<','<<method<<','<<h<<','<<sweeps<<','<<count<<','<<E0<<','<<L0<<','<<energy(s)<<','
                    <<max_dE/E0<<','<<max_dL<<','<<std::sqrt(sum_position_error/(count+1))<<','<<max_position_error<<','
                    <<phase_difference<<','<<max_phase_error<<','<<(c.tangential?"orbital_angle":"harmonic_phase")<<','<<path<<','
                    <<max_rC<<','<<max_rx<<','<<max_kinematic<<','<<max_force<<','<<max_chord_dE<<','<<refinement<<','<<ref_energy_error<<'\n';
                ++runs;
            }
        }
    }
    const double h=1./60;
    const Step first=step(cases[0].initial,h,1,true);
    const double z0=omega*h/2,q0=.25;
    const double exact_q=q0*(1-z0*z0)/(1+z0*z0);
    const double exact_v=-q0*omega*2*z0/(1+z0*z0);
    const double one_step_error=std::max(std::abs(first.next.x.x-(1+exact_q)),std::abs(first.next.v.x-exact_v));
    const Step baseline=step(cases[0].initial,h,1,false);
    const double baseline_denominator=1+omega*omega*h*h;
    const double baseline_error=std::max(std::abs(baseline.next.x.x-(1+q0/baseline_denominator)),
        std::abs(baseline.next.v.x+q0*omega*omega*h/baseline_denominator));
    checks<<"harmonic_baseline_one_step_formula,"<<baseline_error<<",1e-12\n"
          <<"harmonic_midpoint_one_step_formula,"<<one_step_error<<",1e-12\n"
          <<"harmonic_energy_relative_10s,"<<maximum_harmonic_energy<<",1e-10\n"
          <<"midpoint_angular_absolute_10s,"<<maximum_angular<<",1e-10\n"
          <<"midpoint_one_sweep_stationarity,"<<maximum_single_update_residual<<",1e-8\n"
          <<"reference_refinement_position,"<<maximum_reference<<",1e-8\n"
          <<"reference_energy_absolute,"<<maximum_reference_energy<<",1e-8\n"
          <<"one_vs_multiple_sweep_position,"<<maximum_sweep_difference<<",1e-8\n"
          <<"harmonic_10s_discrete_phase_formula,"<<maximum_phase_formula_error<<",1e-8\n";
    require(baseline_error<1e-12,"baseline one-step formula mismatch");
    require(one_step_error<1e-12,"midpoint one-step formula mismatch");
    require(maximum_harmonic_energy<1e-10,"harmonic midpoint energy conservation failed");
    require(maximum_angular<1e-10,"midpoint central-force angular conservation failed");
    require(maximum_reference<1e-8&&maximum_reference_energy<1e-8,"reference not converged");
    require(maximum_sweep_difference<1e-8,"unexpected multi-sweep change for one radial spring");
    require(maximum_single_update_residual<1e-8,"one-sweep midpoint stationarity failure");
    require(maximum_phase_formula_error<1e-8,"harmonic discrete phase formula mismatch");
    std::cout<<"runs="<<runs<<" output="<<out<<"\n"
             <<"harmonic_energy_relative="<<maximum_harmonic_energy<<" angular_absolute="<<maximum_angular
             <<" one_update_stationarity="<<maximum_single_update_residual<<" sweep_difference="<<maximum_sweep_difference<<"\n";
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}
