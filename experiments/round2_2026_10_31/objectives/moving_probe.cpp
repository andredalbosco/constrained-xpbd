#include "objectives.hpp"
#include <fstream>
#include <iomanip>
#include <iostream>
using namespace objectives2;
int main(){
 std::ofstream out("moving_angles.csv"),rod("rod_kinematics.csv");out<<std::setprecision(17);rod<<std::setprecision(17);
 out<<"variant,speed,angle_degrees,hz,h_acceleration_over_speed,converged,prediction_valid_nonzero_speed,acceleration_error_to_physics,acceleration_error_to_predicted_branch,predicted_branch_error_to_physics,velocity_change_x,velocity_change_y,balance,stationarity\n";
 rod<<"case,variant,hz,converged,position_error,velocity_error,energy_error,length_error,tangency,midpoint_position_defect,reference_midpoint_position_defect,kinetic_change,balance,stationarity\n";
 for(auto&q:variants())if(q.id=="original_midpoint_b1"||q.id=="original_midpoint_h2"||q.id=="original_endpoint_b1"||q.id=="quarter_b1"||q.id=="midpotential_b1")for(double speed:{0.,1e-4,.1,.8,2.})for(double degrees:{0.,30.,60.,90.,120.,180.})for(int hz:{80,320,1280,5120}){
  auto c=make_case("radial_orbit");double t=degrees*std::acos(-1.)/180;c.initial.v<<-speed*std::cos(t),speed*std::sin(t);double h=1./hz,E=energy(c.model,c.initial),mass=c.model.mass[0];V a=force(c.model,c.initial.x)/mass;auto r=step(c.model,c.initial,h,q,E,96);V b=(r.state.v-c.initial.v)/h,predicted=a;bool nonzero=speed>0;
  if(q.objective=="original"&&q.kin=="soft_midpoint"&&nonzero)predicted=2*a-a.dot(c.initial.v)/c.initial.v.squaredNorm()*c.initial.v;
  out<<q.id<<','<<speed<<','<<degrees<<','<<hz<<','<<(nonzero?h*a.norm()/speed:-1)<<','<<r.converged<<','<<nonzero<<','<<(b-a).norm()<<','<<(b-predicted).norm()<<','<<(predicted-a).norm()<<','<<b[0]<<','<<b[1]<<','<<r.balance<<','<<r.stationarity<<'\n';
 }
 for(auto name:{"rod_release","rod_kick"})for(auto&q:variants())if(q.id=="original_midpoint_b1"||q.id=="quarter_b1"||q.id=="midpotential_b1"||q.id=="projection_verlet_E"||q.id=="ls_dg_E")for(int hz:{40,80,160,320,640,1280}){
  auto c=make_case(name);double h=1./hz,E=energy(c.model,c.initial)-minimum_energy(c.model);auto r=step(c.model,c.initial,h,q,E,96);auto ref=reference_step(c.model,c.initial,h,256);auto defect=[&](const State&s){return (s.x-c.initial.x-.5*h*(s.v+c.initial.v)).norm();};double K0=.5*c.model.mass[0]*c.initial.v.squaredNorm(),K=.5*c.model.mass[0]*r.state.v.squaredNorm();
  rod<<name<<','<<q.id<<','<<hz<<','<<r.converged<<','<<(r.state.x-ref.x).norm()<<','<<(r.state.v-ref.v).norm()<<','<<std::abs(energy(c.model,r.state)-energy(c.model,c.initial))/E<<','<<std::abs(r.state.x.norm()-1)<<','<<std::abs(r.state.x.dot(r.state.v)/r.state.x.norm())<<','<<defect(r.state)<<','<<defect(ref)<<','<<K-K0<<','<<r.balance<<','<<r.stationarity<<'\n';
 }
 std::cout<<"600 velocity-angle probes and60 rod probes complete\n";
}
