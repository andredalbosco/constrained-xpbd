#define main parent_study_main
#include "study.cpp"
#undef main

struct Holdout {
 double rms=0,max_energy=0,max_L=0,last_L=0,last_angle=0,max_material=0,path=0,min_r=1e30,max_r=0;
 long steps=0,capped=0,skipped=0;bool complete=true;
};
Holdout angular_neutral_run(const Case& c,int substeps,const std::vector<State>& ref) {
 Holdout s;const double h=1./(60*substeps),m=c.model.mass,E0=energy(c.model,c.initial),L0=angular(c.initial,m);State state=c.initial;
 double angle=std::atan2(state.x.y,state.x.x),refangle=angle,last=angle,reflast=angle;
 for(size_t i=1;i<ref.size();++i) {
  const Context ctx=context(c.model,state,h,0,E0);const Iterate z{ctx.predictor};const Rows r=evaluate(ctx,z);
  if(!r.valid){s.complete=false;break;}
  const double dl0=-r.rC/r.A;const V2 dx0=r.n*(dl0/m);const double H0=r.rE+dot(r.g,dx0);
  const V2 nold=state.x/norm(state.x);const double gn=dot(r.g,nold);V2 extra{};
  if(std::abs(H0)<=energy_tolerance(ctx,r)||std::abs(gn)<=gradient_tolerance(r)){++s.skipped;}
  else {
   // Diagnostic only: project the energy direction onto zero-angular-impulse
   // displacements. This changes the original multiplier system.
   extra=nold*(-H0/gn);double blend=.5;
   const double cap=.1*gradient_policies::motion_scale(ctx),size=norm(extra);
   if(blend*size>cap){blend=cap/size;++s.capped;}
   extra=extra*blend;
  }
  const State next{ctx.predictor+dx0+extra,(ctx.predictor+dx0+extra-state.x)/h};
  const double E=energy(c.model,next);
  if(!finite(next.x)||!finite(next.v)||!std::isfinite(E)||std::abs(E-E0)>1e12*scale(c)){s.complete=false;break;}
  ++s.steps;s.rms+=dot(next.x-ref[i].x,next.x-ref[i].x);s.max_energy=std::max(s.max_energy,std::abs(E-E0));
  s.last_L=angular(next,m);s.max_L=std::max(s.max_L,std::abs(s.last_L-L0));
  const double a=std::atan2(next.x.y,next.x.x),ar=std::atan2(ref[i].x.y,ref[i].x.x);
  angle+=std::remainder(a-last,2*std::acos(-1.));refangle+=std::remainder(ar-reflast,2*std::acos(-1.));last=a;reflast=ar;s.last_angle=angle-refangle;
  const double radius=norm(next.x);s.min_r=std::min(s.min_r,radius);s.max_r=std::max(s.max_r,radius);s.path+=norm(next.x-state.x);
  s.max_material=std::max(s.max_material,std::abs(radius-c.model.length+r.at*dl0));state=next;
 }
 s.rms=std::sqrt(s.rms/std::max(1L,s.steps));return s;
}
int main() {
 std::cout<<std::setprecision(17);
 std::cout<<"case,speed_fraction,method,substeps,h,rms_position_error,max_energy_percent,max_angular_error,final_L,final_angle_error,min_radius,max_radius,cap_count,skip_count,max_material_residual,status\n";
 Model m;m.mass=1;m.length=1;m.alpha=.01;m.acceleration={};const double radius=1.1,vcircle=std::sqrt(radius*(radius-m.length)/(m.mass*m.alpha));
 for(double fraction:{.1,.8/vcircle,.5,.75,1.,1.25}) {
  Case c{"eccentricity_family",m,{{radius,0},{0,fraction*vcircle}}};
  for(int substeps:{4,8,16,32,128}) {
   const double h=1./(60*substeps);const auto ref=reference(c,h,600*substeps,16);
   for(const auto& cfg:configurations()) {
    if(cfg.name!="baseline"&&cfg.name!="LE_repair_cap")continue;
    const Stats s=run(c,cfg,substeps,10,ref,nullptr);
    std::cout<<c.name<<','<<fraction<<','<<cfg.name<<','<<substeps<<','<<h<<','<<s.rms_position<<','<<100*s.energy_error/scale(c)<<','<<s.max_torque_error<<",nan,nan,"<<s.min_length<<','<<s.max_length<<','<<s.diagnostic.capped<<','<<s.diagnostic.repair_skips<<','<<s.max_rC<<','<<s.status<<'\n';
   }
   const Holdout s=angular_neutral_run(c,substeps,ref);
   std::cout<<c.name<<','<<fraction<<",angular_neutral_extra,"<<substeps<<','<<h<<','<<s.rms<<','<<100*s.max_energy/scale(c)<<','<<s.max_L<<','<<s.last_L<<','<<s.last_angle<<','<<s.min_r<<','<<s.max_r<<','<<s.capped<<','<<s.skipped<<','<<s.max_material<<','<<(s.complete?"completed":"aborted")<<'\n';
  }
 }
}
