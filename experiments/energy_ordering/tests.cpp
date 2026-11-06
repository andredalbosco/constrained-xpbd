#include "solver.hpp"
#include "scenes.hpp"
#include <Eigen/Dense>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>

using namespace energy_ordering;

namespace {
void require(bool condition,const std::string& message) {
    if(!condition)throw std::runtime_error(message);
}
void near(double value,double expected,double tolerance,const std::string& message) {
    if(!std::isfinite(value)||!std::isfinite(expected)||std::abs(value-expected)>tolerance)
        throw std::runtime_error(message+"; actual="+std::to_string(value)+", expected="+std::to_string(expected));
}
void near(V2 value,V2 expected,double tolerance,const std::string& message) {
    near(norm(value-expected),0,tolerance,message);
}
Eigen::Vector2d eigen(V2 value){return {value.x,value.y};}
V2 vec(const Eigen::Vector2d& value){return {value.x(),value.y()};}

// Independent physical definitions, without calling evaluate/energy/elastic.
double physical_energy(const Model& m,const State& s) {
    const double extension=std::hypot(s.x.x,s.x.y)-m.length;
    return .5*m.mass*(s.v.x*s.v.x+s.v.y*s.v.y)
        +(m.pendulum?0:extension*extension/(2*m.alpha))
        -m.mass*(m.acceleration.x*s.x.x+m.acceleration.y*s.x.y);
}
double physical_defect(const Context& c,V2 x) {
    return physical_energy(c.model,{x,(x-c.old.x)/c.h})-c.target;
}

struct IndependentRows {
    Eigen::Vector2d n,g,R;
    double rc,re,at,bt;
};
IndependentRows independent_rows(const Context& c,const Iterate& z) {
    const auto& m=c.model;
    const Eigen::Vector2d x=eigen(z.x),old=eigen(c.old.x),a=eigen(m.acceleration);
    const double radius=x.norm(),extension=radius-m.length;
    IndependentRows r;
    r.n=x/radius;
    r.g=m.mass*(x-old)/(c.h*c.h)-m.mass*a;
    if(!m.pendulum)r.g+=r.n*(extension/m.alpha);
    r.at=m.pendulum?0:m.alpha/(c.h*c.h);r.bt=c.beta/(c.h*c.h);
    r.rc=extension+r.at*z.lambda;r.re=physical_defect(c,z.x)+r.bt*z.eta;
    const Eigen::Vector2d predictor=old+c.h*eigen(c.old.v)+c.h*c.h*a;
    r.R=m.mass*(x-predictor)-r.n*z.lambda-r.g*z.eta;
    return r;
}

// Solve the uneliminated four-dimensional mixed equations independently.
Eigen::Vector4d mixed_oracle(const Context& c,const Iterate& z,bool retain_rx) {
    const auto r=independent_rows(c,z);
    Eigen::Matrix4d matrix=Eigen::Matrix4d::Zero();
    matrix.topLeftCorner<2,2>()=c.model.mass*Eigen::Matrix2d::Identity();
    matrix.block<2,1>(0,2)=-r.n;matrix.block<2,1>(0,3)=-r.g;
    matrix.block<1,2>(2,0)=r.n.transpose();matrix.block<1,2>(3,0)=r.g.transpose();
    matrix(2,2)=r.at;matrix(3,3)=r.bt;
    Eigen::Vector4d rhs;
    rhs.head<2>()=retain_rx?(-r.R).eval():Eigen::Vector2d::Zero().eval();
    rhs(2)=-r.rc;rhs(3)=-r.re;
    const auto lu=matrix.fullPivLu();require(lu.rank()==4,"mixed oracle must have full rank");
    const Eigen::Vector4d answer=lu.solve(rhs);
    near((matrix*answer-rhs).norm(),0,2e-12,"independent mixed solve residual");
    return answer;
}

void row_derivatives() {
    for(bool pendulum:{false,true}) {
        Model m;m.pendulum=pendulum;m.mass=1.7;m.length=1.2;m.alpha=.035;m.acceleration={.2,-.3};
        const State old{{1.1,.7},{.15,-.12}};
        const Context c=context(m,old,.08,.4,physical_energy(m,old)-.03);
        const Iterate z{{1.06,.76},-.014,.003};const Rows r=evaluate(c,z);
        const auto expected=independent_rows(c,z);
        near(r.C,std::hypot(z.x.x,z.x.y)-m.length,1e-14,"material definition");
        near(r.H,physical_defect(c,z.x),1e-13,"physical energy definition");
        near(r.R,vec(expected.R),1e-14,"position residual definition");
        near(r.rC,expected.rc,1e-13,"material compliance scaling");
        near(r.rE,expected.re,1e-13,"energy compliance scaling");
        for(V2 direction:{V2{.6,.8},V2{-.8,.6}}) {
            const double eps=1e-6;
            const double derivative=(physical_defect(c,z.x+direction*eps)-physical_defect(c,z.x-direction*eps))/(2*eps);
            near(dot(r.g,direction),derivative,2e-8,"energy gradient finite difference");
        }
    }
}

void baseline_spring() {
    double max_x_error=0,max_v_error=0;
    for(double h:{1./240,1./60}) {
        Model m;State old;cxpbd::Model production=cxpbd::make_spring();
        cxpbd::Settings settings;settings.dt=h;settings.iterations=20;settings.mode=cxpbd::Mode::XPBD;
        Options options;options.method=Method::Baseline;
        const int steps=static_cast<int>(std::lround(10/h));
        for(int k=1;k<=steps;++k) {
            const Context c=context(m,old,h,0,physical_energy(m,old));
            Iterate z{c.predictor};Counters count;
            for(int i=0;i<20;++i)require(sweep(c,z,options,count),"baseline spring update failed");
            old={z.x,(z.x-c.old.x)/h};
            const double amplitude=.25*std::pow(1+100*h*h,-.5*k),phase=k*std::atan(10*h);
            max_x_error=std::max(max_x_error,std::abs(old.x.x-(1+amplitude*std::cos(phase))));
            max_v_error=std::max(max_v_error,std::abs(old.v.x+10*amplitude*std::sin(phase)));
            cxpbd::step(production,settings);
            near(old.x,{production.particles[1].x.x(),production.particles[1].x.y()},2e-12,"production XPBD spring position parity");
            near(old.v,{production.particles[1].v.x(),production.particles[1].v.y()},2e-11,"production XPBD spring velocity parity");
        }
        const double energy_ratio=physical_energy(m,old)/3.125;
        near(energy_ratio,std::pow(1+100*h*h,-steps),2e-12,"analytic backward-Euler spring energy decay");
    }
    near(max_x_error,0,2e-12,"analytic backward-Euler position");
    near(max_v_error,0,2e-11,"analytic backward-Euler velocity");
}

void baseline_pendulum() {
    Model m;m.pendulum=true;m.alpha=0;m.acceleration={0,-9.81};
    State old{{std::sin(.65),-std::cos(.65)},{}};
    cxpbd::Model production=cxpbd::make_pendulum();
    cxpbd::Settings settings;settings.dt=1./240;settings.iterations=20;
    Options options;options.method=Method::Baseline;
    for(int k=0;k<2400;++k) {
        const Context c=context(m,old,settings.dt,0,physical_energy(m,old));
        const V2 expected=c.predictor/norm(c.predictor);
        Iterate z{c.predictor};Counters count;
        for(int i=0;i<20;++i)require(sweep(c,z,options,count),"baseline pendulum update failed");
        old={z.x,(z.x-c.old.x)/c.h};cxpbd::step(production,settings);
        near(z.x,expected,1e-14,"analytic radial projection");
        near(norm(z.x),1,1e-14,"baseline pendulum length");
        near(old.x,{production.particles[1].x.x(),production.particles[1].x.y()},2e-10,"production XPBD pendulum position parity");
        near(old.v,{production.particles[1].v.x(),production.particles[1].v.y()},2e-9,"production XPBD pendulum velocity parity");
    }
}

Context exact_rest_spring(double beta=0) {
    Model m;m.alpha=1./64;const State old{{1.25,0},{}};
    return context(m,old,1./8,beta,2);
}
void resting_joint_freeze() {
    const Context c=exact_rest_spring();
    for(bool retain:{false,true}) {
        Options options;options.method=Method::Joint;options.retain_rx=retain;
        Iterate z{c.predictor};Counters count;
        require(sweep(c,z,options,count),"resting joint failed");
        near(z.x,c.old.x,0,"joint rest position must remain frozen");
        near(z.lambda,-.25,0,"joint rest material multiplier");
        near(z.eta,1./64,0,"joint rest energy multiplier");
        const Rows r=evaluate(c,z);
        near(r.H,0,0,"frozen rest preserves physical energy");
        near(r.rC,0,0,"frozen rest material residual");near(r.rE,0,0,"frozen rest energy residual");
        near(norm(r.R),0,0,"frozen rest position residual");
        require(count.singular==0&&count.incompatible==0,"spring rest joint block must be nonsingular");
        // A mathematically consistent root can still suppress the physical oscillation.
        require(std::abs((c.old.x.x-1)*64)>1,"frozen spring must retain a nonzero physical force");
    }
}

void stationary_energy_gradient() {
    const Context c=exact_rest_spring();Iterate z{c.predictor};Counters count;
    require(material_update(c,z,count),"resting spring material update failed");
    near(z.x,{1.125,0},0,"exact spring backward-Euler state");
    Rows r=evaluate(c,z);near(norm(r.g),0,0,"material-first energy gradient must vanish exactly");
    near(r.H,-1,0,"material-first energy defect must remain nonzero");
    const Iterate before=z;
    require(energy_update(c,z,count),"singular row should be explicitly skipped");
    near(z.x,before.x,0,"zero-gradient hard row changed position");near(z.eta,0,0,"zero-gradient hard row changed dual");
    require(count.singular==1&&count.incompatible==1,"incompatible zero-gradient energy row not reported");

    // Softening makes the same scalar row finite but does not restore physical E.
    const Context soft=exact_rest_spring(1./64);Counters soft_count;z=before;
    require(energy_update(soft,z,soft_count),"soft zero-gradient row failed");
    near(z.x,before.x,0,"soft zero-gradient row moved the spring");near(z.eta,1,0,"soft energy multiplier");
    r=evaluate(soft,z);
    near(r.rE,0,0,"soft energy equation residual");near(r.rC,0,0,"soft material equation residual");near(norm(r.R),0,0,"soft position equation residual");
    near(r.H,-1,0,"soft equation must not be confused with hard conservation");
    require(soft_count.singular==0&&soft_count.incompatible==0,"finite compliance incorrectly reported singular");
}

void frozen_linear_system() {
    Model m;m.mass=1.7;m.length=1.2;m.alpha=.035;m.acceleration={.2,-.3};
    const State old{{1.1,.7},{.15,-.12}};
    const Context c=context(m,old,.08,.4,physical_energy(m,old)-.03);
    const Iterate initial{{1.06,.76},-.014,.003};
    const auto r=independent_rows(c,initial);
    const double a=1/m.mass+r.at,b=r.n.dot(r.g)/m.mass,d=r.g.squaredNorm()/m.mass+r.bt;
    const double contraction=b*b/(a*d);
    require(contraction>1e-4&&contraction<.9,"frozen fixture must have nontrivial but convergent coupling");
    for(bool retain:{false,true}) {
        const Eigen::Vector4d expected=mixed_oracle(c,initial,retain);
        for(Method method:{Method::Joint,Method::Frozen}) {
            Options options;options.method=method;options.retain_rx=retain;options.inner=200;
            Iterate actual=initial;Counters count;
            require(frozen_update(c,actual,options,count),"frozen linearized update failed");
            near(actual.x-initial.x,vec(expected.head<2>()),3e-13,"mixed-oracle displacement");
            near(actual.lambda-initial.lambda,expected(2),3e-13,"mixed-oracle lambda increment");
            near(actual.eta-initial.eta,expected(3),3e-13,"mixed-oracle eta increment");
            require(count.evaluations==1,"frozen inner iterations unexpectedly relinearized");
        }

        // One GS pass solves the lower triangular Schur matrix. Jacobi would
        // discard its off-diagonal entry and use a different energy RHS.
        Eigen::Matrix2d schur;schur<<a,b,b,d;
        const Eigen::Vector2d base=retain?(-r.R/m.mass).eval():Eigen::Vector2d::Zero().eval();
        const Eigen::Vector2d rhs(-r.rc-r.n.dot(base),-r.re-r.g.dot(base));
        const Eigen::Vector2d gs=schur.triangularView<Eigen::Lower>().solve(rhs);
        const Eigen::Vector2d jacobi(rhs(0)/a,rhs(1)/d);
        require((gs-jacobi).norm()>1e-6,"one-pass fixture must distinguish GS from Jacobi");
        Options one;one.method=Method::Frozen;one.inner=1;one.retain_rx=retain;
        Iterate result=initial;Counters count;
        require(frozen_update(c,result,one,count),"one-pass frozen GS failed");
        near(result.lambda-initial.lambda,gs(0),3e-13,"one-pass material increment");
        near(result.eta-initial.eta,gs(1),3e-13,"one-pass energy must use updated material RHS");
        near(result.x-initial.x,vec(base+(r.n*gs(0)+r.g*gs(1))/m.mass),3e-13,"one-pass base applied exactly once");
    }

    Options frozen;frozen.method=Method::Frozen;frozen.inner=1;
    Options immediate;immediate.method=Method::Immediate;
    Iterate a_iterate=initial,b_iterate=initial;Counters a_count,b_count;
    require(sweep(c,a_iterate,frozen,a_count)&&sweep(c,b_iterate,immediate,b_count),"schedule distinction fixture failed");
    require(norm(a_iterate.x-b_iterate.x)>1e-6,"immediate relinearization unexpectedly matched frozen one-pass GS");
}

void resting_pendulum_energy_singularity() {
    Model m;m.pendulum=true;m.alpha=0;m.acceleration={0,-9.81};
    const State old{{std::sin(.65),-std::cos(.65)},{}};
    const Context c=context(m,old,1./240,0,physical_energy(m,old));
    const Iterate predictor{c.predictor};const Rows r=evaluate(c,predictor);
    near(norm(r.g),0,0,"gravity-centered resting predictor gradient must be exactly zero");
    near(r.H,-.5*c.h*c.h*9.81*9.81,2e-14,"resting gravity predictor energy defect");
    require(r.H<0,"resting predictor must have a nonzero negative energy defect");
    Iterate result=predictor;Counters count;
    require(energy_update(c,result,count),"zero predictor energy row should be explicitly skipped");
    require(count.singular==1&&count.incompatible==1,"zero predictor energy row incompatibility not reported");
    near(result.x,predictor.x,0,"zero predictor gradient generated a spurious giant step");
    near(result.eta,0,0,"zero predictor gradient changed energy multiplier");
    Options joint;joint.method=Method::Joint;Counters joint_count;
    result=predictor;
    require(frozen_update(c,result,joint,joint_count),"incompatible joint predictor should be explicitly skipped");
    require(joint_count.singular==1&&joint_count.incompatible==1,"joint predictor incompatibility not reported");
    near(result.x,predictor.x,0,"incompatible joint predictor changed position");
}

void hard_pendulum_position_invariance() {
    Model m;m.pendulum=true;m.alpha=0;m.mass=1.3;m.acceleration={0,-9.81};
    const State old{{std::sin(.65),-std::cos(.65)},{1.2*std::cos(.65),1.2*std::sin(.65)}};
    const Context c=context(m,old,.03,0,physical_energy(m,old)+.01);
    const Iterate initial{c.predictor+V2{.0004,-.0002},.0007,-.0003};
    const auto r=independent_rows(c,initial);
    Eigen::Matrix2d constraint;
    constraint.row(0)=r.n.transpose();constraint.row(1)=r.g.transpose();
    require(std::abs(constraint.determinant())>.1*r.g.norm(),"hard pendulum fixture must have independent rows");
    const Eigen::Vector2d expected=constraint.fullPivLu().solve(Eigen::Vector2d(-r.rc,-r.re));
    Iterate result[2];int index=0;
    for(bool retain:{false,true}) {
        result[index]=initial;Options options;options.method=Method::Joint;options.retain_rx=retain;Counters count;
        require(frozen_update(c,result[index],options,count),"hard pendulum block failed");
        const auto& z=result[index];
        near(z.x-initial.x,vec(expected),1e-13,"fullrank hard constraint displacement must be independent of top residual");
        const Eigen::Vector2d frozen_R=r.R+m.mass*eigen(z.x-initial.x)
            -r.n*(z.lambda-initial.lambda)-r.g*(z.eta-initial.eta);
        near(frozen_R.norm(),retain?0:r.R.norm(),2e-13,"frozen top-equation correction");
        const Eigen::Vector4d oracle=mixed_oracle(c,initial,retain);
        near(z.lambda-initial.lambda,oracle(2),2e-13,"hard pendulum dual from mixed oracle");
        near(z.eta-initial.eta,oracle(3),2e-13,"hard pendulum energy dual from mixed oracle");
        ++index;
    }
    require(std::abs(result[0].lambda-result[1].lambda)+std::abs(result[0].eta-result[1].eta)>1e-6,
            "retaining top residual must change this fixture's duals");
}

void weak_soft_feedback() {
    Model m;const State old{{1.25,0},{.4,0}};
    const Context c=context(m,old,.02,1e20,physical_energy(m,old));
    Iterate baseline{c.predictor},joint=baseline;Counters baseline_count,joint_count;
    require(material_update(c,baseline,baseline_count),"weak-feedback baseline failed");
    Options options;options.method=Method::Joint;
    require(frozen_update(c,joint,options,joint_count),"weak-feedback joint failed");
    near(joint.x,baseline.x,1e-13,"weak energy feedback must approach XPBD");
    near(joint.lambda,baseline.lambda,1e-13,"weak energy feedback material multiplier");
    require(std::abs(physical_defect(c,joint.x))>1e-4,"weak-feedback fixture must not accidentally preserve physical energy");
}

void hard_pendulum_feasibility() {
    Model m;m.pendulum=true;m.alpha=0;
    State old{{1,0},{0,1}};
    require(pendulum_feasible(context(m,old,1,0,.5)),"free rotating pendulum with hv<2L must be feasible");
    require(!pendulum_feasible(context(m,old,3,0,.5)),"free rotating pendulum with hv>2L must be infeasible");
    require(!pendulum_feasible(context(m,old,.1,0,-1)),"negative kinetic target must be infeasible without gravity");
}
}

int main() {
    try {
        struct Test { const char* name; void (*run)(); };
        for(const Test test:{Test{"physical rows and gradients",row_derivatives},
            {"baseline spring analytic and production parity",baseline_spring},
            {"baseline pendulum radial and production parity",baseline_pendulum},
            {"hard joint exact resting-spring freeze",resting_joint_freeze},
            {"hard singular and soft nonconserving stationary energy rows",stationary_energy_gradient},
            {"frozen GS and joint versus independent mixed system",frozen_linear_system},
            {"centered resting-pendulum predictor singularity",resting_pendulum_energy_singularity},
            {"hard fullrank pendulum position invariance and dual change",hard_pendulum_position_invariance},
            {"weak soft feedback limit",weak_soft_feedback},
            {"hard pendulum feasibility",hard_pendulum_feasibility}}) {
            test.run();std::cout<<"[pass] "<<test.name<<'\n';
        }
    }catch(const std::exception& error){std::cerr<<"[fail] "<<error.what()<<'\n';return 1;}
}
