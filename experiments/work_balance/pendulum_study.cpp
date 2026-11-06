#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

struct V {
    double x = 0.0, y = 0.0;
    V operator+(V b) const { return {x+b.x,y+b.y}; }
    V operator-(V b) const { return {x-b.x,y-b.y}; }
    V operator*(double a) const { return {a*x,a*y}; }
    V operator/(double a) const { return {x/a,y/a}; }
    double dot(V b) const { return x*b.x+y*b.y; }
    double norm2() const { return dot(*this); }
    double norm() const { return std::sqrt(norm2()); }
    double cross(V b) const { return x*b.y-y*b.x; }
};

struct State { V q,v; };
struct Angle { double theta,omega; };
State state_from_angle(Angle a);
struct Metrics {
    double work_identity = 0, radial_work = 0, gravity_mismatch = 0;
    double projection_loss = 0, impulse_residual = 0;
    double normalized_denominator = 1;
    bool negative_denominator = false;
};

double sinc(double value) {
    if(std::abs(value)<1e-4) {
        const double square=value*value;
        return 1-square/6+square*square/120-square*square*square/5040;
    }
    return std::sin(value)/value;
}

double sinc_derivative(double value) {
    if(std::abs(value)<1e-4) {
        const double square=value*value;
        return value*(-1.0/3+square/30-square*square/840);
    }
    return (value*std::cos(value)-std::sin(value))/(value*value);
}

State angular_discrete_gradient(State old,double h,V gravity) {
    const double theta=std::atan2(old.q.x,-old.q.y);
    const double omega=old.q.cross(old.v);
    const double acceleration=-gravity.y;
    if(acceleration==0) return state_from_angle({theta+h*omega,omega});
    const double bound=.5*h*h*acceleration;
    double low=h*omega-bound,high=h*omega+bound;
    double delta=std::clamp(h*omega-bound*std::sin(theta),low,high);
    // The bracket follows from |sin(midpoint)*sinc(delta/2)|<=1.
    // At the three tested small steps F is strictly monotone.
    for(int iteration=0;iteration<80;++iteration) {
        const double t=delta/2,mid=theta+t;
        const double value=delta-h*omega+bound*std::sin(mid)*sinc(t);
        if(std::abs(value)<2e-16*(1+std::abs(h*omega))) break;
        if(value>0) high=delta; else low=delta;
        const double derivative=1+bound*.5*(std::cos(mid)*sinc(t)
                                     +std::sin(mid)*sinc_derivative(t));
        const double candidate=delta-value/derivative;
        const double next=(candidate>low && candidate<high)?candidate:.5*(low+high);
        if(next==delta || high-low<4e-16*(1+std::abs(delta))) { delta=next; break; }
        delta=next;
    }
    // Momentum equation avoids cancellation in 2*delta/h - omega at turning points.
    const double next_omega=omega-h*acceleration*std::sin(theta+delta/2)*sinc(delta/2);
    return state_from_angle({theta+delta,next_omega});
}

double energy(State a, V gravity) { return .5*a.v.norm2()-gravity.dot(a.q); }
State state_from_angle(Angle a) {
    return {{std::sin(a.theta),-std::cos(a.theta)},
            {a.omega*std::cos(a.theta),a.omega*std::sin(a.theta)}};
}

Angle rk4(Angle a,double h,double gravity) {
    const auto rhs = [gravity](Angle p) { return Angle{p.omega,-gravity*std::sin(p.theta)}; };
    const auto add = [](Angle a,Angle b,double s) { return Angle{a.theta+s*b.theta,a.omega+s*b.omega}; };
    const Angle k1=rhs(a),k2=rhs(add(a,k1,h/2)),k3=rhs(add(a,k2,h/2)),k4=rhs(add(a,k3,h));
    return {a.theta+h*(k1.theta+2*k2.theta+2*k3.theta+k4.theta)/6,
            a.omega+h*(k1.omega+2*k2.omega+2*k3.omega+k4.omega)/6};
}

Angle reference_step(Angle a,double h,double gravity,double target_step=1.0/20000) {
    if (gravity == 0) return {a.theta+h*a.omega,a.omega};
    const int n=static_cast<int>(std::ceil(h/target_step));
    for(int i=0;i<n;++i) a=rk4(a,h/n,gravity);
    return a;
}

State advance(State old,double h,V gravity,const std::string& method,Metrics& metrics) {
    State next;
    if(method=="angular_dg") {
        next=angular_discrete_gradient(old,h,gravity);
    } else if(method=="xpbd") {
        const V predicted=old.q+old.v*h+gravity*(h*h);
        if(predicted.norm()<1e-14) throw std::runtime_error("collapsed XPBD prediction");
        next.q=predicted/predicted.norm();
        next.v=(next.q-old.q)/h;
    } else {
        // Unit radius: the exact midpoint length equation reduces to a reflection.
        const V s=old.q*2+old.v*h+gravity*(.5*h*h);
        const double q_dot_s=old.q.dot(s);
        metrics.normalized_denominator=std::abs(q_dot_s)/(old.q.norm()*s.norm());
        metrics.negative_denominator=q_dot_s<0;
        if(s.norm()<1e-14) throw std::runtime_error("midpoint continuum: s=0");
        if(std::abs(q_dot_s)<1e-13*old.q.norm()*s.norm())
            throw std::runtime_error("midpoint antipodal point: no finite radial multiplier");
        next.q=s*(2*q_dot_s/s.norm2())-old.q;
        next.v=(next.q-old.q)*(2/h)-old.v;
        const V average_q=(next.q+old.q)*.5;
        const V impulse=next.v-old.v-gravity*h;
        const double coefficient=impulse.dot(average_q)/average_q.norm2();
        metrics.impulse_residual=(impulse-average_q*coefficient).norm();
        if(method=="midpoint_project") {
            const double radial=next.q.dot(next.v)/next.q.norm2();
            const V projected=next.v-next.q*radial;
            metrics.projection_loss=.5*radial*radial*next.q.norm2();
            next.v=projected;
        }
    }
    const V average_v=(old.v+next.v)*.5;
    const V total_impulse=next.v-old.v;
    const V radial_impulse=total_impulse-gravity*h;
    const double kinetic_change=.5*(next.v.norm2()-old.v.norm2());
    metrics.work_identity=std::abs(kinetic_change-average_v.dot(total_impulse));
    metrics.radial_work=average_v.dot(radial_impulse);
    metrics.gravity_mismatch=gravity.dot(average_v*h-(next.q-old.q));
    // NaNs can otherwise disappear in std::max-based summaries and evade >-tolerance checks.
    const double finite_values[] = {next.q.x,next.q.y,next.v.x,next.v.y,
        next.q.norm2(),next.v.norm2(),energy(next,gravity),metrics.work_identity,
        metrics.radial_work,metrics.gravity_mismatch,metrics.projection_loss,
        metrics.impulse_residual,metrics.normalized_denominator};
    for(double value:finite_values)
        if(!std::isfinite(value))
            throw std::runtime_error("non-finite pendulum state or diagnostic");
    return next;
}

void study(const std::string& scenario,const std::string& method,double h,
           std::ofstream& summary,std::ofstream& trajectory) {
    const double acceleration=scenario=="gravity"?9.81:0;
    const V gravity{0,-acceleration};
    Angle reference{.65,acceleration==0?2.0:0.0};
    State model=state_from_angle(reference);
    const double initial_energy=energy(model,gravity);
    const double initial_angular=model.q.cross(model.v);
    const int steps=static_cast<int>(std::round(10/h));
    double max_energy=0,max_length=0,max_radial=0,max_position=0,max_velocity=0,max_angular=0;
    double work_identity=0,radial_work=0,gravity_mismatch=0,projection_loss=0,impulse_residual=0;
    double min_denominator=1;
    double max_midpoint_kinematics=0,max_tangency_identity=0;
    int negative_denominators=0,completed=0;
    std::string status="ok";
    for(int i=1;i<=steps;++i) {
        Metrics metrics;
        const State old=model;
        try { model=advance(model,h,gravity,method,metrics); }
        catch(const std::exception& error) { status=error.what(); break; }
        completed=i;
        reference=reference_step(reference,h,acceleration);
        const State truth=state_from_angle(reference);
        const double de=energy(model,gravity)-initial_energy;
        const double radial=model.q.dot(model.v)/model.q.norm();
        const double angular=model.q.cross(model.v);
        max_energy=std::max(max_energy,std::abs(de));
        max_length=std::max(max_length,std::abs(model.q.norm()-1));
        max_radial=std::max(max_radial,std::abs(radial));
        max_position=std::max(max_position,(model.q-truth.q).norm());
        max_velocity=std::max(max_velocity,(model.v-truth.v).norm());
        max_angular=std::max(max_angular,std::abs(angular-initial_angular));
        work_identity=std::max(work_identity,metrics.work_identity);
        radial_work+=metrics.radial_work;
        gravity_mismatch+=metrics.gravity_mismatch;
        projection_loss+=metrics.projection_loss;
        impulse_residual=std::max(impulse_residual,metrics.impulse_residual);
        min_denominator=std::min(min_denominator,metrics.normalized_denominator);
        negative_denominators+=metrics.negative_denominator;
        max_midpoint_kinematics=std::max(max_midpoint_kinematics,
            (model.q-old.q-(model.v+old.v)*(h/2)).norm());
        const double kinetic_change=.5*(model.v.norm2()-old.v.norm2());
        max_tangency_identity=std::max(max_tangency_identity,
            std::abs(model.q.dot(model.v)+old.q.dot(old.v)-h*kinetic_change/2));
        trajectory<<scenario<<','<<method<<','<<h<<','<<i*h<<','
            <<model.q.x<<','<<model.q.y<<','<<model.v.x<<','<<model.v.y<<','
            <<de<<','<<radial<<','<<angular<<','<<truth.q.x<<','<<truth.q.y<<','
            <<metrics.radial_work<<','<<metrics.gravity_mismatch<<','<<metrics.projection_loss<<'\n';
    }
    const State truth=state_from_angle(reference);
    if(h<=1.0/15) {
        if(completed!=steps || max_length>1e-9 || work_identity>1e-10)
            throw std::runtime_error("small-step run failed length/work checks");
        if((method=="midpoint" || method=="angular_dg") && max_energy>1e-8)
            throw std::runtime_error("conservative reference exceeded energy tolerance");
        if(method=="angular_dg" && max_radial>1e-12)
            throw std::runtime_error("angular reference violated endpoint tangency");
        if(method=="midpoint_project" &&
            std::abs(energy(model,gravity)-initial_energy+projection_loss)>1e-8)
            throw std::runtime_error("projected midpoint loss differs from removed radial kinetic energy");
        if(scenario=="free_rotation" && method=="midpoint" && max_angular>1e-8)
            throw std::runtime_error("free midpoint rotation changed angular momentum");
    }
    summary<<scenario<<','<<method<<','<<h<<','<<completed*h<<','<<status<<','
        <<energy(model,gravity)-initial_energy<<','<<max_energy<<','<<max_length<<','
        <<max_radial<<','<<(model.q-truth.q).norm()<<','<<max_position<<','<<max_velocity<<','
        <<max_angular<<','<<work_identity<<','<<radial_work<<','<<gravity_mismatch<<','
        <<projection_loss<<','<<impulse_residual<<','<<min_denominator<<','<<negative_denominators<<','
        <<max_midpoint_kinematics<<','<<max_tangency_identity<<'\n';
    std::cout<<scenario<<' '<<method<<" h="<<h<<" dE="<<energy(model,gravity)-initial_energy
        <<" radial="<<max_radial<<" qerr="<<max_position<<" status="<<status<<'\n';
}

void singularity_probes(std::ofstream& out) {
    out<<"initial_speed,h_factor,h,status,qx,qy,vx,vy,radial_multiplier,normalized_denominator\n";
    const V gravity{0,-9.81};
    const double critical=2/std::sqrt(9.81);
    for(double speed:{0.0,1.0}) for(double factor:{.5,.99,.999999,1.0,1.000001,1.01,2.0}) {
        const State old{{0,1},{speed,0}};
        const double h=critical*factor;
        Metrics metrics;
        try {
            const State next=advance(old,h,gravity,"midpoint",metrics);
            const V average_q=(next.q+old.q)*.5;
            const double coefficient=(next.v-old.v-gravity*h).dot(average_q)/average_q.norm2();
            out<<speed<<','<<factor<<','<<h<<",ok,"<<next.q.x<<','<<next.q.y<<','
                <<next.v.x<<','<<next.v.y<<','<<coefficient<<','<<metrics.normalized_denominator<<'\n';
        } catch(const std::exception& error) {
            out<<speed<<','<<factor<<','<<h<<','<<error.what()<<",nan,nan,nan,nan,nan,nan\n";
        }
    }
}

int main() {
    std::ofstream summary("pendulum_summary.csv"),trajectory("pendulum_trajectories.csv"),
        singularities("pendulum_singularities.csv");
    if(!summary || !trajectory || !singularities)
        throw std::runtime_error("cannot open pendulum result files");
    summary<<std::setprecision(17); trajectory<<std::setprecision(17);
    singularities<<std::setprecision(17); std::cout<<std::setprecision(10);
    summary<<"scenario,method,h,duration,status,final_energy_drift,max_energy_drift,max_length_error,"
        "max_radial_velocity,final_position_error,max_position_error,max_velocity_error,"
        "max_angular_change,max_kinetic_work_identity_error,sum_constraint_impulse_work,"
        "sum_gravity_work_mismatch,sum_projection_loss,max_radial_impulse_residual,"
        "minimum_normalized_denominator,negative_denominator_steps,"
        "max_cartesian_midpoint_drift_residual,max_midpoint_tangency_identity_residual\n";
    trajectory<<"scenario,method,h,time,qx,qy,vx,vy,energy_drift,radial_velocity,angular_momentum,"
        "reference_qx,reference_qy,constraint_impulse_work,gravity_work_mismatch,projection_loss\n";
    for(const std::string scenario:{"free_rotation","gravity"})
        for(double h:{1.0/240,1.0/60,1.0/15,.2,.5,1.0,2.0}) {
            for(const std::string method:{"xpbd","midpoint","midpoint_project"})
                study(scenario,method,h,summary,trajectory);
            if(h<=1.0/15) study(scenario,"angular_dg",h,summary,trajectory);
        }
    singularity_probes(singularities);
    const Angle a=reference_step({.65,0},10,9.81,1.0/20000);
    const Angle b=reference_step({.65,0},10,9.81,1.0/40000);
    if(std::abs(a.theta-b.theta)>1e-10 || std::abs(a.omega-b.omega)>1e-10)
        throw std::runtime_error("RK4 reference step halving did not converge");
    std::cout<<"RK4 reference step halving: angle="<<std::abs(a.theta-b.theta)
        <<" omega="<<std::abs(a.omega-b.omega)<<'\n';
    summary.flush(); trajectory.flush(); singularities.flush();
    if(!summary || !trajectory || !singularities)
        throw std::runtime_error("failed writing pendulum result files");
}
