#include "xpbd.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

using namespace cxpbd;

namespace {

void require(bool condition,const std::string& message) {
    if(!condition)throw std::runtime_error(message);
}
void near(double actual,double expected,double tolerance,const std::string& message) {
    if(!std::isfinite(actual)||std::abs(actual-expected)>tolerance)
        throw std::runtime_error(message+": actual="+std::to_string(actual)+", expected="+std::to_string(expected));
}
double U(const Model& model,std::size_t j) {
    const auto& s=model.springs[j];
    if(s.compliance==0)return 0;
    const double c=(model.particles[s.a].x-model.particles[s.b].x).norm()-s.rest_length;
    return c*c/(2*s.compliance);
}
double E(const Model& model) {
    double sum=0;
    for(const auto& p:model.particles)sum+=p.mass*(.5*p.v.squaredNorm()-model.gravity.dot(p.x));
    for(std::size_t j=0;j<model.springs.size();++j)sum+=U(model,j);
    return sum;
}
double R(const Model& model) {
    double sum=0;
    for(const auto& p:model.particles)sum+=p.reservoir;
    for(const auto& s:model.springs)sum+=s.reservoir;
    return sum;
}
double distance(const Model& a,const Model& b) {
    double sum=0;
    for(std::size_t i=0;i<a.particles.size();++i)sum+=(a.particles[i].x-b.particles[i].x).squaredNorm();
    return std::sqrt(sum);
}

struct Frozen {
    std::vector<Vec3> edge,particle;
    explicit Frozen(const Model& old):edge(old.springs.size(),Vec3::Zero()),particle(old.particles.size(),Vec3::Zero()) {
        for(std::size_t j=0;j<old.springs.size();++j) {
            const auto& s=old.springs[j];
            const Vec3 d=old.particles[s.a].x-old.particles[s.b].x;
            if(d.norm()<=1e-14*s.rest_length)continue;
            const double magnitude=old.previous_dt>0
                ?s.previous_lambda/(old.previous_dt*old.previous_dt)
                :(s.compliance>0?-(d.norm()-s.rest_length)/s.compliance:0);
            edge[j]=magnitude*d/d.norm();
            particle[s.a]+=edge[j];particle[s.b]-=edge[j];
        }
    }
};

// Independent accepted-state accounting; no production energy/balance helpers.
void audit(const Model& old,const Model& now,double h,const Diagnostics& diagnostic) {
    Frozen force(old);
    double scale=1+std::abs(E(old))+std::abs(R(old));
    for(std::size_t i=0;i<old.particles.size();++i) {
        const auto& a=old.particles[i];const auto& b=now.particles[i];
        require(b.x.allFinite()&&b.v.allFinite()&&std::isfinite(b.reservoir),"nonfinite particle");
        const Vec3 dx=b.x-a.x;
        near((b.v-(a.pinned?Vec3::Zero().eval():(dx/h).eval())).norm(),0,1e-11,"velocity reconstruction");
        if(a.pinned)near(dx.norm(),0,1e-14,"pin moved");
        const double defect=.5*a.mass*(b.v.squaredNorm()-a.v.squaredNorm())
            -(force.particle[i]+a.mass*old.gravity).dot(dx);
        near(defect+b.reservoir-a.reservoir,0,2e-10*scale,"particle work account");
    }
    for(std::size_t j=0;j<old.springs.size();++j) {
        const auto& s=old.springs[j];
        const Vec3 dx=(now.particles[s.a].x-old.particles[s.a].x)
                     -(now.particles[s.b].x-old.particles[s.b].x);
        near(U(now,j)-U(old,j)+force.edge[j].dot(dx)+now.springs[j].reservoir-s.reservoir,
             0,2e-10*scale,"spring work account");
        require(std::isfinite(now.springs[j].lambda)&&std::isfinite(now.springs[j].reservoir),"nonfinite spring");
    }
    near(E(now)+R(now),E(old)+R(old),2e-10*scale,"augmented work account");
    near(diagnostic.mechanical_energy,E(now),2e-10*scale,"reported mechanical energy");
    near(diagnostic.reservoir_energy,R(now),2e-10*scale,"reported reservoir energy");
    near(diagnostic.energy_balance_error,E(now)+R(now)-E(old)-R(old),
         2e-10*scale,"reported augmented account change");
    const double values[]={diagnostic.mechanical_energy,diagnostic.reservoir_energy,diagnostic.reservoir_norm,
        diagnostic.energy_balance_error,diagnostic.max_balance_residual,diagnostic.max_feedback_residual,
        diagnostic.max_material_residual,diagnostic.position_residual_norm};
    for(double value:values)require(std::isfinite(value),"nonfinite diagnostic");
}

Model pair_model() {
    Model model;model.gravity=Vec3(.1,-.7,.05);
    Particle a,b;a.mass=1.2;b.mass=2.4;
    a.x=Vec3(-.6,.1,0);b.x=Vec3(.55,-.2,.1);
    a.v=Vec3(.2,.1,-.05);b.v=Vec3(-.08,.16,.03);
    a.reservoir=.08;b.reservoir=-.05;
    model.particles={a,b};
    Spring s;s.a=0;s.b=1;s.rest_length=1.05;s.compliance=.04;s.reservoir=.03;
    model.springs={s};return model;
}
Model shared_model() {
    Model model=pair_model();
    Particle c;c.mass=.8;c.x=Vec3(.1,1.05,-.12);c.v=Vec3(.12,-.08,.06);c.reservoir=.045;
    model.particles.push_back(c);
    Spring s;s.a=0;s.b=2;s.rest_length=1.12;s.compliance=.06;s.reservoir=-.025;
    model.springs.push_back(s);
    model.previous_dt=.03;
    model.springs[0].previous_lambda=-.001;
    model.springs[1].previous_lambda=.0008;
    return model;
}
Settings settings() {
    Settings s;s.mode=Mode::Active;s.active_order=ActiveOrder::Sequential;
    s.dt=.02;s.iterations=1;s.kinetic_reservoir_compliance=.7;s.elastic_reservoir_compliance=.3;
    return s;
}

enum class Wrong { None, ResetMu, CurrentEdgeForce, GlobalKinetic };
struct Oracle {
    Model model;
    std::vector<double> mu,nu;
};

// Literal scalar row equations, without production update or line-search code.
// Only one prescribed sweep is needed to check the actual visitation schedule.
Oracle one_sweep(const Model& old,const Settings& s,Wrong wrong=Wrong::None) {
    Oracle out{old,std::vector<double>(old.particles.size(),0),std::vector<double>(old.springs.size(),0)};
    const Frozen force(old);
    const double h=s.dt,h2=h*h;
    const double beta=s.kinetic_reservoir_compliance/h2,gamma=s.elastic_reservoir_compliance/h2;
    for(auto& p:out.model.particles)if(!p.pinned)p.x+=h*p.v+h2*old.gravity;
    for(auto& spring:out.model.springs)spring.lambda=0;
    const auto kinetic=[&](int i,Vec3 f) {
        auto& p=out.model.particles[i];const auto& before=old.particles[i];
        const double w=p.pinned?0:1/p.mass;
        const Vec3 dx=p.x-before.x;
        const Vec3 gradient=p.mass*dx/h2-f-p.mass*old.gravity;
        const double b=.5*p.mass*(dx.squaredNorm()/h2-before.v.squaredNorm())
                       -(f+p.mass*old.gravity).dot(dx);
        if(wrong==Wrong::ResetMu)out.mu[i]=0;
        const double increment=(before.reservoir-b-beta*out.mu[i])/(w*gradient.squaredNorm()+beta);
        out.mu[i]+=increment;p.x+=w*gradient*increment;
    };
    std::vector<bool> visited(old.particles.size(),false);
    for(std::size_t j=0;j<old.springs.size();++j) {
        auto& spring=out.model.springs[j];auto& a=out.model.particles[spring.a];auto& b=out.model.particles[spring.b];
        const double wa=a.pinned?0:1/a.mass,wb=b.pinned?0:1/b.mass;
        Vec3 edge=a.x-b.x;Vec3 normal=edge/edge.norm();
        const double alpha=spring.compliance/h2;
        const double dl=(-(edge.norm()-spring.rest_length)-alpha*spring.lambda)/(wa+wb+alpha);
        spring.lambda+=dl;a.x+=wa*normal*dl;b.x-=wb*normal*dl;
        edge=a.x-b.x;normal=edge/edge.norm();
        Vec3 gradient=force.edge[j];
        if(spring.compliance>0)gradient+=(edge.norm()-spring.rest_length)*normal/spring.compliance;
        const Vec3 dx=(a.x-old.particles[spring.a].x)-(b.x-old.particles[spring.b].x);
        const double d=U(out.model,j)-U(old,j)+force.edge[j].dot(dx);
        const double dn=(old.springs[j].reservoir-d-gamma*out.nu[j])/((wa+wb)*gradient.squaredNorm()+gamma);
        out.nu[j]+=dn;a.x+=wa*gradient*dn;b.x-=wb*gradient*dn;
        if(wrong!=Wrong::GlobalKinetic) {
            kinetic(spring.a,wrong==Wrong::CurrentEdgeForce?force.edge[j]:force.particle[spring.a]);
            kinetic(spring.b,wrong==Wrong::CurrentEdgeForce?-force.edge[j]:force.particle[spring.b]);
            visited[spring.a]=visited[spring.b]=true;
        }
    }
    for(std::size_t i=0;i<old.particles.size();++i)
        if(!visited[i])kinetic(static_cast<int>(i),force.particle[i]);
    for(std::size_t i=0;i<old.particles.size();++i) {
        auto& p=out.model.particles[i];const auto& before=old.particles[i];
        const Vec3 dx=p.x-before.x;p.v=p.pinned?Vec3::Zero().eval():(dx/h).eval();
        const double balance=.5*p.mass*(p.v.squaredNorm()-before.v.squaredNorm())
            -(force.particle[i]+p.mass*old.gravity).dot(dx);
        p.reservoir=before.reservoir-balance;
    }
    for(std::size_t j=0;j<old.springs.size();++j) {
        auto& spring=out.model.springs[j];
        const Vec3 dx=(out.model.particles[spring.a].x-old.particles[spring.a].x)
                     -(out.model.particles[spring.b].x-old.particles[spring.b].x);
        spring.reservoir=old.springs[j].reservoir-(U(out.model,j)-U(old,j)+force.edge[j].dot(dx));
    }
    return out;
}

void compare_oracle(const Model& old,const Settings& s) {
    Model actual=old;const Oracle expected=one_sweep(old,s);
    const Diagnostics diagnostic=step(actual,s);
    require(diagnostic.backtracks==0&&diagnostic.rejected_updates==0,"oracle fixture took a damped/rejected step");
    near(distance(actual,expected.model),0,2e-13,"ordered scalar update positions");
    for(std::size_t i=0;i<actual.particles.size();++i) {
        near((actual.particles[i].v-expected.model.particles[i].v).norm(),0,2e-11,"ordered velocity");
        near(actual.particles[i].reservoir,expected.model.particles[i].reservoir,2e-11,"ordered particle account");
    }
    for(std::size_t j=0;j<actual.springs.size();++j) {
        near(actual.springs[j].lambda,expected.model.springs[j].lambda,2e-13,"material multiplier visitation");
        near(actual.springs[j].reservoir,expected.model.springs[j].reservoir,2e-11,"ordered spring account");
    }
    double feedback=0;
    const double h2=s.dt*s.dt;
    for(std::size_t i=0;i<actual.particles.size();++i)
        feedback=std::max(feedback,std::abs(expected.model.particles[i].reservoir
            -(s.kinetic_reservoir_compliance/h2)*expected.mu[i]));
    for(std::size_t j=0;j<actual.springs.size();++j)
        feedback=std::max(feedback,std::abs(expected.model.springs[j].reservoir
            -(s.elastic_reservoir_compliance/h2)*expected.nu[j]));
    near(diagnostic.max_feedback_residual,feedback,2e-11,"accepted-state dual feedback residual");
    audit(old,actual,s.dt,diagnostic);
}

void ordered_pair() {
    compare_oracle(pair_model(),settings());
}
void shared_mu_and_incident_work() {
    const Model old=shared_model();const Settings s=settings();
    const Model expected=one_sweep(old,s).model;
    for(Wrong wrong:{Wrong::ResetMu,Wrong::CurrentEdgeForce,Wrong::GlobalKinetic})
        require(distance(expected,one_sweep(old,s,wrong).model)>1e-7,"shared fixture does not distinguish an incorrect schedule");
    compare_oracle(old,s);
}
void isolated_and_fixed_particles() {
    Model isolated;Particle p;p.mass=1.7;p.x=Vec3(.3,-.2,.1);p.v=Vec3(.4,.1,-.2);p.reservoir=.08;
    isolated.particles={p};isolated.gravity=Vec3(.1,-.7,.05);
    compare_oracle(isolated,settings());
    Model one_pin=pair_model();one_pin.particles[0].pinned=true;one_pin.particles[0].v.setZero();
    compare_oracle(one_pin,settings());
    Model fixed=shared_model();
    for(auto& particle:fixed.particles){particle.pinned=true;particle.v.setZero();}
    for(ActiveOrder order:{ActiveOrder::Paired,ActiveOrder::PairedLocal,ActiveOrder::Sequential}) {
        Model actual=fixed;Settings s=settings();s.active_order=order;
        const Diagnostics d=step(actual,s);
        near(distance(actual,fixed),0,0,"fixed geometry moved");
        near(d.max_material_residual,0,2e-13,"fixed material multiplier");
        near(d.max_feedback_residual,0,2e-13,"fixed feedback multipliers");
        audit(fixed,actual,s.dt,d);
    }
}
void inactive_order_parity() {
    require(Settings{}.active_order==ActiveOrder::Paired,"default active order changed");
    for(Mode mode:{Mode::XPBD,Mode::Passive}) {
        Model reference=shared_model(),local=reference,sequential=reference;
        Settings s=settings();s.mode=mode;s.dt=1.0/240;s.iterations=4;
        for(int n=0;n<40;++n) {
            s.active_order=ActiveOrder::Paired;step(reference,s);
            s.active_order=ActiveOrder::PairedLocal;step(local,s);
            s.active_order=ActiveOrder::Sequential;step(sequential,s);
            for(std::size_t i=0;i<reference.particles.size();++i) {
                for(const Model* candidate:{&local,&sequential}) {
                    require((candidate->particles[i].x.array()==reference.particles[i].x.array()).all(),"inactive order changed positions");
                    require((candidate->particles[i].v.array()==reference.particles[i].v.array()).all(),"inactive order changed velocities");
                    require(candidate->particles[i].reservoir==reference.particles[i].reservoir,"inactive order changed account");
                }
            }
            for(std::size_t j=0;j<reference.springs.size();++j) {
                for(const Model* candidate:{&local,&sequential}) {
                    require(candidate->springs[j].lambda==reference.springs[j].lambda,"inactive order changed material multiplier");
                    require(candidate->springs[j].previous_lambda==reference.springs[j].previous_lambda,"inactive order changed force history");
                    require(candidate->springs[j].reservoir==reference.springs[j].reservoir,"inactive order changed spring account");
                }
            }
        }
    }
}
void weak_feedback_limit() {
    const Model initial=shared_model();
    Settings base=settings();base.mode=Mode::XPBD;base.iterations=3;
    Model xpbd=initial;step(xpbd,base);
    for(ActiveOrder order:{ActiveOrder::Paired,ActiveOrder::PairedLocal,ActiveOrder::Sequential}) {
        Model active=initial;Settings s=base;s.mode=Mode::Active;s.active_order=order;
        s.kinetic_reservoir_compliance=s.elastic_reservoir_compliance=1e16;
        const Diagnostics d=step(active,s);
        near(distance(active,xpbd),0,1e-12,"weak feedback did not approach XPBD");
        audit(initial,active,s.dt,d);
    }
}
void separate_compliances() {
    const Model initial=pair_model();Settings base=settings();
    Model xpbd=initial;Settings x=base;x.mode=Mode::XPBD;step(xpbd,x);
    for(bool kinetic:{false,true}) {
        Settings strong=base,weak=base;
        if(kinetic) {
            strong.elastic_reservoir_compliance=weak.elastic_reservoir_compliance=1e16;
            strong.kinetic_reservoir_compliance=.02;weak.kinetic_reservoir_compliance=2;
        } else {
            strong.kinetic_reservoir_compliance=weak.kinetic_reservoir_compliance=1e16;
            strong.elastic_reservoir_compliance=.02;weak.elastic_reservoir_compliance=2;
        }
        compare_oracle(initial,strong);compare_oracle(initial,weak);
        Model a=initial,b=initial;step(a,strong);step(b,weak);
        require(distance(a,xpbd)>1.1*distance(b,xpbd),"isolated compliance channel did not weaken its one-sweep correction");
        require(distance(a,b)>1e-8,"separate compliance channel was inactive");
    }
}
void local_work_accounts_all_orders() {
    for(ActiveOrder order:{ActiveOrder::Paired,ActiveOrder::PairedLocal,ActiveOrder::Sequential}) {
        Model model=shared_model();Settings s=settings();s.active_order=order;s.iterations=3;
        for(int n=0;n<40;++n) {
            const Model before=model;const Diagnostics d=step(model,s);audit(before,model,s.dt,d);
        }
    }
}

} // namespace

int main() {
    int failures=0;
    const auto run=[&](const char* name,void(*test)()) {
        try{test();std::cout<<"[pass] "<<name<<'\n';}
        catch(const std::exception& e){++failures;std::cerr<<"[fail] "<<name<<": "<<e.what()<<'\n';}
    };
    run("independent one-edge lambda/nu/mu_a/mu_b update",ordered_pair);
    run("shared mu and all-incident frozen work",shared_mu_and_incident_work);
    run("isolated particles and pinned multipliers",isolated_and_fixed_particles);
    run("baseline and passive order invariance",inactive_order_parity);
    run("weak feedback limit",weak_feedback_limit);
    run("separate one-sweep compliance responses",separate_compliances);
    run("independent local work accounts in all orders",local_work_accounts_all_orders);
    return failures?1:0;
}
