#include "frozen_schedules.hpp"
#include <Eigen/Dense>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

using namespace energy_ordering;

namespace {
void require(bool condition,const std::string& message) {
    if(!condition) throw std::runtime_error(message);
}
void near(double actual,double expected,double tolerance,const std::string& message) {
    if(!std::isfinite(actual)||!std::isfinite(expected)||std::abs(actual-expected)>tolerance)
        throw std::runtime_error(message+"; actual="+std::to_string(actual)+", expected="+std::to_string(expected));
}
void near(V2 actual,V2 expected,double tolerance,const std::string& message) {
    near(norm(actual-expected),0,tolerance,message);
}
Eigen::Vector2d eigen(V2 x) { return {x.x,x.y}; }
V2 vec(const Eigen::Vector2d& x) { return {x.x(),x.y()}; }

// These physical definitions do not use the implementation's row evaluator.
double physical_energy(const Model& m,const State& s) {
    const double extension=std::hypot(s.x.x,s.x.y)-m.length;
    const double elastic=m.pendulum?0:extension*extension/(2*m.alpha);
    return .5*m.mass*(s.v.x*s.v.x+s.v.y*s.v.y)+elastic
        -m.mass*(m.acceleration.x*s.x.x+m.acceleration.y*s.x.y);
}
struct PhysicalRows {
    Eigen::Vector2d n,g,R;
    double rC,rE,at,bt;
};
PhysicalRows physical_rows(const Context& c,const Iterate& z) {
    const auto& m=c.model;
    const auto x=eigen(z.x),old=eigen(c.old.x),acc=eigen(m.acceleration);
    PhysicalRows r;
    const double extension=x.norm()-m.length;
    r.n=x/x.norm();
    r.g=m.mass*(x-old)/(c.h*c.h)-m.mass*acc;
    if(!m.pendulum)r.g+=r.n*(extension/m.alpha);
    r.at=m.pendulum?0:m.alpha/(c.h*c.h);r.bt=c.beta/(c.h*c.h);
    r.rC=extension+r.at*z.lambda;
    r.rE=physical_energy(m,{z.x,(z.x-c.old.x)/c.h})-c.target+r.bt*z.eta;
    const Eigen::Vector2d predictor=old+c.h*eigen(c.old.v)+c.h*c.h*acc;
    r.R=m.mass*(x-predictor)-r.n*z.lambda-r.g*z.eta;
    return r;
}
Eigen::Vector4d mixed_oracle(const Context& c,const Iterate& z,bool retain) {
    const auto r=physical_rows(c,z);
    Eigen::Matrix4d matrix=Eigen::Matrix4d::Zero();
    matrix.topLeftCorner<2,2>()=c.model.mass*Eigen::Matrix2d::Identity();
    matrix.block<2,1>(0,2)=-r.n;matrix.block<2,1>(0,3)=-r.g;
    matrix.block<1,2>(2,0)=r.n.transpose();matrix.block<1,2>(3,0)=r.g.transpose();
    matrix(2,2)=r.at;matrix(3,3)=r.bt;
    Eigen::Vector4d rhs;
    rhs.head<2>()=retain?Eigen::Vector2d(-r.R):Eigen::Vector2d::Zero();
    rhs(2)=-r.rC;rhs(3)=-r.rE;
    const auto lu=matrix.fullPivLu();
    require(lu.rank()==4,"independent mixed fixture must have full rank");
    const Eigen::Vector4d answer=lu.solve(rhs);
    near((matrix*answer-rhs).norm(),0,2e-12,"independent mixed oracle residual");
    return answer;
}
struct Oracle {
    Iterate result;
    Eigen::Vector2d increments;
    double contraction;
};
Oracle coordinate_oracle(const Context& c,const Iterate& initial,bool retain,
                         const std::vector<int>& order,int passes) {
    const auto r=physical_rows(c,initial);
    Eigen::Matrix<double,2,2> columns;
    columns.col(0)=r.n;columns.col(1)=r.g;
    Eigen::Matrix2d matrix=columns.transpose()*columns/c.model.mass;
    matrix(0,0)+=r.at;matrix(1,1)+=r.bt;
    const Eigen::Vector2d base=retain?Eigen::Vector2d(-r.R/c.model.mass):Eigen::Vector2d::Zero();
    const Eigen::Vector2d rhs=Eigen::Vector2d(-r.rC,-r.rE)-columns.transpose()*base;
    Eigen::Vector2d increments=Eigen::Vector2d::Zero();
    // Independent implementation applies each *change* to the pending position.
    // Repeated visits must not add the entire accumulated multiplier again.
    Eigen::Vector2d position=eigen(initial.x)+base;
    for(int pass=0;pass<passes;++pass)for(int row:order) {
        const double correction=(rhs(row)-matrix.row(row).dot(increments))/matrix(row,row);
        increments(row)+=correction;
        position+=columns.col(row)*(correction/c.model.mass);
    }
    Oracle answer;
    answer.result={vec(position),initial.lambda+increments(0),initial.eta+increments(1)};
    answer.increments=increments;
    answer.contraction=matrix(0,1)*matrix(0,1)/(matrix(0,0)*matrix(1,1));
    near((position-(eigen(initial.x)+base+columns*increments/c.model.mass)).norm(),0,2e-13,
         "independent immediate and deferred formulas");
    return answer;
}

// Jacobi advances the complete dual vector from a single old residual vector.
// This independently checks that neither row sees the other's new increment.
Oracle jacobi_oracle(const Context& c,const Iterate& initial,bool retain,int passes) {
    const auto r=physical_rows(c,initial);
    Eigen::Matrix2d columns;
    columns.col(0)=r.n;columns.col(1)=r.g;
    Eigen::Matrix2d matrix=columns.transpose()*columns/c.model.mass;
    matrix(0,0)+=r.at;matrix(1,1)+=r.bt;
    const Eigen::Vector2d base=retain?Eigen::Vector2d(-r.R/c.model.mass):Eigen::Vector2d::Zero();
    const Eigen::Vector2d rhs=Eigen::Vector2d(-r.rC,-r.rE)-columns.transpose()*base;
    const Eigen::Matrix2d diagonal_inverse=matrix.diagonal().cwiseInverse().asDiagonal();
    Eigen::Vector2d increments=Eigen::Vector2d::Zero();
    for(int pass=0;pass<passes;++pass) {
        const Eigen::Vector2d correction=diagonal_inverse*(rhs-matrix*increments);
        increments+=correction;
    }
    Oracle answer;
    answer.result={vec(eigen(initial.x)+base+columns*increments/c.model.mass),
                   initial.lambda+increments(0),initial.eta+increments(1)};
    answer.increments=increments;
    answer.contraction=std::abs(matrix(0,1))/std::sqrt(matrix(0,0)*matrix(1,1));
    return answer;
}

Context fixture() {
    Model m;m.mass=1.7;m.length=1.2;m.alpha=.035;m.acceleration={.2,-.3};
    const State old{{1.1,.7},{.15,-.12}};
    return context(m,old,.08,.4,physical_energy(m,old)-.03);
}
const Iterate fixture_initial{{1.06,.76},-.014,.003};

struct Order { Schedule schedule; std::vector<int> rows; };
const std::vector<Order> orders{
    {Schedule::LE,{0,1}},{Schedule::EL,{1,0}},
    {Schedule::LEL,{0,1,0}},{Schedule::ELE,{1,0,1}}
};

void compare(Iterate actual,Iterate expected,double tolerance,const std::string& label) {
    near(actual.x,expected.x,tolerance,label+" position");
    near(actual.lambda,expected.lambda,tolerance,label+" material multiplier");
    near(actual.eta,expected.eta,tolerance,label+" energy multiplier");
}

void ordered_frozen_solves() {
    const Context c=fixture();
    for(bool retain:{false,true})for(int passes:{1,3,9})for(const auto& order:orders) {
        const auto oracle=coordinate_oracle(c,fixture_initial,retain,order.rows,passes);
        require(oracle.contraction>1e-4&&oracle.contraction<.9,
                "ordering fixture must have nontrivial convergent coupling");
        for(Application application:{Application::Deferred,Application::Immediate}) {
            Iterate actual=fixture_initial;Counters count;
            require(schedule_update(c,actual,order.schedule,retain,count,passes,application),
                    "regular ordered frozen solve failed");
            const std::string label=std::string(name(order.schedule))+" passes="+std::to_string(passes)
                +(retain?" retain rx":" drop rx");
            compare(actual,oracle.result,5e-13,label);
            require(count.evaluations==1,label+" reevaluated the frozen gradients");
            require(count.scalar_solves==static_cast<long>(order.rows.size())*passes,
                    label+" executed the wrong row order/count");
            require(count.singular==0&&count.incompatible==0,label+" falsely reported a singular row");
        }
    }
    // Reversal and the two symmetric visits must actually change this fixture.
    std::vector<Iterate> results;
    for(const auto& order:orders) {
        Iterate actual=fixture_initial;Counters count;
        require(schedule_update(c,actual,order.schedule,false,count),"order distinction failed");
        results.push_back(actual);
    }
    for(size_t i=0;i<results.size();++i)for(size_t j=i+1;j<results.size();++j)
        require(norm(results[i].x-results[j].x)>1e-7,"different orderings collapsed to the same one-pass result");
}

void equivalence_and_convergence() {
    const Context c=fixture();
    for(bool retain:{false,true}) {
        const Eigen::Vector4d full=mixed_oracle(c,fixture_initial,retain);
        const Iterate expected{fixture_initial.x+vec(full.head<2>()),
                               fixture_initial.lambda+full(2),fixture_initial.eta+full(3)};
        Iterate joint=fixture_initial;Counters joint_count;
        require(schedule_update(c,joint,Schedule::Joint,retain,joint_count),"joint frozen solve failed");
        compare(joint,expected,5e-13,"joint versus independent uneliminated system");
        for(const auto& order:orders)for(Application application:{Application::Deferred,Application::Immediate}) {
            Iterate actual=fixture_initial;Counters count;
            require(schedule_update(c,actual,order.schedule,retain,count,200,application),
                    "many-pass GS solve failed");
            compare(actual,expected,8e-13,std::string(name(order.schedule))+" converged versus full mixed system");
            require(count.evaluations==1,"inner passes relinearized");
        }
        // New LE must preserve the preexisting frozen-L/E implementation exactly
        // up to harmless different position summation order.
        for(int passes:{1,4,20}) {
            Iterate actual=fixture_initial,old=fixture_initial;Counters count,old_count;
            Options options;options.method=Method::Frozen;options.retain_rx=retain;options.inner=passes;
            require(schedule_update(c,actual,Schedule::LE,retain,count,passes),"new LE failed");
            require(frozen_update(c,old,options,old_count),"existing frozen LE failed");
            compare(actual,old,2e-13,"existing frozen LE parity");
        }
    }
    Iterate frozen=fixture_initial,relinearized=fixture_initial;Counters f_count,r_count;
    require(schedule_update(c,frozen,Schedule::LE,false,f_count),"frozen fixture failed");
    Options old;old.method=Method::Immediate;
    require(sweep(c,relinearized,old,r_count),"relinearized fixture failed");
    require(norm(frozen.x-relinearized.x)>1e-6,"fixture must expose reevaluation versus frozen ordering");
    require(f_count.evaluations==1&&r_count.evaluations==2,"reevaluation count distinction lost");
}

void frozen_equation_residuals() {
    const Context c=fixture();const auto r=physical_rows(c,fixture_initial);
    const auto mass=c.model.mass;
    require(r.R.norm()>.01,"stationarity fixture must exercise nonzero base displacement");
    for(bool retain:{false,true})for(const auto& order:orders) {
        Iterate actual=fixture_initial;Counters count;
        require(schedule_update(c,actual,order.schedule,retain,count,5,Application::Immediate),
                "retained residual update failed");
        const auto dx=eigen(actual.x-fixture_initial.x);
        const double dl=actual.lambda-fixture_initial.lambda,de=actual.eta-fixture_initial.eta;
        const Eigen::Vector2d top=mass*dx-r.n*dl-r.g*de;
        near((top-(retain?Eigen::Vector2d(-r.R):Eigen::Vector2d::Zero())).norm(),0,5e-13,
             "position stationarity base must be applied exactly once");
        const int last=order.rows.back();
        const double last_residual=last==0?r.rC+r.n.dot(dx)+r.at*dl:r.rE+r.g.dot(dx)+r.bt*de;
        near(last_residual,0,2e-12,"last visited frozen row must be satisfied");
    }
}

void jacobi_solves() {
    const Context c=fixture();const auto r=physical_rows(c,fixture_initial);
    for(bool retain:{false,true}) {
        for(int passes:{1,3,9}) {
            const auto oracle=jacobi_oracle(c,fixture_initial,retain,passes);
            require(oracle.contraction>1e-3&&oracle.contraction<1,
                    "Jacobi fixture must have nontrivial convergent coupling");
            for(Application app:{Application::Deferred,Application::Immediate}) {
                Iterate actual=fixture_initial;Counters count;
                require(schedule_update(c,actual,Schedule::Jacobi,retain,count,passes,app),
                        "frozen Jacobi failed");
                compare(actual,oracle.result,5e-13,"Jacobi must use both previous increments");
                require(count.evaluations==1,"Jacobi inner passes reevaluated gradients");
                require(count.scalar_solves==2*passes,"Jacobi scalar solve count");
                require(count.singular==0&&count.incompatible==0,"regular Jacobi classified singular");
                const Eigen::Vector2d dx=eigen(actual.x-fixture_initial.x);
                const double dl=actual.lambda-fixture_initial.lambda,de=actual.eta-fixture_initial.eta;
                const Eigen::Vector2d top=c.model.mass*dx-r.n*dl-r.g*de;
                near((top-(retain?Eigen::Vector2d(-r.R):Eigen::Vector2d::Zero())).norm(),0,5e-13,
                     "Jacobi must apply retained stationarity displacement exactly once");
            }
        }
        const Eigen::Vector4d full=mixed_oracle(c,fixture_initial,retain);
        const Iterate expected{fixture_initial.x+vec(full.head<2>()),
                               fixture_initial.lambda+full(2),fixture_initial.eta+full(3)};
        for(Application app:{Application::Deferred,Application::Immediate}) {
            Iterate actual=fixture_initial;Counters count;
            require(schedule_update(c,actual,Schedule::Jacobi,retain,count,200,app),
                    "converged Jacobi failed");
            compare(actual,expected,8e-13,"Jacobi convergence to uneliminated mixed-system oracle");
            require(count.evaluations==1,"converged Jacobi relinearized");
        }
        Iterate jacobi=fixture_initial;Counters count;
        require(schedule_update(c,jacobi,Schedule::Jacobi,retain,count),"Jacobi first pass failed");
        for(const auto& order:orders) {
            Iterate gs=fixture_initial;Counters gs_count;
            require(schedule_update(c,gs,order.schedule,retain,gs_count),"GS distinction fixture failed");
            require(norm(jacobi.x-gs.x)>1e-7,"Jacobi fixture must differ from each sequential ordering");
        }
    }
}

void jacobi_singular_rows() {
    Model m;m.pendulum=true;m.alpha=0;m.acceleration={0,-9.81};
    const State old{{std::sin(.65),-std::cos(.65)},{}};
    const Context hard=context(m,old,1./60,0,physical_energy(m,old));
    const Iterate start{hard.predictor};const Rows r=evaluate(hard,start);
    require(norm(r.g)==0&&r.H<0,"Jacobi singular fixture must have an incompatible zero energy gradient");
    for(bool retain:{false,true})for(int passes:{1,3,9})
        for(Application app:{Application::Deferred,Application::Immediate}) {
        Iterate actual=start;Counters count;
        require(schedule_update(hard,actual,Schedule::Jacobi,retain,count,passes,app),
                "Jacobi singular row policy must skip explicitly");
        require(count.singular==passes&&count.incompatible==passes,"Jacobi must report every unavailable energy visit");
        require(count.rank_deficient_blocks==1&&count.evaluations==1,"Jacobi rank classification must happen once per frozen solve");
        require(count.scalar_solves==passes,"Jacobi singular energy visits must not count as solves");
        near(actual.eta,0,0,"Jacobi hard zero-gradient row must preserve eta");
        near(actual.x,start.x-r.n*(r.rC/(m.mass*r.A)),2e-15,"Jacobi must still execute the material row");

        Context soft=hard;soft.beta=.1;
        actual=start;count={};
        require(schedule_update(soft,actual,Schedule::Jacobi,retain,count,passes,app),
                "Jacobi compliant zero-gradient row must solve");
        require(count.singular==0&&count.incompatible==0&&count.scalar_solves==2*passes,
                "Jacobi compliant zero-gradient classification");
        near(actual.eta,-r.H/(soft.beta/(soft.h*soft.h)),1e-15,"Jacobi compliant zero-gradient multiplier");
        near(actual.x,start.x-r.n*(r.rC/(m.mass*r.A)),2e-15,
             "Jacobi compliant zero-gradient multiplier must not add motion");
    }
}

void singular_rows() {
    Model m;m.pendulum=true;m.alpha=0;m.acceleration={0,-9.81};
    const State old{{std::sin(.65),-std::cos(.65)},{}};
    const Context hard=context(m,old,1./60,0,physical_energy(m,old));
    const Iterate start{hard.predictor};const Rows r=evaluate(hard,start);
    require(norm(r.g)==0&&r.H<0,"resting pendulum must have zero energy gradient and negative defect");
    for(const auto& order:orders)for(Application app:{Application::Deferred,Application::Immediate}) {
        Iterate actual=start;Counters count;
        require(schedule_update(hard,actual,order.schedule,false,count,1,app),"singular row should skip explicitly");
        const long energy_visits=order.schedule==Schedule::ELE?2:1;
        require(count.singular==energy_visits&&count.incompatible==energy_visits,
                "zero-gradient hard energy visits must each be reported incompatible");
        require(count.rank_deficient_blocks==1,"rank-deficient frozen system not reported");
        near(actual.x,start.x-r.n*(r.rC/(m.mass*r.A)),2e-15,"singular GS still executes material row");
        near(actual.eta,0,0,"singular hard row must not invent an energy multiplier");
    }
    Iterate actual=start;Counters count;
    require(schedule_update(hard,actual,Schedule::Joint,false,count),"singular joint policy should return a reported skip");
    compare(actual,start,0,"incompatible joint block must skip entirely");
    require(count.singular==1&&count.incompatible==1,"joint incompatible skip missing counters");

    Context soft=hard;soft.beta=.1;
    for(const auto& order:orders) {
        actual=start;count={};
        require(schedule_update(soft,actual,order.schedule,false,count),"compliant zero-gradient row must solve");
        require(count.singular==0&&count.incompatible==0,"compliance should remove zero-gradient row singularity");
        near(actual.eta,-r.H/(soft.beta/(soft.h*soft.h)),1e-15,"compliant zero-gradient dual");
        near(actual.x,start.x-r.n*(r.rC/(m.mass*r.A)),2e-15,"compliant zero-gradient dual must have no position contribution");
    }
    // Redundant but compatible hard energy row: choose delta eta=0 and solve L.
    Context compatible=hard;
    compatible.target=physical_energy(m,{start.x,(start.x-old.x)/hard.h});
    actual=start;count={};
    require(schedule_update(compatible,actual,Schedule::Joint,false,count),"compatible redundant block failed");
    require(count.singular==1&&count.incompatible==0,"compatible redundant row classification");
    near(actual.eta,0,0,"compatible redundant row representative");
    near(actual.x,start.x-r.n*(r.rC/(m.mass*r.A)),2e-15,"compatible block should still correct material");
}

void baseline_parity() {
    const Context c=fixture();Iterate actual=fixture_initial,expected=fixture_initial;
    Counters count,expected_count;
    require(schedule_update(c,actual,Schedule::Baseline,false,count),"new baseline failed");
    require(material_update(c,expected,expected_count),"existing baseline failed");
    compare(actual,expected,0,"existing baseline parity");
}

void resting_spring_branches() {
    Model m;m.alpha=1./64;
    const State old{{1.25,0},{}};
    const Context c=context(m,old,1./8,0,2);
    struct Answer { Schedule order; Iterate result; double energy_defect; };
    const std::vector<Answer> answers{
        {Schedule::Joint,{{1.25,0},-.25,1./64},0},
        {Schedule::Jacobi,{{1.125,0},-.125,0},-1},
        {Schedule::LE,{{1.25,0},-.125,1./128},0},
        {Schedule::EL,{{1.125,0},-.125,0},-1},
        {Schedule::LEL,{{1.1875,0},-.1875,1./128},-.75},
        {Schedule::ELE,{{1.25,0},-.125,1./128},0}
    };
    for(const auto& answer:answers) {
        Iterate actual{c.predictor};Counters count;
        require(schedule_update(c,actual,answer.order,false,count),"resting scalar spring update failed");
        compare(actual,answer.result,0,std::string(name(answer.order))+" analytic resting spring branch");
        near(physical_energy(m,{actual.x,(actual.x-old.x)/c.h})-c.target,
             answer.energy_defect,0,"analytic physical energy after one frozen pass");
    }
}

} // namespace

int main() {
    try {
        baseline_parity();
        ordered_frozen_solves();
        equivalence_and_convergence();
        frozen_equation_residuals();
        jacobi_solves();
        jacobi_singular_rows();
        singular_rows();
        resting_spring_branches();
        std::cout<<"PASS: frozen Jacobi and LE, EL, LEL, ELE orderings; immediate/deferred equivalence; mixed-system oracle; "
                    "retained stationarity once; joint convergence; singular policies; baseline/old LE parity\n";
        return 0;
    } catch(const std::exception& error) {
        std::cerr<<"FAIL: "<<error.what()<<'\n';return 1;
    }
}
