#include "independent.hpp"
#include <iostream>
#include <iomanip>
using namespace compose;
int main(){
    using namespace independent_detail;
    double worst_gradient=0,worst_hessian=0,worst_balance=0;int checks=0;
    const std::vector<std::string> works={"start","end","force_average","separate_average","dg_start","dg_end","dg_average","frozen_power","kinetic_start","kinetic_end"};
    for(auto name:{"scalar_mixed","radial_orbit","rod_kick","chain_mixed"})for(auto work:works)for(auto kinematics:{"hard_endpoint","soft_endpoint","soft_midpoint"}){
        auto c=make_case(name);Config cfg;cfg.work=work;cfg.kinematics=kinematics;cfg.treatment="penalized";cfg.beta=2;cfg.tangency=c.model.rigid;cfg.curvature="exact";
        double E=std::max(1.,energy(c.model,c.initial)-minimum_energy(c.model));Problem p(c.model,c.initial,.03,cfg,E);V z=p.initial();
        for(int k=0;k<z.size();++k)z[k]+=.025*std::sin(k+1.3);
        Evaluation e=p.evaluate(z);const double eps=2e-5;
        for(int k=0;k<z.size();++k){V a=z,b=z;a[k]+=eps;b[k]-=eps;auto ea=p.evaluate(a),eb=p.evaluate(b);
            double numeric=(ea.objective.a-eb.objective.a)/(2*eps);worst_gradient=std::max(worst_gradient,std::abs(numeric-e.objective.g[k])/(1+std::abs(numeric)));
            V nh=(ea.objective.g-eb.objective.g)/(2*eps);worst_hessian=std::max(worst_hessian,(nh-e.objective.H.col(k)).norm()/(1+nh.norm()));
            for(int r=0;r<int(e.rows.size());++r){double ng=(ea.rows[r].value.a-eb.rows[r].value.a)/(2*eps);worst_gradient=std::max(worst_gradient,std::abs(ng-e.rows[r].value.g[k])/(1+std::abs(ng)));V hh=(ea.rows[r].value.g-eb.rows[r].value.g)/(2*eps);worst_hessian=std::max(worst_hessian,(hh-e.rows[r].value.H.col(k)).norm()/(1+hh.norm()));}
            ++checks;
        }
    }
    // The DG-average estimate must exactly exchange quadratic elastic energy,
    // including rotating/stretching two-body edges and affine scalar springs.
    double dg_defect=0;
    for(auto name:{"scalar_mixed","radial_orbit","chain_mixed"}){
        auto c=make_case(name);Config cfg;cfg.work="dg_average";cfg.treatment="prescribed";double E=std::max(1.,energy(c.model,c.initial)-minimum_energy(c.model));Problem p(c.model,c.initial,.03,cfg,E);V z=p.initial();for(int k=0;k<p.d;++k)z[k]+=.3*std::cos(k+1.);auto e=p.evaluate(z);dg_defect=std::max(dg_defect,e.estimate_defect);for(auto&r:e.rows)if(r.kind==Storage||r.kind==GravityStorage){std::cerr<<"identity row retained\n";return 2;}
    }
    // The global soft-velocity formulation should be able to solve a release
    // and a moving state without relying on the endpoint energy identity.
    for(auto name:{"scalar_release","scalar_mixed","radial_orbit","rod_kick","chain_mixed"}){
        auto c=make_case(name);Config cfg;cfg.work="dg_average";cfg.treatment="global";cfg.beta=2;cfg.iterations=30;cfg.curvature="exact";cfg.tangency=c.model.rigid;double E=std::max(1.,energy(c.model,c.initial)-minimum_energy(c.model));auto r=independent_step(c.model,c.initial,.02,cfg,E);
        std::cout<<name<<" valid="<<r.valid<<" balance="<<r.balance<<" stationarity="<<r.rx<<" iterations="<<r.iterations<<" visits="<<r.visits<<" evaluations="<<r.evaluations<<"\n";
        if(!r.valid)return 3;worst_balance=std::max(worst_balance,r.balance);
        if(r.balance>2e-7||r.rx>1e-5){std::cerr<<"global solve failed\n";return 4;}
        if(std::string(name)=="scalar_release"&&(r.state.x-c.initial.x).norm()<1e-7){std::cerr<<"release unexpectedly frozen\n";return 5;}
    }
    // Identical free/global parameterizations must give identical steps.
    auto c=make_case("chain_mixed");Config a;a.beta=2;a.iterations=15;Config b=a;b.treatment="free";b.schedule="constraint_block";double E=energy(c.model,c.initial)-minimum_energy(c.model);auto ra=independent_step(c.model,c.initial,.02,a,E),rb=independent_step(c.model,c.initial,.02,b,E);if((ra.state.x-rb.state.x).norm()+(ra.state.v-rb.state.v).norm()>1e-12)return 6;
    // A shared particle really appears in both edge blocks; one work correction
    // enters its kinetic row and the supplying edge's storage row oppositely.
    Config local;local.treatment="penalized";local.schedule="constraint_block";Problem lp(c.model,c.initial,.02,local,E);V lz=lp.initial();auto le=lp.evaluate(lz);auto lb=groups(le,lp);int shared_visits=0;for(auto&block:lb)for(int r:block)if(le.rows[r].kind==Kinetic&&le.rows[r].index==0)++shared_visits;if(shared_visits!=2)return 8;
    lz[2*lp.d+lp.work_parent[1]]+=.25;auto lt=lp.evaluate(lz);for(int r=0;r<int(le.rows.size());++r){double change=lt.rows[r].value.a-le.rows[r].value.a;double expect=(le.rows[r].kind==Kinetic&&le.rows[r].index==0)?-.25:((le.rows[r].kind==Storage&&le.rows[r].index==1)?.25:0.);if(std::abs(change-expect)>1e-12)return 9;}
    // Finite-budget block solves may retain residuals: test finite execution and
    // visit semantics, without treating those residuals as solved physics.
    int local_runs=0;double local_worst_balance=0;
    for(auto name:{"scalar_mixed","rod_kick","chain_mixed"})for(auto treatment:{"prescribed","penalized"})for(auto schedule:{"joint","constraint_block","separate","symmetric"})for(auto work:{"end","dg_average","kinetic_start"}){
        auto cc=make_case(name);Config cfg;cfg.treatment=treatment;cfg.schedule=schedule;cfg.work=work;cfg.beta=2;cfg.iterations=4;cfg.tangency=cc.model.rigid;double scale=std::max(1.,energy(cc.model,cc.initial)-minimum_energy(cc.model));auto r=independent_step(cc.model,cc.initial,.02,cfg,scale);if(!r.valid||!std::isfinite(r.balance)||!std::isfinite(r.rx)){std::cerr<<"local schedule fatal: "<<name<<" "<<treatment<<" "<<schedule<<" "<<work<<" "<<r.reason<<"\n";return 10;}++local_runs;local_worst_balance=std::max(local_worst_balance,r.balance);
    }
    std::cout<<"finite_local_runs="<<local_runs<<" maximum_unsolved_local_balance="<<local_worst_balance<<" shared_vertex_block_visits="<<shared_visits<<"\n";
    std::cout<<std::setprecision(12)<<"derivative_checks="<<checks<<" worst_relative_gradient="<<worst_gradient<<" worst_relative_hessian="<<worst_hessian<<" dg_storage_defect="<<dg_defect<<" worst_solved_balance="<<worst_balance<<"\n";
    if(worst_gradient>2e-6||worst_hessian>2e-5||dg_defect>1e-12)return 7;
    return 0;
}
