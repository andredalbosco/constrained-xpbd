#include "dual_gs.hpp"
#include "problem.hpp"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <set>
#include <sstream>
#include <stdexcept>

using dual_gs::Problem;
using dual_gs::Kind;
using compose::V;
using compose::Mat;

namespace {
struct Audit {
    std::ofstream csv{"validation.csv"};
    int checks=0;
    double worst_gradient=0,worst_hessian=0,worst_objective=0,worst_update=0,worst_frozen=0;
    Audit(){csv<<"check,case,error,tolerance,passed\n"<<std::setprecision(16);}
    void check(bool pass,const std::string&name,const std::string&scene,double error=0,double tolerance=0){
        ++checks;csv<<name<<','<<scene<<','<<error<<','<<tolerance<<','<<pass<<'\n';csv.flush();
        if(!pass)throw std::runtime_error(name+" / "+scene+" error="+std::to_string(error));
    }
    void near(double error,double tolerance,const std::string&name,const std::string&scene){check(std::isfinite(error)&&error<=tolerance,name,scene,error,tolerance);}
};
double relative(const V&a,const V&b){return (a-b).norm()/(1+std::max(a.norm(),b.norm()));}
double relative(const Mat&a,const Mat&b){return (a-b).norm()/(1+std::max(a.norm(),b.norm()));}
V dense(const dual_gs::RowEval&e,int n){V g=V::Zero(n);for(const auto&a:e.grad)g[a.index]+=a.value;return g;}
V pick(const V&v,const std::vector<int>&s){V q(s.size());for(int i=0;i<int(s.size());++i)q[i]=v[s[i]];return q;}
Mat pick(const Mat&A,const std::vector<int>&s){Mat q(s.size(),s.size());for(int i=0;i<int(s.size());++i)for(int j=0;j<int(s.size());++j)q(i,j)=A(s[i],s[j]);return q;}
Mat metric_oracle(const Problem&p){
    Mat A=Mat::Zero(p.n,p.n);
    for(int i=0;i<p.d;++i){double c=p.model.mass[i/p.model.dim]/p.mref;
        A(i,i)=c*(1+4/p.beta);A(p.d+i,p.d+i)=c/p.beta;
        A(i,p.d+i)=A(p.d+i,i)=-2*c/p.beta;
    }
    for(int j=2*p.d;j<p.n;++j)A(j,j)=p.gamma;
    return A;
}
V initial_oracle(const Problem&p){
    V z=V::Zero(p.n);
    for(int i=0;i<p.d;++i){double s=p.old.v[i]+.25*p.h*(p.aold[i]+p.model.gravity[i%p.model.dim]);z[i]=s/p.vs;z[p.d+i]=(2*s-p.old.v[i])/p.vs;}
    z.tail(p.wcount)=p.what;return z;
}
Mat gradients(const Problem&p,const V&z){Mat G(p.n,p.rows.size());for(int a=0;a<int(p.rows.size());++a)G.col(a)=dense(p.eval_row(a,z),p.n);return G;}
V values(const Problem&p,const V&z){V b(p.rows.size());for(int a=0;a<int(p.rows.size());++a)b[a]=p.eval_row(a,z).value;return b;}
V compliances(const Problem&p){V D(p.rows.size());for(int a=0;a<int(p.rows.size());++a)D[a]=p.compliance(a);return D;}
V stationarity_oracle(const Problem&p,const V&z,const V&y){return metric_oracle(p)*(z-initial_oracle(p))-gradients(p,z)*y;}
std::vector<int> expected_support(const Problem&p,int row){
    const auto&r=p.rows[row];std::set<int> indices;int dim=p.model.dim;
    auto particle=[&](int i){if(i<0||i>=p.model.count)return;for(int k=0;k<dim;++k){indices.insert(i*dim+k);indices.insert(p.d+i*dim+k);}};
    auto work=[&](int j){if(j>=0)indices.insert(j);};
    if(r.kind==Kind::Kinetic){int i=r.index;particle(i);work(p.work_child[i]);if(i+1<p.model.count)work(p.work_parent[i+1]);work(p.work_gravity[i]);}
    else if(r.kind==Kind::GravityStorage){particle(r.index);work(p.work_gravity[r.index]);}
    else{int j=r.index;if(!(r.kind==Kind::Storage&&p.model.rigid)){particle(j);particle(j-1);}if(r.kind==Kind::Storage){work(p.work_child[j]);work(p.work_parent[j]);}}
    return {indices.begin(),indices.end()};
}
int hour_kind(Kind k){
    if(k==Kind::Kinetic)return int(compose::independent_detail::Kinetic);
    if(k==Kind::Storage)return int(compose::independent_detail::Storage);
    if(k==Kind::GravityStorage)return int(compose::independent_detail::GravityStorage);
    if(k==Kind::Length)return int(compose::independent_detail::Position);
    if(k==Kind::Tangency)return int(compose::independent_detail::Tangency);
    return -1;
}
void derivatives(Audit&a,const compose::Case&c,double h){
    double E=std::max(1e-8,compose::energy(c.model,c.initial)-compose::minimum_energy(c.model));
    Problem p(c.model,c.initial,h,E);hour_local::Problem old(c.model,c.initial,h,E);
    const std::string scene=c.name+"_h"+std::to_string(h);Mat A=metric_oracle(p);
    a.near(relative(p.dense_metric(),A),1e-13,"base_metric",scene);
    a.near(relative(p.initial(),initial_oracle(p)),1e-13,"zero_dual_base_target",scene);
    a.near(p.base_gradient(p.initial()).norm(),1e-11,"base_target_stationarity",scene);
    a.check(p.n==old.p.n,"unchanged_primal_dimension",scene);
    int expected_rows=p.model.count*(p.model.rigid?3:3)+(compose::has_gravity(p.model)?p.model.count:0)+(p.model.rigid?p.model.count:0);
    a.check(int(p.rows.size())==expected_rows,"lifted_row_count",scene);
    for(int i=0;i<p.model.count;++i){int count=0;for(const auto&r:p.rows)count+=r.kind==Kind::Kinetic&&r.index==i;a.check(count==1,"single_shared_kinetic_multiplier",scene);}
    for(int row=0;row<int(p.rows.size());++row){auto s=p.support(row),expected=expected_support(p,row);std::sort(s.begin(),s.end());a.check(s==expected,"structural_metric_closure",scene+"_r"+std::to_string(row));}
    for(int sample=0;sample<3;++sample){
        V z=p.initial();for(int j=0;j<p.n;++j)z[j]+=.023*sample*std::sin(.7*(j+1)+.3);
        V y(p.rows.size());for(int j=0;j<y.size();++j)y[j]=.03*std::cos(.9*(j+1));
        auto oracle=old.oracle(z);double objective_error=std::abs(p.objective(z)-oracle.objective.a)/(1+std::abs(oracle.objective.a));a.worst_objective=std::max(a.worst_objective,objective_error);a.near(objective_error,1e-11,"objective_value_oracle",scene);
        Mat G=gradients(p,z);V base=A*(z-initial_oracle(p));a.near(relative(p.base_gradient(z),base),1e-12,"base_gradient",scene);a.near(relative(p.stationarity(z,y),base-G*y),1e-11,"whole_stationarity",scene);
        double lifted_value=p.base_objective(z);V eliminated_gradient=base;Mat eliminated_H=A;
        V sum_gradient=V::Zero(p.n);double sum_value=0;
        for(int row=0;row<int(p.rows.size());++row){
            auto r=p.eval_row(row,z);V g=G.col(row),fd(p.n);Mat H(p.n,p.n);
            for(int j=0;j<p.n;++j){double eps=2e-6;V zp=z,zm=z;zp[j]+=eps;zm[j]-=eps;auto rp=p.eval_row(row,zp),rm=p.eval_row(row,zm);fd[j]=(rp.value-rm.value)/(2*eps);H.col(j)=(dense(rp,p.n)-dense(rm,p.n))/(2*eps);}
            double ge=relative(g,fd);a.worst_gradient=std::max(a.worst_gradient,ge);a.near(ge,3e-7,"sparse_row_gradient_FD",scene+"_r"+std::to_string(row));
            std::set<int> support(p.support(row).begin(),p.support(row).end());for(const auto&e:r.grad)a.check(support.count(e.index)!=0,"gradient_inside_support",scene);
            auto local=p.stationarity_local(row,z,y);a.near(relative(local,pick(V(base-G*y),p.support(row))),1e-11,"local_stationarity_closure",scene+"_r"+std::to_string(row));
            V local_rhs=pick(V(base-G*y),p.support(row));V inv_expected=pick(A,p.support(row)).ldlt().solve(local_rhs);a.near(relative(p.apply_inv_local(p.support(row),local_rhs),inv_expected),1e-12,"local_metric_inverse",scene);
            if(p.rows[row].kind==Kind::Material){
                double alpha=p.compliance(row),lambda=-r.value/alpha;lifted_value-=lambda*r.value+.5*alpha*lambda*lambda;
                eliminated_gradient-=lambda*g;eliminated_H+=(g*g.transpose()+r.value*H)/alpha;
                double expected_alpha=4*p.mref/(p.h*p.h*p.model.stiffness[p.rows[row].index]);a.near(std::abs(alpha-expected_alpha)/(1+expected_alpha),1e-13,"material_compliance",scene);
            }else{
                int wanted=hour_kind(p.rows[row].kind),found=-1;for(int j=0;j<int(oracle.rows.size());++j)if(int(oracle.rows[j].kind)==wanted&&oracle.rows[j].index==p.rows[row].index)found=j;
                a.check(found>=0,"independent_row_match",scene);
                const auto&o=oracle.rows[found].value;a.near(std::abs(r.value-o.a)/(1+std::abs(o.a)),1e-11,"row_value_AD_oracle",scene);a.near(relative(g,o.g),1e-11,"row_gradient_AD_oracle",scene);
                double he=relative(H,o.H);a.worst_hessian=std::max(a.worst_hessian,he);a.near(he,5e-6,"row_hessian_AD_oracle",scene);
                if(p.rows[row].kind==Kind::Kinetic||p.rows[row].kind==Kind::Storage||p.rows[row].kind==Kind::GravityStorage){sum_value+=r.value;sum_gradient+=g;}
            }
        }
        a.near(std::abs(lifted_value-oracle.objective.a)/(1+std::abs(oracle.objective.a)),1e-11,"material_lift_value_elimination",scene);
        a.near(relative(eliminated_gradient,oracle.objective.g),2e-10,"material_lift_gradient_elimination",scene);
        double hh=relative(eliminated_H,oracle.objective.H);a.worst_hessian=std::max(a.worst_hessian,hh);a.near(hh,5e-6,"material_lift_hessian_elimination",scene);
        auto state=p.state(z);double de=(compose::energy(p.model,state)-compose::energy(p.model,p.old))/E;a.near(std::abs(sum_value-de)/(1+std::abs(de)),1e-11,"sum_balances_total_energy",scene);
        for(int j=2*p.d;j<p.n;++j)a.near(std::abs(sum_gradient[j]),1e-14,"work_cancellation",scene);
    }
}
void local_visits(Audit&a,const compose::Case&c){
    double E=std::max(1e-8,compose::energy(c.model,c.initial)-compose::minimum_energy(c.model));Problem p(c.model,c.initial,1./120,E);Mat A=metric_oracle(p);
    V z=p.initial(),y(p.rows.size());for(int j=0;j<p.n;++j)z[j]+=.017*std::sin(.8*j+.2);for(int j=0;j<y.size();++j)y[j]=.012*std::cos(.7*j+.6);
    for(int row=0;row<int(p.rows.size());++row)for(std::string method:{"refreshed_drop","refreshed_local_rx"}){
        V actual_z=z,actual_y=y;dual_gs::Counters counts;auto support=p.support(row);int n=support.size();auto r=p.eval_row(row,z);V g=dense(r,p.n),gs=pick(g,support);double residual=r.value+p.compliance(row)*y[row];V rz=stationarity_oracle(p,z,y),delta;
        if(method=="refreshed_local_rx"){
            Mat K=Mat::Zero(n+1,n+1);K.topLeftCorner(n,n)=pick(A,support);K.block(0,n,n,1)=-gs;K.block(n,0,1,n)=gs.transpose();K(n,n)=p.compliance(row);V rhs(n+1);rhs.head(n)=-pick(rz,support);rhs[n]=-residual;delta=K.fullPivLu().solve(rhs);
        }else{Mat AS=pick(A,support);V q=AS.ldlt().solve(gs);delta=V::Zero(n+1);delta[n]=-residual/(gs.dot(q)+p.compliance(row));delta.head(n)=q*delta[n];}
        V expected_z=z,expected_y=y;for(int j=0;j<n;++j)expected_z[support[j]]+=delta[j];expected_y[row]+=delta[n];
        std::string failure=dual_gs::visit_refreshed(p,row,actual_z,actual_y,method,counts);a.check(failure.empty(),"local_visit_finite",c.name+"_"+method+"_r"+std::to_string(row));
        double error=std::max(relative(actual_z,expected_z),relative(actual_y,expected_y));a.worst_update=std::max(a.worst_update,error);a.near(error,5e-11,"restricted_KKT_oracle",c.name+"_"+method+"_r"+std::to_string(row));
        std::set<int> touched(support.begin(),support.end());for(int j=0;j<p.n;++j)if(!touched.count(j))a.check(actual_z[j]==z[j],"untouched_primal_coordinate",c.name);
        for(int j=0;j<y.size();++j)if(j!=row)a.check(actual_y[j]==y[j],"single_dual_coordinate",c.name);
        a.check(counts.visits==1,"single_visit_count",c.name);
        if(method=="refreshed_drop"){
            V drift=rz-(gradients(p,actual_z)-gradients(p,z))*actual_y;V new_r=stationarity_oracle(p,actual_z,actual_y);a.near(relative(new_r,drift),5e-11,"dropped_rx_drift_identity",c.name);
        }else{
            V local_linear=pick(rz,support)+pick(A,support)*delta.head(n)-gs*delta[n];a.near(local_linear.norm()/(1+pick(rz,support).norm()),5e-11,"restricted_stationarity_restored",c.name);
            double local_row=residual+gs.dot(delta.head(n))+p.compliance(row)*delta[n];a.near(std::abs(local_row)/(1+std::abs(residual)),5e-11,"restricted_linear_row_restored",c.name);
        }
    }
}
void frozen_sweeps(Audit&a,const compose::Case&c){
    double E=std::max(1e-8,compose::energy(c.model,c.initial)-compose::minimum_energy(c.model));Problem p(c.model,c.initial,1./120,E);Mat A=metric_oracle(p),inv=A.inverse();V D=compliances(p);
    for(std::string order:{"material_storage_kinetic","kinetic_storage_material","storage_kinetic_material","edge","edge_storage_first","edge_reverse"}){
        auto schedule=p.schedule_indices(order);a.check(!schedule.empty(),"nonempty_schedule",c.name+"_"+order);for(int r:schedule)a.check(r>=0&&r<int(p.rows.size()),"schedule_row_identity",c.name);
        a.check(schedule==p.schedule_indices(order),"deterministic_schedule",c.name+"_"+order);
        if(c.model.count>=3&&order.find("edge")==0){int repeats=0;for(int r:schedule)repeats+=p.rows[r].kind==Kind::Kinetic&&p.rows[r].index==1;a.check(repeats==2,"middle_mu_revisited_not_duplicated",c.name+"_"+order);}
        V z=initial_oracle(p),y=V::Zero(p.rows.size());
        for(int outer=0;outer<2;++outer){
            V z0=z,y0=y;Mat G=gradients(p,z0),S=G.transpose()*inv*G;S.diagonal()+=D;V r=values(p,z0)+D.cwiseProduct(y0),rz=A*(z0-initial_oracle(p))-G*y0,t=-r+G.transpose()*inv*rz,dy=V::Zero(y.size());
            for(int row:schedule){double denominator=S(row,row);a.check(denominator>0,"positive_schur_diagonal",c.name);dy[row]+=(t[row]-S.row(row).dot(dy))/denominator;}
            z=z0-inv*rz+inv*G*dy;y=y0+dy;
            auto actual=dual_gs::solve(p,"frozen_full_rx",order,outer+1);a.check(actual.failure.empty(),"frozen_sweep_finite",c.name+"_"+order);
            double error=std::max(relative(actual.z,z),relative(actual.y,y));a.worst_frozen=std::max(a.worst_frozen,error);a.near(error,2e-10,"explicit_Schur_GS",c.name+"_"+order+"_s"+std::to_string(outer+1));
            V restored=rz+A*(z-z0)-G*(y-y0);a.near(restored.norm()/(1+rz.norm()),1e-10,"frozen_stationarity_identity",c.name);
        }
    }
}
}
int main(){try{
    Eigen::setNbThreads(1);Audit audit;std::vector<compose::Case> cases;
    for(const std::string name:{"scalar_release","radial_orbit","rod_kick","elastic_mixed","chain_mixed"})cases.push_back(compose::make_case(name));
    cases.push_back(hour_local::chain_case(3));
    for(const auto&c:cases){for(double h:{1./60,1./120})derivatives(audit,c,h);local_visits(audit,c);frozen_sweeps(audit,c);}
    std::cout<<std::setprecision(16)<<"PASS checks="<<audit.checks<<" max_gradient="<<audit.worst_gradient<<" max_hessian="<<audit.worst_hessian<<" max_objective="<<audit.worst_objective<<" max_local_update="<<audit.worst_update<<" max_frozen_schur="<<audit.worst_frozen<<'\n';return 0;
}catch(const std::exception&e){std::cerr<<"VALIDATION FAILED: "<<e.what()<<'\n';return 1;}}
