#include "oracle.hpp"
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>

using namespace compose;
namespace fs=std::filesystem;
void require(bool ok,const std::string&message){if(!ok)throw std::runtime_error(message);}

double angle_energy_error(const Model&m,const State&old,const State&s){
    const double L=m.length[0],g=-m.gravity[1];
    const double q=std::atan2(s.x[0],-s.x[1]),q0=std::atan2(old.x[0],-old.x[1]);
    const double w=(s.x[0]*s.v[1]-s.x[1]*s.v[0])/(L*L);
    const double w0=(old.x[0]*old.v[1]-old.x[1]*old.v[0])/(L*L);
    return std::abs(w*w-w0*w0-2*g/L*(std::cos(q)-std::cos(q0)));
}

// Reconstruct the lifted stationarity directly from the physical spring terms,
// independently of the sparse Problem or multiplier-GS implementation.
std::pair<double,double> lift_check(const hour_local::Problem&p,const dual_gs::OracleResult&r){
    auto e=p.evaluate(r.z);const auto&m=p.p.m;const int N=m.count,d=p.p.d,dim=m.dim;
    V gradient=e.objective.g;double material_error=0;
    if(!m.rigid)for(int j=0;j<N;++j){
        const V delta=edge(m,r.state.x,j),normal=edge_direction(m,delta);
        const double C=edge_length(m,delta)-m.length[j];
        const double alpha=4*m.mass.mean()/(p.p.h*p.p.h*m.stiffness[j]);
        material_error=std::max(material_error,std::abs(C/p.p.xs+alpha*r.y_lifted[j]));
        const V quarter_gradient=p.p.xs*m.stiffness[j]*C/(4*p.p.E)*normal;
        gradient.segment(j*dim,dim)-=quarter_gradient+r.y_lifted[j]*normal;
        if(j)gradient.segment((j-1)*dim,dim)+=quarter_gradient+r.y_lifted[j]*normal;
    }
    for(const auto&row:e.rows){int index=-1;
        if(row.kind==independent_detail::Position)index=row.index;
        else if(row.kind==independent_detail::Kinetic)index=N+row.index;
        else if(row.kind==independent_detail::Storage)index=2*N+row.index;
        else if(row.kind==independent_detail::GravityStorage)index=3*N+row.index;
        else if(row.kind==independent_detail::Tangency)index=3*N+(has_gravity(m)?N:0)+row.index;
        require(index>=0,"unknown_lift_row");gradient-=r.y_lifted[index]*row.value.g;
    }
    (void)d;
    return {material_error,gradient.norm()};
}

int main(int argc,char**argv){try{
    fs::path out=argc>1?argv[1]:"oracle_results";fs::create_directories(out);
    std::ofstream cases(out/"oracle_cases.csv"),rank(out/"oracle_rigid_rank.csv"),derivatives(out/"oracle_derivatives.csv");
    cases<<std::setprecision(16)<<"case,h,converged,reason,objective,primal,stationarity,residual,material_residual,lifted_stationarity,rank,rows,regular,tangent_dimension,curvature,local_minimum_certified,angle_energy_error,angle_displacement,expected_angle_branch,evaluations,iterations,fallbacks,milliseconds\n";
    rank<<std::setprecision(16)<<"state,converged,rows,rank,minimum_scalar_denominator,primal,stationarity,regular,local_minimum_certified\n";
    derivatives<<std::setprecision(16)<<"case,direction,relative_KKT_derivative_error\n";
    std::vector<Case> all;for(const auto&name:case_names())all.push_back(make_case(name));
    all.push_back(hour_local::chain_case(3));
    int failures=0,uncertified=0,runs=0;double max_lift=0,max_fd=0;
    for(const auto&c:all)for(double h:{.02,.08,.2}){
        const double E=std::max(1e-8,energy(c.model,c.initial)-minimum_energy(c.model));
        hour_local::Problem p(c.model,c.initial,h,E);
        const auto r=dual_gs::solve_oracle(c,h,E);
        ++runs;if(!r.converged)++failures;if(!r.local_minimum_certified)++uncertified;
        const auto lifted=lift_check(p,r);max_lift=std::max(max_lift,lifted.second);
        double ae=std::numeric_limits<double>::quiet_NaN(),dq=ae;int expected_branch=-1;
        if(c.model.rigid){
            ae=angle_energy_error(c.model,c.initial,r.state);
            dq=std::remainder(std::atan2(r.state.x[0],-r.state.x[1])-std::atan2(c.initial.x[0],-c.initial.x[1]),2*M_PI);
            if(r.converged){require(ae<1e-7,"rigid_angle_energy_curve");
                if(c.name=="rod_release")expected_branch=dq<0;
                if(c.name=="rod_kick")expected_branch=dq>0;}
        }
        cases<<c.name<<','<<h<<','<<r.converged<<','<<r.reason<<','<<r.objective<<','<<r.primal<<','<<r.stationarity<<','<<r.residual<<','<<lifted.first<<','<<lifted.second<<','<<r.jacobian_rank<<','<<r.constraint_count<<','<<r.regular<<','<<r.tangent_dimension<<','<<r.curvature<<','<<r.local_minimum_certified<<','<<ae<<','<<dq<<','<<expected_branch<<','<<r.evaluations<<','<<r.iterations<<','<<r.fallbacks<<','<<r.milliseconds<<'\n';
        require(lifted.first<1e-10,"material_lift_mapping");
        if(r.converged)require(lifted.second<2e-9,"lifted_stationarity_mapping");
    }
    for(const std::string label:{"bottom_rest","top_rest","nonvertical_rest","near_bottom","zero_gravity_rest"}){
        auto c=make_case("rod_release");double q=0,w=0;
        if(label=="top_rest")q=M_PI;
        if(label=="nonvertical_rest")q=.3;
        if(label=="near_bottom"){q=1e-5;w=1e-5;}
        if(label=="zero_gravity_rest")c.model.gravity.setZero();
        c.initial.x<<std::sin(q),-std::cos(q);c.initial.v<<w*std::cos(q),w*std::sin(q);
        hour_local::Problem p(c.model,c.initial,.02,1.);
        const auto initial=p.evaluate(p.initial());const Mat J=independent_detail::jacobian(initial);
        const Mat B=initial.objective.H.inverse();const Mat S=J*B*J.transpose();
        const auto r=dual_gs::solve_oracle(c,.02,1.);
        const int expected_rows=label=="zero_gravity_rest"?4:5;
        require(int(initial.rows.size())==expected_rows,"unexpected_or_duplicate_rigid_row");
        int lengths=0,tangencies=0;for(const auto&row:initial.rows){lengths+=row.kind==independent_detail::Position;tangencies+=row.kind==independent_detail::Tangency;}
        require(lengths==1&&tangencies==1,"duplicate_geometry_row");
        require(S.diagonal().minCoeff()>1e-6,"scalar_denominators_not_positive");
        const bool singular=label=="bottom_rest"||label=="top_rest"||label=="zero_gravity_rest";
        require(r.converged,"equilibrium_or_turning_oracle_failed_"+label);
        require(r.jacobian_rank==expected_rows-(singular?1:0),"rigid_rank_classification_"+label);
        if(singular)require(!r.local_minimum_certified,"singular_constraint_qualification_overclaimed");
        rank<<label<<','<<r.converged<<','<<r.constraint_count<<','<<r.jacobian_rank<<','<<S.diagonal().minCoeff()<<','<<r.primal<<','<<r.stationarity<<','<<r.regular<<','<<r.local_minimum_certified<<'\n';
    }
    for(const std::string name:{"scalar_mixed","rod_kick","chain_mixed"}){
        auto c=make_case(name);double E=std::max(1e-8,energy(c.model,c.initial)-minimum_energy(c.model));
        hour_local::Problem p(c.model,c.initial,.08,E);V z=p.initial();
        for(int i=0;i<z.size();++i)z[i]+=.003*std::sin(i+1.);
        auto e=p.evaluate(z);V y(e.rows.size());for(int i=0;i<y.size();++i)y[i]=.07*std::cos(i+1.);
        const Mat A=dual_gs::oracle_detail::newton_matrix(e,y);
        for(int probe=0;probe<3;++probe){V direction(z.size()+y.size());for(int i=0;i<direction.size();++i)direction[i]=std::sin((i+1.)*(probe+1.3));direction.normalize();const double eps=1e-6;
            V za=z+eps*direction.head(z.size()),zb=z-eps*direction.head(z.size()),ya=y+eps*direction.tail(y.size()),yb=y-eps*direction.tail(y.size());
            V fd=(dual_gs::oracle_detail::residual(p.evaluate(za),ya)-dual_gs::oracle_detail::residual(p.evaluate(zb),yb))/(2*eps);
            const double error=(fd-A*direction).norm()/(1+(A*direction).norm());max_fd=std::max(max_fd,error);
            derivatives<<name<<','<<probe<<','<<error<<'\n';require(error<2e-7,"minus_sign_exact_KKT_derivative");
        }
    }
    std::cout<<"Oracle runs="<<runs<<" converged="<<runs-failures<<" local_minima="<<runs-uncertified<<" max_lifted_stationarity="<<max_lift<<" max_KKT_derivative_error="<<max_fd<<'\n';
    // A legitimate failed solve or saddle is a retained reference outcome.
    // Consumers must gate root comparisons on converged/curvature flags.
    return 0;
}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}
