#include "dual_gs.hpp"
#include "oracle.hpp"
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>

using compose::V;
using compose::Mat;

namespace {
double full_residual(const dual_gs::Residuals&r){
    return std::max({r.material,r.kinetic,r.storage,r.geometry,r.tangency,r.stationarity});
}
V sweep(const dual_gs::Problem&p,const V&w,const std::string&method,const std::string&schedule){
    V z=w.head(p.n),y=w.tail(p.rows.size());
    dual_gs::Counters counts;
    for(int row:p.schedule_indices(schedule)){
        const auto reason=dual_gs::visit_refreshed(p,row,z,y,method,counts,1.);
        if(!reason.empty())throw std::runtime_error(reason);
    }
    V out(w.size());out<<z,y;
    if(!out.allFinite())throw std::runtime_error("nonfinite_map");
    return out;
}
Mat difference(const dual_gs::Problem&p,const V&w,const std::string&method,const std::string&schedule,double eps){
    Mat J(w.size(),w.size());
    for(int j=0;j<w.size();++j){
        V a=w,b=w;a[j]+=eps;b[j]-=eps;
        J.col(j)=(sweep(p,a,method,schedule)-sweep(p,b,method,schedule))/(2*eps);
    }
    return J;
}
// Independent derivative of the dropped map at a root, from STABILITY_MATH.md.
Mat drop_derivative(const dual_gs::Problem&p,const V&w,const std::string&schedule){
    const int q=p.rows.size(),n=p.n,total=n+q;
    Mat result=Mat::Identity(total,total);
    const Mat A=p.dense_metric();
    for(int row:p.schedule_indices(schedule)){
        const auto evaluation=p.eval_row(row,w.head(n));
        V g=V::Zero(n);for(const auto&e:evaluation.grad)g[e.index]+=e.value;
        V direction=A.ldlt().solve(g),u=V::Zero(total),v=V::Zero(total);
        u.head(n)=direction;u[n+row]=1;v.head(n)=g;v[n+row]=p.compliance(row);
        result=(Mat::Identity(total,total)-u*v.transpose()/(g.dot(direction)+p.compliance(row)))*result;
    }
    return result;
}
Mat local_derivative(const dual_gs::Problem&p,const V&w,const Mat&H,const std::string&schedule){
    const int q=p.rows.size(),n=p.n,total=n+q;
    Mat result=Mat::Identity(total,total),R(n,total);
    R.leftCols(n)=H;R.rightCols(q)=-p.dense_gradients(w.head(n));
    const Mat Ainv=p.dense_metric().inverse();
    for(int row:p.schedule_indices(schedule)){
        const auto evaluation=p.eval_row(row,w.head(n));
        V g=V::Zero(n);for(const auto&e:evaluation.grad)g[e.index]+=e.value;
        Mat B=Mat::Zero(n,n);
        for(int i:p.support(row))for(int j:p.support(row))B(i,j)=Ainv(i,j);
        V direction=B*g,u=V::Zero(total),v=V::Zero(total);
        u.head(n)=direction;u[n+row]=1;v.head(n)=g;v[n+row]=p.compliance(row);
        Mat embedded=Mat::Zero(total,n);embedded.topRows(n)=B;
        const Mat visit=Mat::Identity(total,total)-embedded*R
            +u*(-v.transpose()+g.transpose()*B*R)/(g.dot(direction)+p.compliance(row));
        result=visit*result;
    }
    return result;
}
struct Spectrum {
    double rho=std::numeric_limits<double>::quiet_NaN(),re=rho,im=rho,norm=rho;
    int unit=0,unstable=0;
};
Spectrum spectrum(const Mat&J){
    Eigen::EigenSolver<Mat> eig(J,false);
    if(eig.info()!=Eigen::Success)throw std::runtime_error("eigen_solver_failed");
    Spectrum s;s.rho=-1;
    for(auto lambda:eig.eigenvalues()){
        if(std::abs(lambda)>s.rho){s.rho=std::abs(lambda);s.re=lambda.real();s.im=lambda.imag();}
        if(std::abs(lambda-std::complex<double>(1.,0.))<1e-5)++s.unit;
        if(std::abs(lambda)>1+1e-5)++s.unstable;
    }
    Eigen::JacobiSVD<Mat> svd(J);s.norm=svd.singularValues()[0];
    return s;
}
void identity_checks(){
    auto c=compose::make_case("scalar_mixed");
    double E=compose::energy(c.model,c.initial)-compose::minimum_energy(c.model);
    dual_gs::Problem p(c.model,c.initial,.08,E);
    auto root=dual_gs::solve_oracle(c,.08,E);
    if(!root.converged)throw std::runtime_error("identity_oracle_failed");
    V w(p.n+p.rows.size());w<<root.z,root.y_lifted;
    for(std::string schedule:{"material_storage_kinetic","storage_kinetic_material","edge_reverse"}){
        auto J=difference(p,w,"refreshed_drop",schedule,1e-6);
        auto exact=drop_derivative(p,w,schedule);
        double error=(J-exact).norm()/(1+exact.norm());
        if(error>1e-7)throw std::runtime_error("drop_derivative_identity_failed");
        auto s=spectrum(J);
        if(s.unit<p.n)throw std::runtime_error("drop_feasible_tangent_identity_failed");
        for(std::string method:{"refreshed_drop","refreshed_local_rx"})
            if((sweep(p,w,method,schedule)-w).norm()>1e-8)throw std::runtime_error("root_not_fixed");
        std::cout<<"identity "<<schedule<<" relative_J_error="<<error<<" unit_modes="<<s.unit<<'\n';
    }
}
}

int main(int argc,char**argv){try{
    Eigen::setNbThreads(1);identity_checks();
    const std::string output=argc>1?argv[1]:"results/stability.csv";
    const auto parent=std::filesystem::path(output).parent_path();
    if(!parent.empty())std::filesystem::create_directories(parent);
    std::ofstream out(output);if(!out)throw std::runtime_error("output_not_writable");
    out<<std::setprecision(17);
    out<<"case,N,h,method,schedule,epsilon,dimension,primal_dimension,rows,sweep_visits,oracle_converged,oracle_reason,oracle_primal,oracle_stationarity,oracle_curvature,root_full_residual,after_full_residual,fixed_point_defect,spectral_radius,leading_real,leading_imag,unit_modes,unstable_modes,matrix_2norm,jacobian_relative_difference,radius_difference,analytic_J_error,analytic_spectral_radius,analytic_matrix_2norm,status\n";
    std::vector<compose::Case> cases{compose::make_case("chain_mixed"),hour_local::chain_case(3),hour_local::chain_case(8)};
    for(const auto&c:cases)for(double h:{.02,.08,.2}){
        const double E=std::max(1e-8,compose::energy(c.model,c.initial)-compose::minimum_energy(c.model));
        dual_gs::Problem p(c.model,c.initial,h,E);auto root=dual_gs::solve_oracle(c,h,E);
        Mat lifted_H;
        if(root.converged){
            hour_local::Problem independent(c.model,c.initial,h,E);
            lifted_H=dual_gs::oracle_detail::lagrangian_hessian(independent.evaluate(root.z),root.y_unlifted);
            for(int row:p.material_rows){
                V g=V::Zero(p.n);for(const auto&e:p.eval_row(row,root.z).grad)g[e.index]+=e.value;
                lifted_H-=g*g.transpose()/p.compliance(row);
            }
        }
        std::cerr<<"oracle "<<c.name<<" h="<<h<<" converged="<<root.converged<<" reason="<<root.reason<<'\n';
        for(std::string method:{"refreshed_drop","refreshed_local_rx"})
        for(std::string schedule:{"material_storage_kinetic","storage_kinetic_material","edge_reverse"}){
            const double nan=std::numeric_limits<double>::quiet_NaN();
            double before=nan,after=nan,defect=nan,agreement=nan,rho_difference=nan,analytic_error[2]={nan,nan};
            Spectrum spec[2],analytic;std::string status="oracle_not_converged";
            if(root.converged){try{
                if(root.z.size()!=p.n||root.y_lifted.size()!=int(p.rows.size()))throw std::runtime_error("root_dimension_mismatch");
                V w(p.n+p.rows.size());w<<root.z,root.y_lifted;
                before=full_residual(p.residuals(root.z,root.y_lifted));
                if(!(before<1e-7))throw std::runtime_error("lifted_root_residual_gate");
                const Mat exact=method=="refreshed_drop"?drop_derivative(p,w,schedule):local_derivative(p,w,lifted_H,schedule);
                analytic=spectrum(exact);
                const V mapped=sweep(p,w,method,schedule);defect=(mapped-w).norm();
                after=full_residual(p.residuals(mapped.head(p.n),mapped.tail(p.rows.size())));
                Mat J[2]={difference(p,w,method,schedule,1e-5),difference(p,w,method,schedule,1e-6)};
                for(int k=0;k<2;++k)spec[k]=spectrum(J[k]);
                agreement=(J[0]-J[1]).norm()/(1+J[1].norm());
                rho_difference=std::abs(spec[0].rho-spec[1].rho);
                for(int k=0;k<2;++k)analytic_error[k]=(J[k]-exact).norm()/(1+exact.norm());
                status="ok";
                if(defect>1e-7)status="fixed_point_defect_too_large";
                else if(agreement>1e-3||rho_difference>1e-3*std::max(1.,spec[1].rho))status="finite_difference_scale_disagreement";
            }catch(const std::exception&e){status=e.what();}}
            for(int k=0;k<2;++k){
                const auto&s=spec[k];
                out<<c.name<<','<<c.model.count<<','<<h<<','<<method<<','<<schedule<<','<<(k?1e-6:1e-5)<<','<<p.n+p.rows.size()<<','<<p.n<<','<<p.rows.size()<<','<<p.schedule_indices(schedule).size()<<','<<root.converged<<','<<root.reason<<','<<root.primal<<','<<root.stationarity<<','<<root.curvature<<','<<before<<','<<after<<','<<defect<<','<<s.rho<<','<<s.re<<','<<s.im<<','<<s.unit<<','<<s.unstable<<','<<s.norm<<','<<agreement<<','<<rho_difference<<','<<analytic_error[k]<<','<<analytic.rho<<','<<analytic.norm<<','<<status<<'\n';
            }
            out.flush();
        }
    }
    return 0;
}catch(const std::exception&e){std::cerr<<"STABILITY FAILED: "<<e.what()<<'\n';return 1;}}
