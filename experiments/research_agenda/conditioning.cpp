#include "cases.hpp"
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
using namespace agenda;
std::array<double,3> eigenvalues(Matrix a,int n){
    for(int it=0;it<80;++it){int p=0,q=1;double off=0,diag=0;for(int i=0;i<n;++i){diag=std::max(diag,std::abs(a[i][i]));for(int j=i+1;j<n;++j)if(std::abs(a[i][j])>off){off=std::abs(a[i][j]);p=i;q=j;}}
        if(off<=2e-16*diag)break;const double theta=.5*std::atan2(2*a[p][q],a[q][q]-a[p][p]),s=std::sin(theta),c=std::cos(theta);const double app=a[p][p],aqq=a[q][q],apq=a[p][q];
        a[p][p]=c*c*app-2*s*c*apq+s*s*aqq;a[q][q]=s*s*app+2*s*c*apq+c*c*aqq;a[p][q]=a[q][p]=0;
        for(int k=0;k<n;++k)if(k!=p&&k!=q){const double kp=a[k][p],kq=a[k][q];a[k][p]=a[p][k]=c*kp-s*kq;a[k][q]=a[q][k]=s*kp+c*kq;}}
    std::array<double,3> e{};for(int i=0;i<n;++i)e[i]=a[i][i];std::sort(e.begin(),e.begin()+n);return e;
}
int main(int argc,char** argv){try{const std::filesystem::path dir=argc>1?argv[1]:"results";std::filesystem::create_directories(dir);std::ofstream out(dir/"conditioning.csv");out<<std::setprecision(17);
    out<<"case,hz,mass_ratio,stiffness_factor,raw_condition,scaled_condition,min_scaled_eigenvalue,gs_sweeps,scaled_residual,scaling_position_difference,dual_gradient_relative_error,max_constitutive_residual,status\n";
    double max_difference=0,max_dual=0;int count=0;
    for(const auto& original:cases())for(int hz:{240,960})for(double ratio:{1.,10.,100.,10000.})for(double stiffness:{1.,100.,10000.}){
        auto cse=original;for(int i=0;i<cse.model.count;++i){if(cse.model.alpha[i]>0)cse.model.alpha[i]/=stiffness;if(i==cse.model.count-1)cse.model.mass[i]*=ratio;}
        const auto c=context(cse.model,cse.initial,1./hz,energy(cse.model,cse.initial));Iterate z{c.predictor};z.lambda[0]=.001;if(cse.model.count==2)z.lambda[1]=-.002;
        const auto r=evaluate(c,z);const int n=cse.model.count+1;Matrix normalized{};std::array<double,3> scale{},rhs{},d{},u{};
        bool zero=false;for(int i=0;i<n;++i){scale[i]=std::sqrt(r.S[i][i]);if(!(scale[i]>0))zero=true;}
        if(zero){out<<cse.name<<','<<hz<<','<<ratio<<','<<stiffness<<",inf,inf,0,0,inf,0,0,"<<std::abs(r.rC[0])<<",zero_energy_gradient\n";++count;continue;}
        for(int i=0;i<n;++i)for(int j=0;j<n;++j)normalized[i][j]=r.S[i][j]/(scale[i]*scale[j]);
        const auto raw_e=eigenvalues(r.S,n),unit_e=eigenvalues(normalized,n);const double raw_condition=raw_e[0]>0?raw_e[n-1]/raw_e[0]:std::numeric_limits<double>::infinity(),scaled_condition=unit_e[0]>0?unit_e[n-1]/unit_e[0]:std::numeric_limits<double>::infinity();
        for(int i=0;i<n;++i)rhs[i]=-(i==n-1?r.H:r.rC[i]);
        Vec reconstructed{};for(int i=0;i<cse.model.count;++i){reconstructed[i]=(z.x[i]-(c.old.x[i]+cse.model.gravity*(c.h*c.h)))*(cse.model.mass[i]/(c.h*c.h));
            for(int j=0;j<cse.model.count;++j)if(cse.model.alpha[j]>0)reconstructed[i]+=r.J[j][i]*(-z.lambda[j]/(c.h*c.h)+r.rC[j]/cse.model.alpha[j]);}
        Vec error{};for(int i=0;i<cse.model.count;++i)error[i]=reconstructed[i]-r.g[i];const double dual_error=mass_norm(cse.model,error)/(1+mass_norm(cse.model,r.g));max_dual=std::max(max_dual,dual_error);
        for(int pass=1;pass<=16;++pass){for(int i=0;i<n;++i){double b=rhs[i],bs=rhs[i]/scale[i];for(int j=0;j<n;++j)if(i!=j){b-=r.S[i][j]*d[j];bs-=normalized[i][j]*u[j];}d[i]=b/r.S[i][i];u[i]=bs;}
            if(pass!=1&&pass!=4&&pass!=16)continue;
            std::array<double,3> sd{};for(int i=0;i<n;++i)sd[i]=u[i]/scale[i];const Vec dx=displacement(cse.model,r,d),sx=displacement(cse.model,r,sd);Vec delta{};for(int i=0;i<cse.model.count;++i)delta[i]=dx[i]-sx[i];const double difference=mass_norm(cse.model,delta)/(1+mass_norm(cse.model,dx));max_difference=std::max(max_difference,difference);
            double rn=0,bn=0;for(int i=0;i<n;++i){double v=-rhs[i]/scale[i];for(int j=0;j<n;++j)v+=normalized[i][j]*u[j];rn+=v*v;bn+=rhs[i]*rhs[i]/(scale[i]*scale[i]);}
            double constitutive=0;for(int j=0;j<cse.model.count;++j)constitutive=std::max(constitutive,std::abs(r.rC[j]));
            out<<cse.name<<','<<hz<<','<<ratio<<','<<stiffness<<','<<raw_condition<<','<<scaled_condition<<','<<unit_e[0]<<','<<pass<<','<<std::sqrt(rn/(1+bn))<<','<<difference<<','<<dual_error<<','<<constitutive<<",valid\n";++count;
        }
    }
    if(max_difference>1e-9||max_dual>1e-10)throw std::runtime_error("conditioning invariance check failed");std::cout<<count<<" conditioning rows; max scaled-position discrepancy "<<max_difference<<"; residual-complete dual gradient error "<<max_dual<<'\n';
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
