#include "common.hpp"
#include "el.hpp"
#include "/home/blacknoir/Projects/constrained-xpbd/experiments/research_agenda/solver.hpp"
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <random>

namespace {
using namespace compose;
using V2=Eigen::Vector2d;
struct Stat {int count=0;double abs=0,normalized=0,min=1e300,max=-1e300;void add(double d,double scale){++count;abs=std::max(abs,std::abs(d));normalized=std::max(normalized,std::abs(d)/std::max(1e-30,scale));min=std::min(min,d);max=std::max(max,d);}};
std::map<std::string,Stat> stats,checks;
void check(std::string id,double actual,double expected,double scale=1){checks[id].add(actual-expected,std::max(scale,1e-12));}
V2 dir(double a){return V2(std::cos(a),std::sin(a));}
V2 tangent(const V2&x){return V2(-x[1],x[0]);}
const std::vector<std::string> laws={"start_spatial","end_spatial","full_force_average","separate_average","dg_start","dg_end","dg_average","frozen_start_power","kinetic_start_force","kinetic_end_force"};
void fixed_work(const std::filesystem::path&out){
    std::mt19937_64 rng(42);std::uniform_real_distribution<double> u(-1,1);auto r=[&](){return u(rng);};
    std::ofstream probes(out/"fixed_work_incidence.csv");probes<<std::setprecision(17)<<"environment,probe,edge,body,law,work,edge_potential_change,edge_total_work,edge_storage_defect\n";
    for(const std::string env:{"scalar","radial","rod","elastic","chain"})for(int sample=0;sample<1000;++sample){
        int count=env=="chain"?2:1;bool scalar=env=="scalar",rod=env=="rod";double h=.001+.049*(r()+1)/2,k=100;
        std::vector<V2> old(count),now(count),vn(count);std::vector<double> qn(count),q(count);
        for(int j=0;j<count;++j){double ro=rod?1:1+.4*r(),rn=rod?1:1+.4*r(),a=scalar?0:3*r(),b=scalar?0:a+2.8*r();V2 eo=ro*dir(a),en=rn*dir(b);old[j]=eo+(j?old[j-1]:V2::Zero());now[j]=en+(j?now[j-1]:V2::Zero());vn[j]=rod?r()*tangent(eo):V2(r(),scalar?0:r());qn[j]=rod?-(1+20*(r()+1)/2):-k*(ro-1);q[j]=rod?-(1+20*(r()+1)/2):-k*(rn-1);}
        for(int j=0;j<count;++j){V2 eo=old[j]-(j?old[j-1]:V2::Zero()),en=now[j]-(j?now[j-1]:V2::Zero());double ro=eo.norm(),rn=en.norm(),du=rod?0:.5*k*((rn-1)*(rn-1)-(ro-1)*(ro-1));V2 no=eo/ro,nn=en/rn,dg=(eo+en)/(ro+rn);
            std::vector<std::vector<double>> w(laws.size(),std::vector<double>(count,0));
            for(int body:{j,j-1}){if(body<0)continue;double sign=body==j?1:-1;V2 dx=now[body]-old[body],fo=sign*qn[j]*no,fe=sign*q[j]*nn;
                w[0][body]=fo.dot(dx);w[1][body]=fe.dot(dx);w[2][body]=.5*(fo+fe).dot(dx);w[3][body]=sign*.25*(qn[j]+q[j])*(no+nn).dot(dx);w[4][body]=sign*qn[j]*dg.dot(dx);w[5][body]=sign*q[j]*dg.dot(dx);w[6][body]=sign*.5*(qn[j]+q[j])*dg.dot(dx);w[7][body]=h*fo.dot(vn[body]);w[8][body]=.5*fo.dot(dx+h*vn[body]);w[9][body]=.5*fe.dot(dx+h*vn[body]);
                V2 perturb(.031,-.047);check("live_old_force_target_sensitivity",fo.dot(dx+perturb)-w[0][body],fo.dot(perturb),1+std::abs(w[0][body]));check("old_force_kinetic_target_sensitivity",.5*fo.dot(dx+perturb+h*vn[body])-w[8][body],.5*fo.dot(perturb),1+std::abs(w[8][body]));check("frozen_power_target_sensitivity",h*fo.dot(vn[body])-w[7][body],0,1+std::abs(w[7][body]));
            }
            for(size_t l=0;l<laws.size();++l){double sum=0;for(double a:w[l])sum+=a;double defect=du+sum,scale=1+std::abs(du)+std::abs(sum);stats[env+","+laws[l]].add(defect,scale);
                for(int body:{j,j-1})if(body>=0)probes<<env<<','<<sample<<','<<j<<','<<body<<','<<laws[l]<<','<<w[l][body]<<','<<du<<','<<sum<<','<<defect<<'\n';
                if(l==6)check(env+"_DG_average_storage",sum,-du,scale);
                if(rod&&(l==3||l==4||l==5||l==6||l==7))check("rod_zero_"+laws[l],sum,0,scale);
            }
            V2 de=en-eo;check(env+"_discrete_geometry",dg.dot(de),rn-ro,1+ro+rn);
            double cosine=no.dot(nn);check(env+"_average_direction",.5*(no+nn).dot(de),(rn-ro)*(1+cosine)/2,1+ro+rn);
            if(rod){check("rod_full_force_average_formula",w[2][j],.5*(q[j]-qn[j])*(1-cosine),1+std::abs(q[j])+std::abs(qn[j]));}
            if(scalar){double dc=rn-ro;check("scalar_start_storage_square",du+w[0][j],.5*k*dc*dc,1+k*dc*dc);check("scalar_end_storage_square",du+w[1][j],-.5*k*dc*dc,1+k*dc*dc);check("scalar_full_vs_separate",w[2][j],w[3][j],1+std::abs(w[2][j]));check("scalar_separate_vs_DG",w[3][j],w[6][j],1+std::abs(w[3][j]));check("scalar_endpoint_vs_DG",w[1][j],w[5][j],1+std::abs(w[1][j]));}
            if(env=="chain"&&j==1){check("chain_incidence_transfer_DG",w[6][0]+w[6][1],.5*(qn[j]+q[j])*dg.dot(de),1+std::abs(du));}
        }
        if(env=="rod"||env=="elastic"||env=="chain")for(int i=0;i<count;++i){V2 grav(0,-9.81),dx=now[i]-old[i];double du=-grav.dot(dx),w=grav.dot(dx);check("gravity_spatial_exact",du+w,0,1+std::abs(du));}
    }
    std::ofstream f(out/"fixed_work_summary.csv");f<<std::setprecision(17)<<"environment,law,count,max_absolute_storage_defect,max_normalized_storage_defect,min_signed_defect,max_signed_defect\n";for(auto&[id,s]:stats)f<<id<<','<<s.count<<','<<s.abs<<','<<s.normalized<<','<<s.min<<','<<s.max<<'\n';
}
void kinetic(const std::filesystem::path&out){
    std::mt19937_64 rng(42);std::uniform_real_distribution<double>u(-1,1);auto rv=[&](){return V2(u(rng),u(rng));};
    std::ofstream f(out/"kinetic_impulse_probes.csv");f<<std::setprecision(17)<<"probe,mass,closed_defect,uncorrected_ledger_after_unlogged_velocity_change,midpoint_dot_momentum_residual,identity_error,AB_work_A,BA_work_A,allocation_difference,cross_term,total_order_difference\n";
    for(int sample=0;sample<1000;++sample){double m=1+.8*u(rng);V2 vo=rv(),ia=rv(),ib=rv(),v=vo+(ia+ib)/m,vm=.5*(v+vo);double dk=.5*m*(v.squaredNorm()-vo.squaredNorm()),w=vm.dot(ia+ib);check("kinetic_closed_momentum",dk,w,1+std::abs(dk)+std::abs(w));
        V2 altered=v+.2*rv(),altered_mid=.5*(altered+vo),res=m*(altered-vo)-ia-ib;double altered_dk=.5*m*(altered.squaredNorm()-vo.squaredNorm()),altered_w=altered_mid.dot(ia+ib),defect=altered_dk-altered_w,prediction=altered_mid.dot(res);check("kinetic_unlogged_correction_identity",defect,prediction,1+std::abs(altered_dk)+std::abs(altered_w));check("kinetic_stale_recorded_work_correction",altered_dk-w,.5*m*(altered.squaredNorm()-v.squaredNorm()),1+std::abs(altered_dk)+std::abs(w));
        V2 va=vo+ia/m,vb=vo+ib/m;double waAB=ia.dot((vo+va)/2),wbAB=ib.dot((va+v)/2),wbBA=ib.dot((vo+vb)/2),waBA=ia.dot((vb+v)/2),cross=ia.dot(ib)/m;
        check("sequential_AB_telescope",waAB+wbAB,dk,1+std::abs(dk));check("sequential_BA_telescope",waBA+wbBA,dk,1+std::abs(dk));check("sequential_order_cross_term",waBA-waAB,cross,1+std::abs(cross));check("global_midpoint_half_cross_allocation",vm.dot(ia)-waAB,.5*cross,1+std::abs(cross));
        f<<sample<<','<<m<<','<<dk-w<<','<<defect<<','<<prediction<<','<<defect-prediction<<','<<waAB<<','<<waBA<<','<<waBA-waAB<<','<<cross<<','<<waAB+wbAB-waBA-wbBA<<'\n';
        // Two-mass chain: edge zero acts on body zero; edge one acts equally
        // and oppositely on the two bodies. Audit each kinetic row separately.
        V2 v0=rv(),v1=rv(),j0=rv(),j1=rv(),end0=v0+(j0-j1)/m,end1=v1+j1;double dk0=.5*m*(end0.squaredNorm()-v0.squaredNorm()),dk1=.5*(end1.squaredNorm()-v1.squaredNorm());
        check("chain_body0_impulse_balance",dk0,(end0+v0).dot(j0-j1)/2,1+std::abs(dk0));check("chain_body1_impulse_balance",dk1,(end1+v1).dot(j1)/2,1+std::abs(dk1));
    }
}
void scalar_analytic(const std::filesystem::path&out){
    std::ofstream f(out/"scalar_release_analytic.csv");f<<std::setprecision(17)<<"h,old_extension,BE_extension,energy_gradient_at_BE,energy_defect_at_BE,frozen_energy_root,moving_energy_root,exact_extension,moving_displacement_over_exact,velocity_over_exact,first_step_frozen_LE_displacement\n";
    for(double h:{1./30,1./60,1./120,1./240,1./480,1./960}){double k=100,m=1,a=.25,s=k*h*h/m,be=a/(1+s),moving=a*(1-s)/(1+s),exact=a*std::cos(std::sqrt(k/m)*h),grad=m*(be-a)/(h*h)+k*be,def=.5*m*(be-a)*(be-a)/(h*h)+.5*k*be*be-.5*k*a*a;
        check("release_BE_gradient_zero",grad,0,1+k*a);check("release_moving_energy_root",.5*m*(moving-a)*(moving-a)/(h*h)+.5*k*moving*moving,.5*k*a*a,1+k*a*a);
        auto c=make_case("scalar_release");Config cfg;cfg.schedule="LE";cfg.treatment="raw";cfg.keep_rx=false;cfg.relinearize=false;cfg.iterations=1;auto result=el_step(c.model,c.initial,h,cfg,energy(c.model,c.initial));double step=result.state.x[0]-c.initial.x[0];check("release_frozen_LE_freeze",step,0,1);
        f<<h<<','<<a<<','<<be<<','<<grad<<','<<def<<','<<a<<','<<moving<<','<<exact<<','<<(moving-a)/(exact-a)<<','<<((moving-a)/h)/(-std::sqrt(k/m)*a*std::sin(std::sqrt(k/m)*h))<<','<<step<<'\n';
    }
}
void el_diagnostic(const std::filesystem::path&out){
    std::ofstream f(out/"EL_diagnostic_short.csv");f<<std::setprecision(17)<<"case,order,treatment,keep_rx,refresh,sweeps,h,steps,valid,singular_visits,capped_sweeps,skipped_sweeps,max_normalized_energy_error,max_rx,max_step_displacement,min_presweep_energy_gradient,first_step_motion,final_motion_from_initial\n";
    for(const auto&name:case_names())for(std::string order:{"EL","LE"})for(std::string treatment:{"raw","bounded"})for(bool keep:{false,true})for(bool refresh:{false,true})for(int sweeps:{1,4,16}){
        auto problem=make_case(name);const auto&m=problem.model;State s=problem.initial;double h=1./120,initialE=energy(m,s),scale=std::max(1e-8,initialE-minimum_energy(m)),maxE=0,maxR=0,maxDx=0,minG=1e300,first=0;int steps=0;bool valid=true;traditional::Diagnostics diagnostic;
        Config cfg;cfg.schedule=order;cfg.treatment=treatment;cfg.keep_rx=keep;cfg.relinearize=refresh;cfg.iterations=sweeps;
        for(int step=0;step<24;++step){auto tm=el_detail::model(m);auto old=el_detail::state(m,s);auto c=traditional::context(tm,old,h,energy(m,s));traditional::Iterate z{c.predictor};for(int iteration=0;iteration<sweeps;++iteration){auto rows=traditional::evaluate(c,z);minG=std::min(minG,std::sqrt(traditional::inner(rows.g,rows.g,m.count)));if(!el_detail::sweep(c,z,cfg,diagnostic)){valid=false;break;}}
            traditional::State ts;ts.x=z.x;for(int i=0;i<m.count;++i)ts.v[i]=(z.x[i]-old.x[i])/h;State next=el_detail::state(m,ts);
            if(!finite(next)||next.x.norm()>1e6||next.v.norm()>1e8||(m.dim==1&&next.x[0]<=0)){valid=false;break;}
            double dx=(next.x-s.x).norm();maxDx=std::max(maxDx,dx);if(step==0)first=dx;s=next;++steps;maxE=std::max(maxE,std::abs(energy(m,s)-initialE)/scale);auto rr=traditional::evaluate(c,z);maxR=std::max(maxR,traditional::mass_norm(tm,traditional::minv(tm,rr.R)));if(!valid)break;
        }
        f<<name<<','<<order<<','<<treatment<<','<<keep<<','<<refresh<<','<<sweeps<<','<<h<<','<<steps<<','<<valid<<','<<diagnostic.singular<<','<<diagnostic.caps<<','<<diagnostic.skips<<','<<maxE<<','<<maxR<<','<<maxDx<<','<<minG<<','<<first<<','<<(s.x-problem.initial.x).norm()<<'\n';
    }
    // Direct compatibility checks against the historical refreshed raw
    // retained-residual implementation, and frozen zero-residual first sweeps.
    for(auto name:case_names())for(bool EL:{false,true}){
        auto p=make_case(name);auto tm=el_detail::model(p.model);auto ts=el_detail::state(p.model,p.initial);auto c=traditional::context(tm,ts,1./120,energy(p.model,p.initial));traditional::Iterate z{c.predictor};Config cfg;cfg.schedule=EL?"EL":"LE";cfg.treatment="raw";cfg.keep_rx=true;cfg.relinearize=true;traditional::Diagnostics dd;el_detail::sweep(c,z,cfg,dd);agenda::Report rr;auto old=agenda::energy_step(c,EL,true,0,0,rr);for(int i=0;i<tm.count;++i)check("retained_refreshed_matches_agenda",traditional::norm(z.x[i]-old.x[i]),0,1);
        traditional::Iterate a{c.predictor},b{c.predictor};cfg.relinearize=false;el_detail::sweep(c,a,cfg,dd);cfg.keep_rx=false;el_detail::sweep(c,b,cfg,dd);for(int i=0;i<tm.count;++i)check("frozen_first_zero_rx_matches_legacy",traditional::norm(a.x[i]-b.x[i]),0,1);
    }
}
void el_release_trace(const std::filesystem::path&out){
    std::ofstream f(out/"EL_scalar_release_visit_trace.csv");f<<std::setprecision(17)<<"order,step,h,old_position,old_velocity,predictor,H_at_predictor,energy_gradient,energy_denominator,energy_delta_eta,raw_energy_position_correction,material_delta_lambda,final_position,final_velocity,final_energy_change\n";
    for(bool EL:{true,false}){auto p=make_case("scalar_release");auto tm=el_detail::model(p.model);auto old=el_detail::state(p.model,p.initial);double h=1./120;for(int step=0;step<24;++step){auto c=traditional::context(tm,old,h,traditional::energy(tm,old));traditional::Iterate z{c.predictor};auto r=traditional::evaluate(c,z);double dl=0,de=0;if(!EL)dl=-r.rC[0]/r.S[0][0];if(r.S[1][1]>0)de=(-r.H-r.S[1][0]*dl)/r.S[1][1];if(EL)dl=(-r.rC[0]-r.S[0][1]*de)/r.S[0][0];traditional::Options o;o.schedule=EL?traditional::Schedule::EL:traditional::Schedule::LE;traditional::Diagnostics d;bool ok=traditional::update(c,z,o,d);traditional::State next;next.x=z.x;next.v[0]=(z.x[0]-old.x[0])/h;f<<(EL?"EL":"LE")<<','<<step<<','<<h<<','<<old.x[0].x<<','<<old.v[0].x<<','<<c.predictor[0].x<<','<<r.H<<','<<r.g[0].x<<','<<r.S[1][1]<<','<<de<<','<<r.g[0].x*de/tm.mass[0]<<','<<dl<<','<<next.x[0].x<<','<<next.v[0].x<<','<<traditional::energy(tm,next)-traditional::energy(tm,old)<<'\n';if(!ok||next.x[0].x<=0)break;old=next;}}
}
}
int main(int argc,char**argv){
    std::filesystem::path out=argc>1?argv[1]:"audits_results";std::filesystem::create_directories(out);fixed_work(out);kinetic(out);scalar_analytic(out);el_diagnostic(out);el_release_trace(out);
    std::ofstream f(out/"analytic_checks.csv");f<<std::setprecision(17)<<"identity,count,max_absolute_error,max_normalized_error,pass\n";bool pass=true;for(auto&[id,s]:checks){bool good=s.normalized<1e-11;pass=pass&&good;f<<id<<','<<s.count<<','<<s.abs<<','<<s.normalized<<','<<good<<'\n';}
    std::cout<<"fixed work environments=5 probes_per_environment=1000 laws=10 seed=42\n"<<"analytic identity groups="<<checks.size()<<" all_pass="<<pass<<" output="<<out<<'\n';return pass?0:1;
}
