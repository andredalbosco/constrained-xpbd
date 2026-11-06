// Reuse the existing physical cases and independently checked RK4 references.
#define main correction_study_main
#include "study.cpp"
#undef main

namespace {
enum class ScalarMethod {Cold,Warm,PredictorOnly,LambdaOnly};
const char* scalar_name(ScalarMethod s) {
    switch(s){case ScalarMethod::Cold:return "cold_xpbd";case ScalarMethod::Warm:return "consistent_warm";
    case ScalarMethod::PredictorOnly:return "old_force_predictor_only";case ScalarMethod::LambdaOnly:return "lambda_only";}
    return "unknown";
}
struct Scalar {double q,v,lambda;};
Scalar scalar_step(Scalar s,double h,ScalarMethod method) {
    const double m=1,k=100,at=1/(k*h*h),A=1/m+at;
    double q=s.q+h*s.v,lambda=0;
    if(method==ScalarMethod::Warm||method==ScalarMethod::PredictorOnly)q+=s.lambda/m;
    if(method==ScalarMethod::Warm||method==ScalarMethod::LambdaOnly)lambda=s.lambda;
    if(method!=ScalarMethod::PredictorOnly) {
        const double dl=-(q+at*lambda)/A;q+=dl/m;lambda+=dl;
    } else lambda=-h*h*k*q;
    return {q,(q-s.q)/h,lambda};
}
void scalar_runs(const std::filesystem::path& dir) {
    std::ofstream out(dir/"scalar.csv");out<<std::setprecision(17);
    out<<"case,method,h,h_omega,steps,requested,energy_ratio_min,energy_ratio_max,energy_ratio_final,max_error_from_BE,modified_energy_relative_error\n";
    for(bool moving:{false,true})for(double h:{1./240,1./480,1./960,.1,.19,.21})
    for(auto method:{ScalarMethod::Cold,ScalarMethod::Warm,ScalarMethod::PredictorOnly,ScalarMethod::LambdaOnly}) {
        Scalar s{.25,moving?.7:0,-h*h*100*.25},be=s;
        const double E0=.5*s.v*s.v+50*s.q*s.q,I0=E0-50*h*s.q*s.v;
        double minE=1,maxE=1,err=0,modified=0;int done=0,steps=std::lround(10/h);
        for(int j=0;j<steps;++j) {
            s=scalar_step(s,h,method);be=scalar_step(be,h,ScalarMethod::Cold);
            const double E=.5*s.v*s.v+50*s.q*s.q;
            if(!std::isfinite(E)||E/E0>1e20)break;
            ++done;minE=std::min(minE,E/E0);maxE=std::max(maxE,E/E0);
            err=std::max(err,std::hypot(s.q-be.q,(s.v-be.v)/10));
            modified=std::max(modified,std::abs((E-50*h*s.q*s.v)-I0)/std::abs(I0));
        }
        if(method==ScalarMethod::Warm&&err>1e-10)throw std::runtime_error("consistent scalar warm start differs from ordinary XPBD");
        if(method==ScalarMethod::PredictorOnly&&h*10<2&&modified>1e-10)throw std::runtime_error("symplectic modified energy check failed");
        out<<(moving?"moving":"release")<<','<<scalar_name(method)<<','<<h<<','<<h*10<<','<<done<<','<<steps<<','<<minE<<','<<maxE<<','<<(.5*s.v*s.v+50*s.q*s.q)/E0<<','<<err<<','<<modified<<'\n';
    }
}
void first_energy_step(const std::filesystem::path& dir) {
    std::ofstream out(dir/"first_energy_step.csv");out<<std::setprecision(17);
    out<<"h,a,q_old,q_warm,energy_residual,energy_gradient,energy_dx,exact_dx_formula,post_energy_ratio\n";
    for(double h:{1./240,1./480,1./960}) {
        const double qold=.25,a=100*h*h,q=(1-a)*qold;
        const double H=.5/(h*h)*(q-qold)*(q-qold)+50*q*q-50*qold*qold;
        const double g=(q-qold)/(h*h)+100*q,dx=-H/g;
        const double expected=-(1-a)*qold/2;
        if(std::abs(dx-expected)>1e-8)throw std::runtime_error("warm energy first-step formula mismatch");
        const double endpoint=q+dx,E=.5/(h*h)*(endpoint-qold)*(endpoint-qold)+50*endpoint*endpoint;
        out<<h<<','<<a<<','<<qold<<','<<q<<','<<H<<','<<g<<','<<dx<<','<<expected<<','<<E/(50*qold*qold)<<'\n';
    }
}
struct WarmConfig {std::string name;bool warm=false,retain=false;se::Options options;};
void physical_runs(const std::filesystem::path& dir) {
    std::ofstream out(dir/"physical.csv");out<<std::setprecision(17);
    out<<"case,method,substeps,h,completed,requested,max_energy_percent,rms_position,max_rx,max_lambda,status\n";
    std::vector<WarmConfig> configs;
    se::Options baseline;baseline.schedule=se::Schedule::Baseline;
    configs.push_back({"cold_baseline",false,false,baseline});
    configs.push_back({"warm_baseline_drop_rx",true,false,baseline});
    configs.push_back({"warm_baseline_retain_rx",true,true,baseline});
    for(auto s:{se::Schedule::EL,se::Schedule::LE,se::Schedule::Joint})for(bool warm:{false,true})
    for(auto p:{se::Policy::Raw,se::Policy::RepairCap}) {
        se::Options o;o.schedule=s;o.policy=p;
        configs.push_back({std::string(warm?"warm_":"cold_")+se::name(s)+"_"+se::name(p),warm,false,o});
    }
    int runs=0,stops=0;
    for(const auto& c:cases())for(int substeps:{4,8,16}) {
        const double h=1./(60*substeps);const int steps=600*substeps;
        const auto ref=reference(c,h,steps,16),fine=reference(c,h,steps,32);
        for(int i=0;i<=steps;++i)if(norm(ref[i].x-fine[i].x)>1e-7)throw std::runtime_error("reference convergence failed");
        for(const auto& cfg:configs) {
            State state=c.initial;const double E0=energy(c.model,state),es=scale(c);
            // Seed the physical elastic force. A hard rod needs an estimated
            // reaction; start it at zero and subsequently carry solved lambda.
            double previous=c.model.pendulum?0:-h*h*(norm(state.x)-c.model.length)/c.model.alpha;
            double Emax=0,x2=0,Rmax=0,lmax=0;int done=0;std::string status="completed";
            for(int j=0;j<steps;++j) {
                const Context ctx=context(c.model,state,h,0,E0);Iterate z{ctx.predictor};
                if(cfg.warm){z.lambda=previous;z.x+=(state.x/norm(state.x))*(previous/c.model.mass);}
                bool ok=true;se::Diagnostics d;
                if(cfg.retain) {
                    const Rows r=evaluate(ctx,z);const V2 b=r.R*(-1/c.model.mass);
                    const double dl=(-r.rC-dot(r.n,b))/r.A;
                    z.x+=b+r.n*(dl/c.model.mass);z.lambda+=dl;ok=r.valid&&finite(z.x);
                } else ok=se::update(ctx,z,cfg.options,j,d);
                const State next{z.x,(z.x-state.x)/h};const double E=energy(c.model,next);
                if(!ok||!finite(next.x)||!finite(next.v)||!std::isfinite(E)){status="undefined";break;}
                if(norm(next.x)>1e6||norm(next.v)>1e6||std::abs(E-E0)>1e12*es||std::abs(z.lambda)>1e18||std::abs(z.eta)>1e18){status="safety_stop";break;}
                ++done;Emax=std::max(Emax,std::abs(E-E0));x2+=dot(next.x-ref[j+1].x,next.x-ref[j+1].x);
                Rmax=std::max(Rmax,norm(evaluate(ctx,z).R)/(h*h));lmax=std::max(lmax,std::abs(z.lambda));
                state=next;previous=z.lambda;
            }
            out<<c.name<<','<<cfg.name<<','<<substeps<<','<<h<<','<<done<<','<<steps<<','<<100*Emax/es<<','<<std::sqrt(x2/std::max(1,done))<<','<<Rmax<<','<<lmax<<','<<status<<'\n';
            if(cfg.options.schedule==se::Schedule::Baseline&&done!=steps)
                throw std::runtime_error("vanilla warm-start validation did not complete");
            ++runs;stops+=done<steps;
        }
    }
    std::cout<<runs<<" physical runs; "<<stops<<" explicit stops\n";
}
}
int main(int argc,char** argv) {
    try {
        const std::filesystem::path dir=argc>1?argv[1]:"warm_start_results";
        std::filesystem::create_directories(dir);
        scalar_runs(dir);first_energy_step(dir);physical_runs(dir);
        std::cout<<"PASS: affine warm-start equivalence, symplectic invariant, first energy increment, reference convergence\n";
        return 0;
    }catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}
