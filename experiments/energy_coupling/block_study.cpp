// Direct energy-row coupling: frozen-system control and turning-point checks.
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>

void check(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

int main() {
    try {
        std::cout << std::setprecision(17)
            << "alpha_tilde,gradient_angle,normalized_condition,gs_spectral_radius,"
               "joint_displacement_norm,joint_linear_residual,gs_10_displacement_norm,gs_10_linear_residual,status\n";
        for (double alpha : {0.0, 1e-6, 0.01, 36.0}) {
            for (double theta : {1.0, 0.1, 0.01, 0.0001, 0.0}) {
                // M=I, J=(1,0), g=(cos(theta),sin(theta)), residuals (.1,-.2).
                const double c=std::cos(theta), s=std::sin(theta), a=1+alpha;
                const double d=c*c+s*s;
                const double schur=s*s+alpha*c*c/a;
                const double rho=std::abs(c)/std::sqrt(a*d);
                const double condition=rho<1 ? (1+rho)/(1-rho) : std::numeric_limits<double>::infinity();
                double dl=0,de=0;
                for(int iteration=0;iteration<10;++iteration) {
                    dl=(-.1-c*de)/a;
                    de=(.2-c*dl)/d;
                }
                const double gs_norm=std::hypot(dl+c*de,s*de);
                const double gs_residual=std::hypot(a*dl+c*de+.1,c*dl+d*de-.2);
                if(schur==0) {
                    std::cout<<alpha<<','<<theta<<','<<condition<<','<<rho*rho
                        <<",nan,nan,"<<gs_norm<<','<<gs_residual<<",singular\n";
                    continue;
                }
                const double energy=(.2+.1*c/a)/schur;
                const double material=(-.1-c*energy)/a;
                const double dx=material+c*energy,dy=s*energy;
                const double defect=std::hypot(dx+alpha*material+.1,c*dx+s*dy-.2);
                check(std::isfinite(defect) && defect<1e-7,"joint frozen-system solve failed");
                std::cout<<alpha<<','<<theta<<','<<condition<<','<<rho*rho<<','
                    <<std::hypot(dx,dy)<<','<<defect<<','<<gs_norm<<','<<gs_residual<<",ok\n";
            }
        }

        // A compliant spring released from rest has a frozen, fully feasible branch.
        const double h=1.0/60, k=100, C=.25, alpha=1/(k*h*h), grad=k*C;
        const double a=1+alpha, b=grad, d=grad*grad;
        const double schur=alpha*d/a;
        const double energy=(b*C/a)/schur;
        const double material=(-C-b*energy)/a;
        const double dx=material+grad*energy;
        check(std::abs(dx)<1e-16 && std::abs(C+alpha*material)<1e-14,
              "resting spring frozen branch not reproduced");
        check(std::abs(energy-h*h)<1e-17,"resting spring energy multiplier is not h squared");

        // Material-first reaches the energy minimum, whose gradient is zero.
        const double projected_C=C/(1+k*h*h);
        const double projected_dx=projected_C-C;
        const double energy_gradient=projected_dx/(h*h)+k*projected_C;
        const double H=.5*projected_dx*projected_dx/(h*h)+.5*k*projected_C*projected_C-.5*k*C*C;
        check(std::abs(energy_gradient)<1e-11 && H<0,"sequential spring singularity not reproduced");
        std::cerr<<"rest spring: joint dx="<<dx<<", eta/h^2="<<energy/(h*h)
            <<"; after material projection gradH="<<energy_gradient<<", H="<<H<<'\n';
        return 0;
    } catch(const std::exception& e) {
        std::cerr<<e.what()<<'\n'; return 1;
    }
}
