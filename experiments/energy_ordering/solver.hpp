#pragma once

#include <algorithm>
#include <cmath>
#include <limits>
#include <string>

// Direct endpoint energy experiments. No reservoir or velocity correction.
namespace energy_ordering {
struct V2 {
    double x = 0, y = 0;
    V2 operator+(V2 b) const { return {x+b.x, y+b.y}; }
    V2 operator-(V2 b) const { return {x-b.x, y-b.y}; }
    V2 operator*(double s) const { return {x*s, y*s}; }
    V2 operator/(double s) const { return {x/s, y/s}; }
    V2& operator+=(V2 b) { x+=b.x; y+=b.y; return *this; }
};
inline double dot(V2 a, V2 b) { return a.x*b.x+a.y*b.y; }
inline double norm(V2 a) { return std::hypot(a.x,a.y); }
inline V2 tangent(V2 a) { return {-a.y,a.x}; }
inline bool finite(V2 a) { return std::isfinite(a.x)&&std::isfinite(a.y); }

struct Model {
    bool pendulum = false;
    double mass = 1, length = 1, alpha = .01;
    V2 acceleration{};
};
struct State { V2 x{1.25,0}, v{}; };
struct Context {
    Model model;
    State old;
    V2 predictor;
    double h = 1./60, beta = 0, target = 0;
};
struct Iterate { V2 x; double lambda = 0, eta = 0; };
enum class Method { Baseline, Joint, Immediate, Frozen, Reverse, Symmetric };
inline const char* name(Method m) {
    switch(m) {
    case Method::Baseline: return "baseline";
    case Method::Joint: return "joint";
    case Method::Immediate: return "immediate";
    case Method::Frozen: return "frozen";
    case Method::Reverse: return "reverse";
    case Method::Symmetric: return "symmetric";
    }
    return "unknown";
}
struct Options {
    Method method = Method::Baseline;
    int sweeps = 5, inner = 1;
    bool retain_rx = false;
};
struct Counters {
    long evaluations = 0, scalar_solves = 0, block_solves = 0;
    long singular = 0, incompatible = 0, rank_deficient_blocks = 0;
};
struct Rows {
    V2 n{}, g{}, R{};
    double C=0,H=0,rC=0,rE=0,A=0,B=0,D=0,schur=0;
    double at=0,bt=0,p=0,q=0,gradient_scale=1;
    bool valid = true;
};
inline double elastic(const Model& m,V2 x) {
    const double c=norm(x)-m.length;
    return m.pendulum ? 0 : c*c/(2*m.alpha);
}
inline double energy(const Model& m,const State& s) {
    return .5*m.mass*dot(s.v,s.v)+elastic(m,s.x)
        -m.mass*dot(m.acceleration,s.x);
}
inline Context context(Model m,State s,double h,double beta,double target) {
    return {m,s,s.x+s.v*h+m.acceleration*(h*h),h,beta,target};
}
inline Rows evaluate(const Context& c,const Iterate& z) {
    Rows r;
    const auto& m=c.model;
    const double h2=c.h*c.h, radius=norm(z.x);
    if (!(radius>1e-14*m.length) || !finite(z.x)) { r.valid=false; return r; }
    r.n=z.x/radius;
    r.C=radius-m.length;
    const V2 u=m.pendulum ? V2{} : r.n*(r.C/m.alpha);
    const V2 inertia=(z.x-c.old.x)*(m.mass/h2);
    // Center gravity before subtraction: the rest predictor has exactly g=0.
    const V2 center=c.old.x+m.acceleration*h2;
    r.g=(z.x-center)*(m.mass/h2)+u;
    r.H=energy(m,{z.x,(z.x-c.old.x)/c.h})-c.target;
    r.at=m.pendulum ? 0 : m.alpha/h2;
    r.bt=c.beta/h2;
    r.rC=r.C+r.at*z.lambda;
    r.rE=r.H+r.bt*z.eta;
    r.R=(z.x-c.predictor)*m.mass-r.n*z.lambda-r.g*z.eta;
    r.p=dot(r.n,r.g);
    r.q=dot(tangent(r.n),r.g);
    r.A=1/m.mass+r.at;
    r.B=r.p/m.mass;
    r.D=dot(r.g,r.g)/m.mass+r.bt;
    // Nonnegative Schur complement; avoids D-B*B/A cancellation.
    r.schur=r.q*r.q/m.mass+r.bt+(r.at/r.A)*r.p*r.p/m.mass;
    r.gradient_scale=1+norm(inertia)+norm(u)+m.mass*norm(m.acceleration);
    r.valid=finite(r.R)&&finite(r.g)&&std::isfinite(r.H)&&std::isfinite(r.schur);
    return r;
}
inline double gradient_tolerance(const Rows& r) {
    return 64*std::numeric_limits<double>::epsilon()*r.gradient_scale;
}
inline bool energy_row_zero(const Rows& r) {
    return r.bt==0 && norm(r.g)<=gradient_tolerance(r);
}
inline double energy_tolerance(const Context& c,const Rows& r) {
    return 128*std::numeric_limits<double>::epsilon()
        *(1+std::abs(c.target)+std::abs(r.H)+std::abs(r.bt));
}

// One joint block or an inner scalar GS iteration of ONE frozen Schur system.
// The top residual base displacement is applied once, after the inner solve.
inline bool frozen_update(const Context& c,Iterate& z,const Options& o,Counters& count) {
    const Rows r=evaluate(c,z); ++count.evaluations;
    if (!r.valid) return false;
    const double m=c.model.mass;
    const V2 base=o.retain_rx ? r.R*(-1/m) : V2{};
    const double fC=-r.rC-dot(r.n,base), fE=-r.rE-dot(r.g,base);
    double dl=0,de=0;
    const bool zero=energy_row_zero(r);
    const bool dependent=r.at==0 && r.bt==0 && std::abs(r.q)<=gradient_tolerance(r);
    if (zero||dependent) ++count.rank_deficient_blocks;
    if (o.method==Method::Joint && (zero||dependent)) {
        ++count.singular;
        // A dependent hard row is compatible only if its reduced RHS vanishes.
        const double reduced=fE-r.B*fC/r.A;
        if (std::abs(reduced)>energy_tolerance(c,r)) {
            ++count.incompatible;
            return true; // Explicit skipped block, never denominator regularization.
        }
        dl=fC/r.A; // Compatible redundant row: choose zero energy increment.
        ++count.scalar_solves;
    } else if (o.method==Method::Joint) {
        ++count.block_solves;
        if (r.at==0 && r.bt==0) {
            // The two hard rows determine both position coordinates, independent
            // of R. Avoid cancellation between R and huge multiplier forces.
            const double normal=-r.rC, transverse=(-r.rE-r.p*normal)/r.q;
            const V2 dx=r.n*normal+tangent(r.n)*transverse;
            de=m*(transverse-dot(tangent(r.n),base))/r.q;
            dl=m*(normal-dot(r.n,base))-r.p*de;
            if (!finite(dx)||!std::isfinite(dl)||!std::isfinite(de)) return false;
            z.x+=dx; z.lambda+=dl; z.eta+=de;
            return true;
        }
        de=(fE-r.B*fC/r.A)/r.schur;
        dl=(fC-r.B*de)/r.A;
    } else {
        for (int i=0;i<o.inner;++i) {
            dl=(fC-r.B*de)/r.A; ++count.scalar_solves;
            if (zero) {
                ++count.singular;
                if (std::abs(fE-r.B*dl)>energy_tolerance(c,r)) ++count.incompatible;
            } else {
                de=(fE-r.B*dl)/r.D; ++count.scalar_solves;
            }
        }
    }
    V2 dx=base+(r.n*dl+r.g*de)/m;
    if (o.method==Method::Joint && !zero && !dependent) {
        // Evaluate the same dx without subtracting large parallel multiplier forces.
        const double normal=dot(r.n,base)+(fC+r.at*r.p*de)/(m*r.A);
        const double transverse=dot(tangent(r.n),base)+r.q*de/m;
        dx=r.n*normal+tangent(r.n)*transverse;
    }
    if (!finite(dx)||!std::isfinite(dl)||!std::isfinite(de)) return false;
    z.x+=dx; z.lambda+=dl; z.eta+=de;
    return true;
}
inline bool material_update(const Context& c,Iterate& z,Counters& count) {
    const Rows r=evaluate(c,z); ++count.evaluations;
    if (!r.valid) return false;
    const double dl=-r.rC/r.A;
    z.x+=r.n*(dl/c.model.mass); z.lambda+=dl; ++count.scalar_solves;
    return finite(z.x)&&std::isfinite(z.lambda);
}
inline bool energy_update(const Context& c,Iterate& z,Counters& count) {
    const Rows r=evaluate(c,z); ++count.evaluations;
    if (!r.valid) return false;
    if (energy_row_zero(r)) {
        ++count.singular;
        if (std::abs(r.rE)>energy_tolerance(c,r)) ++count.incompatible;
        return true;
    }
    const double de=-r.rE/r.D;
    z.x+=r.g*(de/c.model.mass); z.eta+=de; ++count.scalar_solves;
    return finite(z.x)&&std::isfinite(z.eta);
}
inline bool sweep(const Context& c,Iterate& z,const Options& o,Counters& count) {
    switch(o.method) {
    case Method::Baseline: return material_update(c,z,count);
    case Method::Joint:
    case Method::Frozen: return frozen_update(c,z,o,count);
    case Method::Immediate:
        return material_update(c,z,count)&&energy_update(c,z,count);
    case Method::Reverse:
        return energy_update(c,z,count)&&material_update(c,z,count);
    case Method::Symmetric:
        return material_update(c,z,count)&&energy_update(c,z,count)&&material_update(c,z,count);
    }
    return false;
}
// Existence of a hard rod/endpoint-energy intersection. Not a local rank test.
inline bool pendulum_feasible(const Context& c) {
    const auto& m=c.model;
    const V2 center=c.old.x+m.acceleration*(c.h*c.h);
    const double R2=2*c.h*c.h/m.mass*(c.target+m.mass*dot(m.acceleration,c.old.x)
        +.5*m.mass*c.h*c.h*dot(m.acceleration,m.acceleration));
    const double tol=1e-12*(1+m.length+norm(center));
    if (R2 < -tol*tol) return false;
    const double R=std::sqrt(std::max(0.,R2)), d=norm(center);
    return d+tol>=std::abs(m.length-R) && d<=m.length+R+tol;
}
} // namespace energy_ordering
