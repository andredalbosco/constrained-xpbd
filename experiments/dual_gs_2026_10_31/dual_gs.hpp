#pragma once
#include "predictions.hpp"
#include <chrono>
#include <limits>
#include <numeric>
#include <string>
#include <vector>

namespace dual_gs {
using compose::V; using compose::Mat; using compose::Model; using compose::State;
enum class Kind { Material, Kinetic, Storage, GravityStorage, Length, Tangency };
inline const char* kind_name(Kind kind);
struct Entry { int index; double value; };
using Sparse = std::vector<Entry>;
struct Row { Kind kind; int index; double D=0; std::vector<int> support; std::vector<int> incident_rows; };
struct RowEval { double value=0; Sparse grad; };
struct Counters { long long visits=0,row_evaluations=0,coordinate_writes=0,stationarity_passes=0; };
struct Residuals { double material=0,kinetic=0,storage=0,geometry=0,tangency=0,stationarity=0,energy=0; };
struct Result : Counters, Residuals {
    V z,y; std::string failure; int completed_sweeps=0; double milliseconds=0; bool converged=false;
};
struct Problem {
    Model model; State old;
    double h,E,mref,vs,xs,beta,gamma; int d,n,wcount;
    V aold,what,old_parts; // what is W_hat/E, in incidence order.
    // Work maps contain GLOBAL z-coordinate indices; absent incidences are -1.
    std::vector<int> work_child,work_parent,work_gravity;
    std::vector<Row> rows;
    std::vector<int> material_rows,kinetic_rows,storage_rows,gravity_rows,length_rows,tangency_rows;
    Problem(const Model&,const State&,double h,double E,double beta=1,double gamma=1);
    V initial() const;
    State state(const V& z) const;
    RowEval eval_row(int rowID,const V& z) const;
    const std::vector<int>& support(int rowID) const { return rows.at(rowID).support; }
    double compliance(int rowID) const { return rows.at(rowID).D; }
    std::vector<int> schedule_indices(const std::string&) const;
    V base_grad_local(const V& z,const std::vector<int>& support) const;
    V apply_inv_local(const std::vector<int>& support,const V& values) const;
    V stationarity_local(int rowID,const V& z,const V& y,Counters* counts=nullptr) const;
    double base_objective(const V& z) const;
    double objective(const V& z) const; // Material-eliminated original objective.
    // Whole-system diagnostics/oracles; never used inside a local row visit.
    V base_gradient(const V& z) const;
    V apply_inv(const V& values) const;
    V stationarity(const V& z,const V& y) const;
    Residuals residuals(const V& z,const V& y) const;
    Mat dense_metric() const;
    Mat dense_gradients(const V& z) const;
};
inline std::string visit_refreshed(const Problem& p,int rowID,V& z,V& y,
    const std::string& method,Counters& counts,double omega=1);
inline Result solve(const Problem& p,const std::string& method,const std::string& schedule,
                    int sweeps,double omega=1,int max_visits=-1);
}

namespace dual_gs {
inline const char* kind_name(Kind k) {
    switch(k){case Kind::Material:return "material";case Kind::Kinetic:return "kinetic";
    case Kind::Storage:return "storage";case Kind::GravityStorage:return "gravity_storage";
    case Kind::Length:return "length";case Kind::Tangency:return "tangency";} return "unknown";
}
inline Problem::Problem(const Model& m,const State& s,double dt,double scale,double b,double g)
 :model(m),old(s),h(dt),E(std::max(scale,1e-10)),mref(m.mass.mean()),beta(b),gamma(g) {
    if(!(h>0&&beta>0&&gamma>0&&mref>0)||!std::isfinite(h+E+beta+gamma))
        throw std::invalid_argument("invalid_metric_or_timestep");
    if(model.count<1||model.dim<1||model.mass.size()!=model.count||
       !(model.mass.array()>0).all())throw std::invalid_argument("invalid_model");
    if(model.rigid&&(model.count!=1||model.dim!=2))
        throw std::invalid_argument("rigid_setup_currently_single_planar_rod");
    d=compose::dof(model);vs=std::sqrt(E/mref);xs=h*vs;wcount=0;
    work_child.resize(model.count);work_parent.assign(model.count,-1);work_gravity.assign(model.count,-1);
    for(int j=0;j<model.count;++j){work_child[j]=2*d+wcount++;if(j)work_parent[j]=2*d+wcount++;}
    if(compose::has_gravity(model))for(int i=0;i<model.count;++i)work_gravity[i]=2*d+wcount++;
    n=2*d+wcount;old_parts=compose::potential_parts(model,old.x);
    aold=local_work::components(model,old,false).acceleration;
    what=local_work::work_prediction(model,old,h,"implicit_midpoint_dg")/E;
    if(what.size()!=wcount)throw std::runtime_error("prediction_incidence_count");
    auto vertex=[&](std::vector<int>& a,int i){if(i<0)return;for(int k=0;k<model.dim;++k){a.push_back(i*model.dim+k);a.push_back(d+i*model.dim+k);}};
    auto add=[&](Kind kind,int i,double D,std::vector<int> support){
        std::sort(support.begin(),support.end());support.erase(std::unique(support.begin(),support.end()),support.end());
        int id=rows.size();rows.push_back({kind,i,D,std::move(support),{}});return id;};
    for(int j=0;j<model.count;++j){std::vector<int>a;vertex(a,j);vertex(a,j-1);
        if(model.rigid)length_rows.push_back(add(Kind::Length,j,0,a));
        else {if(!(model.stiffness[j]>0))throw std::invalid_argument("nonpositive_compliant_stiffness");
            material_rows.push_back(add(Kind::Material,j,4*E/(model.stiffness[j]*xs*xs),a));}}
    for(int i=0;i<model.count;++i){std::vector<int>a;vertex(a,i);a.push_back(work_child[i]);
        if(i+1<model.count)a.push_back(work_parent[i+1]);if(work_gravity[i]>=0)a.push_back(work_gravity[i]);
        kinetic_rows.push_back(add(Kind::Kinetic,i,0,a));}
    for(int j=0;j<model.count;++j){std::vector<int>a;if(!model.rigid){vertex(a,j);vertex(a,j-1);}a.push_back(work_child[j]);if(j)a.push_back(work_parent[j]);
        storage_rows.push_back(add(Kind::Storage,j,0,a));}
    if(compose::has_gravity(model))for(int i=0;i<model.count;++i){std::vector<int>a;vertex(a,i);a.push_back(work_gravity[i]);gravity_rows.push_back(add(Kind::GravityStorage,i,0,a));}
    if(model.rigid)for(int j=0;j<model.count;++j){std::vector<int>a;vertex(a,j);vertex(a,j-1);tangency_rows.push_back(add(Kind::Tangency,j,0,a));}
    // Structural incidence is built once. Per-visit gathers never scan all rows.
    std::vector<std::vector<int>> coordinate_rows(n);
    for(int r=0;r<int(rows.size());++r)for(int k:rows[r].support)coordinate_rows[k].push_back(r);
    for(auto& row:rows){for(int k:row.support)row.incident_rows.insert(row.incident_rows.end(),coordinate_rows[k].begin(),coordinate_rows[k].end());
        std::sort(row.incident_rows.begin(),row.incident_rows.end());row.incident_rows.erase(std::unique(row.incident_rows.begin(),row.incident_rows.end()),row.incident_rows.end());}
}
inline V Problem::initial() const {
    V z=V::Zero(n);
    for(int i=0;i<d;++i){double s=old.v[i]+.25*h*(aold[i]+model.gravity[i%model.dim]);z[i]=s/vs;z[d+i]=(2*s-old.v[i])/vs;}
    z.tail(wcount)=what;return z;
}
inline State Problem::state(const V& z) const { return {old.x+xs*z.head(d),vs*z.segment(d,d)}; }
inline RowEval Problem::eval_row(int id,const V& z) const {
    const auto& row=rows.at(id);const int j=row.index,dim=model.dim;RowEval e;
    auto add=[&](int index,double v){e.grad.push_back({index,v});};
    auto position=[&](int i)->V{return old.x.segment(i*dim,dim)+xs*z.segment(i*dim,dim);};
    if(row.kind==Kind::Storage&&model.rigid){
        // The zero rod potential is symbolically independent of positions.
        e.value=z[work_child[j]];add(work_child[j],1);
        if(j){e.value+=z[work_parent[j]];add(work_parent[j],1);}return e;
    }
    if(row.kind==Kind::Kinetic){e.value=.5*model.mass[j]*(vs*vs*z.segment(d+j*dim,dim).squaredNorm()-old.v.segment(j*dim,dim).squaredNorm())/E;
        for(int k=0;k<dim;++k)add(d+j*dim+k,model.mass[j]/mref*z[d+j*dim+k]);
        for(int w:{work_child[j],j+1<model.count?work_parent[j+1]:-1,work_gravity[j]})if(w>=0){e.value-=z[w];add(w,-1);}
    }else if(row.kind==Kind::GravityStorage){V gradient=-model.mass[j]*model.gravity;
        // Exact linear increment, avoiding subtraction of gravitational offsets.
        e.value=xs/E*gradient.dot(z.segment(j*dim,dim))+z[work_gravity[j]];
        for(int k=0;k<dim;++k)add(j*dim+k,xs/E*gradient[k]);add(work_gravity[j],1);
    }else{
        V relative=position(j);if(j)relative-=position(j-1);
        double length=compose::edge_length(model,relative);V normal=compose::edge_direction(model,relative);
        double C=length-model.length[j];
        if(row.kind==Kind::Material||row.kind==Kind::Length){e.value=C/xs;
            for(int k=0;k<dim;++k){add(j*dim+k,normal[k]);if(j)add((j-1)*dim+k,-normal[k]);}
        }else if(row.kind==Kind::Storage){
            double U=model.rigid?0:.5*model.stiffness[j]*C*C;e.value=(U-old_parts[j])/E+z[work_child[j]];
            if(!model.rigid)for(int k=0;k<dim;++k){double q=xs/E*model.stiffness[j]*C*normal[k];add(j*dim+k,q);if(j)add((j-1)*dim+k,-q);}
            add(work_child[j],1);if(j){e.value+=z[work_parent[j]];add(work_parent[j],1);}
        }else if(row.kind==Kind::Tangency){V relative_u=z.segment(d+j*dim,dim);if(j)relative_u-=z.segment(d+(j-1)*dim,dim);
            e.value=normal.dot(relative_u);V q=dim==1?V::Zero(1):V(xs/length*(relative_u-normal*e.value));
            for(int k=0;k<dim;++k){add(j*dim+k,q[k]);add(d+j*dim+k,normal[k]);if(j){add((j-1)*dim+k,-q[k]);add(d+(j-1)*dim+k,-normal[k]);}}
        }else throw std::runtime_error("unknown_row_kind");
    }
    std::sort(e.grad.begin(),e.grad.end(),[](const Entry&a,const Entry&b){return a.index<b.index;});return e;
}
inline V Problem::base_grad_local(const V& z,const std::vector<int>& a) const {
    V r(a.size());for(int t=0;t<int(a.size());++t){int k=a[t];
        if(k>=2*d){r[t]=gamma*(z[k]-what[k-2*d]);continue;}
        int i=k%d;double ratio=model.mass[i/model.dim]/mref;
        double uold=old.v[i]/vs,position=z[i]-uold-.25*h*aold[i]/vs,velocity=z[d+i]+uold-2*z[i];
        r[t]=k<d?ratio*(position-2*velocity/beta)-.25*model.mass[i/model.dim]*model.gravity[i%model.dim]*xs/E:ratio*velocity/beta;
    }return r;
}
inline V Problem::apply_inv_local(const std::vector<int>& a,const V& r) const {
    if(int(a.size())!=r.size())throw std::invalid_argument("local_metric_size");
    V out=V::Zero(a.size());
    for(int t=0;t<int(a.size());++t){int k=a[t];if(k>=2*d){out[t]=r[t]/gamma;continue;}if(k>=d){
        if(!std::binary_search(a.begin(),a.end(),k-d))throw std::runtime_error("support_not_closed_under_metric");continue;}
        auto at=std::lower_bound(a.begin(),a.end(),d+k);if(at==a.end()||*at!=d+k)throw std::runtime_error("support_not_closed_under_metric");
        int partner=at-a.begin();double factor=mref/model.mass[k/model.dim];out[t]=factor*(r[t]+2*r[partner]);out[partner]=factor*(2*r[t]+(beta+4)*r[partner]);}
    return out;
}
inline V Problem::stationarity_local(int id,const V& z,const V& y,Counters* counts) const {
    const auto& row=rows.at(id);const auto& a=row.support;V residual=base_grad_local(z,a);
    for(int r:row.incident_rows){if(y[r]==0)continue;auto e=eval_row(r,z);if(counts)++counts->row_evaluations;
        for(const auto&g:e.grad){auto at=std::lower_bound(a.begin(),a.end(),g.index);if(at!=a.end()&&*at==g.index)residual[at-a.begin()]-=y[r]*g.value;}}
    return residual;
}
inline std::vector<int> Problem::schedule_indices(const std::string& name) const {
    const auto& material=model.rigid?length_rows:material_rows;std::vector<int> result;
    auto append=[&](const std::vector<int>& a){result.insert(result.end(),a.begin(),a.end());};
    auto M=[&](){append(material);};auto S=[&](){append(storage_rows);append(gravity_rows);};auto K=[&](){append(kinetic_rows);};
    if(name=="material_storage_kinetic"){M();S();K();}
    else if(name=="kinetic_storage_material"){K();S();M();}
    else if(name=="storage_kinetic_material"){S();K();M();}
    else if(name=="edge"||name=="edge_storage_first"||name=="edge_reverse"){
        for(int t=0;t<model.count;++t){int j=name=="edge_reverse"?model.count-1-t:t;if(name=="edge_storage_first"){result.push_back(storage_rows[j]);result.push_back(material[j]);}
            else {result.push_back(material[j]);result.push_back(storage_rows[j]);}
            if(!gravity_rows.empty())result.push_back(gravity_rows[j]);result.push_back(kinetic_rows[j]);if(j)result.push_back(kinetic_rows[j-1]);
            if(!tangency_rows.empty())result.push_back(tangency_rows[j]);}
    }else throw std::invalid_argument("unknown_schedule_"+name);
    if(name=="material_storage_kinetic"||name=="kinetic_storage_material"||name=="storage_kinetic_material")append(tangency_rows);
    return result;
}
inline double Problem::base_objective(const V& z) const {
    double value=0;for(int i=0;i<d;++i){double ratio=model.mass[i/model.dim]/mref,uold=old.v[i]/vs;
        double position=z[i]-uold-.25*h*aold[i]/vs,velocity=z[d+i]+uold-2*z[i];
        value+=.5*ratio*(position*position+velocity*velocity/beta)-.25/E*model.mass[i/model.dim]*model.gravity[i%model.dim]*(old.x[i]+xs*z[i]);}
    value+=.5*gamma*(z.tail(wcount)-what).squaredNorm();return value;
}
inline double Problem::objective(const V& z) const {
    double value=base_objective(z);if(!model.rigid)for(int r:material_rows){double c=eval_row(r,z).value;value+=.5*c*c/rows[r].D;}return value;
}
inline V Problem::base_gradient(const V& z) const {std::vector<int>a(n);std::iota(a.begin(),a.end(),0);return base_grad_local(z,a);}
inline V Problem::apply_inv(const V& v) const {std::vector<int>a(n);std::iota(a.begin(),a.end(),0);return apply_inv_local(a,v);}
inline V Problem::stationarity(const V& z,const V& y) const {
    V r=base_gradient(z);for(int j=0;j<int(rows.size());++j){auto e=eval_row(j,z);for(auto g:e.grad)r[g.index]-=y[j]*g.value;}return r;
}
inline Residuals Problem::residuals(const V& z,const V& y) const {
    Residuals out;for(int j=0;j<int(rows.size());++j){double r=std::abs(eval_row(j,z).value+rows[j].D*y[j]);
        switch(rows[j].kind){case Kind::Material:out.material=std::max(out.material,r);break;case Kind::Kinetic:out.kinetic=std::max(out.kinetic,r);break;
        case Kind::Storage:case Kind::GravityStorage:out.storage=std::max(out.storage,r);break;case Kind::Length:out.geometry=std::max(out.geometry,r);break;case Kind::Tangency:out.tangency=std::max(out.tangency,r);break;}}
    out.stationarity=stationarity(z,y).norm();out.energy=std::abs(compose::energy(model,state(z))-compose::energy(model,old))/E;return out;
}
inline Mat Problem::dense_metric() const {
    Mat A=Mat::Zero(n,n);for(int i=0;i<d;++i){double ratio=model.mass[i/model.dim]/mref;
        A(i,i)=ratio*(1+4/beta);A(i,d+i)=A(d+i,i)=-2*ratio/beta;A(d+i,d+i)=ratio/beta;}
    for(int i=2*d;i<n;++i)A(i,i)=gamma;return A;
}
inline Mat Problem::dense_gradients(const V& z) const {Mat G=Mat::Zero(n,rows.size());for(int r=0;r<int(rows.size());++r)for(auto g:eval_row(r,z).grad)G(g.index,r)+=g.value;return G;}
inline V restricted_gradient(const RowEval& e,const std::vector<int>& a) {
    V g=V::Zero(a.size());for(auto entry:e.grad){auto at=std::lower_bound(a.begin(),a.end(),entry.index);
        if(at==a.end()||*at!=entry.index)throw std::runtime_error("gradient_outside_structural_support");g[at-a.begin()]+=entry.value;}return g;
}
inline std::string visit_refreshed(const Problem& p,int id,V& z,V& y,const std::string& method,Counters& counts,double omega) {
    if(method!="refreshed_drop"&&method!="refreshed_local_rx")return "unknown_refreshed_method";
    auto row=p.eval_row(id,z);++counts.row_evaluations;const auto& a=p.rows[id].support;
    V gradient=restricted_gradient(row,a),direction=p.apply_inv_local(a,gradient),restoration=V::Zero(a.size());
    if(method=="refreshed_local_rx")restoration=p.apply_inv_local(a,p.stationarity_local(id,z,y,&counts));
    double denominator=gradient.dot(direction)+p.rows[id].D;
    double rhs=-(row.value+p.rows[id].D*y[id])+gradient.dot(restoration),delta=0;
    if(!std::isfinite(denominator+rhs))return "nonfinite_row_model";
    if(denominator>0)delta=rhs/denominator;
    else if(denominator<0||std::abs(rhs)>1e-12)return "zero_gradient_incompatible";
    y[id]+=omega*delta;for(int k=0;k<int(a.size());++k)z[a[k]]+=omega*(direction[k]*delta-restoration[k]);
    ++counts.visits;counts.coordinate_writes+=a.size()+1;return {};
}
inline Result solve(const Problem& p,const std::string& method,const std::string& schedule,int sweeps,double omega,int max_visits) {
    auto start=std::chrono::steady_clock::now();Result out;out.z=p.initial();out.y=V::Zero(p.rows.size());
    auto order=p.schedule_indices(schedule);
    if(sweeps<0||!(omega>0)||!std::isfinite(omega))throw std::invalid_argument("invalid_budget_or_relaxation");
    if(method!="refreshed_drop"&&method!="refreshed_local_rx"&&method!="frozen_full_rx")throw std::invalid_argument("unknown_method_"+method);
    auto safety=[&](){if(!out.z.allFinite()||!out.y.allFinite())return std::string("nonfinite_iterate");if(out.z.lpNorm<Eigen::Infinity>()>1e8)return std::string("coordinate_limit");return std::string();};
    auto local_safety=[&](int row){if(!std::isfinite(out.y[row]))return std::string("nonfinite_iterate");
        for(int k:p.rows[row].support){if(!std::isfinite(out.z[k]))return std::string("nonfinite_iterate");if(std::abs(out.z[k])>1e8)return std::string("coordinate_limit");}return std::string();};
    out.failure=safety();
    try{
        for(int sweep=0;sweep<sweeps&&out.failure.empty();++sweep){if(max_visits>=0&&out.visits>=max_visits)break;
            V origin;std::vector<RowEval> frozen;
            if(method=="frozen_full_rx"){
                origin=out.z;V rz=p.base_gradient(origin);frozen.reserve(p.rows.size());
                for(int r=0;r<int(p.rows.size());++r){frozen.push_back(p.eval_row(r,origin));++out.row_evaluations;for(auto g:frozen.back().grad)rz[g.index]-=out.y[r]*g.value;}
                // Exactly once per frozen outer sweep; omega damps row pushes only.
                out.z-=p.apply_inv(rz);++out.stationarity_passes;out.coordinate_writes+=p.n;
                out.failure=safety();if(!out.failure.empty())break;
            }
            bool complete=true;
            for(int r:order){if(max_visits>=0&&out.visits>=max_visits){complete=false;break;}
                if(method=="frozen_full_rx"){
                    const auto& row=frozen[r];const auto& a=p.rows[r].support;V gradient=restricted_gradient(row,a),direction=p.apply_inv_local(a,gradient);
                    double value=row.value+p.rows[r].D*out.y[r];for(auto g:row.grad)value+=g.value*(out.z[g.index]-origin[g.index]);
                    double denominator=gradient.dot(direction)+p.rows[r].D,delta=0;
                    if(!std::isfinite(denominator+value)){out.failure="nonfinite_row_model";complete=false;break;}
                    if(denominator>0)delta=-value/denominator;
                    else if(denominator<0||std::abs(value)>1e-12){out.failure="zero_gradient_incompatible";complete=false;break;}
                    out.y[r]+=omega*delta;for(int k=0;k<int(a.size());++k)out.z[a[k]]+=omega*direction[k]*delta;
                    ++out.visits;out.coordinate_writes+=a.size()+1;
                }else out.failure=visit_refreshed(p,r,out.z,out.y,method,out,omega);
                if(out.failure.empty())out.failure=local_safety(r);if(!out.failure.empty()){complete=false;break;}
            }
            if(complete)++out.completed_sweeps;if(!out.failure.empty()||!complete)break;
        }
    }catch(const std::exception& e){out.failure=std::string("evaluation_failure:")+e.what();}
    out.milliseconds=std::chrono::duration<double,std::milli>(std::chrono::steady_clock::now()-start).count();
    // Final diagnostics are outside the numerical work counters and timing.
    auto invalidate_residuals=[&](){double inf=std::numeric_limits<double>::infinity();out.material=out.kinetic=out.storage=out.geometry=out.tangency=out.stationarity=out.energy=inf;};
    if(out.z.allFinite()&&out.y.allFinite())try{static_cast<Residuals&>(out)=p.residuals(out.z,out.y);
        out.converged=std::max({out.material,out.kinetic,out.storage,out.geometry,out.tangency})<1e-8&&out.stationarity<1e-7;
    }catch(const std::exception&e){invalidate_residuals();if(out.failure.empty())out.failure=std::string("diagnostic_failure:")+e.what();}
    else invalidate_residuals();
    return out;
}
}
