"""Single-particle controls: exact nonlinear energy visits and constrained PD.

Every endpoint method reconstructs v=(x-x_old)/h. PD local/global uses exact
physical energy unless the explicitly named surrogate variant is selected.
No final energy repair, velocity rescale, or failure fallback is hidden here.
"""
from dataclasses import dataclass
import math
import numpy as np
from scipy.optimize import minimize_scalar

@dataclass
class Result:
    x: np.ndarray
    v: np.ndarray
    converged: bool = True
    reason: str = ''
    post_E: float = 0.
    residual: float = 0.
    nfev: int = 0
    ties: int = 0
    stationarity: float = math.nan


def potential(x, m, k, l, g):
    return .5*k*(np.linalg.norm(x)-l)**2-m*np.dot(g,x)


def energy(x,v,m,k,l,g):
    return .5*m*np.dot(v,v)+potential(x,m,k,l,g)


def endpoint_energy(x, old, h, m,k,l,g):
    return energy(x,(x-old)/h,m,k,l,g)


def gradient(x,old,h,m,k,l,g):
    r=np.linalg.norm(x)
    if r<1e-14: raise ValueError('collapsed_spring')
    return m*(x-(old+h*h*g))/h**2+k*(r-l)*x/r


class EnergyCurve:
    """Exact physical endpoint-energy contour for one radial spring in 2D.

    Analytic radial roots, smoothly joined using theta=theta_max sin(phi).
    Supports the component not enclosing the origin (all declared test cases).
    Degenerate/origin-enclosing cases return explicit failures.
    """
    def __init__(self,old,h,m,k,l,g,target):
        self.old,self.h,self.m,self.k,self.l,self.g,self.target=old,h,m,k,l,g,target
        self.A=m/h**2+k
        b=m*old/h**2+m*g
        self.B=float(np.linalg.norm(b))
        if self.B<1e-14: raise ValueError('zero_energy_axis')
        self.axis=b/self.B; self.tan=np.array([-self.axis[1],self.axis[0]])
        self.rmin=(self.B+k*l)/self.A
        center=self.rmin*self.axis
        self.gap=target-endpoint_energy(center,old,h,m,k,l,g)
        tol=2e-11*max(1.,abs(target),abs(potential(old,m,k,l,g)))
        if self.gap < -tol: raise ValueError('no_energy_surface')
        self.gap=max(0.,self.gap)
        # Difference of squared coefficients, with cancellation-free angular width.
        b0=self.B+k*l
        discr=b0*b0-2*self.A*self.gap
        if discr<=0: raise ValueError('energy_contour_encloses_origin')
        diff=2*self.A*self.gap/(b0+math.sqrt(discr))
        self.cosmax=1-diff/self.B
        if not (-1 < self.cosmax <= 1): raise ValueError('unsupported_energy_contour')
        self.tmax=2*math.asin(math.sqrt(max(0.,diff/(2*self.B))))
        self.nfev=0

    def point(self,phi):
        self.nfev+=1
        theta=self.tmax*math.sin(phi)
        ct=math.cos(theta)
        cdiff=2*math.sin((self.tmax+theta)/2)*math.sin((self.tmax-theta)/2)
        gaptheta=self.B*cdiff*(self.B*(ct+self.cosmax)+2*self.k*self.l)/(2*self.A)
        rad=(self.B*ct+self.k*self.l)/self.A
        rad+=math.copysign(math.sqrt(max(0.,2*gaptheta/self.A)),math.cos(phi))
        return rad*(ct*self.axis+math.sin(theta)*self.tan)

    def closest(self,z,tie,grid=32):
        if self.gap==0: return self.point(0),1
        step=2*math.pi/grid
        phis=np.arange(grid)*step
        pts=np.array([self.point(float(t)) for t in phis])
        values=np.sum((pts-z)**2,axis=1)
        candidates=[]
        for j in range(grid):
            if values[j]<=values[j-1] and values[j]<=values[(j+1)%grid]:
                f=lambda t: float(np.sum((self.point(t)-z)**2))
                opt=minimize_scalar(f,bounds=(phis[j]-step,phis[j]+step),method='bounded',options={'xatol':2e-12})
                candidates.append(self.point(float(opt.x)))
                candidates.append(pts[j])
        # Feasible tie seed is essential at rest, where the objective is flat.
        if abs(endpoint_energy(tie,self.old,self.h,self.m,self.k,self.l,self.g)-self.target)<2e-10*max(1,abs(self.target)):
            candidates.append(tie.copy())
        if not candidates: raise ValueError('contour_search_failed')
        vals=np.array([np.sum((p-z)**2) for p in candidates])
        tol=1e-12*max(1e-8,float(np.max(values)),self.h*self.h)
        ids=np.flatnonzero(vals <= vals.min()+tol)
        chosen=min(ids,key=lambda j: np.sum((candidates[j]-tie)**2))
        return candidates[chosen],max(0,len(ids)-1)


def circle_project(z,c,r,tie):
    d=z-c; length=np.linalg.norm(d)
    degenerate=length<2e-12*max(1.,np.linalg.norm(c))
    if degenerate:
        d=tie-c; length=np.linalg.norm(d)
        if length<1e-14: d=np.array([1.,0.]); length=1.
    return c+r*d/length,int(degenerate)


def rod_energy_circle(old,h,m,g,target):
    c=old+h*h*g
    rsq=2*h*h/m*(target+m*np.dot(g,old)+.5*m*h*h*np.dot(g,g))
    if rsq < -1e-12*h*h*max(1,abs(target)): raise ValueError('no_energy_surface')
    return c,math.sqrt(max(0,rsq))


def rod_joint(old,v,h,m,l,g,target):
    c,R=rod_energy_circle(old,h,m,g,target)
    d=np.linalg.norm(c)
    if d<1e-14: raise ValueError('concentric_rod_energy_circles')
    axis=c/d; perp=np.array([-axis[1],axis[0]])
    # a = (l²+d²-R²)/(2d); evaluate short transverse height via Heron.
    a=(l*l+d*d-R*R)/(2*d)
    terms=np.array([l+R+d,-l+R+d,l-R+d,l+R-d])
    tol=2e-12*max(l,R,d)
    if np.any(terms < -tol): raise ValueError('no_rod_energy_intersection')
    b=math.sqrt(max(0,float(np.prod(np.maximum(terms,0)))))/(2*d)
    candidates=[a*axis+b*perp,a*axis-b*perp]
    y=old+h*v+h*h*g
    vals=[np.sum((p-y)**2) for p in candidates]
    tolerance=1e-12*max(h*h,1e-8,max(vals))
    ids=[j for j in range(2) if vals[j]<=min(vals)+tolerance]
    j=min(ids,key=lambda j: np.sum((candidates[j]-old)**2))
    x=candidates[j]
    # Return exact old point on a rest tie to prevent roundoff choosing motion.
    if np.linalg.norm(v)<1e-13 and abs(energy(old,v,m,0,l,g)-target)<1e-12:
        x=old.copy()
    return x,len(ids)-1


def step(x,v,h,m,k,l,g,target,method,sweeps=1,rigid=False,grid=32):
    old=np.array(x,dtype=float); v=np.array(v,dtype=float); g=np.array(g,dtype=float)
    z=old+h*v+h*h*g
    curve=None; ties=0; calls=0; post_E=math.nan
    try:
        if method=='pd_joint' and rigid:
            z,ties=rod_joint(old,v,h,m,l,g,target)
        elif method in ('pd_true','pd_inertia'):
            if rigid: raise ValueError('use_pd_joint_for_rigid')
            curve=EnergyCurve(old,h,m,k,l,g,target)
            for _ in range(sweeps):
                p=l*z/np.linalg.norm(z)
                weight=k if method=='pd_true' else 0.
                center=(m/h**2*(old+h*v+h*h*g)+weight*p)/(m/h**2+weight)
                z,nt=curve.closest(center,old,grid); ties+=nt
            post_E=endpoint_energy(z,old,h,m,k,l,g)-target
        elif method=='pd_surrogate':
            if rigid: raise ValueError('surrogate_not_defined_for_rigid')
            A=m/h**2+k
            for _ in range(sweeps):
                p=l*z/np.linalg.norm(z)
                c=(m/h**2*old+m*g+k*p)/A
                zmin=(m/h**2*(old+h*v+h*h*g)+k*p)/A
                emin=.5*m*np.sum(((c-old)/h)**2)+.5*k*np.sum((c-p)**2)-m*np.dot(g,c)
                radius2=2*(target-emin)/A
                if radius2 < -1e-12: raise ValueError('no_surrogate_energy_surface')
                z,nt=circle_project(zmin,c,math.sqrt(max(0,radius2)),old); ties+=nt
            post_E=endpoint_energy(z,old,h,m,k,l,g)-target
        elif method in ('raw_EL','nonlinear_EL','xpbd','pd'):
            lam=0.; eta=0.
            if method=='nonlinear_EL' and not rigid: curve=EnergyCurve(old,h,m,k,l,g,target)
            for _ in range(sweeps):
                if method=='pd':
                    if rigid: z=l*z/np.linalg.norm(z)
                    else:
                        p=l*z/np.linalg.norm(z)
                        z=(m/h**2*(old+h*v+h*h*g)+k*p)/(m/h**2+k)
                    continue
                if method=='raw_EL':
                    grad=gradient(z,old,h,m,k,l,g)
                    e=endpoint_energy(z,old,h,m,k,l,g)-target
                    gg=np.dot(grad,grad)/m
                    grad_scale=1+np.linalg.norm(m*(z-old)/h**2)+abs(k*(np.linalg.norm(z)-l))+m*np.linalg.norm(g)
                    if np.linalg.norm(grad)>64*np.finfo(float).eps*grad_scale:
                        de=-e/gg; z=z+grad/m*de; eta+=de
                    calls+=1
                elif method=='nonlinear_EL':
                    if rigid:
                        c,R=rod_energy_circle(old,h,m,g,target)
                        z,nt=circle_project(z,c,R,old); ties+=nt
                    else:
                        z,nt=curve.closest(z,old,grid); ties+=nt
                    calls+=1
                post_E=endpoint_energy(z,old,h,m,k,l,g)-target
                r=np.linalg.norm(z)
                if r<1e-13: raise ValueError('collapsed_material_row')
                n=z/r; at=0 if rigid else 1/(k*h*h)
                dl=(-(r-l)-at*lam)/(1/m+at)
                z=z+n*dl/m;lam+=dl
        else: raise ValueError('unknown_method')
        vv=(z-old)/h
        res=energy(z,vv,m,k,l,g)-target
        if not np.all(np.isfinite(z)) or not np.all(np.isfinite(vv)): raise ValueError('nonfinite')
        # Numerical feasibility for exact final constraints is verified separately.
        return Result(z,vv,True,'',post_E,res,calls+(curve.nfev if curve else 0),ties)
    except (ValueError,FloatingPointError) as exc:
        return Result(old,v,False,str(exc),post_E,math.inf,calls,ties)
