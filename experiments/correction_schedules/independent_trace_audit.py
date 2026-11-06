import csv,json,math,collections
from pathlib import Path
root=Path('/tmp/xpbd-correction-schedules')
groups=collections.defaultdict(list)
for r in csv.DictReader(open(root/'selected/traces.csv')):
    groups[(r['case'],r['method'])].append(r)
def cross(a,b):return a[0]*b[1]-a[1]*b[0]
def dot(a,b):return a[0]*b[0]+a[1]*b[1]
def add(a,b):return [a[0]+b[0],a[1]+b[1]]
def mul(a,s):return [a[0]*s,a[1]*s]
def phase_path(rows,ref=False):
    out=[];last=None;value=0
    for r in rows:
        now=math.atan2(float(r['ref_y' if ref else 'y']),float(r['ref_x' if ref else 'x']))
        if last is None:value=now
        else:value+=(now-last+math.pi)%(2*math.pi)-math.pi
        out.append(value);last=now
    return out
out={}
for (case,method),rows in groups.items():
    if case not in ['spring_orbit','coupled_motion'] or method not in ['baseline','LE_repair_cap','joint_repair_cap']:continue
    m=float(rows[0]['mass']);h=float(rows[0]['h']);k=1/float(rows[0]['alpha']);L=float(rows[0]['length']);a=[0,-float(rows[0]['gravity'])]
    sum_base=sum_extra=0;closure=0;ratio=0;fracbad=0;capbad=0
    for old,new in zip(rows,rows[1:]):
        x=[float(old['x']),float(old['y'])];v=[float(old['vx']),float(old['vy'])]
        y=add(add(x,mul(v,h)),mul(a,h*h));r=math.hypot(*y);n=mul(y,1/r);C=r-L;u=mul(n,k*C)
        g=add(mul(v,m/h),u);A=1/m+1/(k*h*h);B=dot(n,g)/m;D=dot(g,g)/m
        dl0=-C/A;dx0=mul(n,dl0/m);H0=float(new['H0']);t=float(new['blend'])
        de=-H0/D if D else 0
        if method=='baseline':d=[0,0]
        elif method=='LE_repair_cap':d=mul(g,t*de/m)
        else:
            schur=D-B*B/A;de=-H0/schur if schur else 0
            d=mul(add(g,mul(n,-B/A)),t*de/m)
        xbase=add(y,dx0)
        db=m/h*cross(x,dx0)+.5*m*h*cross(add(x,mul(xbase,-1)),a)
        dd=m/h*cross(x,d)-.5*m*h*cross(d,a)
        measured=float(new['angular_momentum'])-float(old['angular_momentum'])-(float(new['torque_impulse'])-float(old['torque_impulse']))
        closure=max(closure,abs(measured-db-dd));sum_base+=db;sum_extra+=dd
        raw=float(new['candidate_change']);repair=-H0*t*raw
        if method!='baseline':
            bound=float(new['fraction'])*H0*H0
            fracbad+=repair < -1e-20 or repair > bound+1e-16*max(1,bound)
    phases=phase_path(rows);refs=phase_path(rows,True)
    last=rows[-1]
    energies=[]
    for r in rows[-481:]:
        x=[float(r['x']),float(r['y'])];v=[float(r['vx']),float(r['vy'])];rad=math.hypot(*x);vr=dot(x,v)/rad;ell=m*cross(x,v)
        energies.append([.5*m*vr*vr,.5*ell*ell/(m*rad*rad),.5*k*(rad-L)**2])
    out[case+'/'+method]={
       'initial_L':float(rows[0]['angular_momentum']),'final_L':float(last['angular_momentum']),
       'final_torque_balance':float(last['angular_momentum'])-float(rows[0]['angular_momentum'])-float(last['torque_impulse']),
       'counterfactual_material_balance_sum':sum_base,'energy_extra_balance_sum':sum_extra,'max_per_step_decomposition_error':closure,
       'final_angle_difference':phases[-1]-refs[-1],
       'last_second_mean_radial_kinetic':sum(e[0] for e in energies)/len(energies),
       'last_second_mean_tangential_kinetic':sum(e[1] for e in energies)/len(energies),
       'last_second_mean_elastic':sum(e[2] for e in energies)/len(energies),
       'fractional_bound_failures':fracbad,
    }
json.dump(out,open(root/'independent_trace_audit.json','w'),indent=2)
print(json.dumps(out,indent=2))
