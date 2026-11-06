#!/usr/bin/env python3
"""Plot gradient timing and rejected energy updates in the endpoint-energy solver."""
import argparse
import collections
import csv
import json
import math
import os
from pathlib import Path
os.environ.setdefault('MPLCONFIGDIR','/tmp/xpbd-policy-matplotlib')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np

CASES=['spring_rest','spring_radial','spring_tangent','pendulum_rest','pendulum_tangent']
CASE_NAMES=['Spring: released from rest','Spring: initial radial velocity','Spring: initial tangential velocity',
            'Pendulum: released from rest','Pendulum: initial tangential velocity']
POLICIES=['fresh','lag_sweep','frozen_step','previous_step']
NAMES={'fresh':'Current gradient','lag_sweep':'Previous sweep gradient','frozen_step':'Predictor gradient',
       'previous_step':'Previous timestep gradient','none':'Ordinary XPBD'}
COLORS={'fresh':'#1769aa','lag_sweep':'#d47715','frozen_step':'#975399','previous_step':'#138765','none':'#8c9298'}
KEY=['case','order','gradient','gamma','retain_rx','h','sweeps']
TEXT={'case','order','gradient','status'}


def read(path):
    with path.open() as stream:rows=list(csv.DictReader(stream))
    for row in rows:
        for field in row:
            if field not in TEXT:
                row[field]=float(row[field])
                if not math.isfinite(row[field]):raise ValueError(f'Nonfinite {field} in {path}')
    return rows


def key(row):return tuple(row[k] for k in KEY)
def close(a,b):return math.isclose(a,b,rel_tol=1e-12,abs_tol=1e-13)
def gamma_name(gamma):return 'no rejection' if gamma<0 else f'γ={gamma:g}'
def variant_name(policy,gamma):return NAMES[policy]+' · '+gamma_name(gamma)
def condition(order,retain):return ('Energy → material (EL)' if order=='EL' else 'Material → energy (LE)')+f" · {'retain' if retain else 'drop'} rₓ"
def slug_gamma(gamma):return 'off' if gamma<0 else str(gamma).replace('.','p')


class Data:
    def __init__(self,root):
        self.rows=read(root/'summary.csv')
        self.summary={key(r):r for r in self.rows}
        if len(self.summary)!=len(self.rows):raise ValueError('Duplicate summary keys')
        self.traces=collections.defaultdict(list)
        self.energy_check=0.
        for r in read(root/'traces.csv'):
            potential=9.81*r['y'] if r['case'].startswith('pendulum') else 50*(math.hypot(r['x'],r['y'])-1)**2
            E=.5*(r['vx']**2+r['vy']**2)+potential
            discrepancy=abs(E-r['energy'])/max(1,abs(E))
            self.energy_check=max(self.energy_check,discrepancy)
            if discrepancy>2e-12:raise ValueError('Endpoint energy reconstruction failed')
            self.traces[key(r)].append(r)
        for settings,trace in self.traces.items():
            s=self.summary[settings]
            if len(trace)!=s['completed_steps']+1:raise ValueError('Trace length mismatch')
            if any(r['step']!=i or not close(r['time'],i*r['h']) for i,r in enumerate(trace)):
                raise ValueError('Trace step/time mismatch')
            error=max(abs(r['energy']-r['target_energy']) for r in trace)
            if not math.isclose(error,s['max_abs_energy_error'],rel_tol=1e-10,abs_tol=1e-10):
                raise ValueError('Trace and summary energy error mismatch')
            if trace[-1]['oversized_skips']>s['oversized_skips']:raise ValueError('Rejected update count decreased')

    def get(self,case,order,policy,gamma,retain,h=1/60,sweeps=5):
        if order=='baseline':policy,gamma,retain='none',-1,0
        rows=[r for r in self.rows if r['case']==case and r['order']==order and r['gradient']==policy
              and r['gamma']==gamma and r['retain_rx']==retain and close(r['h'],h) and r['sweeps']==sweeps]
        if len(rows)!=1:raise ValueError(f'Expected one row: {case,order,policy,gamma,retain,h,sweeps}, got {len(rows)}')
        s=rows[0]
        return s,self.traces.get(key(s),[])

    def base(self,case):return self.get(case,'baseline','none',-1,0)


def values(trace,reference=False):
    fields=['ref_energy','ref_x','ref_y'] if reference else ['energy','x','y']
    return [[r[f] for r in trace] for f in fields]


def limits(data,case,order,variants,retain):
    ref=data.base(case)[1]
    rv=values(ref,True)
    result=[]
    for col in range(3):
        extent=max(1,max(abs(x) for x in rv[col]))
        normal=list(rv[col])+values(ref)[col]
        for policy,gamma in variants:
            tr=data.get(case,order,policy,gamma,retain)[1]
            v=values(tr)[col]
            extreme=(max(abs(x-tr[0]['target_energy']) for x in v)>10*tr[0]['Escale']) if col==0 else max(abs(x) for x in v)>3*extent
            if not extreme:normal+=v
        lo,hi=min(normal),max(normal)
        pad=max(.07*(hi-lo),.05*extent)
        result.append({'scale':extent,'range':(lo-pad,hi+pad)})
    return result


def setup(ax,base_limit,trace_groups,col):
    vals=[v for trace in trace_groups for v in values(trace)[col]]
    first=trace_groups[0][0];extent=base_limit['scale']
    lo,hi=min(vals),max(vals)
    extreme=(max(abs(v-first['target_energy']) for v in vals)>10*first['Escale']) if col==0 else max(abs(lo),abs(hi))>3*extent
    if extreme:
        lo=lo*1.15 if lo < -extent else lo-.05*extent
        hi=hi*1.15 if hi > extent else hi+.05*extent
        ax.set_yscale('symlog',linthresh=extent)
        locator=ticker.SymmetricalLogLocator(base=10,linthresh=extent);locator.set_params(numticks=5)
        ticks=[v for v in locator.tick_values(lo,hi) if v==0 or abs(v)>=extent]
        visible=[v for v in ticks if lo<=v<=hi]
        if len(visible)<2:
            ticks=sorted(set(visible+[extent*v for v in [-4,-2,-1,0,1,2,4] if lo<=extent*v<=hi]))
            if len(ticks)<2:ticks=list(np.linspace(lo,hi,3))
            ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda v,_:f'{v:g}'))
        ax.yaxis.set_major_locator(ticker.FixedLocator(ticks));ax.yaxis.set_minor_locator(ticker.NullLocator())
        ax.set_ylim(lo,hi)
        ax.text(.98,.96,'symlog',transform=ax.transAxes,ha='right',va='top',fontsize=7,color='#a54032')
    else:ax.set_ylim(*base_limit['range'])
    ax.set_xlim(0,10);ax.grid(True,alpha=.18,linewidth=.7);ax.tick_params(labelsize=7)
    ax.set_ylabel(['Total mechanical E [J]','Position x [m]','Position y [m]'][col],fontsize=8)


def reference(axes,base):
    times=[r['time'] for r in base]
    for ax,v in zip(axes,values(base,True)):ax.plot(times,v,color='#333a43',ls='--',lw=1.0,alpha=.8,label='RK4 reference')
    for ax,v in zip(axes,values(base)):ax.plot(times,v,color='#b7bcc3',lw=1.0,label='Ordinary XPBD')
    axes[0].axhline(base[0]['target_energy'],color='#1d232a',ls=':',lw=.8,label='Initial energy')


def draw(axes,summary,trace,policy,gamma,color=None,label=None):
    color=color or COLORS[policy]
    times=[r['time'] for r in trace]
    for ax,v in zip(axes,values(trace)):
        ax.plot(times,v,color=color,lw=1.35,label=label or variant_name(policy,gamma))
        if summary['completed_steps']<summary['requested_steps']:
            ax.plot(times[-1],v[-1],'x',color=color,ms=7,mew=1.4)
            ax.axvline(times[-1],color=color,lw=.6,ls=':',alpha=.5)


def finish(fig,out,name,title,handles,labels):
    fig.suptitle(title,fontsize=14,y=.995)
    fig.legend(handles,labels,loc='upper center',bbox_to_anchor=(.5,.958),ncol=min(len(labels),4),frameon=False,fontsize=8)
    fig.text(.5,.014,'Original endpoint-energy formulation: v = (x − xₙ)/h. h = 1/60 s; five sweeps. All three residuals use current analytical values.\n'
             'Oversized energy corrections are rejected, never clipped. Bounded panels share ranges; extreme panels use marked symlog ranges. Crosses mark guarded stops.\n'
             'Predictor gradient is frozen for a timestep; previous-sweep and previous-timestep gradients are historical snapshots. Position curves expose freezing and phase drift.',
             ha='center',fontsize=7.3)
    fig.tight_layout(rect=(.005,.067,.995,.925),h_pad=1.65,w_pad=1.5)
    fig.savefig(out/(name+'.png'),dpi=165,facecolor='white');plt.close(fig)
    return name+'.png'


def policy_plot(data,out,order,policy,gamma,retain):
    fig,axes=plt.subplots(5,3,figsize=(13.5,13),squeeze=False)
    variants=[(p,1) for p in POLICIES]+[('fresh',4)]
    for row,case in enumerate(CASES):
        summary,tr=data.get(case,order,policy,gamma,retain);base=data.base(case)[1]
        reference(axes[row],base);draw(axes[row],summary,tr,policy,gamma)
        lim=limits(data,case,order,variants,retain)
        for col,ax in enumerate(axes[row]):
            setup(ax,lim[col],[base,tr],col)
            if row==4:ax.set_xlabel('Time [s]')
        axes[row,0].set_title(CASE_NAMES[row],loc='left',weight='bold',fontsize=10,pad=8)
        axes[row,1].set_title(f'Rejected {int(summary["oversized_skips"])}/{int(summary["energy_attempts"])} energy updates · immobile {int(summary["frozen_steps"])}/{int(summary["completed_steps"])} steps',loc='left',fontsize=7.5,pad=8)
        status='complete' if summary['completed_steps']==summary['requested_steps'] else f'stopped {summary["completed_steps"]*summary["h"]:.3g}s'
        axes[row,2].set_title(f'Path {summary["path_ratio"]:.3g}× RK4 · reversals {int(summary["velocity_reversals"])}/{int(summary["reference_reversals"])} · {status}',loc='left',fontsize=7.5,pad=8)
    handles,labels=axes[0,0].get_legend_handles_labels()
    name=f'{policy}_{order}_gamma{slug_gamma(gamma)}_{"retain" if retain else "drop"}_rx'
    return finish(fig,out,name,NAMES[policy]+' · '+gamma_name(gamma)+'\n'+condition(order,retain),handles,labels)


def comparison(data,out,order,variants,retain,name,title):
    fig,axes=plt.subplots(5,3,figsize=(13.5,13),squeeze=False)
    palette=['#1769aa','#d47715','#138765','#975399','#c44a5a']
    uniform=len({p for p,g in variants})==1
    for row,case in enumerate(CASES):
        base=data.base(case)[1];reference(axes[row],base)
        groups=[base]
        for index,(policy,gamma) in enumerate(variants):
            s,tr=data.get(case,order,policy,gamma,retain);groups.append(tr)
            draw(axes[row],s,tr,policy,gamma,color=palette[index] if uniform else COLORS[policy],
                 label=gamma_name(gamma) if uniform else NAMES[policy])
        lim=limits(data,case,order,variants,retain)
        for col,ax in enumerate(axes[row]):
            setup(ax,lim[col],groups,col)
            if row==4:ax.set_xlabel('Time [s]')
        axes[row,0].set_title(CASE_NAMES[row],loc='left',weight='bold',fontsize=10,pad=8)
    handles,labels=axes[0,0].get_legend_handles_labels()
    return finish(fig,out,name,title+'\n'+condition(order,retain),handles,labels)


def sweep_plot(data,out,retain,h):
    fig,axes=plt.subplots(5,3,figsize=(13.5,12.5),squeeze=False)
    variants=[('fresh',1),('fresh',4),('lag_sweep',1),('frozen_step',1),('previous_step',1)]
    for row,case in enumerate(CASES):
        for policy,gamma in variants:
            summaries=[data.get(case,'EL',policy,gamma,retain,h,s)[0] for s in [1,5,20]]
            color='#bd424f' if policy=='fresh' and gamma==4 else COLORS[policy]
            fields=[[max(s['max_abs_energy_error']/s['Escale'],1e-15) for s in summaries],
                    [max(s['rms_position_error'],1e-15) for s in summaries],
                    [s['oversized_skips']/max(s['energy_attempts'],1) for s in summaries]]
            for col,ax in enumerate(axes[row]):
                accepted=[v if s['completed_steps']==s['requested_steps'] else math.nan for v,s in zip(fields[col],summaries)]
                ax.plot([1,5,20],accepted,'o-',color=color,lw=1.1,ms=4,label=variant_name(policy,gamma))
                for sweeps,v,s in zip([1,5,20],fields[col],summaries):
                    if s['completed_steps']<s['requested_steps']:
                        ax.plot(sweeps,v,'x',color=color,ms=8)
                        ax.annotate(f'{s["completed_steps"]*h:.2g}s',(sweeps,v),fontsize=6,color=color)
        axes[row,0].set_title(CASE_NAMES[row],loc='left',weight='bold',fontsize=10,pad=8)
        for col,ax in enumerate(axes[row]):
            ax.set_xscale('log');ax.set_xticks([1,5,20],['1','5','20']);ax.grid(True,alpha=.18)
            ax.set_ylabel(['Max |E−E₀| / Escale','RMS position error [m]','Rejected fraction'][col],fontsize=8)
            if col<2:ax.set_yscale('log')
            else:ax.set_ylim(-.04,1.04)
            ax.tick_params(labelsize=7)
            if row==4:ax.set_xlabel('Sweeps per timestep')
    fig.suptitle(f'Iteration count: original endpoint-energy solver\n{condition("EL",retain)} · h=1/{round(1/h)} s',fontsize=14,y=.995)
    handles,labels=axes[0,0].get_legend_handles_labels()
    fig.legend(handles,labels,loc='upper center',bbox_to_anchor=(.5,.953),ncol=3,frameon=False,fontsize=8)
    fig.text(.5,.015,'Lines connect completed ten-second runs only. Crosses show guarded prefixes, annotated with completed duration.\n'
             'Higher sweep counts increase solve attempts; these are equal-duration, not equal-cost, comparisons. Error axes are logarithmic; display floor 10⁻¹⁵.',ha='center',fontsize=8)
    fig.tight_layout(rect=(.005,.065,.995,.925),h_pad=1.6,w_pad=1.6)
    filename=f'sweeps_EL_h{round(1/h)}_{"retain" if retain else "drop"}_rx.png'
    fig.savefig(out/filename,dpi=165,facecolor='white');plt.close(fig)
    return filename


def compact_springs(data,out):
    fig,axes=plt.subplots(2,2,figsize=(12,7.2),squeeze=False)
    for row,case in enumerate(['spring_rest','spring_radial']):
        baseline,base=data.base(case);times=[r['time'] for r in base]
        target=base[0]['target_energy']
        axes[row,0].plot(times,[r['ref_energy']/target for r in base],color='#333a43',ls=':',lw=1.1,label='RK4 reference')
        axes[row,1].plot(times,[r['ref_x'] for r in base],color='#333a43',ls=':',lw=1.2,label='RK4 reference')
        for col,field in enumerate(['energy','x']):
            axes[row,col].plot(times,[r[field]/target if col==0 else r[field] for r in base],
                               color='#aeb4bb',lw=1.15,label='Ordinary XPBD')
        summaries=[]
        for policy,gamma in [('fresh',4),('frozen_step',1)]:
            s,tr=data.get(case,'EL',policy,gamma,0);summaries.append(s)
            for col,field in enumerate(['energy','x']):
                axes[row,col].plot([r['time'] for r in tr],[r[field]/target if col==0 else r[field] for r in tr],
                                   color=COLORS[policy],ls='-' if policy=='fresh' else '--',lw=1.4,
                                   label='Current gradient · γ=4' if policy=='fresh' else 'Predictor gradient · γ=1')
                if s['completed_steps']<s['requested_steps']:
                    value=tr[-1][field]/target if col==0 else tr[-1][field]
                    axes[row,col].plot(tr[-1]['time'],value,'x',color=COLORS[policy],ms=7)
        axes[row,0].set_title(CASE_NAMES[row],loc='left',weight='bold',fontsize=11,pad=10)
        fresh,frozen=summaries
        axes[row,1].set_title(f'Current: {int(fresh["oversized_skips"])} rejected, {int(fresh["frozen_steps"])} immobile steps; '
                              f'predictor: {int(frozen["oversized_skips"])} rejected, {int(frozen["frozen_steps"])} immobile steps',
                              loc='left',fontsize=8,pad=10)
        for col,ax in enumerate(axes[row]):
            ax.set_xlim(0,10);ax.grid(True,alpha=.18);ax.tick_params(labelsize=8)
            ax.set_ylabel('Mechanical energy / initial energy' if col==0 else 'Position x [m]',fontsize=9)
            if row==1:ax.set_xlabel('Time [s]')
        axes[row,0].set_ylim(-.05,1.45)
        axes[row,1].set_ylim(.69,1.32)
    handles,labels=axes[0,0].get_legend_handles_labels()
    fig.suptitle('Spring: rejected energy corrections and a frozen gradient\nOriginal endpoint-energy formulation · energy → material · drop rₓ',fontsize=14,y=.995)
    fig.legend(handles,labels,loc='upper center',bbox_to_anchor=(.5,.912),ncol=4,frameon=False,fontsize=8.5)
    fig.text(.5,.018,'h = 1/60 s; five sweeps; ten seconds. Current gradient uses γ=4; predictor gradient stays fixed within each timestep and uses γ=1.\n'
             'γ sets the maximum allowed energy-correction displacement relative to the timestep motion scale. Oversized corrections are rejected, not clipped.\n'
             'Velocity remains v = (x − xₙ)/h. Immobile means displacement below 10⁻¹⁰ m during a timestep; all four trajectories complete ten seconds.',
             ha='center',fontsize=8)
    fig.tight_layout(rect=(.005,.10,.995,.88),h_pad=1.8,w_pad=1.8)
    filename='spring_skip_vs_old_gradient.png'
    fig.savefig(out/filename,dpi=180,facecolor='white');plt.close(fig)
    return filename


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('results',type=Path)
    parser.add_argument('--output',type=Path,help='Default: RESULTS/plots')
    args=parser.parse_args()
    if args.output is None:args.output=args.results/'plots'
    args.output.mkdir(parents=True,exist_ok=True)
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.spines.top':False,'axes.spines.right':False,
                         'axes.edgecolor':'#c0c4c8','text.color':'#252c33','axes.labelcolor':'#333d46'})
    data=Data(args.results);figures=[]
    for retain in [0,1]:
        state='retain' if retain else 'drop'
        for policy in POLICIES:figures.append(policy_plot(data,args.output,'EL',policy,1,retain))
        figures.append(policy_plot(data,args.output,'EL','fresh',4,retain))
        for gamma in [-1,1]:
            figures.append(comparison(data,args.output,'EL',[(p,gamma) for p in POLICIES],retain,
                f'policies_EL_gamma{slug_gamma(gamma)}_{state}_rx','Gradient timing · '+gamma_name(gamma)))
        figures.append(comparison(data,args.output,'LE',[(p,1) for p in POLICIES],retain,
            f'policies_LE_gamma1_{state}_rx','Gradient timing · γ=1'))
        figures.append(comparison(data,args.output,'EL',[('fresh',g) for g in [-1,.25,1,4,16]],retain,
            f'fresh_gamma_EL_{state}_rx','Current gradient · rejection threshold'))
        for h in [1/60,1/240]:figures.append(sweep_plot(data,args.output,retain,h))
    figures.append(compact_springs(data,args.output))
    representative=[r for r in data.rows if close(r['h'],1/60) and r['sweeps']==5]
    with (args.output/'representative.csv').open('w',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(representative[0]));writer.writeheader();writer.writerows(representative)
    audit={'summary_runs':len(data.rows),'trajectory_groups':len(data.traces),
           'max_scaled_endpoint_energy_reconstruction_difference':data.energy_check,'figures':figures,
           'all_time_panels':'Actual mechanical energy, Cartesian x(t), Cartesian y(t). Endpoint velocity is (x-x_n)/h.',
           'gradient_definitions':{'fresh':'Current analytical gradient at the start of the sweep.',
             'lag_sweep':'Gradient from the start of the preceding sweep; initialized from predictor.',
             'frozen_step':'Gradient at the timestep predictor; frozen throughout that timestep.',
             'previous_step':'Gradient at the end of the preceding timestep; first timestep falls back to fresh.'},
           'rejection':'Reject only the energy increment when its displacement exceeds gamma*(h*|v_n|+h^2*|f_n|/m); negative gamma disables rejection.',
           'physical_stationarity':'Reported r_x uses the current analytical energy gradient, including stale-gradient experiments.'}
    (args.output/'plot_audit.json').write_text(json.dumps(audit,indent=2)+'\n')
    print(json.dumps({k:v for k,v in audit.items() if k!='figures'},indent=2))
    print(f'Generated {len(figures)} PNG figures in {args.output}')


if __name__=='__main__':main()
