#!/usr/bin/env python3
"""Plot the independent circular-speed holdout; keep altered directions separate."""
import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
os.environ.setdefault('MPLCONFIGDIR','/tmp/xpbd-eccentricity-matplotlib')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('csv',type=Path);p.add_argument('--output',type=Path,default=Path('results/correction_schedules/plots'))
    a=p.parse_args();a.output.mkdir(parents=True,exist_ok=True)
    with a.csv.open(newline='') as f:rows=list(csv.DictReader(f))
    text={'case','method','status'}
    for r in rows:
        for k in r:
            if k not in text:r[k]=float(r[k])
        for k in ['speed_fraction','substeps','h','rms_position_error','max_energy_percent','max_angular_error']:
            assert math.isfinite(r[k])
        assert math.isclose(r['h'],1/(60*r['substeps']))
    fractions=sorted({r['speed_fraction'] for r in rows});subs=sorted({int(r['substeps']) for r in rows})
    methods=['baseline','LE_repair_cap','angular_neutral_extra']
    assert len(rows)==len(fractions)*len(subs)*len(methods)==90
    assert len({(r['speed_fraction'],r['substeps'],r['method']) for r in rows})==len(rows)
    colors={'baseline':'#8c98a2','LE_repair_cap':'#168270','angular_neutral_extra':'#a34880'}
    labels={'baseline':'Ordinary XPBD','LE_repair_cap':'LE, fractional repair ½ + cap 0.1',
            'angular_neutral_extra':'Projected energy direction: altered equations'}
    def metric(r,i):
        if i==0:return r['rms_position_error']
        if i==1:return r['max_energy_percent']
        ell0=1.1*math.sqrt(11)*r['speed_fraction']
        return 100*r['max_angular_error']/ell0
    titles=['Position RMS error (m)','Maximum energy error (%)','Maximum angular-momentum error (%)']
    for diagnostic in [False,True]:
        plotted=['LE_repair_cap','angular_neutral_extra'] if diagnostic else ['baseline','LE_repair_cap']
        fig,axs=plt.subplots(1,3,figsize=(19,7.0))
        fig.subplots_adjust(left=.065,right=.985,top=.70,bottom=.22,wspace=.28)
        title='Direction-projection diagnostic: this changes the multiplier system' if diagnostic else 'Circular motion succeeds; radial-dominated orbits remain difficult'
        fig.suptitle(title,x=.045,y=.97,ha='left',fontsize=22,weight='bold')
        fig.text(.045,.915,'480 Hz · one frozen solve per substep · 10 seconds · initial radius 1.1 m, zero radial speed · no gravity',fontsize=11)
        fig.legend(handles=[Line2D([],[],color=colors[m],lw=2,label=labels[m]) for m in plotted],
                   loc='upper left',bbox_to_anchor=(.045,.875),ncol=2,frameon=False,fontsize=11)
        for i,ax in enumerate(axs):
            for m in plotted:
                rr=sorted([r for r in rows if r['method']==m and r['substeps']==8],key=lambda r:r['speed_fraction'])
                assert all(r['status']=='completed' for r in rr)
                ax.plot([r['speed_fraction'] for r in rr],[metric(r,i) for r in rr],color=colors[m],lw=2,marker='o',ms=5)
            ax.axvline(1,color='#51606c',ls=':',lw=1);ax.axvline(.8/math.sqrt(11),color='#6d7c88',ls=':',lw=1)
            ax.set_yscale('log');ax.set_title(titles[i],loc='left',fontsize=12)
            ax.set_xticks(fractions,['0.10','0.241','0.50','0.75','1.00','1.25'])
            ax.set_xlabel('Initial tangential speed / circular speed',fontsize=10);ax.grid(alpha=.2,which='both')
            ax.text(1,.97,'circle',ha='center',va='top',transform=ax.get_xaxis_transform(),fontsize=9,color='#51606c')
        caption='Fraction 0.241 is the original v = 0.8 m/s orbit; fraction 1 is the true circular orbit with v = √11 m/s. Every marked point is a separate completed run.'
        if diagnostic:
            caption+='\nThe projected direction enforces zero angular impulse for the extra correction only. It changes the original dual system and is not a tested schedule alternative.'
        else:
            caption+='\nEnergy normalization is E₀ (the potential minimum is zero); angular error is max|L−L₀|/|L₀|. Low energy error alone does not imply correct orbital motion.'
        fig.text(.045,.035,caption,fontsize=10,color='#465561',va='bottom')
        name='eccentricity_projected_direction_diagnostic.png' if diagnostic else 'eccentricity_holdout_480hz.png'
        fig.savefig(a.output/name,dpi=160,facecolor='white');plt.close(fig)
    audit={'source':str(a.csv),'sha256':hashlib.sha256(a.csv.read_bytes()).hexdigest(),
           'rows':len(rows),'speed_fractions':fractions,'substeps':subs,
           'missing_final_L_and_angle':'NaN fields for ordinary and LE methods are not plotted or imputed.',
           'angular_normalization':'100 max_abs_angular_error / (1.1*sqrt(11)*speed_fraction)',
           'scope':'Single spring mass, no gravity; altered-direction diagnostic is separate from unchanged-equation schedules.'}
    (a.output/'eccentricity_plot_audit.json').write_text(json.dumps(audit,indent=2)+'\n')
    print(json.dumps(audit,indent=2))


if __name__=='__main__':main()
