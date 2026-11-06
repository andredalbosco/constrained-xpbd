#!/usr/bin/env python3
"""Audit and graph EL/LE with frozen or refreshed rows."""
import argparse, csv, json, math, os
from collections import defaultdict
from pathlib import Path
os.environ.setdefault('MPLCONFIGDIR','/tmp/xpbd-traditional-matplotlib')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import SymmetricalLogLocator
from matplotlib.lines import Line2D
import numpy as np
from plot_results import energy as physical_energy
from plot_el import read

CASE_NAMES={'spring_release':'Spring released from rest','pendulum_rest':'Rigid pendulum released from rest',
            'spring_orbit':'Spring with tangential velocity','pendulum_kick':'Rigid pendulum with tangential velocity'}
COLOR={'xpbd':'#87929d','EL':'#ce6b20','LE':'#167b9b'}

def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('results',type=Path);p.add_argument('--output',type=Path,required=True)
 a=p.parse_args();a.output.mkdir(parents=True,exist_ok=True)
 models={r['case']:r for r in read(a.results/'cases.csv')}
 summary={(r['case'],r['method'],int(r['hz'])):r for r in read(a.results/'summary.csv')}
 groups=defaultdict(list)
 for r in read(a.results/'traces.csv'):groups[r['case'],r['method']].append(r)
 audit={'runs':len(summary),'trace_groups':len(groups),'max_scaled_energy_error':0.,'max_rms_error':0.}
 for (c,m),tr in groups.items():
  s=summary[c,m,480];assert s['completed']==s['steps']==len(tr)-1,(c,m,s['status'])
  assert all(r['step']==i and math.isclose(r['time'],i/480,abs_tol=1e-14) for i,r in enumerate(tr))
  err=max(abs(physical_energy(models[c],r)-r['energy'])/s['energy_scale'] for r in tr)
  rms=math.sqrt(sum((r['x1']-r['ref_x1'])**2+(r['y1']-r['ref_y1'])**2 for r in tr[1:])/(len(tr)-1))
  emax=max(abs(r['energy']-s['initial_energy'])/s['energy_scale'] for r in tr[1:])
  assert err<1e-9 and math.isclose(rms,s['rms_position'],abs_tol=1e-12,rel_tol=1e-12)
  assert math.isclose(emax,s['max_energy_rel'],abs_tol=1e-12,rel_tol=1e-12)
  audit['max_scaled_energy_error']=max(audit['max_scaled_energy_error'],err)
  audit['max_rms_error']=max(audit['max_rms_error'],abs(rms-s['rms_position']))
 plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False})
 def draw(ax,c,m,metric='energy',linewidth=1.2):
  tr=groups[c,m];s=summary[c,m,480];t=np.array([r['time'] for r in tr]);color=COLOR['xpbd' if m=='xpbd' else m[:2]]
  if metric=='energy':v=np.array([100*abs(r['energy']-s['initial_energy'])/s['energy_scale'] for r in tr])
  elif metric=='error':v=np.array([math.hypot(r['x1']-r['ref_x1'],r['y1']-r['ref_y1']) for r in tr])
  elif c.startswith('spring'):v=np.array([r['x1']-(models[c]['l1'] if c=='spring_release' else 0) for r in tr])
  else:v=np.array([math.atan2(r['x1'],-r['y1']) for r in tr])
  ax.plot(t,v,color=color,lw=linewidth,ls='--' if '_refresh_' in m else '-')
 def reference(ax,c):
  tr=groups[c,'xpbd'];t=[r['time'] for r in tr]
  if c.startswith('spring'):v=[r['ref_x1']-(models[c]['l1'] if c=='spring_release' else 0) for r in tr]
  else:v=[math.atan2(r['ref_x1'],-r['ref_y1']) for r in tr]
  ax.plot(t,v,color='#171717',ls=':',lw=1.2)
 def energy_axis(ax):
  ax.set_yscale('symlog',linthresh=1e-4);loc=SymmetricalLogLocator(base=10,linthresh=1e-4);loc.set_params(numticks=7)
  ax.yaxis.set_major_locator(loc);ax.set_ylim(bottom=0);ax.set_ylabel('Absolute energy error (%) · symlog')
 def decorate(fig,axes,title,rebuilt=False):
  fig.subplots_adjust(top=.80,bottom=.17,left=.08,right=.98,hspace=.48,wspace=.22)
  fig.suptitle(title,x=.08,y=.975,ha='left',fontsize=19,fontweight='bold')
  fig.text(.08,.933,'480 Hz · 10 s · one pass per small step · all multipliers reset each step',fontsize=11)
  handles=[Line2D([],[],color=COLOR['xpbd'],lw=2,label='Ordinary XPBD'),Line2D([],[],color=COLOR['EL'],lw=2,label='EL'),Line2D([],[],color=COLOR['LE'],lw=2,label='LE'),Line2D([],[],color='black',ls=':',label='Reference motion')]
  if rebuilt:handles += [Line2D([],[],color='#555',ls='-',label='Solid: frozen rows'),Line2D([],[],color='#555',ls='--',label='Dashed: rebuild next row')]
  fig.legend(handles=handles,loc='upper left',bbox_to_anchor=(.075,.91),frameon=False,ncol=3 if rebuilt else 4)
  for ax in axes.flat:ax.set_xlim(0,10);ax.set_xlabel('Time (s)');ax.grid(alpha=.2)
 def note(c,policy,refresh=False):
  names=['EL_'+('refresh_' if refresh else '')+policy,'LE_'+('refresh_' if refresh else '')+policy]
  return ';  '.join(m.split('_')[0]+f": max E {100*summary[c,m,480]['max_energy_rel']:.3g}%, position RMS {100*summary[c,m,480]['rms_position']:.3g} cm" for m in names)
 for c in CASE_NAMES:
  for rebuilt in (False,True):
   fig,axs=plt.subplots(2,2,figsize=(14,9))
   decorate(fig,axs,CASE_NAMES[c]+(': rebuild after each visit' if rebuilt else ': EL versus LE'),rebuilt)
   for row,policy in enumerate(['raw','repair_cap']):
    methods=['xpbd','EL_'+policy,'LE_'+policy]
    if rebuilt:methods+=['EL_refresh_'+policy,'LE_refresh_'+policy]
    for m in methods:draw(axs[row,0],c,m);draw(axs[row,1],c,m,metric='position')
    reference(axs[row,1],c);energy_axis(axs[row,0])
    title='Raw updates' if row==0 else 'Bounded repair + cap'
    axs[row,0].set_title(title+' — energy',loc='left',fontweight='bold');axs[row,1].set_title(title+' — motion',loc='left',fontweight='bold')
    axs[row,1].set_ylabel('x − rest length (m)' if c=='spring_release' else 'x coordinate (m)' if c=='spring_orbit' else 'Angle from vertical (rad)')
   text=['Frozen raw: '+note(c,'raw'),'Frozen bounded: '+note(c,'repair_cap')]
   if rebuilt:text=['Rebuilt raw: '+note(c,'raw',True),'Rebuilt bounded: '+note(c,'repair_cap',True)]
   text+=['Energy error uses initial energy above minimum potential. Energy panels have separate vertical ranges; spikes are retained.',
          'Rebuilt visits still solve linearized rows, drop curvature and r_x; bounded policy limits repair to 50% and extra motion to 10%.']
   fig.text(.08,.035,'\n'.join(text),fontsize=8.5,color='#333')
   name=c+('_rebuild_comparison.png' if rebuilt else '_EL_LE.png');fig.savefig(a.output/name,dpi=160,facecolor='white');plt.close(fig)
 # Compact figure for the source review: both cases, bounded results, all row variants.
 fig,axs=plt.subplots(2,2,figsize=(14,9));decorate(fig,axs,'Reconstruct positions, then rebuild the next row',True)
 for row,c in enumerate(['spring_release','pendulum_rest']):
  for m in ['xpbd','EL_repair_cap','LE_repair_cap','EL_refresh_repair_cap','LE_refresh_repair_cap']:
   draw(axs[row,0],c,m);draw(axs[row,1],c,m,'error')
  axs[row,1].axhline(0,color='#171717',ls=':',lw=1.2)
  energy_axis(axs[row,0]);axs[row,1].set_ylabel('Position error vs reference (m)');axs[row,1].set_ylim(bottom=0)
  axs[row,0].set_title(CASE_NAMES[c]+' — bounded energy repair',loc='left',fontweight='bold')
  axs[row,1].set_title('Motion accuracy',loc='left',fontweight='bold')
 fig.text(.08,.035,'Spring: rebuilding EL leaves the result unchanged; rebuilding LE improves bounded RMS from 0.959 cm to 0.500 cm.\n'
          'Pendulum: rebuilding EL leaves the result almost unchanged; rebuilding LE changes bounded RMS from 0.491 cm to 0.748 cm.\n'
          'Every visit is a linearized solve; this comparison retains the same whole-proposal filter and drops curvature and r_x.\n'
          'This figure contains bounded methods. Full raw-update plots preserve the spring energy spikes and frozen LE branch.',fontsize=9,color='#333')
 fig.savefig(a.output/'reconstruction_overview.png',dpi=160,facecolor='white');plt.close(fig)
 (a.output/'audit.json').write_text(json.dumps(audit,indent=2)+'\n');print(json.dumps(audit,indent=2))
 print('Figures:',len(list(a.output.glob('*.png'))))
if __name__=='__main__':main()
