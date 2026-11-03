"""Plot completed direct fields for the cleanup full filter without touching EM."""
from pathlib import Path
from datetime import datetime,timezone
import json,re
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from geometry.taper.full_filter import circuit_values,ladder

CASE='full_filter_10_regular_3_tapers_each_side_T2_cleanup_2026-11-01'
OUT=Path('output_data/taper_exploration/full_filter')/CASE

def main():
 text=(OUT/'em.log').read_text(encoding='utf-8-sig');rows=[]
 for line in text.splitlines():
  values=line.split()
  if len(values)!=9:continue
  try:v=list(map(float,values))
  except ValueError:continue
  if 4.65<=v[0]<=5.4:rows.append(v)
 a=np.array(sorted(rows,key=lambda v:v[0]));assert len(a)
 selection=json.loads((OUT/'selected_taper3_source.json').read_text(encoding='utf-8-sig'))
 f=np.linspace(4.65e9,5.4e9,7501)
 fig,axes=plt.subplots(2,1,figsize=(15,10),sharex=True)
 for label,nominal,color in [('Nominal target circuit',True,'#d65a31'),('Target-bounded extracted circuit',False,'#408359')]:
  ls,cs,gs=circuit_values(selection,nominal=nominal);m=abs(ladder(f,ls,cs,gs)[:,1,0])
  axes[0].plot(f/1e9,m,'--',lw=1.8,color=color,label=label)
  axes[1].plot(f/1e9,20*np.log10(np.maximum(m,1e-8)),'--',lw=1.8,color=color,label=label)
 for ax,db in zip(axes,[False,True]):
  ax.plot(a[:,0],20*np.log10(np.maximum(a[:,3],1e-12)) if db else a[:,3],'o',color='#235d95',ms=7,label=f'Direct Sonnet calculations ({len(a)} points)',zorder=5)
  ax.grid(alpha=.2);ax.tick_params(labelsize=14);ax.legend(fontsize=12,loc='lower left' if db else 'lower center')
 axes[0].set_ylim(0,1.08);axes[1].set_ylim(-85,2)
 axes[0].set_ylabel(r'$|S_{21}|$',fontsize=23);axes[1].set_ylabel('S21 (dB)',fontsize=20);axes[1].set_xlabel('Frequency (GHz)',fontsize=20)
 axes[0].set_title('Full Filter: Completed Direct Sonnet Calculations',fontsize=23)
 axes[0].text(.02,.97,'10 regular cells + three tapers per side\nAdaptive sweep still running; dots are direct EM results\nDashed curves are circuit predictions',transform=axes[0].transAxes,va='top',fontsize=12,bbox=dict(fc='white',ec='gray',boxstyle='round,pad=.4'))
 fig.tight_layout()
 for ext in ['png','pdf']:fig.savefig(OUT/f'full_filter_partial_direct_samples.{ext}',dpi=170)
 plt.close(fig)
 report=dict(status='PARTIAL: native adaptive sweep active; no converged dense EM curve yet',updated_utc=datetime.now(timezone.utc).isoformat(),completed_direct_points=len(a),current_frequency_GHz=float(re.findall(r'Frequency:\s+([\d.]+) GHZ',text)[-1]),samples=[dict(f_GHz=v[0],S21_magnitude=v[3],S21_dB=float(20*np.log10(v[3]))) for v in a],plot=str((OUT/'full_filter_partial_direct_samples.png').resolve()))
 (OUT/'full_filter_partial_direct_samples.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
 print(json.dumps({k:v for k,v in report.items() if k!='samples'}))
if __name__=='__main__':main()
