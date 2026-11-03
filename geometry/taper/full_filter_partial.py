"""Plot completed direct EM points while ABS is still running. No EM launch."""
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from geometry.taper.full_filter import OUT,circuit_values,ladder


def main():
    from pathlib import Path
    plan=json.loads(Path('geometry/taper/full_filter_work_plan.json').read_text(encoding='utf-8-sig'))
    rows=plan['completed_direct_samples']
    selection=json.loads((OUT/'selected_taper3_source.json').read_text(encoding='utf-8-sig'))
    f=np.linspace(4.65e9,5.4e9,7501)
    fig,axes=plt.subplots(2,1,figsize=(14,10),sharex=True)
    comparisons={}
    for label,nominal in [('Nominal target circuit',True),('Extracted cell circuit',False)]:
        ls,cs,gs=circuit_values(selection,nominal=nominal)
        mag=abs(ladder(f,ls,cs,gs)[:,1,0])
        axes[0].plot(f/1e9,mag,'--',lw=1.6,label=label)
        axes[1].plot(f/1e9,20*np.log10(np.maximum(mag,1e-8)),'--',lw=1.6,label=label)
        point_mag=abs(ladder(np.array([r['f_GHz'] for r in rows])*1e9,ls,cs,gs)[:,1,0])
        comparisons[label]=[float(x) for x in point_mag]
    axes[0].plot([r['f_GHz'] for r in rows],[r['S21_magnitude'] for r in rows],'o',ms=8,color='black',label='Completed direct Sonnet points')
    axes[1].plot([r['f_GHz'] for r in rows],[r['S21_dB'] for r in rows],'o',ms=8,color='black',label='Completed direct Sonnet points')
    axes[0].set_title('Full Filter: 16 Direct Sonnet Samples (Dense Curve Unverified)',fontsize=21)
    axes[0].set_ylabel(r'$|S_{21}|$',fontsize=22);axes[0].set_ylim(0,1.1)
    axes[1].set_ylabel('S21 (dB)',fontsize=20);axes[1].set_xlabel('Frequency (GHz)',fontsize=20);axes[1].set_ylim(-65,2)
    for ax in axes:ax.tick_params(labelsize=14);ax.grid(alpha=.2);ax.legend(fontsize=12)
    fig.tight_layout();fig.savefig(OUT/'full_filter_partial_direct_samples.png',dpi=170);plt.close(fig)
    report=dict(status='PARTIAL: direct samples valid; cache-only dense interpolation rejected',direct_samples=rows,circuit_prediction_magnitudes_at_these_frequencies=comparisons,warning='Native samples are discrete points, not a converged full-filter curve.')
    (OUT/'full_filter_partial_direct_samples.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(report,indent=2))

if __name__=='__main__':main()
