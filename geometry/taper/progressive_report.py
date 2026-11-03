"""Readable matched-control evidence and a compact machine-readable index."""
import json
import numpy as np
from geometry.taper.progressive_fixtures import OUT,NATIVE


def main():
    case=OUT/'07_plain_T2_direct_regular'
    items=[(case/'junction_check','Physical Wire: Original LC Model'),
           (case/'mutual_diagnostic','Physical Wire: Coupled-L Diagnostic'),
           (case/'paired_ideal/twoC','Ideal Inductor: Original LC Model')]
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(1,3,figsize=(18,6),sharex=True,sharey=True)
    for ax,(folder,title) in zip(axes,items):
        data=np.load(folder/'response.npz');r=json.loads((folder/'fit.json').read_text())
        actual=data['actual'];model=data['model']
        if actual.ndim==3:actual=actual[:,1,0];model=model[:,1,0]
        ax.plot(data['f']/1e9,abs(actual),color='#b22222',lw=2.6,label='EM response')
        ax.plot(data['f']/1e9,abs(model),'k--',lw=2.1,label='Complex fit')
        text=f"Complex RMS = {r['complex_RMS']:.5f}"
        if 'mutual_nH' in r:text+=f"\nEffective M = {r['mutual_nH']:.4f} nH"
        ax.text(.03,.96,text,transform=ax.transAxes,va='top',fontsize=13,bbox=dict(boxstyle='round,pad=.4',fc='white',ec='gray'))
        ax.set_title(title,fontsize=16,pad=13);ax.set_xlabel('Frequency (GHz)',fontsize=16)
        ax.tick_params(labelsize=13);ax.grid(alpha=.25);ax.legend(fontsize=11,loc='upper right')
        ax.set_ylim(0,1.22)
    axes[0].set_ylabel(r'$|S_{21}|$',fontsize=22)
    fig.suptitle('Matched 0.733 nH T2 Controls: Same Fixed Downstream Circuit',fontsize=22)
    fig.tight_layout();path=OUT/'matched_controls_and_mutual.png';fig.savefig(path,dpi=190);plt.close(fig)
    rows=[]
    for p in OUT.rglob('target_review.json'):
        r=json.loads(p.read_text())
        rows.append(dict(case=str(p.parent.relative_to(OUT)),L_nH=r['L_nH'][0],C0_T2_fF=r['C_fF'][0],
                         input_Cg_fF=r['G_fF'][0],complex_RMS=r['complex_RMS'],
                         target_errors_percent=r['target_errors_percent'],active_bounds=r['active_bounds'],
                         project=str(NATIVE/p.parent.relative_to(OUT)/'generated_ncell_3.sonx'),accepted=False))
    summary=dict(final_T2_design=None,geometry_trials=rows,matched_control_plot=str(path),
                 scope='No design is accepted merely for fitting a curve; component targets and physical validation are required')
    (OUT/'review_index.json').write_text(json.dumps(summary,indent=2)+'\n')
    return summary


if __name__=='__main__':main()
