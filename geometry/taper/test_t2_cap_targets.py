"""Circuit-only counterfactual: set the three T2 capacitances to their targets.

Keep every other extracted circuit value, all inductances, and the fitted
fixture delay unchanged. Compare against the fully nominal target circuit.
"""
from pathlib import Path
from hashlib import sha256
import json
import numpy as np
from scipy.optimize import minimize_scalar
from scipy.signal import find_peaks
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from circuit_common import read_snp
from geometry.taper.progressive_fitting import transmission
from geometry.taper.staged_fitting import circuit
from geometry.taper.taper2_grid4 import OLD_FIT
from geometry.specifications import TAPER_2,TAPER_1,UNIT_CELL


def main():
    response=OLD_FIT.parent.parent
    fit_path=response/'extraction_three_caps_adopted_reference/fit.json'
    fit=json.loads(fit_path.read_text());source=Path(fit['provenance']['source'])
    if sha256(source.read_bytes()).hexdigest()!=fit['provenance']['sha256']:raise ValueError('Data changed')
    net=read_snp(source);f=net.f;observed=net.s[:,1,0]
    ls=list(fit['L_nH']);cs=list(fit['C_fF']);gs=list(fit['G_fF']);delay=fit['delay_ps']
    cs[0]=TAPER_2.C0_fF;gs[0]=TAPER_2.Cg_input_fF;gs[1]=TAPER_1.Cg_input_fF
    hypothetical=transmission(f,ls,cs,gs,delay)
    # Independently check the established ABCD implementation and index mapping.
    check=circuit(f,ls,cs,gs,left_ps=delay/2,right_ps=delay/2)[:,1,0]
    if np.max(abs(check-hypothetical))>1e-10:raise ValueError('Circuit implementations disagree')
    targets=[TAPER_2,TAPER_1,UNIT_CELL,UNIT_CELL]
    nominal=transmission(f,[t.L0_nH for t in targets],[t.C0_fF for t in targets],
                         [t.Cg_input_fF for t in targets],0)
    def score(curve):
        objective=lambda t:float(np.mean(abs(curve-nominal*np.exp(-2j*np.pi*f*t*1e-12))**2))
        grid=np.linspace(-100,100,801);i=int(np.argmin([objective(t) for t in grid]))
        q=minimize_scalar(objective,bounds=(max(-100,grid[i]-.5),min(100,grid[i]+.5)),method='bounded')
        return dict(complex_RMS=float(np.sqrt(q.fun)),target_fixture_delay_ps=float(q.x))
    new_score=score(hypothetical);old_score=score(observed)
    same_delay=nominal*np.exp(-2j*np.pi*f*delay*1e-12)
    current=transmission(f,fit['L_nH'],fit['C_fF'],fit['G_fF'],delay)
    result=dict(scope='Circuit prediction only; not a simulated or achieved geometry',
        source_fit=str(fit_path),source_fit_sha256=sha256(fit_path.read_bytes()).hexdigest(),
        changes=[dict(name='T2 C0',old=fit['C_fF'][0],new=cs[0],unit='fF'),
                 dict(name='Input-to-T2 Cg',old=fit['G_fF'][0],new=gs[0],unit='fF'),
                 dict(name='T2-to-T1 Cg',old=fit['G_fF'][1],new=gs[1],unit='fF')],
        L_nH=ls,C_fF=cs,G_fF=gs,unchanged_model_fixture_delay_ps=delay,
        target_comparison=new_score,current_Sonnet_target_comparison=old_score,
        current_fitted_model_target_comparison=score(current),
        complex_RMS_with_identical_fixture_delay=float(np.sqrt(np.mean(abs(hypothetical-same_delay)**2))),
        magnitude_RMS=float(np.sqrt(np.mean((abs(hypothetical)-abs(nominal))**2))),
        cutoff=.025,frequency_range_GHz=[float(f[0]/1e9),float(f[-1]/1e9)],samples=len(f),
        target_peaks_GHz=(f[find_peaks(abs(nominal),prominence=.05)[0]]/1e9).tolist(),
        hypothetical_peaks_GHz=(f[find_peaks(abs(hypothetical),prominence=.05)[0]]/1e9).tolist(),
        ABCD_implementation_max_difference=float(np.max(abs(check-hypothetical))))
    out=response/'hypothetical_T2_capacitances_at_targets';out.mkdir(exist_ok=True)
    (out/'comparison.json').write_text(json.dumps(result,indent=2)+'\n')
    np.savez(out/'comparison.npz',f=f,target=nominal,hypothetical=hypothetical,current_Sonnet=observed)
    fig,ax=plt.subplots(figsize=(14,8))
    ax.plot(f/1e9,abs(nominal),'k',lw=2.6,label='Specified target circuit')
    ax.plot(f/1e9,abs(hypothetical),'--',color='#843ca0',lw=2.6,label='Three T2 capacitances set to targets')
    ax.plot(f/1e9,abs(observed),color='#888',lw=1.8,label='Current original Sonnet design')
    text=f'Hypothetical target RMS: {new_score["complex_RMS"]:.5f}\nCurrent Sonnet target RMS: {old_score["complex_RMS"]:.5f}\nRequired: 0.025 or below\n\nT2 C0 = 149 fF\nInput-to-T2 Cg = 153 fF\nT2-to-T1 Cg = 15.3 fF'
    ax.text(.025,.97,text,transform=ax.transAxes,va='top',fontsize=14,bbox=dict(boxstyle='round,pad=.5',fc='white',ec='gray'))
    ax.set(xlim=(4.65,5.4),ylim=(0,max(abs(hypothetical).max(),abs(nominal).max(),abs(observed).max())*1.45),
        xlabel='Frequency (GHz)',ylabel=r'$|S_{21}|$',title='Circuit Prediction: T2 Capacitances at Their Targets')
    ax.xaxis.label.set_size(20);ax.yaxis.label.set_size(23);ax.title.set_size(22);ax.tick_params(labelsize=16)
    ax.grid(alpha=.2);ax.legend(fontsize=12,loc='upper right')
    fig.text(.5,.015,'All other extracted L/C values and model delay unchanged. Only target-reference delay is matched for the complex RMS.',ha='center',fontsize=11)
    fig.tight_layout(rect=[0,.035,1,1])
    for ext in ('png','pdf'):fig.savefig(out/f's21_comparison.{ext}',dpi=190)
    plt.close(fig);print(json.dumps(result,indent=2),flush=True)


if __name__=='__main__':main()
