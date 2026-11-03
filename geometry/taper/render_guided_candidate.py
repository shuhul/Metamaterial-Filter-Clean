"""Readable magnitude plot with fixed-target comparison and fit heuristics."""
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from geometry.specifications import TAPER_2,TAPER_1
from geometry.taper.capacitance_guided_iteration import OUT


def render(name):
    folder=OUT/name;fit=json.loads((folder/'extraction/extraction.json').read_text())
    score=json.loads((folder/'idc/target_comparison.json').read_text())
    if score['source_sha256']!=fit['provenance']['source_sha256']:
        raise ValueError('Fit and target comparison must use exactly the same EM response')
    d=np.load(folder/'idc/target_comparison.npz');extra=np.load(folder/'extraction/extraction.npz')
    if not np.array_equal(d['f'],extra['f']):raise ValueError('Frequency grids differ')
    values=fit['values'];l=fit['provenance']['fixed_L_nH'][0]
    rows=[]
    for label,value,target,unit in [('$L_{T2}$',l,TAPER_2.L0_nH,'nH'),
            ('$C_{0,T2}$',values['C0_T2_fF'],TAPER_2.C0_fF,'fF'),
            ('$C_{g,T2}$ (input)',values['Cg_input_T2_fF'],TAPER_2.Cg_input_fF,'fF'),
            ('$C_{g,T1}$ (T2–T1)',values['Cg_T2_T1_fF'],TAPER_1.Cg_input_fF,'fF')]:
        rows.append(f'{label} = {value:.4f} {unit} ({100*(value/target-1):+.2f}%)')
    rows+=['Capacitances: effective fit estimates; not unique',
        f"Free T1 C0 = {values['C0_T1_fF']:.2f} fF; junction Cg = {values['Cg_T1_unit_fF']:.3f} fF",
        f"Target complex RMS = {score['complex_RMS']:.5f} (limit 0.05)",
        f"Circuit-fit S21 RMS = {fit['complex_RMS']['S21']:.5f}"]
    fig,ax=plt.subplots(figsize=(15,9))
    ax.plot(d['f']/1e9,abs(d['actual']),color='#2768aa',lw=3,label='Sonnet EM')
    ax.plot(d['f']/1e9,abs(d['target']),color='#cf6b2e',lw=2.6,ls='--',label='Specified Target Circuit')
    ax.plot(extra['f']/1e9,abs(extra['model'][:,1]),color='#202a35',lw=1.8,ls=':',label='Effective Circuit Fit')
    ax.text(.025,.97,'\n'.join(rows),transform=ax.transAxes,va='top',fontsize=16,
        bbox=dict(boxstyle='round,pad=.6',facecolor='white',edgecolor='#aaaaaa',alpha=.97))
    ax.set(xlim=(4.65,5.4),ylim=(0,1.6),title='Taper 2 + Taper 1 Transmission')
    ax.set_xlabel('Frequency (GHz)',fontsize=23);ax.set_ylabel(r'$|S_{21}|$',fontsize=26)
    ax.title.set_fontsize(27);ax.tick_params(labelsize=19);ax.grid(alpha=.22);ax.legend(fontsize=17,loc='upper right')
    fig.tight_layout()
    for suffix in ('png','pdf'):fig.savefig(folder/f's21_target_with_parameters.{suffix}',dpi=180)
    plt.close(fig)
    print(folder/'s21_target_with_parameters.png')


if __name__=='__main__':
    import sys
    render(sys.argv[1])
