"""Minimum-L investigation: cached ideal loads first, physical meanders second.

No IDC, unchanged T1/downstream circuit. These are diagnostic fixtures, not
new taper targets. Run from the repository with python -m geometry.taper.minimum_inductance.
"""
import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'
os.environ['OMP_NUM_THREADS'] = '1'
import json
from pathlib import Path
import numpy as np
import skrf
from scipy.signal import find_peaks
from circuit_common import read_snp
from geometry.taper.compact_stages import ROOT
from geometry.taper.ideal_wire_control import load_inductor
from geometry.taper.compact_fit import fit_direct
from geometry.taper.run_bridged_stage import digest
from geometry.taper.response_cases import W

OUT = ROOT / 'output_data/taper_exploration/minimum_inductance_2026-10-23'


def half_power_features(f, s21):
    """Absolute half-power crossings; never confuse half-prominence with FWHM.

    A neighboring peak can change prominence even when its local width is
    unchanged. Return None if the band does not contain both crossings.
    """
    power=abs(s21)**2
    peaks,_=find_peaks(abs(s21),prominence=.025)
    rows=[]
    for i in peaks:
        level=power[i]/2
        left=np.flatnonzero(power[:i] <= level)
        right=np.flatnonzero(power[i+1:] <= level)
        width=None
        if len(left) and len(right):
            j=left[-1];k=i+1+right[0]
            low=f[j]+(level-power[j])/(power[j+1]-power[j])*(f[j+1]-f[j])
            high=f[k-1]+(level-power[k-1])/(power[k]-power[k-1])*(f[k]-f[k-1])
            # If the half-power region contains another maximum, it is a
            # merged region, not an independently resolved peak linewidth.
            if not any(j < other < k and other != i for other in peaks):
                width=float((high-low)/1e6)
        rows.append(dict(GHz=float(f[i]/1e9),magnitude=float(abs(s21[i])),power_FWHM_MHz=width))
    return rows


def assess(r, arrays=None):
    """Preregistered curve agreement, separate from geometry/parameter targets."""
    a, b = r['peaks'], r['model_peaks']
    if arrays is not None:
        a=half_power_features(arrays['f'],arrays['actual'][:,1,0])
        b=half_power_features(arrays['f'],arrays['model'][:,1,0])
    pairs = (len(a) == len(b) and len(a) >= 2
             and all(p['power_FWHM_MHz'] is not None for p in a+b))
    errors = []
    if pairs:
        for p, q in zip(a, b):
            errors.append(dict(
                height_relative=abs(p['magnitude']-q['magnitude'])/p['magnitude'],
                width_relative=abs(p['power_FWHM_MHz']-q['power_FWHM_MHz'])/p['power_FWHM_MHz'],
                frequency_in_linewidths=abs(p['GHz']-q['GHz'])*1000/p['power_FWHM_MHz']))
    passed = (r['S21_complex_RMSE'] < .01 and r['S21_max_complex_error'] < .05
              and not any(r['active_bounds']) and pairs
              and all(e['height_relative'] < .1 and e['width_relative'] < .1
                      and e['frequency_in_linewidths'] < .1 for e in errors))
    return dict(curve_agreement_pass=bool(passed), peak_errors=errors,
                observed_half_power_peaks=a,model_half_power_peaks=b,
                criteria=dict(complex_RMS_max=.01, complex_pointwise_max=.05,
                              peak_height_relative_max=.1, peak_width_relative_max=.1,
                              peak_frequency_linewidth_max=.1, no_active_bounds=True),
                scope='Curve fit only, 4.65-5.4 GHz. Not final T2 target acceptance.')


def ideal_sweep(values=(5, 4, 3, 2, 1.5, 1, .74)):
    source=ROOT/'sonnet_files/taper_exploration/ideal_wire_diagnostic/wire_terminal_control.snp'
    control=json.loads((ROOT/'output_data/taper_exploration/ideal_wire_diagnostic/open_control.json').read_text())
    if not control['control_passed']:
        raise ValueError('Open control must pass before ideal loading')
    n=W['read_network'](source,4)
    rows=[]
    for inductance in values:
        folder=OUT/'ideal'/f'L{inductance:g}nH'
        folder.mkdir(parents=True,exist_ok=True)
        marker=folder/'completed.json'
        provenance=dict(operator_sha256=digest(source), L_nH=inductance)
        if marker.exists():
            saved=json.loads(marker.read_text())
            if saved['provenance'] != provenance:
                raise ValueError('Changed ideal operator')
            saved['assessment']=assess(saved['fit'],np.load(folder/'fit/response.npz'))
            marker.write_text(json.dumps(saved,indent=2)+'\n')
            rows.append(saved)
            continue
        s=load_inductor(n.s,n.f,inductance)
        path=folder/'ideal_loaded.s2p'
        skrf.Network(frequency=skrf.Frequency.from_f(n.f,unit='hz'),s=s,z0=50).write_touchstone(str(path))
        r=fit_direct(path,folder/'fit',inductance,data_label=f'EM + ideal {inductance:g} nH (no IDC)')
        row=dict(provenance=provenance, assessment=assess(r,np.load(folder/'fit/response.npz')), fit=r)
        marker.write_text(json.dumps(row,indent=2)+'\n')
        rows.append(row)
        print('IDEAL SUMMARY',inductance,row['assessment']['curve_agreement_pass'],r['S21_complex_RMSE'],flush=True)
    (OUT/'ideal_sweep.json').write_text(json.dumps(rows,indent=2)+'\n')
    return rows


def review():
    """Regenerate readable plots/tables from existing results; no EM calls."""
    import csv
    import matplotlib.pyplot as plt
    completed=sorted((OUT/'ideal').glob('*/completed.json'))
    rows=[json.loads(p.read_text()) for p in completed]
    rows.sort(key=lambda r:r['provenance']['L_nH'])
    fig,axes=plt.subplots(2,2,figsize=(15,10),sharex=True)
    for ax,l in zip(axes.flat,(.74,1.5,3,5)):
        folder=OUT/'ideal'/f'L{l:g}nH'
        a=np.load(folder/'fit/response.npz');r=json.loads((folder/'completed.json').read_text())['fit']
        ax.plot(a['f']/1e9,abs(a['actual'][:,1,0]),color='#b22222',lw=2.5,label='EM + ideal inductor')
        ax.plot(a['f']/1e9,abs(a['model'][:,1,0]),'k--',lw=2,label='Fixed downstream model + fit')
        ax.set_title(f'Ideal L = {l:g} nH; complex RMS = {r["S21_complex_RMSE"]:.4f}',fontsize=16)
        ax.set_xlabel('Frequency (GHz)',fontsize=15);ax.set_ylabel(r'$|S_{21}|$',fontsize=18)
        ax.tick_params(labelsize=13);ax.grid(alpha=.25);ax.legend(fontsize=10)
    fig.suptitle('No-IDC Test: Ideal-Inductance Sweep',fontsize=23)
    fig.tight_layout();fig.savefig(OUT/'ideal_comparison.png',dpi=180);plt.close(fig)
    fig,ax=plt.subplots(figsize=(10,6))
    for key,label in [('S21_complex_RMSE','Complex RMS error'),('S21_max_complex_error','Maximum complex error')]:
        ax.plot([r['provenance']['L_nH'] for r in rows],[r['fit'][key] for r in rows],'o-',label=label,lw=2)
    ax.set(xlabel='Ideal inductance (nH)',ylabel='Absolute S21 error',title='Fit Error vs. Ideal Inductance (No IDC)')
    ax.tick_params(labelsize=14);ax.legend(fontsize=13);ax.grid(alpha=.25)
    fig.tight_layout();fig.savefig(OUT/'ideal_error_vs_L.png',dpi=180);plt.close(fig)
    with (OUT/'ideal_summary.csv').open('w',newline='') as stream:
        fields=['L_nH','complex_RMS','complex_max','bare_C_fF','side_Cg_fF','strict_pass']
        writer=csv.DictWriter(stream,fieldnames=fields);writer.writeheader()
        for r in rows:
            f=r['fit'];writer.writerow(dict(L_nH=r['provenance']['L_nH'],complex_RMS=f['S21_complex_RMSE'],
                complex_max=f['S21_max_complex_error'],bare_C_fF=f['Cparasitic_bare_T2_fF'],
                side_Cg_fF=f['Cg_T2_T1_fF'],strict_pass=r['assessment']['curve_agreement_pass']))


if __name__ == '__main__':
    ideal_sweep()
    review()
