"""Corrected capacitance-only fits of the reviewed simple T2 layouts.

All inductances and regular/end-cell capacitances remain fixed. T1 C0 and
both adjacent couplings are free. With an IDC, T2 C0 and input Cg are free
as well. Without the IDC, T2 C0 stays zero. No hidden release of unit values.
"""
import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'
os.environ['OMP_NUM_THREADS'] = '1'
from pathlib import Path
import json
import numpy as np
from geometry.taper.progressive_fitting import fit
from geometry.taper.staged_fitting import LT1, LR
from geometry.specifications import specification

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT/'output_data/taper_exploration/corrected_iteration_2026-10-23'
REGULAR_C = 188.8781163861302
REGULAR_G = 8.241620220231225


def review(source, folder, inductance, *, idc=True, output_end=True,data_label='Sonnet simulation'):
    """Fit the explicit physical topology using one consistent fixed reference."""
    ls = [inductance, LT1, LR] + ([LR] if output_end else [])
    cs = [100 if idc else 0, 137.235, REGULAR_C] + ([REGULAR_C] if output_end else [])
    gs = [459 if idc else None, 20.771, 8.247] + ([REGULAR_G] if output_end else [])
    free = ([('C', 0), ('G', 0)] if idc else []) + [('C', 1), ('G', 1), ('G', 2)]
    seeds = [[100,459,100,56.7,8.2865],[50,300,100,56.7,8.2865],
             [150,700,100,56.7,8.2865]] if idc else [[100,56.7,8.2865]]
    result = fit(Path(source), Path(folder), ls, cs, gs, free,
                 title='Taper 2 Transmission: Free Taper Capacitances',capacitance_seeds=seeds,label=data_label)
    summary = dict(source=str(source),data_label=data_label, topology='T2 -> T1 -> regular' + (' -> output end' if output_end else ''),
                   specification=specification(), fixed_inductances_nH=ls,
                   fixed_regular_C_fF=REGULAR_C, fixed_regular_output_Cg_fF=REGULAR_G if output_end else None,
                   T2_shunt_C_fixed_zero=not idc, free_capacitances=free,
                   fitted=dict(L_T2_nH=inductance, C0_T2_fF=result['C_fF'][0], Cg_T2_fF=result['G_fF'][0],
                               C0_T1_fF=result['C_fF'][1], Cg_T1_fF=result['G_fF'][1],
                               Cg_T1_regular_fF=result['G_fF'][2]),
                   complex_RMS=result['complex_RMS'], complex_max=result['complex_max'],
                   curve_pass=result['curve_pass'],curve_acceptance=result['curve_acceptance'],
                   active_bounds=result['active_bounds'], delay_ps=result['delay_ps'],
                   accepted=False, limitation='Curve fitting does not alone establish unique physical capacitor values.')
    render(folder, result, summary, idc=idc, output_end=output_end)
    (Path(folder)/'review.json').write_text(json.dumps(summary, indent=2)+'\n')
    print('CORRECTED', Path(folder).name, json.dumps(summary['fitted']), result['complex_RMS'], flush=True)
    return summary


def render(folder, result, summary, *, idc, output_end):
    """Readable magnitude plot plus a separate comparison with the specification."""
    import matplotlib.pyplot as plt
    from geometry.taper.staged_fitting import circuit
    folder=Path(folder);data=np.load(folder/'response.npz');f=data['f']
    actual=data['actual'][:,1,0];model=data['model'][:,1,0]
    values=summary['fitted']
    targets=dict(L_T2_nH=.74,C0_T2_fF=100,Cg_T2_fF=459,Cg_T1_fF=56.7)
    errors={k:100*(values[k]/v-1) for k,v in targets.items() if values[k] is not None}
    summary['target_errors_percent']=errors
    summary['passes_four_parameter_targets']=bool(idc and all(abs(v)<=2 for v in errors.values()))
    target=circuit(f,[.74,6.06,4.9194]+([4.9194] if output_end else []),
                   [100 if idc else 0,100,190.4]+([190.4] if output_end else []),
                   [459 if idc else None,56.7,8.2865]+([8.2865] if output_end else []),
                   left_ps=result['delay_ps'])[:,1,0]
    summary['target_magnitude_RMS']=float(np.sqrt(np.mean((abs(actual)-abs(target))**2)))
    summary['target_magnitude_max']=float(np.max(abs(abs(actual)-abs(target))))
    for overlay in (False,True):
        fig,ax=plt.subplots(figsize=(13,7.5))
        ax.plot(f/1e9,abs(actual),color='#b22222',lw=2.8,label=summary.get('data_label','Sonnet simulation'))
        ax.plot(f/1e9,abs(model),'k--',lw=2.2,label='Complex S21 fit')
        if overlay:ax.plot(f/1e9,abs(target),color='#1565c0',ls=':',lw=2.2,label='Target circuit')
        rows=[f"$L_{{T2}}$ = {values['L_T2_nH']:.5f} nH ({errors['L_T2_nH']:+.2f}%)"]
        if idc:
            rows += [f"$C_{{0,T2}}$ = {values['C0_T2_fF']:.2f} fF ({errors['C0_T2_fF']:+.2f}%)",
                     f"$C_{{g,T2}}$ = {values['Cg_T2_fF']:.2f} fF ({errors['Cg_T2_fF']:+.2f}%)"]
        else:rows+=['No IDC; $C_{0,T2}$ fixed at zero']
        rows += [f"$C_{{g,T1}}$ = {values['Cg_T1_fF']:.3f} fF ({errors['Cg_T1_fF']:+.2f}%)",
                 f"Complex S21 RMS = {result['complex_RMS']:.5f}",
                 f"Free $C_{{0,T1}}$ = {values['C0_T1_fF']:.2f} fF; $C_{{g,unit}}$ = {values['Cg_T1_regular_fF']:.3f} fF",
                 'All L and regular-cell parameters fixed']
        ax.text(.025,.97,'\n'.join(rows),transform=ax.transAxes,va='top',fontsize=14,
                bbox=dict(boxstyle='round,pad=.55',fc='white',ec='#999999',alpha=.95))
        ax.set(xlabel='Frequency (GHz)',ylabel=r'$|S_{21}|$',ylim=(0,1.5),title='Taper 2 Transmission')
        ax.xaxis.label.set_size(19);ax.yaxis.label.set_size(23);ax.title.set_size(24)
        ax.tick_params(labelsize=16);ax.legend(fontsize=13,loc='upper right');ax.grid(alpha=.22)
        fig.tight_layout();fig.savefig(folder/('s21_target_comparison.png' if overlay else 's21_fit_linear.png'),dpi=190)
        if not overlay:fig.savefig(folder/'s21_fit.png',dpi=190)
        plt.close(fig)


def cached():
    OUT.mkdir(parents=True, exist_ok=True)
    cases = []
    base = ROOT/'output_data/taper_exploration/compact_centered'
    for name in ('plain_w500_wire52_idc26_r584_t24_gap6',
                 'plain_w500_wire52_idc26_span174_l440',
                 'plain_w500_wire52_idc30_span172_l504'):
        old = json.loads((base/name/'idc/fit/fit.json').read_text())
        cases.append(review(old['source'], OUT/name, old['L_T2_nH']))
    base = ROOT/'output_data/taper_exploration/progressive_extraction_2026-10-23'
    for name in ('trial_w700_wire68_idc36_span266_gap2_left600',
                 'trial_w848_wire128_idc48_span400_gap2_left600'):
        old = json.loads((base/name/'target_review.json').read_text())
        source = ROOT/'sonnet_files/taper_exploration/progressive_extraction_2026-10-23'/name/'generated_ncell_3.snp'
        cases.append(review(source, OUT/name, old['L_nH'][0], output_end=False))
    (OUT/'cached_summary.json').write_text(json.dumps(cases, indent=2)+'\n')


if __name__ == '__main__':
    cached()
