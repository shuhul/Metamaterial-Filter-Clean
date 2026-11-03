"""Saved-response staged bounds study; no EM jobs and no reference updates."""
import json
import time
from geometry.taper.taper3 import OUT, NATIVE, REFERENCE_FIT, require_reference, digest
from geometry.taper.progressive_fitting import fit

CASE = 'regular_w80_n72_dl476_dr500_padL460_R24_top2_i84ni33_stub200_FIXED'

def main():
    require_reference()
    root = OUT / CASE
    destination = root / 'response/staged_bounds'
    destination.mkdir(parents=True, exist_ok=True)
    snp = NATIVE / CASE / 'response/generated_ncell_4.snp'
    ls = [json.loads((root/'inductance.json').read_text())['L0_nH']] + REFERENCE_FIT['L_nH']
    cs = [251.] + REFERENCE_FIT['C_fF']
    gs = [261.] + REFERENCE_FIT['G_fF']
    strict = json.loads((root/'response/fit_fixed_T2junction/fit.json').read_text())
    rows = []
    started = time.perf_counter()
    previous = strict
    for percent in (1, 2, 3):
        for count, name in ((1, 'T2'), (2, 'T2_T1'), (3, 'T2_T1_regular')):
            free = [('C', 0), ('G', 0)] + [(kind, i) for i in range(1, count+1) for kind in ('C', 'G')]
            bounds = [(0., 1000.), (.1, 3000.)]
            for kind, i in free[2:]:
                ref = (cs if kind == 'C' else gs)[i]
                bounds.append((ref*(1-percent/100), ref*(1+percent/100)))
            seeds = [[r['C_fF'][i] if k == 'C' else r['G_fF'][i] for k,i in free] for r in (strict, previous)]
            folder = destination / f'{percent}pct' / name
            tick = time.perf_counter()
            result = fit(snp, folder, ls, cs, gs, free, max_nfev=600,
                         capacitance_bounds=bounds, capacitance_seeds=seeds,
                         title=f'T3 Transmission: Older Capacitances Within {percent}% ({name})')
            changes = {f'{k}{i}': 100*((result['C_fF'] if k=='C' else result['G_fF'])[i]/(cs if k=='C' else gs)[i]-1) for k,i in free[2:]}
            c_error = 100*(result['C_fF'][0]/251-1)
            g_error = 100*(result['G_fF'][0]/261-1)
            new_inside = not any(result['active_bounds'][i] for i in (0, 1, len(free)))
            row = dict(bound_percent=percent, scope=name, complex_RMS=result['complex_RMS'],
                       T3_C0_fF=result['C_fF'][0], T3_input_Cg_fF=result['G_fF'][0],
                       T3_C0_target_error_percent=c_error, T3_input_target_error_percent=g_error,
                       older_changes_percent=changes, active_bounds=result['active_bounds'],
                       fit_pass=bool(result['optimizer_success'] and result['complex_RMS']<=.025),
                       components_pass=bool(abs(c_error)<=5 and abs(g_error)<=5),
                       new_parameters_interior=new_inside, elapsed_seconds=time.perf_counter()-tick,
                       fit_path=str(folder/'fit.json'), plot_path=str(folder/'s21_fit.png'))
            row['combined_pass'] = row['fit_pass'] and row['components_pass'] and new_inside
            rows.append(row)
            previous = result
            summary = dict(case=CASE, source_response=str(snp), source_sha256=digest(snp),
                           baseline_complex_RMS=strict['complex_RMS'], all_inductances_fixed_nH=ls,
                           reference_C_fF=cs, reference_G_fF=gs,
                           policy='Bounds always measured from selected T2 reference. G1 (T3-T2 junction) belongs to T2 and is bounded. Output-end C4,G4 remain fixed. T3 C0/input Cg and total fixture delay free. No EM simulations.',
                           elapsed_seconds=time.perf_counter()-started, results=rows)
            (destination/'summary.json').write_text(json.dumps(summary, indent=2)+'\n')
            print('STAGED_RESULT', json.dumps(row), flush=True)
    lines = ['T3 staged bounded-capacitance fits', summary['policy'], '',
             'Bounds | Scope | Complex RMS | T3 C0 (fF) | T3 input Cg (fF) | Fit pass | Components pass']
    for row in rows:
        lines.append(f"{row['bound_percent']}% | {row['scope']} | {row['complex_RMS']:.6f} | {row['T3_C0_fF']:.3f} | {row['T3_input_Cg_fF']:.3f} | {row['fit_pass']} | {row['components_pass']}")
    (destination/'SUMMARY.txt').write_text('\n'.join(lines)+'\n')
    print('\n'.join(lines), flush=True)

if __name__ == '__main__':
    main()
