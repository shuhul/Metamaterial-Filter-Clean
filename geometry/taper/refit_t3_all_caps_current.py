"""Refit saved final T3 response; all Cs bounded around immutable current fit."""
from pathlib import Path
import json,time
from geometry.taper.progressive_fitting import fit

def main():
    selected=json.loads(Path('geometry/taper/selected_taper3.json').read_text(encoding='utf-8-sig'))
    anchor_path=Path(selected['result']['fit_path'])
    anchor=json.loads(anchor_path.read_text(encoding='utf-8-sig'))
    snp=Path(selected['result']['project']).with_suffix('.snp')
    out=anchor_path.parent.parent/'all_caps_current_anchor_1_2_3pct'
    out.mkdir(exist_ok=True)
    (out/'anchor_fit.json').write_text(json.dumps(anchor,indent=2)+'\n',encoding='utf-8')
    free=[(k,i) for i in range(len(anchor['C_fF'])) for k in ('C','G') if anchor['C_fF' if k=='C' else 'G_fF'][i] is not None]
    values=[anchor['C_fF' if k=='C' else 'G_fF'][i] for k,i in free]
    report=dict(anchor_fit=str(anchor_path),source_snp=str(snp),baseline_complex_RMS=anchor['complex_RMS'],policy='Every C0 and Cg free, including regular and output-end and T3-T2 junction. Bounds centered on same immutable current final T3 fit for all three tests, not cumulative. All inductances fixed; complex S21 and phase delay fitted. Final selection unchanged.',tests=[])
    previous=None
    for pct in (1,2,3):
        started=time.monotonic()
        seeds=[values]
        if previous is not None:seeds.append([previous['C_fF' if k=='C' else 'G_fF'][i] for k,i in free])
        result=fit(snp,out/f'{pct}pct',anchor['L_nH'],anchor['C_fF'],anchor['G_fF'],free,capacitance_bounds=[(v*(1-pct/100),v*(1+pct/100)) for v in values],capacitance_seeds=seeds,max_nfev=700,title=f'T3: All Capacitances Within {pct}% of Current Fit')
        report['tests'].append(dict(bound_percent=pct,complex_RMS=result['complex_RMS'],improvement_percent=100*(1-result['complex_RMS']/anchor['complex_RMS']),optimizer_success=result['optimizer_success'],active_bounds=result['active_bounds'],C_fF=result['C_fF'],G_fF=result['G_fF'],delay_ps=result['delay_ps'],elapsed_seconds=time.monotonic()-started,fit_path=str(out/f'{pct}pct'/'fit.json')))
        (out/'comparison.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
        print('COMPARISON',json.dumps(report['tests'][-1]),flush=True)
        previous=result
    print('REPORT',out/'comparison.json',flush=True)
if __name__=='__main__':main()
