"""Compare historical T3 response fits under identical anchored bounds. No EM."""
import json
import time
from geometry.taper.taper3 import OUT, NATIVE, REFERENCE_FIT, require_reference, digest
from geometry.taper.progressive_fitting import fit

DEST = OUT / 'historical_staged_bounds'

def main():
    require_reference()
    DEST.mkdir(exist_ok=True)
    cases=[]
    for folder in OUT.iterdir():
        snp=NATIVE/folder.name/'response/generated_ncell_4.snp'
        meta=snp.parent/'geometry.json'
        cal=folder/'inductance.json'
        if not (snp.exists() and meta.exists() and cal.exists()): continue
        m=json.loads(meta.read_text());p=m['t3_settings']
        if p.get('pad_right') is not None: continue
        baseline=folder/'response/fit_fixed_T2junction/fit.json'
        if not baseline.exists():baseline=folder/'response/fit/fit.json'
        if not baseline.exists():continue
        cases.append((folder,snp,p,baseline))
    cases.sort(key=lambda q:(q[2]['body']['cap']['t_top'],q[0].name))
    rows=[]
    started=time.perf_counter()
    for folder,snp,p,baseline in cases:
        ls=[json.loads((folder/'inductance.json').read_text())['L0_nH']]+REFERENCE_FIT['L_nH']
        cs=[251.]+REFERENCE_FIT['C_fF'];gs=[261.]+REFERENCE_FIT['G_fF']
        initial=json.loads(baseline.read_text());previous=initial
        for percent in (1,2,3):
            for count,name in ((1,'T2'),(2,'T2_T1'),(3,'T2_T1_regular'),(4,'T2_T1_regular_output_end')):
                free=[('C',0),('G',0)]+[(k,i) for i in range(1,count+1) for k in ('C','G')]
                bounds=[(0.,1000.),(.1,3000.)]+[((cs if k=='C' else gs)[i]*(1-percent/100),(cs if k=='C' else gs)[i]*(1+percent/100)) for k,i in free[2:]]
                seeds=[[(r['C_fF'] if k=='C' else r['G_fF'])[i] for k,i in free] for r in (initial,previous)]
                out=folder/'response/historical_staged_bounds'/f'{percent}pct'/name
                r=fit(snp,out,ls,cs,gs,free,capacitance_bounds=bounds,capacitance_seeds=seeds,max_nfev=600,title=f'T3 Transmission: {percent}% Bounds ({name})')
                previous=r
                ce=100*(r['C_fF'][0]/251-1);ge=100*(r['G_fF'][0]/261-1)
                inside=not any(r['active_bounds'][i] for i in (0,1,len(free)))
                row=dict(case=folder.name,top_thickness_um=p['body']['cap']['t_top'],pad=p['body']['pad'],d_top_um=p.get('d_top',8),
                         junction_geometry=dict(gap=p['d_t2'],fingers=p['n_t2']),bound_percent=percent,scope=name,
                         RMS=r['complex_RMS'],T3_C0_fF=r['C_fF'][0],T3_Cg_fF=r['G_fF'][0],
                         C0_error_percent=ce,Cg_error_percent=ge,L_nH=ls[0],L_error_percent=100*(ls[0]/2.77-1),
                         fit_pass=bool(r['optimizer_success'] and r['complex_RMS']<=.025),components_pass=abs(ce)<=5 and abs(ge)<=5 and abs(ls[0]/2.77-1)<=.02,
                         active_bounds=r['active_bounds'],new_interior=inside,fit_path=str(out/'fit.json'),plot_path=str(out/'s21_fit.png'),
                         native_project=str(snp.with_suffix('.sonx')),response_sha256=digest(snp),reference_L_fixed=ls)
                row['combined_pass']=row['fit_pass'] and row['components_pass'] and inside
                rows.append(row)
                (DEST/'results.json').write_text(json.dumps(dict(policy='T3 C0/input Cg/delay free; old Cs bounded against immutable selected T2. All L fixed. G1 is bounded with T2. Historical junction geometry can differ; recorded explicitly. Symmetric pad settings only. No EM.',total_cases=len(cases),completed_fits=len(rows),elapsed_seconds=time.perf_counter()-started,results=rows),indent=2)+'\n')
        best=min([r for r in rows if r['case']==folder.name],key=lambda r:r['RMS'])
        print('CASE_DONE',json.dumps(best),flush=True)
    bests=[min([r for r in rows if r['case']==c[0].name],key=lambda r:r['RMS']) for c in cases]
    bests.sort(key=lambda r:r['RMS'])
    (DEST/'best_by_case.json').write_text(json.dumps(bests,indent=2)+'\n')
    lines=['Historical T3 bounded-fit comparison','All inductances fixed; bounds anchored to selected T2, not accumulated. Best RMS per geometry; small RMS alone does not establish uniqueness or target compliance.','', 'RMS | C0 fF | Input Cg fF | Bound % | Scope | Both target and fit pass | Case']
    for r in bests:lines.append(f"{r['RMS']:.6f} | {r['T3_C0_fF']:.3f} | {r['T3_Cg_fF']:.3f} | {r['bound_percent']} | {r['scope']} | {r['combined_pass']} | {r['case']}")
    (DEST/'SUMMARY.txt').write_text('\n'.join(lines)+'\n')
    print('FINISHED',len(cases),len(rows),time.perf_counter()-started,flush=True)

if __name__=='__main__':main()
