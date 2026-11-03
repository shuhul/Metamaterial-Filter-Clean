"""Recorded width search, using one independently calibrated EM case at a time."""
import json
from shutil import copyfile
from geometry.taper.minimum_inductance_meander import measure,response,prepare,NATIVE
from geometry.taper.minimum_inductance_review import OUT,review_width,PROFILES


def run_case(width):
    if (OUT/'PAUSE_REQUESTED').exists():raise RuntimeError('Search pause requested; no new solver launched')
    cal=measure(width)
    if cal['R_squared']<.9999:
        raise ValueError('Independent L calibration insufficiently linear')
    pair=cal['adjacent_pair_L0_nH']
    if (max(pair)-min(pair))/cal['L0_nH']>.02:
        raise ValueError('Calibration spread exceeds 2%; do not infer an L threshold')
    if (OUT/'PAUSE_REQUESTED').exists():raise RuntimeError('Search pause requested after L calibration')
    response(width,cal['L0_nH'],fit_legacy=False)
    r=review_width(width)
    report()
    return r


def report():
    rows=[json.loads(p.read_text()) for p in OUT.glob('meander_w*/review/summary.json')]
    rows.sort(key=lambda r:r['L_nH'])
    passes=[r for r in rows if r['primary_pass']]
    best=min(passes,key=lambda r:r['L_nH']) if passes else None
    lower=[r for r in rows if best and r['width_um']<best['width_um'] and not r['primary_pass']]
    failed=max(lower,key=lambda r:r['width_um']) if lower else None
    final=bool(best and (best['width_um']==16 or (failed and best['width_um']-failed['width_um']==4)))
    result=dict(primary_profile='full_T1_staged_calibration',tested_widths=[r['width_um'] for r in rows],
                lowest_verified_passing_L_nH=best['L_nH'] if best else None,
                lowest_verified_passing_width_um=best['width_um'] if best else None,
                adjacent_lower_failed_width_um=failed['width_um'] if failed else None,
                adjacent_lower_failed_L_nH=failed['L_nH'] if failed else None,
                grid_bracket_complete=final,final_T2_target_design=False,
                scope='Lowest verified passing width in the declared 58-run, 4 um trace family; not a universal physical L limit',
                cases=rows)
    if best:
        native=NATIVE/f"meander_w{best['width_um']}"/'direct/generated_ncell_3.sonx'
        plot=OUT/f"meander_w{best['width_um']}"/'review/full_T1_staged_calibration/s21_fit.png'
        result.update(project=str(native),plot=str(plot))
        copyfile(plot,OUT/'lowest_verified_s21.png')
    (OUT/'search_status.json').write_text(json.dumps(result,indent=2)+'\n')
    import matplotlib.pyplot as plt
    fig,ax=plt.subplots(figsize=(11,6))
    for profile in PROFILES:
        valid=[r for r in rows if profile in r['profiles']]
        ax.plot([r['L_nH'] for r in valid],[r['profiles'][profile]['complex_RMS'] for r in valid],
                'o-',lw=2,label=profile)
    ax.axhline(.01,color='gray',ls='--',label='RMS limit (other gates also apply)')
    ax.set_xlabel('Independently measured T2 inductance (nH)',fontsize=16)
    ax.set_ylabel('Complex S21 RMS error',fontsize=16);ax.set_title('Width-Only Meander Search',fontsize=22)
    ax.tick_params(labelsize=13);ax.legend(fontsize=10);ax.grid(alpha=.25)
    fig.tight_layout();fig.savefig(OUT/'fit_error_vs_inductance.png',dpi=190);plt.close(fig)
    print('SEARCH STATUS',json.dumps({k:v for k,v in result.items() if k!='cases'}),flush=True)
    return result


def next_width(result):
    if result['grid_bracket_complete']:return None
    hi=result['lowest_verified_passing_width_um'];lo=result['adjacent_lower_failed_width_um']
    if hi is None:raise ValueError('Establish a passing high-L reference before narrowing')
    if lo is None:return 16 if hi>16 else None
    cases={r['width_um']:r for r in result['cases']}
    low_fit=cases[lo]['profiles']['full_T1_staged_calibration']
    high_fit=cases[hi]['profiles']['full_T1_staged_calibration']
    low_error=low_fit['complex_RMS'];high_error=high_fit['complex_RMS']
    proposed=(hi+lo)/2
    if low_error>.01>high_error and low_error>high_error:
        proposed=lo+(hi-lo)*(low_error-.01)/(low_error-high_error)
        # Avoid a long sequence of tiny moves if the response is curved.
        proposed=max(lo+.2*(hi-lo),min(hi-.2*(hi-lo),proposed))
    width=4*round(proposed/4)
    width=max(lo+4,min(hi-4,width))
    if not lo<width<hi:raise ValueError('No distinct grid point inside bracket')
    return width


def finish_bracket(max_new_cases=6):
    """Bounded interpolated bracket refinement; stop on any exception."""
    for _ in range(max_new_cases):
        result=report();width=next_width(result)
        if width is None:return result
        if (OUT/'PAUSE_REQUESTED').exists():raise RuntimeError('Search pause requested')
        prepare(width,'direct')
        run_case(width)
    result=report()
    if not result['grid_bracket_complete']:raise RuntimeError('Bounded search budget reached; inspect results before any further solve')
    return result


if __name__=='__main__':
    import sys
    if len(sys.argv)>1 and sys.argv[1]=='auto':finish_bracket()
    elif len(sys.argv)>1:run_case(int(sys.argv[1]))
    else:report()
