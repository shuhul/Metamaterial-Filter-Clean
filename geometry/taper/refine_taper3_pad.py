"""Saved-data refinement and target-based selection, with immutable old-C anchors."""
import json
import time
from geometry.taper.taper3 import OUT
from geometry.taper.taper3_pad_tuning import bounded_fit

CASE='w130_t6g2n68_dl342_dr340_pad184x594_i54_j28n40_pad_tune'
BOUNDS=[1,1.5,1.6,2,3]


def main():
    root=OUT/CASE
    # Wait for the existing serial EM/controller, not for any approval.
    started=time.monotonic()
    while not (root/'pad_tune_result.json').exists():
        status=root/'response/run_status.json'
        if status.exists() and json.loads(status.read_text()).get('status')=='failed':raise RuntimeError('Native response failed; inspect existing run before refining')
        if time.monotonic()-started>4*3600:raise TimeoutError('No completed response after four hours')
        time.sleep(3)
    l=json.loads((root/'inductance.json').read_text())['L0_nH']
    rows=[]
    for bound in BOUNDS:
        r=bounded_fit(CASE,l,bound,refinement=True)
        r['max_new_cap_error_percent']=max(abs(r[k]) for k in ('C0_error_percent','Cg_input_error_percent','junction_error_percent'))
        r['within_4_percent']=r['max_new_cap_error_percent']<=4
        rows.append(r)
    accepted=[r for r in rows if r['accepted']]
    qualified4=[r for r in accepted if r['within_4_percent']]
    winner=min(qualified4 or accepted,key=lambda r:(r['max_new_cap_error_percent'],r['older_bound_percent'],r['fit_RMS'])) if accepted else None
    report=dict(case=CASE,policy='All L fixed. Older Cs anchored to original selected T2. Qualify RMS<=.025 and new Cs within5%; prefer4% and rank by maximum new-capacitance target error, not by lowest residual. No new EM for refinement.',all_fits=rows,recommended=winner)
    (root/'refinement_selection.json').write_text(json.dumps(report,indent=2)+'\n')
    print('REFINEMENT_RECOMMENDATION',json.dumps(winner),flush=True)

if __name__=='__main__':main()
