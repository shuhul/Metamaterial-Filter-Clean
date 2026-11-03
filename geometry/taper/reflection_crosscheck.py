"""Use already calculated complex reflections to test S21-ambiguous capacitors.

Only one shared output reference delay is adjusted. Capacitor candidates come
from the S21 profiles; S11 remains a separate prediction check.
"""
import json
import numpy as np
from geometry.taper.clearance_extraction import OUT
from geometry.taper.corrected_iteration import REGULAR_C,REGULAR_G
from geometry.taper.staged_fitting import circuit
from geometry.taper.assess_corrected import best_target_delay


def run(case='frame700_wire88_left214'):
    root=OUT/case;folder=root/'staged_joint'
    r=json.loads((folder/'results.json').read_text());v=r['fitted']
    values=[v[k] for k in ('C0_T2_fF','Cg_input_T2_fF','C0_T1_fF','Cg_T2_T1_fF','Cg_T1_regular_fF','delay_bare_ps','delay_IDC_ps')]
    candidates=[('best_S21',values)]+[(f"C0_fixed_{q['fixed_value_fF']:.3f}",q['all_values']) for q in r['profiles'][0]['rows']]
    data=[]
    for stage in ('direct','idc'):
        with np.load(root/stage/'fit/response.npz') as a:data.append({k:a[k] for k in ('f','actual')})
    rows=[]
    for name,x in candidates:
        base=[]
        for j,a in enumerate(data):
            base.append(circuit(a['f'],r['fixed_inductances_nH'][j],[x[0] if j else 0,x[2],REGULAR_C,REGULAR_C],
                         [x[1] if j else None,x[3],x[4],REGULAR_G]))
        # S22 sees only twice the common output delay when port 1 is matched.
        freq=np.concatenate([2*a['f'] for a in data])
        right,_=best_target_delay(freq,np.concatenate([s[:,1,1] for s in base]),
                                     np.concatenate([a['actual'][:,1,1] for a in data]))
        s22=[];s11=[]
        for j,a in enumerate(data):
            w=2*np.pi*a['f']
            s22.append(float(np.sqrt(np.mean(abs(base[j][:,1,1]*np.exp(-2j*w*right*1e-12)-a['actual'][:,1,1])**2))))
            s11.append(float(np.sqrt(np.mean(abs(base[j][:,0,0]*np.exp(-2j*w*(x[5+j]-right)*1e-12)-a['actual'][:,0,0])**2))))
        rows.append(dict(candidate=name,capacitors_fF=x[:5],output_delay_ps=right,
                         S22_complex_RMS=s22,S11_prediction_complex_RMS=s11))
    result=dict(candidates=rows,scope='S21 capacitor candidates tested against complex S22; only output phase fitted. S11 not fitted. No new EM.')
    (folder/'reflection_crosscheck.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2),flush=True)
    return result


if __name__=='__main__':run()
