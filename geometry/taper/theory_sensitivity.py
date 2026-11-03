"""Cheap circuit-only sensitivity guide; does not alter EM or the target."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
import json
import numpy as np
from geometry.specifications import TAPER_2, TAPER_1, UNIT_CELL
from geometry.taper.progressive_fitting import transmission
from geometry.taper.response_diagnostics import peak_features
from geometry.taper.plain_gap_study import OUT


def main():
    f=np.linspace(4.65e9,5.4e9,7501)
    ls=[TAPER_2.L0_nH,TAPER_1.L0_nH,UNIT_CELL.L0_nH,UNIT_CELL.L0_nH]
    cs=[TAPER_2.C0_fF,TAPER_1.C0_fF,UNIT_CELL.C0_fF,UNIT_CELL.C0_fF]
    gs=[TAPER_2.Cg_input_fF,TAPER_1.Cg_input_fF,UNIT_CELL.Cg_input_fF,UNIT_CELL.Cg_input_fF]
    base=peak_features(f,transmission(f,ls,cs,gs,0));rows=[]
    for name,kind,index in [('T2 C0','C',0),('T1 C0','C',1),('input IDC Cg','G',0),('T2-T1 Cg','G',1)]:
        cp=cs.copy();gp=gs.copy();(cp if kind=='C' else gp)[index]*=1.1
        peaks=peak_features(f,transmission(f,ls,cp,gp,0))
        row=dict(parameter=name,increase_percent=10,peaks=peaks)
        if len(peaks)==len(base)==2:
            row['changes_vs_target']=[dict(frequency_MHz=(new['frequency_GHz']-old['frequency_GHz'])*1000,
                linewidth_percent=(new['power_FWHM_MHz']/old['power_FWHM_MHz']-1)*100)
                for new,old in zip(peaks,base)]
        rows.append(row)
    result=dict(scope='Circuit-only diagnostic. Target acceptance and every EM geometry remain unchanged. '
        'These are electrical-parameter sensitivities, not geometry derivatives.',baseline_peaks=base,perturbations=rows)
    OUT.mkdir(parents=True,exist_ok=True)
    (OUT/'theory_sensitivity.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
