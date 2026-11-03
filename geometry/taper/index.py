"""Refresh the compact, machine-readable taper trial index from saved records."""
from pathlib import Path
import json
ROOT=Path(__file__).resolve().parents[2]

def refresh():
    trials=[]
    for path in sorted((ROOT/'output_data/taper_exploration').glob('*/trial.json')):
        r=json.loads(path.read_text())
        f=r['fits'].get('complex_with_delay') or r['fits']['complex']
        if 'L_t2_nH' in f:
            primary=r.get('network_fit',f)
            trials.append({'name':r['name'],'stage':'taper_2',
                           'primary_fit':r.get('primary_fit',primary['objective']),
                           'L0_nH':primary['L_t2_nH'],'C0_fF':primary['C0_t2_fF'],
                           'Cg_input_fF':primary['Cg_t2_fF'],'Cg_junction_fF':primary['Cg_t1_fF'],
                           'Cg_to_regular_fF':primary['Cg_fF'],
                           's21_complex_rmse':primary['s21_complex_rmse'],
                           'errors_percent':r['errors_percent'],'passes_2_percent':r['passes_2_percent'],
                           'record':str(path.relative_to(ROOT))})
            continue
        trials.append({'name':r['name'],'primary_fit':r['primary_fit'],'L0_nH':f['L_taper_nH'],'C0_fF':f['C_taper_fF'],
                       'Cg_input_fF':f['Cg_input_fF'],'Cg_junction_fF':f['Cg_junction_fF'],
                       'C_regular_fF':f['C_regular_fF'],'Cg_output_fF':f['Cg_output_fF'],
                       's21_complex_rmse':f['s21_complex_rmse'],'errors_percent':r['errors_percent'],
                       'passes_2_percent':r['passes_2_percent'],'record':str(path.relative_to(ROOT))})
    p=ROOT/'output_data/taper_exploration/summary.json';p.write_text(json.dumps({'trials':trials},indent=2)+'\n')
    print(json.dumps(trials,indent=2))
    return trials

if __name__=='__main__':refresh()
