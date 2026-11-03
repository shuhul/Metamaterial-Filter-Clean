"""Monitor one targeted field, then re-test full ABS from cache without new fields."""
from pathlib import Path
import json,time,os
import numpy as np
from geometry.taper.full_filter import OUT
from geometry.taper.full_filter_cache_probe import main as cache_probe
from circuit_common import read_snp

ROOT=OUT/'solver_efficiency_review/seed_4p8375'
def followup(root=ROOT):
    plan_path=Path('geometry/taper/full_filter_work_plan.json')
    plan=json.loads(plan_path.read_text(encoding='utf-8-sig'));plan.update(full_followup_pid=os.getpid(),full_followup_root=str(root));plan_path.write_text(json.dumps(plan,indent=2)+'\n',encoding='utf-8')
    while True:
        try:status=json.loads((root/'run_status.json').read_text(encoding='utf-8-sig'))
        except (OSError,json.JSONDecodeError):time.sleep(1);continue
        if status['status']=='complete':break
        if status['status']=='failed':raise RuntimeError(status.get('reason','Targeted native field failed'))
        try:
            resources=json.loads((root/'em.log.resources.json').read_text(encoding='utf-8-sig'));sample=resources['samples'][-1]
            print(json.dumps(dict(phase='one_targeted_field',preserved_fields=15,current_frequency_GHz=4.8375,elapsed_min=round(sample['elapsed_seconds']/60,1),active_cores=round(sample['child_core_equivalents'],2),tail=(root/'em.log').read_text(encoding='utf-8-sig')[-350:])),flush=True)
        except (OSError,json.JSONDecodeError,IndexError):print('Waiting for the next resource sample; native solver guard remains active.',flush=True)
        time.sleep(50)
    print('Targeted field completed; starting zero-new-field cache tests.',flush=True)
    report=cache_probe(cache_source=status['cache_directory'])
    plan=json.loads(plan_path.read_text(encoding='utf-8-sig'))
    samples={round(v['f_GHz']*1e9):v for v in plan.get('completed_direct_samples',[])}
    net=read_snp(Path(status['project']).with_suffix('.snp'))
    for k,f in enumerate(net.f):
        samples[round(float(f))]=dict(f_GHz=float(f/1e9),S11_mag=float(abs(net.s[k,0,0])),S21_magnitude=float(abs(net.s[k,1,0])),S21_phase_deg=float(np.angle(net.s[k,1,0],deg=True)),S21_dB=float(20*np.log10(max(abs(net.s[k,1,0]),1e-15))))
    plan.update(full_state='cache_retest_complete_after_seed',full_phase='targeted_seed_and_zero_field_retests_complete',full_completed_field_frequencies=16,completed_direct_samples=list(sorted(samples.values(),key=lambda v:v['f_GHz'])),full_latest_cache_retest=str(Path(report['tests'][0]['log']).parent.parent/'results.json'),full_seed_result=status,full_last_field_em_pid=plan.get('full_em_pid'),full_em_pid=None,full_current_frequency_GHz=None)
    plan_path.write_text(json.dumps(plan,indent=2)+'\n',encoding='utf-8')
    print('FULL_SEED_FOLLOWUP_COMPLETE',plan['full_latest_cache_retest'],flush=True)

if __name__=='__main__':followup()
