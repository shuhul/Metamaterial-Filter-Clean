"""Two isolated four-thread Sonnet trials, with serial retry on failure.

Call run([(case, settings, reason), ...], label) from Python. Only the private
worker entry point uses an argument. Accepted results are finalized by the
parent after both workers finish, avoiding competing selection-file writes.
"""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
from pathlib import Path
from dataclasses import asdict
from datetime import datetime,timezone
import ctypes,json,shutil,subprocess,sys,time
from geometry.taper import taper2_plain as w


def available_memory():
    class Memory(ctypes.Structure):
        _fields_=[('length',ctypes.c_ulong),('load',ctypes.c_ulong)]+[(n,ctypes.c_ulonglong) for n in
            ('total','available','page_total','page_available','virtual_total','virtual_available','extended')]
    m=Memory();m.length=ctypes.sizeof(m)
    if os.name!='nt' or not ctypes.windll.kernel32.GlobalMemoryStatusEx(ctypes.byref(m)):return None
    return m.available


def worker(path):
    p=Path(path);job=json.loads(p.read_text());s=w.from_dict(job['settings']);name=job['case']
    measured=w.measure(name,s)
    if job.get('inductance_only'):
        result=dict(case=name,stage='inductance_only',L_nH=measured['L0_nH'],
            L_error_percent=measured['target_error_percent'],
            inductance_pass=bool(abs(measured['target_error_percent'])<=2))
        p.with_suffix('.result.json').write_text(json.dumps(result,indent=2)+'\n')
        print(json.dumps(result,indent=2),flush=True)
        return
    if abs(measured['target_error_percent'])>2:raise ValueError('Independent L failed; skip response')
    score,fitted=w.response(name,s,measured['L0_nH'])
    result=dict(case=name,L_nH=measured['L0_nH'],L_error_percent=measured['target_error_percent'],
        target_RMS=score['complex_RMS'],target_pass=score['pass'],fit_RMS=fitted['complex_RMS'],
        fit_pass=fitted['curve_pass'])
    p.with_suffix('.result.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2),flush=True)


def run(cases,label,*,inductance_only_cases=()):
    if len(cases)!=2 or len({q[0] for q in cases})!=2:raise ValueError('Exactly two independent cases required')
    if not set(inductance_only_cases).issubset({q[0] for q in cases}):
        raise ValueError('Unknown inductance-only case')
    folder=w.ROOT/'investigations/parallel_t2_2026-10-28'/label
    folder.mkdir(parents=True,exist_ok=True)
    if (folder/'report.json').exists():raise ValueError('Choose a new batch label')
    memory=available_memory();limit=2 if memory is not None and memory>6*1024**3 else 1
    report=dict(started_utc=datetime.now(timezone.utc).isoformat(),available_RAM_bytes=memory,
        max_jobs=limit,threads_per_solver=4,maximum_solver_threads=4*limit,jobs=[],
        fallback_to_serial=False,acceptance='Fixed nominal complex S21 RMS<=0.025; no changed targets')
    groups={};paths=[]
    for index,(name,s,reason) in enumerate(cases):
        _,out,_=w.prepare(name,s);_,_,cal=w.prepare(name,s,isolation=True)
        groups.setdefault(json.dumps(cal['unit'],sort_keys=True),[]).append((name,s))
        job=dict(case=name,settings=asdict(s),reason=reason,
                 inductance_only=name in inductance_only_cases)
        path=folder/f'job{index+1}.json';path.write_text(json.dumps(job,indent=2)+'\n');paths.append(path)
        (out/'proposal.json').write_text(json.dumps(job,indent=2)+'\n')
    # A shared central L cache must be complete before concurrent readers use it.
    for group in groups.values():
        if len(group)>1:w.measure(*group[0])
    def start(path,slots,suffix=''):
        log=path.with_suffix('.worker'+suffix+'.log').open('w')
        env=os.environ.copy();env['SONNET_MAX_JOBS']=str(slots)
        began=time.monotonic()
        proc=subprocess.Popen([sys.executable,'-B','-m','geometry.taper.parallel_trials',str(path)],
            cwd=w.ROOT,stdout=log,stderr=subprocess.STDOUT,env=env)
        return proc,log,began
    progress=w.ROOT/'geometry/taper/taper2_progress.json'
    state=json.loads(progress.read_text())
    state.update(new_sims_running=True,active_cases=[q[0] for q in cases],
        active_geometry_cases={name:asdict(s) for name,s,_ in cases},
        active_project_cases={name:str(w.ROOT/'sonnet_files/taper_exploration'/w.CAMPAIGN/name/
                                     'response/generated_ncell_3.sonx') for name,_,_ in cases},
        active_parallel_batch=str(folder),max_concurrent_solver_jobs=limit)
    progress.write_text(json.dumps(state,indent=2)+'\n')
    begin=time.monotonic();running=[]
    for path in paths:
        proc,log,began=start(path,limit);running.append((path,proc,log,began))
        if limit==1:proc.wait()
    for path,proc,log,began in running:
        code=proc.wait();log.close()
        row=dict(manifest=str(path),pid=proc.pid,returncode=code,elapsed_to_collection_seconds=time.monotonic()-began)
        if code:
            report['fallback_to_serial']=True
            out=w.ROOT/'output_data/taper_exploration'/w.CAMPAIGN/json.loads(path.read_text())['case']
            shutil.copytree(out,folder/(path.stem+'_failed_attempt'),dirs_exist_ok=True)
        report['jobs'].append(row)
    # Retry only after both first attempts have stopped; never overlap a serial retry.
    for row,path in zip(report['jobs'],paths):
        if row['returncode']:
            proc,log,began=start(path,1,'.serial_retry');row['serial_retry_returncode']=proc.wait()
            row['serial_retry_seconds']=time.monotonic()-began;log.close()
        result=path.with_suffix('.result.json')
        if result.exists():row['result']=json.loads(result.read_text())
    report.update(total_worker_wall_seconds=time.monotonic()-begin,ended_utc=datetime.now(timezone.utc).isoformat())
    (folder/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    passing=[r['result'] for r in report['jobs'] if r.get('result',{}).get('target_pass') and r['result']['fit_pass']]
    if passing:
        best=min(passing,key=lambda r:r['target_RMS']);w.finalize(best['case'])
    from geometry.taper.taper2_cleared_tuning import ledger
    ledger()
    state=json.loads(progress.read_text());state.update(new_sims_running=False,active_cases=[],
        active_geometry_cases={},active_project_cases={})
    progress.write_text(json.dumps(state,indent=2)+'\n')
    print(json.dumps(report,indent=2),flush=True)
    return report


if __name__=='__main__':
    if len(sys.argv)!=2:raise SystemExit('Worker invoked by run(cases,label); edit trial settings in the caller.')
    worker(sys.argv[1])
