"""Basic checks that reject the observed corrupt dense Sonnet cache export.

This is a saved-response check only. It launches no solver and never stops one.
"""
from pathlib import Path
import json,re
import numpy as np
from circuit_common import read_snp

def check(project,log,report_path):
    project=Path(project);log=Path(log);report_path=Path(report_path)
    net=read_snp(project.with_suffix('.snp'));text=log.read_text(encoding='utf-8-sig',errors='replace')
    issues=[]
    if not np.all(np.isfinite(net.s)) or not np.all(np.isfinite(net.f)):issues.append('Nonfinite exported values')
    if not np.all(np.diff(net.f)>0):issues.append('Frequency samples not strictly increasing')
    if 'Em simulation completed' not in text or 'Em errors detected: 0' not in text:issues.append('Native successful completion not confirmed')
    rows=[]
    number=r'[-+\d.eE]+'
    for line in text.splitlines():
        if re.fullmatch(r'\s*'+number+(r'\s+'+number)*8+r'\s*',line):
            values=list(map(float,line.split()));f=values[0]*1e9
            if net.f[0]-1<=f<=net.f[-1]+1:
                k=int(np.argmin(abs(net.f-f)))
                item=dict(f_GHz=values[0],nearest_export_distance_Hz=float(abs(net.f[k]-f)))
                if item['nearest_export_distance_Hz']<1:
                    direct=values[3]*np.exp(1j*np.deg2rad(values[4]))
                    item['S21_complex_export_difference']=float(abs(net.s[k,1,0]-direct))
                    if item['S21_complex_export_difference']>.0005:issues.append(f'Direct/export disagreement at {values[0]} GHz')
                rows.append(item)
    norm=np.max(abs(net.s),axis=(1,2));bad=np.flatnonzero((norm[1:-1]>.05)&(norm[:-2]<1e-6)&(norm[2:]<1e-6))+1
    if len(bad):issues.append('Isolated finite samples surrounded by near-zero complete S matrices: observed invalid interpolation pattern')
    snr=re.findall(r'ABS: SNR = (\S+) dB \((\S+) dB Goal\)',text)
    if snr:
        last,goal=map(float,snr[-1])
        if np.isfinite(last) and np.isfinite(goal) and last<goal:issues.append('Last recorded ABS accuracy check did not meet its goal')
    report=dict(basic_integrity_pass=not issues,issues=issues,source_project=str(project),native_log=str(log),exported_samples=len(net.f),frequency_range_GHz=[float(net.f[0]/1e9),float(net.f[-1]/1e9)],completed_direct_rows=rows,last_ABS_SNR_and_goal_dB=snr[-1] if snr else None,isolated_bad_sample_frequencies_GHz=(net.f[bad]/1e9).tolist(),scope='Basic native completion/direct-data/known-interpolation-artifact checks, not independent mesh convergence or parameter uniqueness.')
    report_path.parent.mkdir(parents=True,exist_ok=True);report_path.write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    if issues:raise ValueError('; '.join(issues))
    return report
