"""Summarize completed native full-filter EM; never launches a solver."""
from pathlib import Path
import json
import re
import numpy as np
from circuit_common import read_snp
from geometry.taper.full_filter import OUT, NATIVE, PROJECT_NAME


def main():
    status=json.loads((OUT/'run_status.json').read_text(encoding='utf-8-sig'))
    if status.get('status')!='complete':
        raise RuntimeError('Full EM is not complete; no final summary generated')
    net=read_snp(NATIVE/Path(PROJECT_NAME).with_suffix('.snp'))
    f=net.f/1e9
    mag=np.abs(net.s[:,1,0])
    db=20*np.log10(np.maximum(mag,1e-12))
    peak=int(np.argmax(mag))
    def main_band(threshold):
        if db[peak]<threshold:
            return None
        a=b=peak
        while a>0 and db[a-1]>=threshold:a-=1
        while b<len(db)-1 and db[b+1]>=threshold:b+=1
        lo=float(f[a]) if a==0 else float(np.interp(threshold,db[a-1:a+1],f[a-1:a+1]))
        hi=float(f[b]) if b==len(db)-1 else float(np.interp(threshold,db[b:b+2][::-1],f[b:b+2][::-1]))
        return dict(low_GHz=lo,high_GHz=hi,width_MHz=(hi-lo)*1000,
                    low_edge_truncated=a==0,high_edge_truncated=b==len(db)-1,
                    definition='Contiguous interval containing the largest transmission maximum')
    window=(f>=4.8)&(f<=5.2)
    power=np.abs(net.s[:,0,0])**2+mag**2
    log=(OUT/'em.log').read_text(encoding='utf-8-sig')
    direct=[float(v) for v in re.findall(r'Frequency:\s+([\d.]+) GHZ',log)]
    report=dict(project=str(NATIVE/PROJECT_NAME),status=status,
        dense_export_samples=len(f),direct_EM_frequency_count=len(direct),direct_EM_frequencies_GHz=direct,
        peak_GHz=float(f[peak]),peak_magnitude=float(mag[peak]),peak_dB=float(db[peak]),
        band_3dB_absolute=main_band(-3),band_1dB_absolute=main_band(-1),
        nominal_4p8_to_5p2_GHz=dict(min_dB=float(db[window].min()),max_dB=float(db[window].max()),
                                  peak_to_peak_dB=float(np.ptp(db[window]))),
        endpoints_dB=[float(db[0]),float(db[-1])],
        max_reciprocity_error=float(np.max(np.abs(net.s[:,0,1]-net.s[:,1,0]))),
        max_S11_S22_difference=float(np.max(np.abs(net.s[:,0,0]-net.s[:,1,1]))),
        input_power_balance_min=float(power.min()),input_power_balance_max=float(power.max()),
        nominal_and_extracted_comparison=json.loads((OUT/'full_filter_analysis.json').read_text(encoding='utf-8-sig')),
        interpretation='Full geometry EM validation; dense ABS interpolation is distinguished from directly solved frequencies. No full-filter electrical-parameter refit.')
    (OUT/'full_filter_summary.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(report,indent=2))
    return report

if __name__=='__main__':main()
