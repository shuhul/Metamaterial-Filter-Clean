"""Produce a compact, readable comparison from saved historical fits."""
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from geometry.taper.taper3_historical_bounds import DEST


def main():
    summary=json.loads((DEST/'results.json').read_text());rows=summary['results']
    best=min(rows,key=lambda r:r['RMS'])
    closest=min([r for r in rows if r['fit_pass']],key=lambda r:max(abs(r['C0_error_percent']),abs(r['Cg_error_percent'])))
    retained=[r for r in rows if r['junction_geometry']==dict(gap=22,fingers=70)]
    retained_best=min(retained,key=lambda r:r['RMS']) if retained else None
    closest_one=min([r for r in rows if r['fit_pass'] and r['bound_percent']==1],key=lambda r:max(abs(r['C0_error_percent']),abs(r['Cg_error_percent'])))
    picks=[('Lowest RMS',best),('Closest T3 capacitances among passing fits',closest),('Closest T3 targets with 1% bounds',closest_one)]
    if retained_best is not None:picks.append(('Lowest RMS with original T3–T2 IDC',retained_best))
    report=['HISTORICAL T3 REFIT RESULTS',f"{summary['total_cases']} geometries; {len(rows)} saved-response fits; no EM simulations.",
            summary['policy'],'All curves use the saved 4.65–5.40 GHz window. Fit error is complex S21 RMS to simulation, not error to target S21.',
            'Low residual alone does not establish parameter uniqueness. Older capacitance bound hits are permitted but recorded.',
            f"Combined fit/component passing fits: {sum(r['combined_pass'] for r in rows)}",'']
    for label,r in picks:
        a=np.load(str(__import__('pathlib').Path(r['fit_path']).parent/'response.npz'))
        fig,(ax,tableax)=plt.subplots(1,2,figsize=(16,7.5),gridspec_kw={'width_ratios':[2.4,1]})
        ax.plot(a['f']/1e9,abs(a['actual'][:,1,0]),lw=2.8,label='Sonnet',color='#2363a0')
        ax.plot(a['f']/1e9,abs(a['model'][:,1,0]),'--',lw=2.4,label='Complex circuit fit',color='#d65a31')
        ax.set_xlabel('Frequency (GHz)',fontsize=18);ax.set_ylabel(r'$|S_{21}|$',fontsize=22)
        ax.set_title('T3 Transmission',fontsize=23);ax.tick_params(labelsize=14);ax.grid(alpha=.2);ax.legend(fontsize=14)
        tableax.axis('off')
        text=(f"Complex S21 RMS: {r['RMS']:.5f}\n\nT3 L: {r['L_nH']:.4f} nH\nTarget: 2.77 ({r['L_error_percent']:+.2f}%)\n\nT3 C0: {r['T3_C0_fF']:.2f} fF\nTarget: 251 ({r['C0_error_percent']:+.2f}%)\n\nT3 input Cg: {r['T3_Cg_fF']:.2f} fF\nTarget: 261 ({r['Cg_error_percent']:+.2f}%)\n\nOlder bounds: ±{r['bound_percent']}%\nScope: {r['scope'].replace('_',' + ')}\nAll inductances fixed")
        tableax.text(.02,.95,text,va='top',fontsize=14,bbox=dict(boxstyle='round,pad=.6',facecolor='#f4f6f9',edgecolor='#aab3c0'))
        fig.suptitle(label,fontsize=17);fig.tight_layout()
        stem={'Closest T3 targets with 1% bounds':'closest_targets_1pct','Lowest RMS':'lowest_rms','Closest T3 capacitances among passing fits':'closest_components','Lowest RMS with original T3–T2 IDC':'original_junction_best'}[label]
        fig.savefig(DEST/(stem+'.png'),dpi=190);plt.close(fig)
        report.extend([label,json.dumps(r,indent=2),''])
    (DEST/'REPORT.txt').write_text('\n'.join(report)+'\n',encoding='utf-8')
    print(json.dumps(dict(picks=picks,combined_passes=sum(r['combined_pass'] for r in rows)),indent=2))

if __name__=='__main__':main()
