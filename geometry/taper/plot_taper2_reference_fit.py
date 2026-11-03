"""Plot the original 2um T2 fit, nominal response and all parameter targets."""
import json
from hashlib import sha256
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from geometry.taper.taper2_grid4 import OLD_FIT
from geometry.taper.progressive_fitting import transmission
from geometry.specifications import TAPER_2,TAPER_1,UNIT_CELL
from circuit_common import read_snp
from pathlib import Path


def main():
    response=OLD_FIT.parent.parent
    out=response/'extraction_three_caps_adopted_reference'
    fitted=json.loads((out/'fit.json').read_text())
    score=json.loads((response/'target_comparison.json').read_text())
    source=Path(fitted['provenance']['source']);digest=sha256(source.read_bytes()).hexdigest()
    if digest!=fitted['provenance']['sha256'] or digest!=score['source_sha256']:
        raise ValueError('Fit, target comparison and saved data differ')
    net=read_snp(source);f=net.f;actual=net.s[:,1,0]
    ls=fitted['L_nH'];cs=fitted['C_fF'];gs=fitted['G_fF']
    model=transmission(f,ls,cs,gs,fitted['delay_ps'])
    targets=[TAPER_2,TAPER_1,UNIT_CELL,UNIT_CELL]
    target=transmission(f,[t.L0_nH for t in targets],[t.C0_fF for t in targets],
        [t.Cg_input_fF for t in targets],score['fixture_delay_ps'])
    rms=float(np.sqrt(np.mean(abs(model-actual)**2)))
    target_rms=float(np.sqrt(np.mean(abs(target-actual)**2)))
    if not np.isclose(rms,fitted['complex_RMS'],atol=1e-10):raise ValueError('Fit RMS mismatch')
    if not np.isclose(target_rms,score['complex_RMS'],atol=1e-10):raise ValueError('Target RMS mismatch')
    names=[('T2 L0','T2 C0','Input-to-T2 Cg'),('T1 L0','T1 C0','T2-to-T1 Cg'),
           ('Regular L0','Regular C0','T1-to-regular Cg'),
           ('Output-end L0','Output-end C0','Regular-to-output Cg')]
    rows=[];records=[]
    free={tuple(q) for q in fitted['free_capacitances']}
    for i,t in enumerate(targets):
        for label,kind,value,nominal,unit in zip(names[i],['L','C','G'],[ls[i],cs[i],gs[i]],
                [t.L0_nH,t.C0_fF,t.Cg_input_fF],['nH','fF','fF']):
            status='Fitted' if (kind,i) in free else 'Fixed'
            change=100*(value/nominal-1)
            rows.append([label,f'{value:.3f} {unit}',f'{nominal:g} {unit}',f'{change:+.2f}%',status])
            records.append(dict(parameter=label,value=value,target=nominal,unit=unit,
                                target_difference_percent=change,status=status))
    fig,(ax,table_ax)=plt.subplots(1,2,figsize=(20,9),gridspec_kw={'width_ratios':[1.22,1]})
    fig.suptitle('Second Taper Cell: Transmission and Parameter Comparison',fontsize=24,y=.98)
    ax.plot(f/1e9,abs(actual),lw=2.8,color='#2768aa',label='Original Sonnet EM (2 um mesh)')
    ax.plot(f/1e9,abs(model),'k--',lw=2.3,label='Three-capacitance complex fit')
    ax.plot(f/1e9,abs(target),':',lw=2.7,color='#db6d28',label='Specified target circuit')
    ax.text(.025,.97,f'Complex S21 fit RMS: {rms:.6f}\nSimulation-to-target RMS: {target_rms:.5f}\nTarget-response limit: 0.025',
        transform=ax.transAxes,va='top',fontsize=14,bbox=dict(boxstyle='round,pad=.5',fc='white',ec='gray'))
    ax.set(xlim=(4.65,5.4),ylim=(0,max(abs(actual).max(),abs(model).max(),abs(target).max())*1.5),
           xlabel='Frequency (GHz)',ylabel=r'$|S_{21}|$')
    ax.xaxis.label.set_size(20);ax.yaxis.label.set_size(23);ax.tick_params(labelsize=15)
    ax.grid(alpha=.2);ax.legend(fontsize=11,loc='upper right')
    table_ax.axis('off');table_ax.set_title('All Circuit Parameters',fontsize=21,pad=15)
    table=table_ax.table(cellText=rows,colLabels=['Parameter','Used value','Target','Difference','Status'],
        cellLoc='left',colLoc='left',colWidths=[.30,.22,.20,.16,.12],bbox=[0,.18,1,.80])
    table.auto_set_font_size(False);table.set_fontsize(11.5)
    for (row,col),cell in table.get_celld().items():
        cell.set_edgecolor('#d6dce1')
        if row==0:cell.set_facecolor('#e0e8ef');cell.set_text_props(weight='bold')
        else:
            group=(row-1)//3
            cell.set_facecolor('#f4f7fa' if group%2==0 else 'white')
            if rows[row-1][-1]=='Fitted':cell.set_text_props(weight='bold')
    table_ax.text(0,.12,'Differences are relative to the nominal electrical targets.\n'
        'Older capacitances are fixed at the adopted +/-1% fit values.\n'
        'All inductances are fixed; fixture delay is fitted.',
        transform=table_ax.transAxes,va='top',fontsize=11.5,linespacing=1.5)
    fig.tight_layout(rect=[0,0,1,.95],w_pad=2)
    for ext in ('png','pdf'):fig.savefig(out/f's21_fit_target_comparison.{ext}',dpi=180)
    plt.close(fig)
    (out/'parameter_target_comparison.json').write_text(json.dumps(dict(parameters=records,
        fit_RMS=rms,target_RMS=target_rms,source_sha256=digest),indent=2)+'\n')
    print(out/'s21_fit_target_comparison.png')


if __name__=='__main__':main()
