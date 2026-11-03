"""Readable target/fit overlays with extracted-value comparisons."""
from pathlib import Path
import json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from geometry.specifications import TAPER_1, TAPER_2, UNIT_CELL


def plot(stage, heuristic='heuristic_S21_calibrated_L'):
    stage=Path(stage)
    metric=json.loads((stage/'target_comparison.json').read_text())
    fit=json.loads((stage/heuristic/'fit.json').read_text())
    with np.load(stage/'target_comparison.npz') as a:
        f=a['f']/1e9;actual=a['actual'];target=a['target']
    with np.load(stage/heuristic/'response.npz') as a:model=a['model'][:,1,0]
    rows=[('L T2',fit['L_nH'][0],TAPER_2.L0_nH,'nH'),
          ('C0 T2',fit['C_fF'][0],TAPER_2.C0_fF,'fF'),
          ('Cg T2 (input)',fit['G_fF'][0],TAPER_2.Cg_input_fF,'fF'),
          ('L T1',fit['L_nH'][1],TAPER_1.L0_nH,'nH'),
          ('C0 T1',fit['C_fF'][1],TAPER_1.C0_fF,'fF'),
          ('Cg T1 (T2 to T1)',fit['G_fF'][1],TAPER_1.Cg_input_fF,'fF'),
          ('Cg (T1 to unit)',fit['G_fF'][2],UNIT_CELL.Cg_input_fF,'fF')]
    for kind in ('target','fit'):
        with plt.rc_context({'font.size':15}):
            fig,(ax,panel)=plt.subplots(1,2,figsize=(19,8),gridspec_kw={'width_ratios':[1.35,1]})
            ax.plot(f,abs(actual),color='#2768aa',lw=3,label='Sonnet EM')
            ax.plot(f,abs(target if kind=='target' else model),'--',color='#db6d28',lw=2.5,
                    label='Specified target circuit' if kind=='target' else 'Fitted equivalent circuit')
            ax.set(xlabel='Frequency (GHz)',ylabel=r'$|S_{21}|$',xlim=(f.min(),f.max()),ylim=(0,1.06))
            ax.xaxis.label.set_size(20);ax.yaxis.label.set_size(23);ax.tick_params(labelsize=15)
            ax.grid(alpha=.2);ax.spines[['top','right']].set_visible(False)
            ax.legend(fontsize=14,loc='upper right');panel.axis('off')
            panel.text(0,.96,'Extracted Values and Targets',fontsize=21,weight='bold')
            panel.text(0,.90,'Heuristic value (target; signed difference)',fontsize=14,color='#555555')
            for i,(label,value,goal,unit) in enumerate(rows):
                y=.81-i*.092
                panel.text(0,y,label,fontsize=16,weight='bold')
                panel.text(0,y-.04,f'{value:.3f} {unit}  (target {goal:g}; {(value/goal-1)*100:+.2f}%)',fontsize=15)
            panel.text(0,.105,f"Target complex RMS = {metric['complex_RMS']:.4f} (limit 0.05)",fontsize=16,weight='bold')
            panel.text(0,.035,f"Fitted complex RMS = {fit['complex_RMS']:.4f} (diagnostic)\n"
                       'L T2: current calibration; L T1: prior reference.\n'
                       'Capacitances guide geometry; they do not decide success.'+
                       ('\nSome fitted parameters reach bounds.' if any(fit['active_bounds']) else ''),
                       fontsize=12,color='#555555',va='top',linespacing=1.5)
            fig.suptitle('Taper 2 + Taper 1: '+('Target Transmission Comparison' if kind=='target' else 'Transmission Fit'),fontsize=26,weight='bold')
            fig.tight_layout(rect=(0,.03,1,.94),w_pad=3)
            name='s21_target_with_parameters' if kind=='target' else 's21_fit_linear'
            for suffix in ('png','pdf'):fig.savefig(stage/f'{name}.{suffix}',dpi=200,facecolor='white',bbox_inches='tight')
            plt.close(fig)
    return stage/'s21_target_with_parameters.png'
