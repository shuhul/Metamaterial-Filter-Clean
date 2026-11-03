"""Small symmetric-pad width tuning around the closest historical T3 design.

Parameter changes only: symmetric pad width, lead count and optional input IDC span.
Topology, top trace, flat ground and downstream assembly retained.
Explicit capacitor L calibration,
then one guarded full EM sweep; 1% old-C bounds, 2% only if RMS > .025.
"""
import json
import time
from dataclasses import replace
from geometry.taper.taper3 import OUT,NATIVE,ROOT,REFERENCE_FIT,from_dict,prepare,measure,digest,require_reference
from geometry.taper.sonnet import configure,simulate
from geometry.taper.progressive_fitting import fit

BASE_CASE='w130_t6g2n71_dl342_dr340_pad220x594_i56_j28n40'
PAD_WIDTHS=[184]
PAD_LEAD_COUNTS={184:68}
INPUT_IDC_SPANS={184:54}


def bounded_fit(case,l0,percent, *, refinement=False,junction_bound=None):
    root=OUT/case/'response'
    ls=[l0]+REFERENCE_FIT['L_nH'];cs=[251.]+REFERENCE_FIT['C_fF'];gs=[261.]+REFERENCE_FIT['G_fF']
    free=[('C',0),('G',0)]+[(k,i) for i in range(1,5) for k in ('C','G')]
    if junction_bound==0:free.remove(('G',1))
    bounds=[(0.,1000.),(.1,3000.)]+[((cs if k=='C' else gs)[i]*(1-percent/100),(cs if k=='C' else gs)[i]*(1+percent/100)) for k,i in free[2:]]
    if junction_bound is not None:
        if not 0<=junction_bound<=percent:raise ValueError('Junction prior must be within overall older bound')
        if junction_bound>0:bounds[free.index(('G',1))]=(gs[1]*(1-junction_bound/100),gs[1]*(1+junction_bound/100))
    fit_folder=f'fit_bounded_{percent}pct'+('' if junction_bound is None else '_jfixed' if junction_bound==0 else f'_j{junction_bound}pct')
    old=json.loads((OUT/BASE_CASE/'response/historical_staged_bounds/1pct/T2_T1_regular_output_end/fit.json').read_text())
    seed=[(old['C_fF'] if k=='C' else old['G_fF'])[i] for k,i in free]
    fitted=fit(NATIVE/case/'response/generated_ncell_4.snp',root/fit_folder,ls,cs,gs,free,capacitance_bounds=bounds,capacitance_seeds=[seed],max_nfev=600,title=f'Third Taper Transmission: {percent}% Older Capacitance Bounds')
    ce=100*(fitted['C_fF'][0]/251-1);ge=100*(fitted['G_fF'][0]/261-1);le=100*(l0/2.77-1)
    result=dict(case=case,refinement_authorized=refinement,older_bound_percent=percent,junction_bound_percent=junction_bound if junction_bound is not None else percent,all_inductances_fixed=True,fit_RMS=fitted['complex_RMS'],C0_fF=fitted['C_fF'][0],Cg_input_fF=fitted['G_fF'][0],C0_error_percent=ce,Cg_input_error_percent=ge,L_nH=l0,L_error_percent=le,
                junction_Cg_fF=fitted['G_fF'][1],junction_error_percent=100*(fitted['G_fF'][1]/153-1),
                fit_pass=bool(fitted['optimizer_success'] and fitted['complex_RMS']<=.025),components_pass=abs(ce)<=5 and abs(ge)<=5 and abs(le)<=2 and abs(fitted['G_fF'][1]/153-1)<=.05,
                new_interior=not any(fitted['active_bounds'][i] for i in (0,1,len(free))),project=str(NATIVE/case/'response/generated_ncell_4.sonx'),fit_path=str(root/fit_folder/'fit.json'))
    result['accepted']=result['fit_pass'] and result['components_pass'] and result['new_interior']
    (root/fit_folder/'result.json').write_text(json.dumps(result,indent=2)+'\n')
    print('PAD_RESULT',json.dumps(result),flush=True)
    plot_result(case,percent,result,fit_folder=fit_folder)
    return result


def plot_result(case,percent,result, *, fit_folder=None):
    import numpy as np
    import matplotlib.pyplot as plt
    from geometry.specifications import TAPER_3,TAPER_2,TAPER_1,UNIT_CELL
    folder=OUT/case/'response'/(fit_folder or f'fit_bounded_{percent}pct')
    r=json.loads((folder/'fit.json').read_text());a=np.load(folder/'response.npz')
    fig,(ax,tab)=plt.subplots(1,2,figsize=(21,9),gridspec_kw={'width_ratios':[1.25,1]})
    ax.plot(a['f']/1e9,abs(a['actual'][:,1,0]),lw=2.8,label='Sonnet',color='#2363a0');ax.plot(a['f']/1e9,abs(a['model'][:,1,0]),'--',lw=2.4,label='Complex circuit fit',color='#d65a31')
    ax.set_xlabel('Frequency (GHz)',fontsize=20);ax.set_ylabel(r'$|S_{21}|$',fontsize=24);ax.tick_params(labelsize=15);ax.grid(alpha=.2);ax.legend(fontsize=14)
    policy=f'All L fixed; older Cs bounded to {percent}%'
    if result.get('junction_bound_percent')==0:policy+='\nJunction Cg held at saved T2 reference'
    ax.text(.03,.97,f"Fit RMS: {r['complex_RMS']:.5f}\n{policy}",transform=ax.transAxes,va='top',fontsize=15,bbox=dict(boxstyle='round,pad=.5',fc='white',ec='gray'))
    ax.set_ylim(0,max(abs(a['actual'][:,1,0]).max(),abs(a['model'][:,1,0]).max())*1.30)
    rows=[]
    for i,(name,target) in enumerate(zip(['T3','T2','T1','Regular','Output-end'],[TAPER_3,TAPER_2,TAPER_1,UNIT_CELL,UNIT_CELL])):
        for label,v,t,u in [('L',r['L_nH'][i],target.L0_nH,'nH'),('C0',r['C_fF'][i],target.C0_fF,'fF'),('input Cg',r['G_fF'][i],target.Cg_input_fF,'fF')]:
            rows.append([name+' '+label,f'{v:.3f} {u}',f'{t:g} {u}',f'{100*(v/t-1):+.2f}%'])
    tab.axis('off');tab.set_title('Parameter Comparison',fontsize=22)
    table=tab.table(cellText=rows,colLabels=['Parameter','Extracted','Target','Difference'],cellLoc='left',colLoc='left',bbox=[0,.12,1,.84]);table.auto_set_font_size(False);table.set_fontsize(12)
    for (i,j),cell in table.get_celld().items():
        cell.set_edgecolor('#d6dce1');cell.set_facecolor('#e0e8ef' if i==0 else '#f4f7fa' if ((i-1)//3)%2==0 else 'white')
        if i<=3:cell.set_text_props(weight='bold')
    tab.text(0,.06,'Older bounds are relative to the saved T2 fit; no accumulated shifts.\nTarget S21 error is diagnostic only.',fontsize=11,va='top')
    fig.suptitle('Third Taper Cell: Transmission and Parameter Comparison',fontsize=24);fig.tight_layout(rect=[0,0,1,.95],w_pad=2)
    for ext in ('png','pdf'):fig.savefig(folder/f's21_fit_parameters.{ext}',dpi=180)
    plt.close(fig)


def run_width(width, *, n=71, d_input=56):
    require_reference()
    meta=json.loads((NATIVE/BASE_CASE/'response/geometry.json').read_text());baseline=from_dict(meta['t3_settings'])
    p=replace(baseline,d_input=d_input,body=replace(baseline.body,pad=replace(baseline.body.pad,w=width),meander=replace(baseline.body.meander,n=n)))
    case=f'w130_t6g2n{n}_dl342_dr340_pad{width}x594_i{d_input}_j28n40_pad_tune'
    root=OUT/case;root.mkdir(parents=True,exist_ok=True)
    (root/'iteration_plan.json').write_text(json.dumps(dict(baseline=BASE_CASE,change=dict(pad_width_old_um=220,pad_width_new_um=width,meander_n_old=71,meander_n_new=n,input_IDC_span_old_um=56,input_IDC_span_new_um=d_input),policy='Symmetric pad width tuning; optional n adjustment to keep L in tolerance. Cap-rail to pad gaps remain 116 um; top trace, meander width/thickness, junction IDC and downstream assembly unchanged. Input IDC span optionally adjusted, retaining finger count/thickness/gaps.'),indent=2)+'\n')
    print('PAD_TUNE_START',case,flush=True)
    measured=measure(case,p,caps_fF=[450,600,750,900]);l0=measured['L0_nH']
    print('PAD_L',l0,flush=True)
    if abs(l0/2.77-1)>.02:raise ValueError('L outside 2%: inspect before full response')
    project,out,meta=prepare(case,p,stage='response')
    status=out/'run_status.json'
    if status.exists() and json.loads(status.read_text()).get('status')=='complete':
        saved=json.loads(status.read_text());assert digest(project)==saved['project_sha256'] and digest(project.with_suffix('.snp'))==saved['snp_sha256']
    else:
        configure(project,4.65,5.4,7500);simulate(project,out)
    result=bounded_fit(case,l0,1)
    if not result['fit_pass']:result=bounded_fit(case,l0,2)
    (root/'pad_tune_result.json').write_text(json.dumps(result,indent=2)+'\n')
    return result


def main():
    for width in PAD_WIDTHS:
        if run_width(width,n=PAD_LEAD_COUNTS.get(width,71),d_input=INPUT_IDC_SPANS.get(width,56))['accepted']:break

if __name__=='__main__':main()
