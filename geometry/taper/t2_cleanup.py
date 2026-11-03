"""Local, reversible T2 cleanup trials on the accepted T3 assembly.

All original projects and selections remain immutable. Each variant has an
independent native project, measured response and target-bounded fit.
"""
import os
os.environ['SONNET_MAX_JOBS']='1'
os.environ['OPENBLAS_NUM_THREADS']='1'
from pathlib import Path
from types import SimpleNamespace
from datetime import datetime,timezone
import json,copy,shutil,xml.etree.ElementTree as ET
import numpy as np
import gdstk
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from geometry.taper import taper3 as t3
from geometry.prepare_sonnet import tag,points,main as import_native
from geometry.taper.audit import check_geometry
from geometry.taper.sonnet import configure,simulate
from geometry.taper.progressive_fitting import fit
from geometry.taper.explicit_inductance import measure
from geometry.specifications import TAPER_3,TAPER_2,TAPER_1,UNIT_CELL

ROOT=t3.ROOT
CAMPAIGN='T2_cleanup_2026-11-01'
NATIVE=ROOT/'sonnet_files/taper_exploration'/CAMPAIGN
OUT=ROOT/'output_data/taper_exploration'/CAMPAIGN
SELECTION=json.loads((ROOT/'geometry/taper/selected_taper3.json').read_text(encoding='utf-8-sig'))
SOURCE=Path(SELECTION['project'])
ANCHOR=json.loads(Path(SELECTION['result']['fit_path']).read_text(encoding='utf-8-sig'))
TARGETS=[TAPER_3,TAPER_2,TAPER_1,UNIT_CELL,UNIT_CELL]

def read(p):return json.loads(Path(p).read_text(encoding='utf-8-sig'))
def save(p,x):
    p=Path(p);p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(x,indent=2)+'\n',encoding='utf-8')
def raw_metal(project):
    return [gdstk.Polygon(np.column_stack((points(p)[:,0],448-points(p)[:,1])),layer=1) for p in ET.parse(project).iter(tag('PlanarPolygon'))]
def box(x0,y0,x1,y1):return gdstk.rectangle((x0,y0),(x1,y1),layer=1)

def edited(source, *, trim=False, center=False, calibration=None):
    meta=read(source.parent/'geometry.json');polys=raw_metal(source)
    windows=[];cuts=[];adds=[];changes={}
    if calibration=='T3':
        c3=meta['cell_centers'][1];c2=None
    elif calibration=='T2':
        c3=None;c2=meta['cell_centers'][1]
    else:c3,c2=meta['cell_centers'][1:3]
    # The 40-finger junction starts at y=-24 and ends at -340 um.
    lower=-24-(40-1)*8-4
    if trim:
        if c3 is not None:
            x=c3+828/2;cuts.append(box(x-8,-598,x,lower))
        if c2 is not None:
            x=c2-564/2;cuts.append(box(x,-598,x+8,lower))
        windows.extend(cuts)
        changes['IDC_rail_bottom_old_new_um']=[-598,lower]
        changes['IDC_fingers_unchanged']=True
    if center and c2 is not None:
        # Move the exact existing native coil, keeping length/n/trace width.
        bounds=[(c2-18,-618,c2+130,-22),(c2-126,-618,c2-102,-218),(c2+208,-618,c2+232,-218)]
        shifts=[-56,-52,-54]
        for b,dx in zip(bounds,shifts):
            roi=box(*b);shape=gdstk.boolean(polys,[roi],'and',precision=.001,layer=1)
            assert shape,('Missing T2 part',b)
            cuts.append(roi);adds.extend(p.copy().translate(dx,0) for p in shape)
            windows.extend([roi,box(b[0]+dx,b[1],b[2]+dx,b[3])])
        changes['T2_coil_shift_um']=-56
        changes['T2_pad_shifts_um']=[-52,-54]
        changes['T2_meander_center_um']=c2
        changes['T2_pad_centers_um']=[c2-166,c2+166]
        changes['pad_to_meander_gap_um']=80
        if calibration=='T2':
            meta['inductance_pin_x']=c2
            meta['unit']['meander']['x_offset_from_cap_center']=0
    result=gdstk.boolean(polys,cuts,'not',precision=.001,layer=1) if cuts else polys
    result=gdstk.boolean(result+adds,[],'or',precision=.001,layer=1)
    delta=gdstk.boolean(polys,result,'xor',precision=.001)
    assert not gdstk.boolean(delta,windows,'not',precision=.001),'Nonlocal geometry change'
    assert len(result)==meta['expected_metal_components'],'Connectivity changed'
    assert all(np.allclose(p.points/2,np.round(p.points/2)) for p in result)
    meta.update(layout_kind='T2_cleanup_'+('cal_'+calibration if calibration else 'response'),cleanup=dict(source=str(source),source_sha256=t3.digest(source),trim=trim,center=center,calibration=calibration,changes=changes,outside_edit_windows_XOR_empty=True))
    return result,meta

def validate(gds,meta):
    cfg=meta['cleanup'];source=Path(cfg['source']);assert t3.digest(source)==cfg['source_sha256']
    expected,check=edited(source,trim=cfg['trim'],center=cfg['center'],calibration=cfg['calibration']);assert check==meta
    lib=gdstk.read_gds(str(gds));top=lib.top_level();assert len(top)==1 and np.isclose(lib.unit,1e-6)
    actual=top[0].get_polygons();assert not gdstk.boolean(actual,expected,'xor',precision=.001)
    return SimpleNamespace(n_cells=meta['parameters']['n_cells'],l_sim=1500,l_top_ground=448),actual

def prepare(case, *, trim=False,center=False,calibration=None):
    source=SOURCE
    if calibration:
        s=SELECTION if calibration=='T3' else read(ROOT/'geometry/taper/selected_taper2_cap5.json')
        source=Path(s['project']).parent.parent/'Lcal/generated_ncell_1.sonx'
    stage='response' if calibration is None else 'Lcal_'+calibration
    native=NATIVE/case/stage;out=OUT/case/stage
    native.mkdir(parents=True,exist_ok=True);out.mkdir(parents=True,exist_ok=True)
    polys,meta=edited(source,trim=trim,center=center,calibration=calibration)
    gds=native/'geometry.gds'
    if gds.exists():assert read(gds.with_suffix('.json'))==meta
    else:
        lib=gdstk.Library(unit=1e-6,precision=1e-9);lib.new_cell('TOP').add(*polys);lib.write_gds(str(gds));save(gds.with_suffix('.json'),meta)
    n=meta['parameters']['n_cells'];project=native/f'generated_ncell_{n}.sonx'
    if not project.exists():import_native(n,gds_path=gds,output_dir=native,geometry_validator=validate)
    check_geometry(project,gds,meta,geometry_validator=validate)
    if calibration is None:
        fig,ax=plt.subplots(figsize=(19,7))
        for p in polys:ax.add_patch(Polygon(p.points,fc='#297b54',ec='none'))
        c3,c2=meta['cell_centers'][1:3]
        ax.set(xlim=(c3-440,c2+340),ylim=(-650,30),aspect='equal',xlabel='x (um)',ylabel='y (um)',title=f'T2 Cleanup: {case}')
        ax.axvline(c2,ls=':',color='#9e4747',alpha=.65,lw=1)
        fig.tight_layout();fig.savefig(native/'T2_closeup.png',dpi=190);plt.close(fig)
    save(out/'geometry_checks.json',dict(project=str(project),source=str(source),source_sha256=t3.digest(source),GDS_native_XOR_empty=True,mesh_um=2,material_stack_ports_control_unchanged=True,**meta['cleanup']['changes']))
    return project,out,meta

def fit_response(project,out,ls):
    cs=[t.C0_fF for t in TARGETS];gs=[t.Cg_input_fF for t in TARGETS]
    free=[(k,i) for i in range(5) for k in ('C','G')]
    vals=[(cs if k=='C' else gs)[i] for k,i in free]
    seed=[ANCHOR['C_fF' if k=='C' else 'G_fF'][i] for k,i in free]
    result=fit(project.with_suffix('.snp'),out/'fit_all_caps_target5',ls,cs,gs,free,capacitance_bounds=[(.95*v,1.05*v) for v in vals],capacitance_seeds=[seed],max_nfev=900,title='T3 Assembly: All Capacitances Within 5% of Targets')
    # Large readable plot with extracted values and target differences.
    a=np.load(out/'fit_all_caps_target5/response.npz')
    fig,(ax,tab)=plt.subplots(1,2,figsize=(21,9),gridspec_kw={'width_ratios':[1.2,1]})
    ax.plot(a['f']/1e9,abs(a['actual'][:,1,0]),lw=2.8,label='Sonnet');ax.plot(a['f']/1e9,abs(a['model'][:,1,0]),'--',lw=2.2,label='Complex circuit fit')
    ax.set_xlabel('Frequency (GHz)',fontsize=20);ax.set_ylabel(r'$|S_{21}|$',fontsize=24);ax.tick_params(labelsize=15);ax.legend(fontsize=14);ax.grid(alpha=.2)
    ax.text(.03,.97,f"Complex S21 RMS: {result['complex_RMS']:.5f}\nAll capacitances within 5% of targets\nAll inductances fixed; phase delay fitted",transform=ax.transAxes,va='top',fontsize=14,bbox=dict(fc='white',ec='gray',boxstyle='round,pad=.5'))
    ax.set_ylim(0,1.32*max(abs(a['actual'][:,1,0]).max(),abs(a['model'][:,1,0]).max()))
    rows=[]
    for i,(name,t) in enumerate(zip(['T3','T2','T1','Regular','Output-end'],TARGETS)):
        for label,v,target,u in [('L',ls[i],t.L0_nH,'nH'),('C0',result['C_fF'][i],t.C0_fF,'fF'),('input Cg',result['G_fF'][i],t.Cg_input_fF,'fF')]:
            rows.append([name+' '+label,f'{v:.3f} {u}',f'{target:g} {u}',f'{100*(v/target-1):+.2f}%'])
    tab.axis('off');table=tab.table(cellText=rows,colLabels=['Parameter','Fit / fixed L','Target','Difference'],bbox=[0,.12,1,.84],cellLoc='left');table.auto_set_font_size(False);table.set_fontsize(12)
    for (i,j),cell in table.get_celld().items():cell.set_facecolor('#e0e8ef' if i==0 else '#f4f7fa' if ((i-1)//3)%2==0 else 'white');cell.set_edgecolor('#d6dce1')
    tab.text(0,.06,'Bounded equivalent-circuit fit; values at limits are not uniquely identified.\nNo nominal-target S21 acceptance requirement.',fontsize=11)
    fig.suptitle('Taper Assembly: Transmission and Parameter Comparison',fontsize=24);fig.tight_layout(rect=[0,0,1,.95])
    for ext in ('png','pdf'):fig.savefig(out/f'final_s21_fit.{ext}',dpi=180)
    plt.close(fig)
    return result

def simulate_case(case,trim,center,ls):
    project,out,meta=prepare(case,trim=trim,center=center)
    status=out/'run_status.json'
    if status.exists() and read(status)['status']=='complete':
        r=read(status);assert t3.digest(project)==r['project_sha256'] and t3.digest(project.with_suffix('.snp'))==r['snp_sha256']
    else:configure(project,4.65,5.4,7500);simulate(project,out)
    fitted=fit_response(project,out,ls)
    return dict(case=case,trim=trim,center=center,project=str(project),fit_path=str(out/'fit_all_caps_target5/fit.json'),fit_RMS=fitted['complex_RMS'],fit_pass=fitted['optimizer_success'] and fitted['complex_RMS']<=.025,C_fF=fitted['C_fF'],G_fF=fitted['G_fF'],L_nH=ls,plot=str(out/'final_s21_fit.png'),geometry=str(project.parent/'T2_closeup.png'))

def calibrate(case,stage,trim,center):
    project,out,meta=prepare(case,trim=trim,center=center,calibration=stage)
    target=TAPER_3 if stage=='T3' else TAPER_2
    caps=[450,600,750,900] if stage=='T3' else [250,350,450,550]
    return measure(project,project.parent.parent/('Lcal_explicit_'+stage),out.parent/('Lcal_explicit_'+stage),name=case+'_'+stage,target_L_nH=target.L0_nH,caps_fF=caps,metadata=meta,samples=1601,range_GHz=(2.4,4.4))

def main():
    OUT.mkdir(parents=True,exist_ok=True)
    save(OUT/'workflow_status.json',dict(state='preparing',started_utc=datetime.now(timezone.utc).isoformat(),launcher_pid=os.getpid(),policy='Serial four-thread solver. Preserve baseline and every trial; do not stop a healthy full ABS run just because it adds frequencies.'))
    for name,trim,center in [('trimmed_IDC',True,False),('centered_T2',False,True),('trimmed_IDC_centered_T2',True,True)]:prepare(name,trim=trim,center=center)
    print('ALL_PREVIEWS_PREPARED',str(NATIVE),flush=True)

if __name__=='__main__':main()
