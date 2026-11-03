"""Full 10-regular-cell filter with mirrored T3/T2/T1 taper sections.



Build from verified native 2 um metal, so mirroring never changes grid tie

rounding. No EM runs on import. prepare() is separate from simulate_full().

"""

import json

from pathlib import Path

from types import SimpleNamespace

import xml.etree.ElementTree as ET

import numpy as np

import gdstk

from geometry.taper import taper3 as t3

from geometry.prepare_sonnet import main as import_native, points, tag

from geometry.taper.audit import check_geometry

from geometry.taper.sonnet import configure, simulate

from geometry.specifications import TAPER_3,TAPER_2,TAPER_1,UNIT_CELL



N_REGULAR=10

CASE='full_filter_10_regular_3_tapers_each_side'

NATIVE=t3.ROOT/'sonnet_files/taper_exploration/full_filter'/CASE

OUT=t3.ROOT/'output_data/taper_exploration/full_filter'/CASE

PROJECT_NAME='full_filter_10_regular_3_tapers_each_side.sonx'





def source_default():

    selected=json.loads((t3.ROOT/'geometry/taper/selected_taper3.json').read_text())

    if not selected.get('T3_accepted'):raise ValueError('Finalize T3 before assembling the full filter')

    return Path(selected['project'])





def build(source_project=None,n_regular=N_REGULAR):
    t3.require_reference()

    if n_regular<1:raise ValueError('Need at least one regular cell')

    source_project=Path(source_project) if source_project else source_default()

    source_meta=json.loads((source_project.parent/'geometry.json').read_text())

    raw=[gdstk.Polygon(np.column_stack((points(p)[:,0],448-points(p)[:,1])),layer=1) for p in ET.parse(source_project).iter(tag('PlanarPolygon'))]

    if any(not np.allclose(p.points/2,np.round(p.points/2)) for p in raw):raise ValueError('Source native metal is not on 2 um grid')

    shift=source_meta['downstream_shift_um']

    first_left=t3.SOURCE_META['cell_lefts'][3]+shift

    unit_width=t3.SOURCE_META['regular']['cap']['w']

    gap=t3.SOURCE_META['gaps'][3]

    pitch=unit_width+gap

    array_width=n_regular*unit_width+(n_regular-1)*gap

    width=2*first_left+array_width

    half=gdstk.boolean(raw,[gdstk.rectangle((0,-1052),(first_left,448))],'and',precision=.001,layer=1)

    mirrored=[gdstk.Polygon(np.column_stack((width-p.points[:,0],p.points[:,1])),layer=1) for p in half]

    source_top=gdstk.read_gds(str(t3.SOURCE_GDS)).top_level()[0]

    reg_ref=next(r for r in source_top.references if np.allclose(r.origin,(t3.SOURCE_META['cell_centers'][3],0)))

    unit=[]

    for p in reg_ref.get_polygons():

        xy=p.points+np.array([shift,0]);xy=np.floor(xy/2+.5)*2;xy[:,0]-=first_left

        unit.append(gdstk.Polygon(xy,layer=1))

    if not np.allclose(gdstk.Cell('UNIT_CHECK').add(*[q.copy() for q in unit]).bounding_box(),((0,-618),(unit_width,-8))):

        raise ValueError('Regular unit extent changed')

    parts=half+mirrored

    parts.extend([gdstk.rectangle((first_left,0),(width-first_left,448),layer=1),gdstk.rectangle((first_left,-1052),(width-first_left,-618),layer=1)])

    centers=[]

    for i in range(n_regular):

        x=first_left+i*pitch;centers.append(x+unit_width/2)

        parts.extend(p.copy().translate(x,0) for p in unit)

        if i<n_regular-1:

            mid=x+unit_width+gap/2

            parts.append(gdstk.rectangle((mid-2,-618),(mid+2,0),layer=1))

    top=gdstk.Cell('TOP');top.add(*parts)

    meta=dict(layout_kind='full_filter_mirrored_tapers',source_project=str(source_project),source_project_sha256=t3.digest(source_project),source_snp_sha256=t3.digest(source_project.with_suffix('.snp')),

              n_regular=n_regular,total_resonators=n_regular+6,resonator_centers_um=source_meta['cell_centers'][1:4]+centers+[width-x for x in source_meta['cell_centers'][3:0:-1]],stage_order=['T3','T2','T1']+['regular']*n_regular+['T1_mirrored','T2_mirrored','T3_mirrored'],

              w_sim=width,first_regular_left_um=first_left,regular_centers_um=centers,regular_width_um=unit_width,regular_gap_um=gap,regular_ground_width_um=4,

              expected_metal_components=3,mirror_policy='Reflect already snapped native taper metal in x only. All regular cells retain the verified regular orientation.',

              parameters=dict(n_cells=n_regular+6,l_sim=1500,l_top_ground=448,w_side=412,t_feed=10,unit=dict(d_top=8,cap=dict(t_top=14))))

    # Verify the retained taper/fixture regions exactly, and every repeat's body.

    checks=dict(total_resonators=n_regular+6,regular_count=n_regular,right_taper_mirror_xor_empty=True,regular_unit_copies_equal=True,grid_um=2,regular_ground_traces_centered=True)

    actual=top.get_polygons()

    if len(gdstk.boolean(actual,[],'or',precision=.001))!=3:raise ValueError('Expected ground/resonators plus two isolated feeds')

    left=gdstk.boolean(actual,[gdstk.rectangle((0,-1052),(first_left,448))],'and',precision=.001)

    if gdstk.boolean(left,half,'xor',precision=.001):raise ValueError('Left taper clipping changed metal')

    right=gdstk.boolean(actual,[gdstk.rectangle((width-first_left,-1052),(width,448))],'and',precision=.001)

    if gdstk.boolean(right,mirrored,'xor',precision=.001):raise ValueError('Right taper is not a true mirror')

    for i in range(n_regular):

        x=first_left+i*pitch

        roi=gdstk.rectangle((x,-618),(x+unit_width,-8))

        cell=gdstk.boolean(actual,[roi],'and',precision=.001)

        expected=[p.copy().translate(x,0) for p in unit]

        if gdstk.boolean(cell,expected,'xor',precision=.001):raise ValueError('Regular cell body changed')

    if any(not np.allclose(p.points/2,np.round(p.points/2)) for p in actual):raise ValueError('Full geometry off grid')

    return top,meta,checks





def validate(path,meta):

    top,expected,checks=build(meta['source_project'],meta['n_regular'])

    if meta!=expected:raise ValueError('Full filter metadata changed')

    lib=gdstk.read_gds(str(path));tops=lib.top_level()

    if len(tops)!=1 or tops[0].name!='TOP' or not np.isclose(lib.unit,1e-6):raise ValueError('GDS units/top changed')

    actual=tops[0].get_polygons()

    if gdstk.boolean(actual,top.get_polygons(),'xor',precision=.001):raise ValueError('Full filter GDS differs from builder')

    return SimpleNamespace(n_cells=meta['total_resonators'],l_sim=1500,l_top_ground=448),actual





def plot_layout(top,meta):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.patches import Polygon
    fig,ax=plt.subplots(figsize=(23,4.7))
    for p in top.get_polygons():ax.add_patch(Polygon(p.points,fc='#297b54',ec='none'))
    labels=['T3','T2','T1']+[f'R{i+1}' for i in range(meta['n_regular'])]+['T1','T2','T3']
    for x,label in zip(meta['resonator_centers_um'],labels):ax.text(x,35,label,ha='center',va='center',fontsize=12,color='white',weight='bold')
    ax.set(xlim=(0,meta['w_sim']),ylim=(-650,75),aspect='equal',xlabel='x (um)',ylabel='y (um)',title='Full Filter Geometry: 10 Regular Cells and Mirrored Taper Sections')
    fig.tight_layout()
    for dest in (OUT,NATIVE):fig.savefig(dest/'full_filter_layout.png',dpi=220)
    plt.close(fig)


def prepare():

    top,meta,checks=build()

    NATIVE.mkdir(parents=True,exist_ok=True);OUT.mkdir(parents=True,exist_ok=True)

    gds=NATIVE/'geometry.gds'

    if gds.exists():

        if json.loads(gds.with_suffix('.json').read_text())!=meta:raise ValueError('Use a fresh full filter CASE for changed geometry')

    else:

        lib=gdstk.Library(unit=1e-6,precision=1e-9);lib.add(top);lib.write_gds(str(gds));gds.with_suffix('.json').write_text(json.dumps(meta,indent=2)+'\n')

    validate(gds,meta)

    project=NATIVE/PROJECT_NAME

    if not project.exists():

        import_native(meta['total_resonators'],gds_path=gds,output_dir=NATIVE,geometry_validator=validate)

        project.write_bytes((NATIVE/f"generated_ncell_{meta['total_resonators']}.sonx").read_bytes())

    check_geometry(project,gds,meta,geometry_validator=validate)

    if list(ET.parse(project).iter(tag('Component'))):raise ValueError('Full filter cannot contain ideal components')

    selected=t3.ROOT/'geometry/taper/selected_taper3.json'

    (OUT/'selected_taper3_source.json').write_bytes(selected.read_bytes())

    (OUT/'geometry_checks.json').write_text(json.dumps(checks,indent=2)+'\n')

    plot_layout(top,meta)
    print('FULL_FILTER_PREPARED',project,meta['w_sim'],flush=True)

    return project,meta





def ladder(f,ls,cs,series_gaps,delay_ps=0):

    if len(series_gaps)!=len(ls)+1:raise ValueError('Need series coupling at both ends and all internal junctions')

    w=2*np.pi*np.asarray(f);a=np.ones(len(w),complex);b=np.zeros(len(w),complex);c=b.copy();d=a.copy()

    for i in range(len(ls)+1):

        z=1/(1j*w*series_gaps[i]*1e-15);b=b+a*z;d=d+c*z

        if i<len(ls):

            y=1j*w*cs[i]*1e-15+1/(1j*w*ls[i]*1e-9);a=a+b*y;c=c+d*y

    den=a+b/50+50*c+d;s=np.zeros((len(w),2,2),complex)

    s[:,1,0]=s[:,0,1]=2/den

    s[:,0,0]=(a+b/50-50*c-d)/den;s[:,1,1]=(-a+b/50-50*c+d)/den

    return s*np.exp(-1j*w*delay_ps*1e-12)[:,None,None]





def circuit_values(selection,n=N_REGULAR,nominal=False):

    if nominal:

        values=[TAPER_3,TAPER_2,TAPER_1,UNIT_CELL]

        ls=[q.L0_nH for q in values];cs=[q.C0_fF for q in values];gs=[q.Cg_input_fF for q in values]+[UNIT_CELL.Cg_input_fF]

    else:

        fit_path=Path(selection['analysis_fit_path']) if selection.get('analysis_fit_path') else Path(selection['project']).parent.parent/'evidence/response'/selection['fit_folder']/'fit.json'

        fit=json.loads(fit_path.read_text(encoding='utf-8-sig'))

        ls=fit['L_nH'][:4];cs=fit['C_fF'][:4];gs=fit['G_fF']

    return ls[:3]+[ls[3]]*n+ls[2::-1],cs[:3]+[cs[3]]*n+cs[2::-1],gs[:4]+[gs[4]]*(n-1)+[gs[3],gs[2],gs[1],gs[0]]





def analyze(project=None):

    from circuit_common import read_snp

    from scipy.optimize import minimize_scalar

    import matplotlib

    matplotlib.use('Agg')

    import matplotlib.pyplot as plt

    project=Path(project) if project else NATIVE/PROJECT_NAME

    from geometry.taper.response_integrity import check as check_response_integrity

    check_response_integrity(project,OUT/'em.log',OUT/'response_integrity.json')

    net=read_snp(project.with_suffix('.snp'));obs=net.s[:,1,0]

    selection=json.loads((OUT/'selected_taper3_source.json').read_text())

    models={};report={}

    for label,nominal in [('Nominal targets',True),('Extracted cell model',False)]:

        ls,cs,gs=circuit_values(selection,nominal=nominal);zero=ladder(net.f,ls,cs,gs)

        objective=lambda d:np.mean(abs(zero[:,1,0]*np.exp(-1j*2*np.pi*net.f*d*1e-12)-obs)**2)
        delays=np.linspace(-200,200,161);seed=delays[int(np.argmin([objective(d) for d in delays]))]
        optimum=minimize_scalar(objective,bounds=(max(-200,seed-2.5),min(200,seed+2.5)),method='bounded')

        model=ladder(net.f,ls,cs,gs,optimum.x);models[label]=model

        report[label]=dict(L_nH=ls,C_fF=cs,series_Cg_fF=gs,total_delay_ps=float(optimum.x),complex_S21_RMS=float(np.sqrt(np.mean(abs(model[:,1,0]-obs)**2))),magnitude_S21_RMS=float(np.sqrt(np.mean((abs(model[:,1,0])-abs(obs))**2))))

    OUT.mkdir(parents=True,exist_ok=True)

    fig,axes=plt.subplots(2,1,figsize=(15,11),sharex=True)

    for ax,db in zip(axes,[False,True]):

        conv=lambda x:20*np.log10(np.maximum(abs(x),1e-8)) if db else abs(x)

        ax.plot(net.f/1e9,conv(obs),lw=2.8,label='Full Sonnet EM',color='#2363a0')

        for label,model in models.items():ax.plot(net.f/1e9,conv(model[:,1,0]),'--',lw=1.8,label=label)

        ax.set_ylabel('S21 (dB)' if db else r'$|S_{21}|$',fontsize=20);ax.tick_params(labelsize=14);ax.grid(alpha=.2);ax.legend(fontsize=12)

    axes[0].set_title('Full Filter Transmission: 10 Regular Cells and Three Tapers per Side',fontsize=21)

    axes[1].set_xlabel('Frequency (GHz)',fontsize=18);axes[1].set_ylim(-80,2)

    fig.tight_layout();fig.savefig(OUT/'full_filter_s21.png',dpi=190);fig.savefig(OUT/'full_filter_s21.pdf');plt.close(fig)

    fig,ax=plt.subplots(figsize=(15,8))
    ax.plot(net.f/1e9,abs(obs),lw=2.8,label='Full Sonnet EM',color='#2363a0')
    ax.plot(net.f/1e9,abs(models['Nominal targets'][:,1,0]),'--',lw=1.8,label='Nominal target circuit',color='#d65a31')
    ax.plot(net.f/1e9,abs(models['Extracted cell model'][:,1,0]),':',lw=2,label='Extracted cell circuit',color='#3e8052')
    ax.set_title('Full Filter Transmission',fontsize=25);ax.set_xlabel('Frequency (GHz)',fontsize=20);ax.set_ylabel(r'$|S_{21}|$',fontsize=25)
    ax.tick_params(labelsize=16);ax.grid(alpha=.2);ax.legend(fontsize=14,loc='upper right')
    ax.text(.02,.98,'10 regular cells; three tapers per side\n2 um mesh; 50 ohm ports\nNo full-filter capacitance refit',transform=ax.transAxes,va='top',fontsize=14,bbox=dict(boxstyle='round,pad=.5',fc='white',ec='gray'))
    ax.set_ylim(0,max(abs(obs).max(),abs(models['Nominal targets'][:,1,0]).max())*1.25)
    fig.tight_layout()
    for ext in ('png','pdf'):fig.savefig(OUT/f'full_filter_s21_magnitude.{ext}',dpi=190)
    plt.close(fig)
    report.update(project=str(project),source_sha256=t3.digest(project.with_suffix('.snp')),frequency_range_GHz=[float(net.f[0]/1e9),float(net.f[-1]/1e9)],frequency_samples=len(net.f),policy='No electrical refit of full filter. Nominal and measured assembly compared; one phase delay per model only. These RMS values describe agreement, not the staged extraction cutoff.')

    (OUT/'full_filter_analysis.json').write_text(json.dumps(report,indent=2)+'\n')

    np.savez(OUT/'full_filter_curves.npz',f=net.f,actual=net.s,nominal=models['Nominal targets'],extracted=models['Extracted cell model'])

    print('FULL_FILTER_ANALYSIS',json.dumps(report),flush=True)

    return report





def preflight(project):
    import re
    import ctypes
    from geometry.taper.sonnet import BIN,guarded_solver,solver_environment
    signature=OUT/'preflight_source.sha256'
    log_path=OUT/'preflight.log'
    resources=OUT/'preflight.log.resources.json'
    cached=(signature.exists() and signature.read_text().strip()==t3.digest(project)
            and log_path.exists() and resources.exists()
            and json.loads(resources.read_text()).get('returncode')==0)
    if not cached:
        signature.write_text(t3.digest(project)+'\n')
        with log_path.open('w') as log:
            run=guarded_solver([str(BIN/'em.exe'),'-v','-N',str(project)],cwd=project.parent,stdout=log,env=solver_environment())
        run.check_returncode()
    log=log_path.read_text()
    match=re.search(r'requires\s+([\d,]+)\s+subsections\s+and\s+([\d,.]+)\s+(MB|GB)',log,re.I)
    if not match:raise ValueError('Inspect preflight log: cannot identify memory estimate')
    class Memory(ctypes.Structure):
        _fields_=[('length',ctypes.c_ulong),('load',ctypes.c_ulong)]+[(k,ctypes.c_ulonglong) for k in ('total','available','total_page','available_page','total_virtual','available_virtual','extended')]
    memory=Memory();memory.length=ctypes.sizeof(memory)
    if not ctypes.windll.kernel32.GlobalMemoryStatusEx(ctypes.byref(memory)):raise OSError('Cannot query RAM')
    estimate=float(match.group(2).replace(',',''))*(1024 if match.group(3).upper()=='GB' else 1)
    available=memory.available/2**20
    report=dict(project_sha256=t3.digest(project),subsections=int(match.group(1).replace(',','')),Sonnet_estimate_MB=estimate,conservative_peak_estimate_MB=1.9*estimate,available_RAM_MB=available,reserve_MB=8192,reused_subsection_check=cached)
    report['memory_pass']=report['conservative_peak_estimate_MB']+8192<available
    (OUT/'preflight.json').write_text(json.dumps(report,indent=2)+'\n')
    print('FULL_FILTER_PREFLIGHT',json.dumps(report),flush=True)
    if not report['memory_pass']:raise RuntimeError('Full-filter estimated peak plus reserve exceeds available RAM; no full solve launched')
    return report



def simulate_full():

    project,meta=prepare()

    marker=OUT/'run_status.json'

    if marker.exists() and json.loads(marker.read_text()).get('status')=='complete':

        r=json.loads(marker.read_text());assert t3.digest(project)==r['project_sha256'] and t3.digest(project.with_suffix('.snp'))==r['snp_sha256']

    else:

        configure(project,4.65,5.4,7500);preflight(project);simulate(project,OUT)

    check_geometry(project,NATIVE/'geometry.gds',meta,geometry_validator=validate)

    return analyze(project)



if __name__=='__main__':

    # Explicit final-step entry point; requires accepted selected_taper3.json.

    simulate_full()

