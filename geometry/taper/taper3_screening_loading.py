"""Check reduced EM loading against saved full EM; never a final native result.

A cheap through-line calibration estimates the added output feed. Dense native
S-parameters are de-embedded and cascaded with the saved downstream LC model.
No sparse-frequency admittance fit, and no new L/C fitted to the full benchmark.
"""
import json
from types import SimpleNamespace
import gdstk
import numpy as np
import skrf
from scipy.optimize import minimize_scalar
from geometry.taper import taper3 as work
from geometry.taper.progressive_fitting import fit
from circuit_common import read_snp


def line(angle, impedance):
    result=np.zeros((len(angle),2,2),complex)
    result[:,0,0]=result[:,1,1]=np.cos(angle)
    result[:,0,1]=1j*impedance*np.sin(angle)
    result[:,1,0]=1j*np.sin(angle)/impedance
    return result


def balanced(f,length_um):
    # The balanced 12um remnants lie inside the 412um side ground regions.
    # Values are the native port-mode results, essentially constant over4-6GHz.
    return line(2*np.pi*f*(length_um*1e-6)*np.sqrt(6.5538)/299792458,56.8874)


def through_geometry(width):
    top=gdstk.Cell('TOP')
    parts=[gdstk.rectangle((0,0),(width,448)),
        gdstk.rectangle((0,-1052),(width,-618)),
        gdstk.rectangle((0,-1052),(412,-24)),
        gdstk.rectangle((width-412,-1052),(width,-24)),
        gdstk.rectangle((0,-18),(width,-8))]
    top.add(*gdstk.boolean(parts,[],'or',precision=.001,layer=1))
    meta=dict(layout_kind='uniform_feed_fixture_calibration',w_sim=width,expected_metal_components=3,
        parameters=dict(n_cells=0,l_sim=1500,l_top_ground=448,w_side=412,t_feed=10,
            unit=dict(d_top=8,cap=dict(t_top=14))))
    return top,meta


def validate(path,meta):
    expected,rebuild=through_geometry(meta['w_sim'])
    top=gdstk.read_gds(str(path)).top_level()
    if meta!=rebuild or len(top)!=1 or top[0].name!='TOP':
        raise ValueError('Through-line identity differs')
    polys=top[0].get_polygons()
    if gdstk.boolean(polys,expected.get_polygons(),'xor',precision=.001):
        raise ValueError('Through-line geometry differs')
    return SimpleNamespace(n_cells=0,l_sim=1500,l_top_ground=448),polys


def calibrate(width):
    native=work.NATIVE/f'fixture_through_w{width:g}'
    out=work.OUT/f'fixture_through_w{width:g}'
    native.mkdir(parents=True,exist_ok=True);out.mkdir(parents=True,exist_ok=True)
    top,meta=through_geometry(width);gds=native/'geometry.gds'
    if not gds.exists():
        lib=gdstk.Library(unit=1e-6,precision=1e-9);lib.add(top);lib.write_gds(str(gds))
        gds.with_suffix('.json').write_text(json.dumps(meta,indent=2)+'\n')
    project=native/'generated_ncell_0.sonx'
    if not project.exists():
        work.import_native(0,gds_path=gds,output_dir=native,geometry_validator=validate)
    work.check_geometry(project,gds,meta,geometry_validator=validate)
    status=out/'run_status.json'
    if status.exists() and json.loads(status.read_text()).get('status')=='complete':
        previous=json.loads(status.read_text())
        if work.digest(project)!=previous['project_sha256'] or work.digest(project.with_suffix('.snp'))!=previous['snp_sha256']:
            raise ValueError('Completed through-line changed')
    else:
        work.configure(project,4,6,10000);work.simulate(project,out)
    net=read_snp(project.with_suffix('.snp'));a=net.a
    if not np.allclose(skrf.Network(frequency=net.frequency,a=a,z0=50).s,net.s,atol=1e-10):
        raise ValueError('S/ABCD roundtrip failed')
    end=balanced(net.f,12)
    middle=np.linalg.inv(end)@a@np.linalg.inv(end)
    angle=np.arccos((middle[:,0,0]+middle[:,1,1])/2)
    impedance=np.sqrt(middle[:,0,1]/middle[:,1,0])
    impedance=np.where(impedance.real<0,-impedance,impedance)
    reconstructed=line(angle,impedance)
    error=float(np.max(abs(reconstructed-middle))/np.max(abs(middle)))
    np.savez(out/'line_model.npz',f=net.f,angle_per_um=angle/(width-824),impedance=impedance)
    report=dict(project=str(project),project_sha256=work.digest(project),snp_sha256=work.digest(project.with_suffix('.snp')),
        width_um=width,balanced_remnant_um_per_end=12,asymmetric_line_um=width-824,
        balanced_mode_ohm=56.8874,balanced_mode_Eeff=6.5538,
        impedance_real_range_ohm=[float(impedance.real.min()),float(impedance.real.max())],
        ABCD_reconstruction_relative_max=error,
        scope='Dense native through-line fixture estimate; loading predictions must be benchmarked against full EM.')
    (out/'calibration.json').write_text(json.dumps(report,indent=2)+'\n')
    print('THROUGH CALIBRATION',json.dumps(report),flush=True)
    return out


def complex_interp(f,source_f,values):
    return np.interp(f,source_f,values.real)+1j*np.interp(f,source_f,values.imag)


def reconnect(case,branches,calibration):
    stage='screen_T3_T2'+('_T1' if branches==3 else '')
    native=work.NATIVE/case/stage;source=work.OUT/case/stage
    out=source/'loaded_prediction';out.mkdir(exist_ok=True)
    meta=json.loads((native/'geometry.json').read_text())
    net=read_snp(native/f'generated_ncell_{branches}.snp')
    model=np.load(calibration/'line_model.npz');f=net.f
    theta=complex_interp(f,model['f'],model['angle_per_um'])
    z=complex_interp(f,model['f'],model['impedance'])
    length=meta['w_sim']-412-meta['clip_boundary_um']
    output_feed=line(theta*length,z)@balanced(f,12)
    upstream=net.a@np.linalg.inv(output_feed)
    ls=[json.loads((work.OUT/case/'inductance.json').read_text())['L0_nH']]+work.REFERENCE_FIT['L_nH']
    cs=[251]+work.REFERENCE_FIT['C_fF'];gs=[261]+work.REFERENCE_FIT['G_fF']
    downstream=np.tile(np.eye(2,dtype=complex),(len(f),1,1));w=2*np.pi*f
    for index in range(branches,len(ls)):
        s=np.tile(np.eye(2,dtype=complex),(len(f),1,1));s[:,0,1]=1/(1j*w*gs[index]*1e-15)
        y=np.tile(np.eye(2,dtype=complex),(len(f),1,1));y[:,1,0]=1j*w*cs[index]*1e-15+1/(1j*w*ls[index]*1e-9)
        downstream=downstream@s@y
    last_right=2*work.SOURCE_META['cell_centers'][-1]-work.SOURCE_META['cell_lefts'][-1]+meta['downstream_shift_um']
    final_length=meta['w_sim']-412-last_right
    predicted=skrf.Network(frequency=net.frequency,a=upstream@downstream@line(theta*final_length,z)@balanced(f,12),z0=50)
    predicted.write_touchstone(filename=str(out/'prediction'),form='ri')
    saved=out/'prediction.snp';saved.write_bytes((out/'prediction.s2p').read_bytes())
    narrow=predicted[(f>=4.65e9)&(f<=5.4e9)]
    narrow.write_touchstone(filename=str(out/'same_band'),form='ri')
    band=out/'same_band.snp';band.write_bytes((out/'same_band.s2p').read_bytes())
    fitted=fit(band,out/'fit',ls,cs,gs,[('C',0),('G',0),('G',1)],
        capacitance_seeds=[[100,261,153],[251,261,153],[400,120,150]],
        title='Reduced EM + Saved Downstream Circuit: Prediction Only')
    actual=read_snp(work.NATIVE/case/'response/generated_ncell_4.snp')
    values=complex_interp(actual.f,f,predicted.s[:,1,0]);observed=actual.s[:,1,0]
    objective=lambda tau:float(np.mean(abs(values*np.exp(-2j*np.pi*actual.f*tau*1e-12)-observed)**2))
    aligned=minimize_scalar(objective,bounds=(-30,30),method='bounded')
    reference=json.loads((work.OUT/case/'response/fit/fit.json').read_text())
    expected=np.array([reference['C_fF'][0],reference['G_fF'][0],reference['G_fF'][1]])
    extracted=np.array([fitted['C_fF'][0],fitted['G_fF'][0],fitted['G_fF'][1]])
    differences=100*(extracted/expected-1)
    report=dict(case=case,branches=branches,prediction_only=True,final_acceptance_permitted=False,
        output_asymmetric_line_removed_um=length,through_calibration=str(calibration),
        raw_full_EM_complex_RMS=float(np.sqrt(objective(0))),
        phase_aligned_full_EM_complex_RMS=float(np.sqrt(aligned.fun)),alignment_delay_ps=float(aligned.x),
        full_EM_magnitude_RMS=float(np.sqrt(np.mean((abs(values)-abs(observed))**2))),
        predicted_fit_values_fF=extracted.tolist(),actual_full_fit_values_fF=expected.tolist(),
        parameter_transfer_difference_percent=differences.tolist(),
        provisional_transfer_pass=bool(np.sqrt(aligned.fun)<=.025 and max(abs(differences))<=5),
        limitations='Downstream circuit and feed deembedding are approximations; no omitted metal mutual fields are restored.')
    (out/'comparison.json').write_text(json.dumps(report,indent=2)+'\n')
    fig,ax=work.plt.subplots(figsize=(12,7))
    ax.plot(actual.f/1e9,abs(observed),lw=2.5,label='Actual full Sonnet EM')
    ax.plot(actual.f/1e9,abs(values),'--',lw=2,label=f'{branches}-branch EM + downstream circuit prediction')
    ax.set(xlabel='Frequency (GHz)',ylabel=r'$|S_{21}|$',title='Reduced-Model Transfer Check: Prediction vs Full EM')
    ax.legend();ax.grid(alpha=.2);fig.tight_layout();fig.savefig(out/'full_EM_comparison.png',dpi=180);work.plt.close(fig)
    print('LOADING CHECK',json.dumps(report),flush=True)
    return report
