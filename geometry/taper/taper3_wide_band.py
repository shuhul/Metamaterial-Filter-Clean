"""Same-metal wider-band T3 diagnostic. Preserves original response and fit."""
import json,shutil
import numpy as np
from geometry.taper import taper3 as work
from geometry.taper.sonnet import configure,simulate
from geometry.taper.progressive_fitting import fit
from geometry.taper.taper3_fit_sensitivity import analyze
from geometry.taper.finalize_taper3 import plot_case
from circuit_common import read_snp

CASE='w130_t6g2n66_dl342_dr340_pad80x590_mc24_i134ni40_j28n16'
STAGE='response_wide_4to5p4'


def run(case=CASE):
    source=work.NATIVE/case/'response'
    native=work.NATIVE/case/STAGE
    out=work.OUT/case/STAGE
    source_fit=json.loads((work.OUT/case/'response/fit/fit.json').read_text())
    manifest=dict(source_project=str(source/'generated_ncell_4.sonx'),
        source_project_sha256=work.digest(source/'generated_ncell_4.sonx'),
        source_snp_sha256=work.digest(source/'generated_ncell_4.snp'),
        range_GHz=[4.,5.4],target_intervals=14000,
        same_geometry_mesh_materials_ports=True,older_values_fixed=True,
        purpose='Test outer-stage constraint and window dependence. Diagnostic, not an accepted design.')
    out.mkdir(parents=True,exist_ok=True)
    if not native.exists():shutil.copytree(source,native)
    marker=out/'manifest.json'
    if marker.exists() and json.loads(marker.read_text())!=manifest:
        raise ValueError('Source or diagnostic configuration changed')
    marker.write_text(json.dumps(manifest,indent=2)+'\n')
    project=native/'generated_ncell_4.sonx'
    meta=json.loads((native/'geometry.json').read_text())
    if work.digest(native/'geometry.gds')!=work.digest(source/'geometry.gds'):
        raise ValueError('Metal geometry changed')
    status=out/'run_status.json'
    if status.exists() and json.loads(status.read_text()).get('status')=='complete':
        previous=json.loads(status.read_text())
        if work.digest(project)!=previous['project_sha256'] or work.digest(project.with_suffix('.snp'))!=previous['snp_sha256']:
            raise ValueError('Completed wider diagnostic changed')
    else:
        configure(project,4.,5.4,14000)
        work.check_geometry(project,native/'geometry.gds',meta,geometry_validator=work.validate)
        simulate(project,out)
    seeds=[[251,261,153],[321,117,145],[100,400,150],[400,120,150],[200,300,150]]
    free=[('C',0),('G',0),('G',1)]
    fitted=fit(project.with_suffix('.snp'),out/'fit',source_fit['L_nH'],
        [251]+work.REFERENCE_FIT['C_fF'],[261]+work.REFERENCE_FIT['G_fF'],free,
        capacitance_seeds=seeds,title='Third Taper: Wider-Band Diagnostic')
    values=[fitted['C_fF'][0],*fitted['G_fF'][:2]]
    errors=(100*(np.array(values)/[251,261,153]-1)).tolist()
    result=dict(case=case,stage=STAGE,values_fF=values,targets_fF=[251,261,153],
        errors_percent=errors,max_cap_error_percent=max(abs(v) for v in errors),
        fit_RMS=fitted['complex_RMS'],components_pass=max(abs(v) for v in errors)<=5,
        fit_pass=fitted['curve_pass'] and fitted['interior_free_parameters'],
        diagnostic_only=True,project=str(project),all_older_parameters_fixed=True,
        measured_L_nH=source_fit['L_nH'][0],T2_reference_sha256=work.REFERENCE_HASH)
    (out/'result.json').write_text(json.dumps(result,indent=2)+'\n')
    network=read_snp(project.with_suffix('.snp'))
    narrow=network[(network.f>=4.65e9-1)&(network.f<=5.4e9+1)]
    narrow.write_touchstone(str(out/'narrow_band'),form='ri')
    nf=fit(out/'narrow_band.s2p',out/'fit_narrow',source_fit['L_nH'],
        [251]+work.REFERENCE_FIT['C_fF'],[261]+work.REFERENCE_FIT['G_fF'],free,
        capacitance_seeds=seeds,title='Third Taper: Narrow Band from Wider EM Sweep')
    old=read_snp(source/'generated_ncell_4.snp')
    newer=network.interpolate(old.frequency)
    comparison=dict(wide_values_fF=values,wide_RMS=fitted['complex_RMS'],
        original_values_fF=[source_fit['C_fF'][0],*source_fit['G_fF'][:2]],original_RMS=source_fit['complex_RMS'],
        narrow_refit_values_fF=[nf['C_fF'][0],*nf['G_fF'][:2]],narrow_refit_RMS=nf['complex_RMS'],
        old_vs_extended_EM_complex_RMS_in_original_band=float(np.sqrt(np.mean(abs(newer.s[:,1,0]-old.s[:,1,0])**2))),
        diagnostic_only=True)
    (out/'window_comparison.json').write_text(json.dumps(comparison,indent=2)+'\n')
    plot_case(case,stage=STAGE)
    analyze(case,stage=STAGE)
    print(json.dumps(comparison,indent=2),flush=True)
    return result


if __name__=='__main__':run()
