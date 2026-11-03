"""Run an explicit taper trial. Settings are Python objects, not CLI flags."""
from dataclasses import asdict
from pathlib import Path
import json
from hashlib import sha256
from geometry.taper.layout import TaperParams,write,validate
from geometry.taper.sonnet import configure,simulate,add_capacitor,solver_environment
from geometry.taper.fitting import fit
from geometry.prepare_sonnet import main as prepare
from extract_generated_inductance import main as extract_l
from geometry.specifications import TAPER_1,UNIT_CELL

ROOT=Path(__file__).resolve().parents[2]


def run_case(name,taper,*,d_input=6,d_junction=64,inductance_reference=None,input_style='end_cell',l_input=None,t_input=24,fit_response=True):
    """Each new name owns its artifacts. Explicit matching L records can be reused."""
    folder=ROOT/'output_data/taper_exploration'/name
    if (folder/'trial.json').exists():raise ValueError('Completed trial exists; use a new name')
    folder.mkdir(parents=True,exist_ok=True)
    if inductance_reference is None:
        inductance_reference,lr=measure_unit(name,taper)
    else:
        inductance_reference=Path(inductance_reference).resolve()
        iso_name=inductance_reference.parent.name
        metadata=json.loads((ROOT/'geometry/taper/output'/iso_name).with_suffix('.json').read_text())
        if metadata['taper']!=asdict(taper):raise ValueError('Inductance source geometry differs from this taper')
        gds=(ROOT/'geometry/taper/output'/iso_name).with_suffix('.gds');validate(gds,metadata)
        lr=json.loads(inductance_reference.read_text())
    gds=write(taper,ROOT/'geometry/taper/output'/name,d_input=d_input,d_junction=d_junction,input_style=input_style,l_input=l_input,t_input=t_input)
    project_dir=ROOT/'sonnet_files/taper_exploration'/name
    prepare(2,gds_path=gds,output_dir=project_dir,geometry_validator=validate)
    project=configure(project_dir/'generated_ncell_2.sonx')
    snp=simulate(project,folder)
    if not fit_response:
        return {'project':str(project),'snp':str(snp),'inductance_reference':str(inductance_reference),'status':'simulated_unfitted'}
    fits=fit(snp,folder,l_taper=lr['L0_nH'],l_regular=4.922723945068658,c_regular=188.8781163861302,free_regular=True)
    targets={'L0_nH':TAPER_1.L0_nH,'C0_fF':TAPER_1.C0_fF,
             'Cg_input_fF':TAPER_1.Cg_input_fF,'Cg_junction_fF':UNIT_CELL.Cg_input_fF}
    primary=fits['complex_with_delay']
    errors={key:100*(value/target-1) for key,value,target in [
        ('L0',lr['L0_nH'],targets['L0_nH']),('C0',primary['C_taper_fF'],targets['C0_fF']),
        ('Cg_input',primary['Cg_input_fF'],targets['Cg_input_fF']),('Cg_junction',primary['Cg_junction_fF'],targets['Cg_junction_fF'])]}
    result={'name':name,'taper':asdict(taper),'d_input':d_input,'d_junction':d_junction,'input_style':input_style,'l_input':l_input,'t_input':t_input,
            'inductance_reference':str(inductance_reference),'inductance_sha256':sha256(Path(inductance_reference).read_bytes()).hexdigest(),
            'primary_fit':'complex S21 with effective fixture delay; both self capacitances and all couplings free',
            'targets':targets,'fits':fits,'errors_percent':errors,'passes_2_percent':primary['success'] and all(abs(v)<=2 for v in errors.values()),
            'gds':str(gds),'project':str(project),'snp':str(snp),
            'gds_sha256':sha256(gds.read_bytes()).hexdigest()}
    (folder/'trial.json').write_text(json.dumps(result,indent=2)+'\n')
    print('TARGET ERRORS:',errors,flush=True)
    return result


def measure_unit(name,taper,*,sweep_range=None,w_top_slot=0,l_top_slot=100):
    iso_name=name+'_isolation';gds=write(taper,ROOT/'geometry/taper/output'/iso_name,isolation=True,
                                       w_top_slot=w_top_slot,l_top_slot=l_top_slot)
    project_dir=ROOT/'sonnet_files/taper_exploration'/iso_name
    prepare(1,gds_path=gds,output_dir=project_dir,geometry_validator=validate)
    project=project_dir/'generated_ncell_1.sonx'
    add_capacitor(project,json.loads(gds.with_suffix('.json').read_text()))
    if sweep_range is not None:
        configure(project,*sweep_range,target=round((sweep_range[1]-sweep_range[0])*1000))
    result_dir=ROOT/'output_data/taper_exploration'/iso_name
    lr=extract_l(project_file=project,output_dir=result_dir,fit_capacitance=False,target_l0_nh=TAPER_1.L0_nH,
                 solver_env=solver_environment())
    lr['target_L0_nH']=TAPER_1.L0_nH;lr['target_difference_percent']=100*(lr['L0_nH']/TAPER_1.L0_nH-1)
    (result_dir/'results.json').write_text(json.dumps(lr,indent=2)+'\n')
    inductance_reference=result_dir/'results.json'
    return inductance_reference,lr
