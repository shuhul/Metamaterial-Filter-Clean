"""Build, simulate and fit an isolated N=4 array with 186 um meander spans."""
from dataclasses import asdict
from datetime import datetime,timezone
from hashlib import sha256
from pathlib import Path
import json
import argparse
import re
import subprocess
import time
import sys
import xml.etree.ElementTree as ET

BASE=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(BASE))
import gdstk
import numpy as np
import constants as cfg
from geometry.meander import MeanderParams
from geometry.unit_cell import UnitCellParams
from geometry.ncell import NCellParams,make_ncell
from geometry.generate_ncell import save_preview
from geometry.prepare_sonnet import main as prepare,tag,write_xml
from circuit_common import read_snp,magnitude_db
from circuit_models import s21_n_cells
from extract_circuit_params_ncell import main as fit


def resource_snapshot():
    """Best-effort Windows snapshot; diagnostic failures do not block a run."""
    command = "Get-CimInstance Win32_PerfFormattedData_PerfProc_Process | Where-Object { $_.Name -notin @('_Total','Idle') } | Sort-Object PercentProcessorTime -Descending | Select-Object -First 6 Name,IDProcess,PercentProcessorTime | ConvertTo-Json -Compress"
    try:
        result=subprocess.run(['powershell.exe','-NoProfile','-Command',command],capture_output=True,text=True,timeout=8)
        return json.loads(result.stdout) if result.returncode==0 else {'error':result.stderr}
    except Exception as error:
        return {'error':str(error)}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--d-cell',type=float,default=60,help='Outer cap-to-cap spacing in um')
    parser.add_argument('--run-name',help='Optional unique label to preserve earlier experiment folders')
    args=parser.parse_args()
    if not np.isfinite(args.d_cell) or args.d_cell<=0 or args.d_cell%2:
        parser.error('d-cell must be a positive multiple of 2 um')
    if args.run_name is not None and not re.fullmatch(r'[A-Za-z0-9_-]+',args.run_name):
        parser.error('run-name must contain only letters, numbers, underscores or hyphens')
    suffix='' if args.d_cell==60 else f'_d{args.d_cell:g}'
    if args.run_name:
        suffix+='_'+args.run_name
    stem=BASE/'geometry'/'output'/('ncell_4_w186'+suffix)
    folder=cfg.SONNET_FILES_DIR/('generated_ncell_4_w186'+suffix)
    output=cfg.OUTPUT_DATA_DIR/('generated_ncell_4_w186'+suffix)
    if folder.exists() or output.exists():
        raise FileExistsError('Refusing to overwrite an existing 186 um array experiment')
    output.mkdir(parents=True)
    status={'state':'running','started_utc':datetime.now(timezone.utc).isoformat()}
    status_path=output/'run_status.json'
    status_path.write_text(json.dumps(status,indent=2)+'\n')
    timings={}
    stage_start=time.perf_counter()
    try:
        # Only the horizontal run span changes. All six physical cells share it.
        params=NCellParams(n_cells=4,d_cell=args.d_cell,unit=UnitCellParams(meander=MeanderParams(w=186)))
        cell=make_ncell(params)
        lib=gdstk.Library(unit=1e-6,precision=1e-9)
        lib.add(cell,*cell.dependencies(True));top=lib.new_cell('TOP');top.add(gdstk.Reference(cell))
        lib.write_gds(str(stem.with_suffix('.gds')))
        save_preview(cell,params,stem.with_suffix('.png'))
        metadata={'parameters':asdict(params),'units':'um','total_cells':params.total_cells,
                  'w_sim':params.w_sim,'l_total':params.l_total,'intended_box':params.box,
                  'cell_centers':params.cell_centers}
        stem.with_suffix('.json').write_text(json.dumps(metadata,indent=2)+'\n')
        prepare(4,gds_path=stem.with_suffix('.gds'),output_dir=folder)
        project=folder/'generated_ncell_4.sonx'
        # Explicit naming avoids native tools retaining a template export basename.
        tree=ET.parse(project)
        for export in tree.findall(f'.//{tag("OutputFiles")}//{tag("File")}'):
            if export.get('Type')=='TOUCH':export.set('Name',project.stem+'.snp')
        write_xml(tree,project)
        validation_path=folder/'validation.json'
        validation=json.loads(validation_path.read_text())
        validation['project_sha256']=sha256(project.read_bytes()).hexdigest()
        validation_path.write_text(json.dumps(validation,indent=2)+'\n')
        timings['prepare_seconds']=time.perf_counter()-stage_start
        (output/'resources_before_solver.json').write_text(json.dumps(resource_snapshot(),indent=2)+'\n')
        stage_start=time.perf_counter()
        print(f'Running four-cell simulation: {project}',flush=True)
        with (folder/'simulation.log').open('w') as log:
            process=subprocess.run([str(cfg.EM_EXE),'-v',str(project)],cwd=folder,
                                   stdout=log,stderr=subprocess.STDOUT,text=True)
        timings['solver_seconds']=time.perf_counter()-stage_start
        (output/'resources_after_solver.json').write_text(json.dumps(resource_snapshot(),indent=2)+'\n')
        stage_start=time.perf_counter()
        process.check_returncode()
        log=(folder/'simulation.log').read_text()
        if 'Em simulation completed' not in log or 'Em errors detected: 0' not in log:
            raise RuntimeError('Solver log does not confirm successful completion')
        single=cfg.OUTPUT_DATA_DIR/'meander_tuning_ref400'/'best_refined'/'results.json'
        result=json.loads(single.read_text())
        if result['meander_w_um']!=186:raise ValueError('Inductance result does not match meander span')
        inductance=result['L0_nH']
        snp=project.with_suffix('.snp')
        fit(['--no-simulate','--snp',str(snp),'--n-cells','4','--l0-nh',str(inductance),
             '--output-dir',str(output)])
        fitted=json.loads((output/'fit_parameters.json').read_text())
        net=read_snp(snp)
        response=s21_n_cells(2*np.pi*net.f,fitted['C0_fF'],fitted['Cg_fF'],inductance*1e-9,4)
        fitted.update(meander_w_um=186,d_cell_um=args.d_cell,
            target_Cg_fF=cfg.TARGET_CG_FF,
            Cg_error_percent=100*(fitted['Cg_fF']/cfg.TARGET_CG_FF-1),
            Cg_within_2_percent=bool(abs(fitted['Cg_fF']/cfg.TARGET_CG_FF-1)<=.02),
            inductance_source=str(single),
            inductance_source_sha256=sha256(single.read_bytes()).hexdigest(),
            project_sha256=sha256(project.read_bytes()).hexdigest(),
            snp_sha256=sha256(snp.read_bytes()).hexdigest(),
            fit_rms_db=float(np.sqrt(np.mean((magnitude_db(response,clip=True)-magnitude_db(net.s[:,1,0],clip=True))**2))),
            note='L0 transferred from single-cell extraction; fitted effective capacitances of the N=4 model.')
        (output/'fit_parameters.json').write_text(json.dumps(fitted,indent=2)+'\n')
        print(json.dumps(fitted,indent=2),flush=True)
    except Exception as error:
        status.update(state='failed',error=str(error));raise
    else:
        status.update(state='complete',fit_sha256=sha256((output/'fit_parameters.json').read_bytes()).hexdigest())
    finally:
        timings['last_stage_seconds']=time.perf_counter()-stage_start
        (output/'stage_timings.json').write_text(json.dumps(timings,indent=2)+'\n')
        status['updated_utc']=datetime.now(timezone.utc).isoformat()
        status_path.write_text(json.dumps(status,indent=2)+'\n')


if __name__=='__main__':
    main()
