"""Edit the settings below, then run this file. No command-line arguments needed."""
from dataclasses import asdict, replace
from datetime import datetime, timezone
from hashlib import sha256
from pathlib import Path
import json
import re
import subprocess
import sys
import time
import xml.etree.ElementTree as ET

BASE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(BASE))
import numpy as np
import gdstk
import constants as cfg
from geometry.capacitor import CapParams
from geometry.meander import MeanderParams
from geometry.pad import PadParams
from geometry.unit_cell import UnitCellParams
from geometry.ncell import NCellParams, make_ncell
from geometry.generate_ncell import save_preview
from geometry.validation import require, validate_gds_metadata
from geometry.prepare_sonnet import main as prepare_sonnet, tag, write_xml, points, SOURCE
from geometry.prepare_inductance import main as prepare_inductance
from extract_generated_inductance import main as extract_inductance
from extract_circuit_params_ncell import main as fit_array
from circuit_common import read_snp, magnitude_db
from circuit_models import s21_n_cells
from geometry.taper.solver_guard import run as guarded_solver

# ======================== EDIT YOUR SETTINGS HERE ========================
RUN_NAME = "my_filter_01"  # New name when generating a changed design.

GENERATE_GEOMETRY = True
PREPARE_SONNET = True
SIMULATE_INDUCTANCE = False  # True: build N=1 and run the four C_add values.
SIMULATE_ARRAY = True       # True: run the N-cell EM analysis once.
FIT_S21 = True              # Fit C0/Cg and save s21_fit.png.
USE_SAVED_INDUCTANCE = True # False forces FIXED_L0_NH when the sweep is off.

# All geometry is in um. n_cells is the number of INTERIOR cells (+2 ends).
PARAMS = NCellParams(
    unit=UnitCellParams(
        cap=CapParams(w=302, l=600, t_top=14, t_side=6),
        meander=MeanderParams(w=186, l=596, l_top=10, n=72, t=4, d=4),
        pad=PadParams(w=24, l=586),
        d_top=8,
        d_pad=12,
    ),
    n_cells=4,
    d_cell=56,
    t_ground=4,
    d_end_sep=14,
    w_side=412,
    l_top_ground=448,
    l_sim=1500,
    t_feed=10,
    d_feed=6,
)

# When the sweep is off: None reuses this RUN_NAME's extracted L0 if available,
# otherwise uses FIXED_L0_NH. Or supply an explicit results.json path.
# Re-extract L0 when changing the meander, material or substantial loading.
INDUCTANCE_RESULT_FILE = None
FIXED_L0_NH = 4.922723945068658

ARRAY_SWEEP_GHZ = (4.5, 5.5)
ARRAY_SWEEP_TARGET = 500
# Covers the added-cap resonances of the current geometry; expand if needed.
INDUCTANCE_SWEEP_GHZ = (4.35, 4.95)
INDUCTANCE_SWEEP_TARGET = 1500
C_ADD_SWEEP_PF = (0.02, 0.03, 0.04, 0.05)  # 20, 30, 40, 50 fF.

# Validated import uses 2 x 2 um mesh and 400 um reference planes.
# Materials/stackup come from the original four-cell Sonnet template.
# ======================= END OF EDITABLE SETTINGS ========================

GEOMETRY_ROOT = BASE / "geometry/output"
PROJECT_ROOT = cfg.SONNET_FILES_DIR
RESULT_ROOT = cfg.OUTPUT_DATA_DIR


def digest(path):
    return sha256(Path(path).read_bytes()).hexdigest()


def save_json(path, data):
    Path(path).write_text(json.dumps(data, indent=2)+"\n", encoding="utf-8")


def check_settings():
    require(bool(re.fullmatch(r"[A-Za-z0-9_-]+", RUN_NAME)),
            "RUN_NAME may contain only letters, numbers, underscores and hyphens")
    for flag in (GENERATE_GEOMETRY, PREPARE_SONNET, SIMULATE_INDUCTANCE, SIMULATE_ARRAY, FIT_S21, USE_SAVED_INDUCTANCE):
        require(isinstance(flag, bool), "Stage flags must be True or False")
    require(any((GENERATE_GEOMETRY, PREPARE_SONNET, SIMULATE_INDUCTANCE, SIMULATE_ARRAY, FIT_S21)),
            "Enable at least one stage")
    require(not (FIT_S21 and not SIMULATE_ARRAY and (GENERATE_GEOMETRY or PREPARE_SONNET)),
            "For saved-data fitting, also set GENERATE_GEOMETRY and PREPARE_SONNET to False")
    for bounds, target in ((ARRAY_SWEEP_GHZ, ARRAY_SWEEP_TARGET),
                           (INDUCTANCE_SWEEP_GHZ, INDUCTANCE_SWEEP_TARGET)):
        require(len(bounds)==2 and np.all(np.isfinite(bounds)) and 0<bounds[0]<bounds[1],
                "Frequency ranges need positive increasing endpoints in GHz")
        require(isinstance(target, int) and not isinstance(target, bool) and target>0,
                "Sweep targets must be positive integers")
    require(len(C_ADD_SWEEP_PF)>=3 and np.all(np.isfinite(C_ADD_SWEEP_PF))
            and min(C_ADD_SWEEP_PF)>0 and np.all(np.diff(C_ADD_SWEEP_PF)>0),
            "C_ADD_SWEEP_PF needs at least three positive increasing values")
    if FIT_S21 and not SIMULATE_INDUCTANCE and (not USE_SAVED_INDUCTANCE or INDUCTANCE_RESULT_FILE is None):
        require(np.isfinite(FIXED_L0_NH) and FIXED_L0_NH>0, "FIXED_L0_NH must be positive")


def generate(params, stem):
    cell = make_ncell(params)
    lib = gdstk.Library(unit=1e-6, precision=1e-9)
    lib.add(cell, *cell.dependencies(True))
    top = lib.new_cell("TOP")
    top.add(gdstk.Reference(cell))
    lib.write_gds(str(stem.with_suffix(".gds")))
    save_preview(cell, params, stem.with_suffix(".png"))
    save_json(stem.with_suffix(".json"), {
        "parameters":asdict(params), "units":"um", "total_cells":params.total_cells,
        "w_sim":params.w_sim, "l_total":params.l_total,
        "intended_box":params.box, "cell_centers":params.cell_centers,
    })
    print(f"Generated: {stem.with_suffix('.gds')}", flush=True)


def check_geometry(params, stem):
    metadata = json.loads(stem.with_suffix(".json").read_text())
    require(metadata["parameters"]==asdict(params),
            "Saved geometry differs from PARAMS. Use a new RUN_NAME and regenerate")
    validate_gds_metadata(stem.with_suffix(".gds"), metadata)


def configure_project(project, bounds, target):
    tree = ET.parse(project)
    active = 0
    for sweep_set in tree.findall(f"{tag('Sweeps')}/{tag('Set')}"):
        if sweep_set.get("On")=="TRUE":
            sweeps=sweep_set.findall(f"{tag('Frequencies')}/{tag('Sweep')}")
            require(len(sweeps)==1, "Expected one active frequency sweep")
            sweep=sweeps[0]
            sweep.set("Start", str(bounds[0])); sweep.set("Stop", str(bounds[1]))
            sweep.set("Target", str(target))
            active+=1
    require(active==1, "Expected one active sweep set")
    for export in tree.findall(f"{tag('OutputFiles')}//{tag('File')}"):
        if export.get("Type")=="TOUCH":
            export.set("Name", project.stem+".snp")
    write_xml(tree, project)
    validation_path=project.parent/"validation.json"
    if validation_path.exists():
        record=json.loads(validation_path.read_text())
        record.update(project_sha256=digest(project), configured_range_GHz=list(bounds), configured_target=target)
        save_json(validation_path, record)


def check_project_geometry(project, params, stem, *, allow_component=False):
    """Catch externally edited/stale native geometry before using it."""
    _, polygons=validate_gds_metadata(stem.with_suffix('.gds'),json.loads(stem.with_suffix('.json').read_text()))
    expected=[]
    for polygon in polygons:
        xy=polygon.points.copy()
        xy[:,1]+=params.l_sim-params.l_top_ground
        xy=np.floor(xy/2+.5)*2
        xy[:,1]=params.l_sim-xy[:,1]
        expected.append(gdstk.Polygon(xy))
    geometry=ET.parse(project).find(tag('Geometry'))
    actual=[gdstk.Polygon(points(p)) for p in geometry.findall('.//'+tag('PlanarPolygon'))]
    require(not gdstk.boolean(actual,expected,'xor',precision=.001),
            'Sonnet polygons differ from configured GDS; prepare a fresh project')
    require(allow_component or not geometry.findall('.//'+tag('Component')),
            'Array project contains unexpected ideal components')
    size=geometry.find(f"{tag('Box')}/{tag('Size')}")
    require(float(size.get('X'))==params.w_sim and float(size.get('Y'))==params.l_sim,
            'Sonnet box differs from configured geometry')


def input_identity(project, stem):
    return {"project_sha256":digest(project), "gds_sha256":digest(stem.with_suffix('.gds')),
            "metadata_sha256":digest(stem.with_suffix('.json'))}


def simulate_array(project, stem, analysis):
    print(f"Running array EM: {project}", flush=True)
    log_path=analysis/"array_simulation.log"
    with log_path.open("w") as log:
        process=guarded_solver([str(cfg.EM_EXE), "-v", str(project)], cwd=project.parent,
                               stdout=log)
    process.check_returncode()
    text=log_path.read_text()
    require("Em simulation completed" in text and "Em errors detected: 0" in text,
            "Sonnet did not confirm successful array analysis; inspect array_simulation.log")
    snp=project.with_suffix('.snp')
    read_snp(snp)
    (analysis/"array_response.snp").write_bytes(snp.read_bytes())
    (analysis/"array_project.sonx").write_bytes(project.read_bytes())
    save_json(project.parent/"simulation_record.json", {
        **input_identity(project,stem), "snp_sha256":digest(snp), "log":str(log_path),
        "completed_utc":datetime.now(timezone.utc).isoformat(),
    })


def check_saved_array(project, stem):
    record=json.loads((project.parent/"simulation_record.json").read_text())
    for key,value in input_identity(project,stem).items():
        require(record[key]==value, "Project/GDS changed since simulation. Run SIMULATE_ARRAY=True")
    require(record["snp_sha256"]==digest(project.with_suffix('.snp')),
            "Saved S-parameter file changed since simulation")


def saved_inductance():
    if not USE_SAVED_INDUCTANCE:
        return float(FIXED_L0_NH), {"type":"user_fixed", "L0_nH":float(FIXED_L0_NH)}
    if INDUCTANCE_RESULT_FILE is None:
        pointer=PROJECT_ROOT/RUN_NAME/'inductance'/'inductance_record.json'
        if not pointer.exists():
            print(f'Using explicit FIXED_L0_NH={FIXED_L0_NH:.9f}; no new inductance sweep',flush=True)
            return float(FIXED_L0_NH), {"type":"user_fixed", "L0_nH":float(FIXED_L0_NH)}
        saved=json.loads(pointer.read_text())
        path=Path(saved['path'])
        require(saved['sha256']==digest(path),'Previously extracted L0 result changed')
    else:
        path=Path(INDUCTANCE_RESULT_FILE)
    result=json.loads(path.read_text())
    require(np.isfinite(result['L0_nH']) and result['L0_nH']>0, "Invalid saved inductance")
    if 'parameters' in result:
        require(result['parameters']==asdict(replace(PARAMS,n_cells=1)),
                "Saved inductance geometry differs. Enable SIMULATE_INDUCTANCE or choose an explicit fixed L0")
        require(result['template_sha256']==digest(SOURCE), "Material/template changed since L0 extraction")
    else:
        require(result.get('meander_w_um')==PARAMS.unit.meander.w,
                "Legacy saved L0 uses a different meander span; extract L0 for this geometry")
        print("Using historical native L0: only its meander span is machine-checked; "
              "other geometry/material transfer assumptions remain yours to assess.", flush=True)
    return float(result['L0_nH']), {"type":"saved", "path":str(path), "sha256":digest(path)}


def main():
    check_settings()
    stamp=datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    stem=GEOMETRY_ROOT/RUN_NAME
    array_folder=PROJECT_ROOT/RUN_NAME/"array"
    project=array_folder/f"generated_ncell_{PARAMS.n_cells}.sonx"
    single_stem=GEOMETRY_ROOT/(RUN_NAME+"_single")
    single_folder=PROJECT_ROOT/RUN_NAME/"inductance"
    single_base=single_folder/"generated_ncell_1.sonx"
    single_project=single_folder/"generated_inductance.sonx"
    single_params=replace(PARAMS,n_cells=1)
    # Validate early, before creating files or launching expensive stages.
    from geometry.validation import validate_grid_params
    validate_grid_params(PARAMS)
    if GENERATE_GEOMETRY:
        for candidate in [stem]+([single_stem] if SIMULATE_INDUCTANCE else []):
            require(not candidate.with_suffix('.gds').exists() and not candidate.with_suffix('.json').exists(),
                    "RUN_NAME already has geometry. Pick a new name, or disable generation for saved-data reuse")
    if PREPARE_SONNET:
        require(not project.exists(), "Array project already exists. Pick a new RUN_NAME or disable preparation")
        if SIMULATE_INDUCTANCE:
            require(not single_base.exists() and not single_project.exists(),
                    "Inductance project exists. Use a new name or disable preparation")
    if not GENERATE_GEOMETRY:
        check_geometry(PARAMS,stem)
        if SIMULATE_INDUCTANCE:check_geometry(single_params,single_stem)
    if not PREPARE_SONNET and (SIMULATE_ARRAY or FIT_S21):
        require(project.exists(), "Array project missing; enable PREPARE_SONNET")
    if not PREPARE_SONNET and (SIMULATE_ARRAY or FIT_S21):
        check_project_geometry(project,PARAMS,stem)
    if not PREPARE_SONNET and SIMULATE_INDUCTANCE:
        require(single_project.exists(),'Single-cell project missing; enable generation and preparation')
        check_project_geometry(single_project,single_params,single_stem,allow_component=True)
    if FIT_S21 and not SIMULATE_ARRAY:
        check_saved_array(project,stem)
    if SIMULATE_INDUCTANCE:
        inductance, provenance=None, None
    elif FIT_S21:
        inductance, provenance=saved_inductance()
    else:
        inductance, provenance=None, None
    result_root=RESULT_ROOT/RUN_NAME
    analysis=result_root if not result_root.exists() else result_root/("analysis_"+stamp)
    analysis.mkdir(parents=True)
    GEOMETRY_ROOT.mkdir(parents=True,exist_ok=True)
    status={"state":"running", "started_utc":datetime.now(timezone.utc).isoformat()}
    config={"run_name":RUN_NAME,"parameters":asdict(PARAMS),"flags":{
        "generate":GENERATE_GEOMETRY,"prepare":PREPARE_SONNET,"simulate_inductance":SIMULATE_INDUCTANCE,
        "simulate_array":SIMULATE_ARRAY,"fit":FIT_S21},
        "array_range_GHz":ARRAY_SWEEP_GHZ,"array_target":ARRAY_SWEEP_TARGET,
        "inductance_range_GHz":INDUCTANCE_SWEEP_GHZ,"inductance_target":INDUCTANCE_SWEEP_TARGET,
        "C_add_pF":C_ADD_SWEEP_PF,"fixed_L0_nH":FIXED_L0_NH,
        "use_saved_inductance":USE_SAVED_INDUCTANCE,
        "saved_L0_file":str(INDUCTANCE_RESULT_FILE) if INDUCTANCE_RESULT_FILE is not None else None}
    save_json(analysis/"configuration.json",config)
    save_json(analysis/"run_status.json",status)
    timings={}
    started=time.perf_counter()
    try:
        if GENERATE_GEOMETRY:
            generate(PARAMS,stem)
            if SIMULATE_INDUCTANCE:generate(single_params,single_stem)
        check_geometry(PARAMS,stem)
        if PREPARE_SONNET:
            prepare_sonnet(PARAMS.n_cells,gds_path=stem.with_suffix('.gds'),output_dir=array_folder)
            configure_project(project,ARRAY_SWEEP_GHZ,ARRAY_SWEEP_TARGET)
            if SIMULATE_INDUCTANCE:
                prepare_sonnet(1,gds_path=single_stem.with_suffix('.gds'),output_dir=single_folder)
                prepare_inductance(gds_path=single_stem.with_suffix('.gds'),project_path=single_base,
                    output_path=single_project,sweep=(*INDUCTANCE_SWEEP_GHZ,INDUCTANCE_SWEEP_TARGET))
        if SIMULATE_INDUCTANCE:
            require(single_project.exists(),"Inductance project missing; enable generation and preparation")
            check_project_geometry(single_project,single_params,single_stem,allow_component=True)
            configure_project(single_project,INDUCTANCE_SWEEP_GHZ,INDUCTANCE_SWEEP_TARGET)
            step=time.perf_counter()
            result=extract_inductance(project_file=single_project,output_dir=analysis/"inductance",
                                     c_add_sweep_pf=C_ADD_SWEEP_PF,fit_capacitance=False)
            result.update(parameters=asdict(single_params),template_sha256=digest(SOURCE))
            result_path=analysis/"inductance/results.json"
            save_json(result_path,result)
            inductance=float(result['L0_nH'])
            provenance={"type":"new_extraction","path":str(result_path.resolve()),"sha256":digest(result_path)}
            save_json(single_folder/"inductance_record.json",provenance)
            timings['inductance_seconds']=time.perf_counter()-step
        if SIMULATE_ARRAY:
            check_project_geometry(project,PARAMS,stem)
            configure_project(project,ARRAY_SWEEP_GHZ,ARRAY_SWEEP_TARGET)
            step=time.perf_counter()
            simulate_array(project,stem,analysis)
            timings['array_seconds']=time.perf_counter()-step
        if FIT_S21:
            check_saved_array(project,stem)
            print(f"Fitting N={PARAMS.n_cells} with fixed L0={inductance:.9f} nH",flush=True)
            step=time.perf_counter()
            fit_array(['--no-simulate','--snp',str(project.with_suffix('.snp')),
                       '--n-cells',str(PARAMS.n_cells),'--l0-nh',str(inductance),'--output-dir',str(analysis)])
            fitted=json.loads((analysis/'fit_parameters.json').read_text())
            net=read_snp(project.with_suffix('.snp'))
            model=s21_n_cells(2*np.pi*net.f,fitted['C0_fF'],fitted['Cg_fF'],inductance*1e-9,PARAMS.n_cells)
            fitted.update(parameters=asdict(PARAMS),inductance_provenance=provenance,
                project_sha256=digest(project),snp_sha256=digest(project.with_suffix('.snp')),
                Cg_error_percent=100*(fitted['Cg_fF']/cfg.TARGET_CG_FF-1),
                fit_rms_db=float(np.sqrt(np.mean((magnitude_db(model,clip=True)-magnitude_db(net.s[:,1,0],clip=True))**2))))
            save_json(analysis/'fit_parameters.json',fitted)
            # Archive the exact data used, preserving evidence across future reruns.
            (analysis/'array_response.snp').write_bytes(project.with_suffix('.snp').read_bytes())
            timings['fit_seconds']=time.perf_counter()-step
            print(f"S21 fit plot: {(analysis/'s21_fit.png').resolve()}",flush=True)
    except Exception as error:
        status.update(state='failed',error=str(error));raise
    else:
        status.update(state='complete')
        save_json(result_root/'latest_run.json',{'analysis_directory':str(analysis.resolve()),'flags':config['flags']})
    finally:
        timings['total_seconds']=time.perf_counter()-started
        status['finished_utc']=datetime.now(timezone.utc).isoformat()
        save_json(analysis/'stage_timings.json',timings)
        save_json(analysis/'run_status.json',status)
    print(f"Results: {analysis.resolve()}",flush=True)


if __name__=='__main__':
    main()
