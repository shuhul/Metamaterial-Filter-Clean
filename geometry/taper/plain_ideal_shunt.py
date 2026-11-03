"""Known ideal C0 addition to test a search direction, not a physical design.

The source GDS and all planar metal are unchanged. This fixture must never be
published as a completed geometry: a successful direction needs a real layout
and an independent ordinary two-port EM validation.
"""
import json,shutil
import xml.etree.ElementTree as ET
from geometry.taper.compact_stages import ROOT,NATIVE,OUT,GDS
from geometry.taper.sonnet import add_capacitor,configure,simulate
from geometry.taper.run_bridged_stage import digest
from geometry.prepare_sonnet import tag,BIN,run as native_command
from geometry.taper.fit_restored_t2 import fit


def run(cap_fF=420,name='plain_w500_wire52_idc30_span172_l504'):
    source=NATIVE/name/'idc/generated_ncell_3.sonx';meta=json.loads((GDS/name/'idc/geometry.json').read_text())
    dest=ROOT/'sonnet_files/taper_exploration/ideal_shunt_diagnostics'/f'{name}_c{cap_fF:g}'
    out=ROOT/'output_data/taper_exploration/ideal_shunt_diagnostics'/f'{name}_c{cap_fF:g}'
    dest.mkdir(parents=True,exist_ok=True);out.mkdir(parents=True,exist_ok=True)
    project=dest/'diagnostic_ideal_shunt.sonx'
    if not project.exists():
        shutil.copyfile(source,project);add_capacitor(project,meta)
        tree=ET.parse(project);g=tree.find(tag('Geometry'))
        var=next(q for q in g.findall(tag('Variable')) if q.get('Name')=='C_add');var.set('Value',str(cap_fF/1000))
        from geometry.prepare_sonnet import write_xml
        write_xml(tree,project);configure(project,4.65,5.4,7500)
        native_command([BIN/'soncmd.exe','-Read',project],dest/'native_read.log')
        intent=dict(diagnostic_only=True,physical_geometry_design=False,source=str(source),source_sha256=digest(source),
                    ideal_parallel_C_added_fF=cap_fF,unit='C_add in pF',T1_metal_unchanged=True)
        (out/'intent.json').write_text(json.dumps(intent,indent=2)+'\n')
        (dest/'DIAGNOSTIC_ONLY.md').write_text('# Ideal capacitor diagnostic only\n\nThis is not a manufacturable or final design. Planar metal matches the plain source. An ideal capacitor is added in Sonnet to test the direction of a possible T2 ground-capacitance change.\n')
    marker=out/'run_status.json'
    if marker.exists() and json.loads(marker.read_text()).get('status')=='complete':
        record=json.loads(marker.read_text())
        if digest(project)!=record['project_sha256'] or digest(project.with_suffix('.snp'))!=record['snp_sha256']:
            raise ValueError('Completed shunt control changed')
        snp=project.with_suffix('.snp')
    else:snp=simulate(project,out)
    ell=json.loads((OUT/name/'idc_Lcal/results.json').read_text())['L0_nH']
    return fit(snp,out/'ideal_ladder_fit',ell,downstream_g=56.63965936392767)


if __name__=='__main__':run()
