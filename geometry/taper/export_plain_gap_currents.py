"""Export native current/charge diagnostics with installed soncmd, no EM runs."""
import json
import subprocess
import csv
import xml.etree.ElementTree as ET
from geometry.prepare_sonnet import BIN
from geometry.taper.plain_gap_study import NATIVE, OUT


def export(kind='jxy',complex_data=True,*,leads=False):
    if kind not in ('jxy','jx','jy','q'):
        raise ValueError('Installed 19.52 exporter requires q for surface charge; unknown types silently fall back to current')
    project=NATIVE/'baseline_currents/generated_ncell_3.sonx'
    out=OUT/'baseline_currents';dest=out/'exports';dest.mkdir(parents=True,exist_ok=True)
    record=json.loads((out/'run_status.json').read_text())
    if record['status']!='complete' or not record['reference_check_pass']:
        raise ValueError('Current diagnostic must complete and match baseline before export')
    root=ET.Element('JXY_Export_Set')
    files=[]
    stem=kind+('_leads' if leads else '_whole')
    for label,f in zip(('lower','middle','upper'),record['frequencies_GHz']):
        file=dest/f'{stem}_{label}.csv';files.append(file)
        entry=ET.SubElement(root,'JXY_Export',Filename=str(file.resolve()),Label=f'{kind} {f} GHz')
        if leads:ET.SubElement(entry,'Region',Style='Line',Vertical='no',Position='1026')
        else:ET.SubElement(entry,'Region',Style='Whole')
        ET.SubElement(entry,'Levels',Range='Some',Start='0',Stop='0')
        ET.SubElement(entry,'Grid',XStep='2' if leads else '4',YStep='2' if leads else '4')
        ET.SubElement(entry,'Measurement',Complex='Yes' if complex_data else 'No',Type=kind)
        drive=ET.SubElement(entry,'Drive')
        for port,voltage in ((1,1),(2,0)):
            ET.SubElement(drive,'JXYPort',Capacitance='0',Inductance='0',Number=str(port),
                Phase='0',Reactance='0',Resistance='50',Voltage=str(voltage))
        locator=ET.SubElement(entry,'Locator')
        ET.SubElement(locator,'Frequency',Value=str(f*1e9))
    control=dest/f'export_{stem}.xml';ET.ElementTree(root).write(control,encoding='utf-8',xml_declaration=True)
    result=subprocess.run([str(BIN/'soncmd.exe'),'-JXYExport',str(control),str(project)],
                          cwd=project.parent,capture_output=True,text=True,timeout=120)
    (dest/f'export_{stem}.log').write_text(result.stdout+result.stderr)
    manifest=dict(type=kind,complex_data=complex_data,returncode=result.returncode,
                  files=[dict(path=str(p),exists=p.exists(),bytes=p.stat().st_size if p.exists() else 0) for p in files],
                  excitation='Port 1: 1 V source with 50 ohm series resistance; port 2: 50 ohm load',
                  note='Current/charge fields are not deembedded; geometry/stack/ports identical to baseline',
                  accepted_for_field_analysis=False,
                  coordinate_validation='Line format failed metal-support validation' if leads else
                      'Pending plot_plain_gap_fields whole-layout metal-support validation')
    headers=[]
    for p in files:
        if p.exists():
            with p.open() as handle:
                rows=csv.reader(handle)
                for _ in range(7):next(rows)
                headers.append(next(rows))
    manifest['measurement_headers']=headers
    manifest['verified_measurement']=bool(len(headers)==3 and all(
        ('Coulomb' in ','.join(h) if kind=='q' else 'Amps/Meter' in ','.join(h)) for h in headers))
    (dest/f'export_{stem}.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print(json.dumps(manifest,indent=2))
    print((result.stdout+result.stderr)[:2500])
    if not manifest['verified_measurement']:raise ValueError('Exported measurement does not match request')
    return manifest


if __name__=='__main__':
    import sys
    export(sys.argv[1] if len(sys.argv)>1 else 'jxy',leads='--leads' in sys.argv)
