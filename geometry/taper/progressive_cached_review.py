"""Reassess existing low-L field data before deciding which geometry to solve."""
import json
from circuit_common import read_snp
from geometry.taper.progressive_fixtures import ROOT,OUT
from geometry.taper.progressive_fitting import fit
from geometry.taper.minimum_inductance_review import PROFILES,LT1,LR


def main():
    p=PROFILES['full_T1_staged_calibration']
    cases=[('w500_top10_wire52_gap6','direct',.7327695472953403),
           ('plain_w700_wire68_gap6','direct',.7404213459916129),
           ('plain_w1000_wire104_gap6','direct',.7428384176222952),
           ('plain_w500_wire52_idc26_r584_t24_gap6','idc',.7341365873386546),
           ('plain_w500_wire52_idc26_r584_t24_gap2','idc',.7333007434145957),
           ('plain_w700_wire68_idc26_r584_t24_gap6','idc',.7405033660195242)]
    rows={}
    for name,stage,l in cases:
        source=ROOT/f'sonnet_files/taper_exploration/compact_centered/{name}/{stage}/generated_ncell_3.snp'
        if not source.exists():continue
        try:
            net=read_snp(source)
        except ValueError:
            print('SKIP invalid/sparse export',source,flush=True);continue
        if len(net.f)<101:
            print('SKIP sparse screen, not a fitted sweep',source,flush=True);continue
        cs=[100 if stage=='idc' else 0,p['ct'],p['cr'],p['ce']]
        gs=[459 if stage=='idc' else None,56.7,p['gj'],p['go']]
        free=[('C',1),('C',2),('C',3),('G',1),('G',2),('G',3)]
        if stage=='idc':free=[('C',0),('G',0)]+free
        rows[name]=fit(source,OUT/'cached_screen'/name,[l,LT1,LR,LR],cs,gs,free,
                       title=f'Cached Low-L Diagnostic: {name}')
    (OUT/'cached_screen/summary.json').write_text(json.dumps(rows,indent=2)+'\n')


if __name__=='__main__':main()
