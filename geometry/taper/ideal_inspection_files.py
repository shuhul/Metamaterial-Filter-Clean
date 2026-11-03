"""Prepare the requested matched no-L / ideal-5nH Sonnet inspection pair.

Never runs a field solve. Cache-only native recombination is a separate action.
"""
import json
import shutil
import xml.etree.ElementTree as ET
import numpy as np
from circuit_common import read_snp
from geometry.taper.compact_stages import ROOT
from geometry.prepare_sonnet import tag,write_xml,BIN,run
from geometry.taper.bridged_fixtures import require_native_locked
from geometry.taper.run_bridged_stage import digest

DEST=ROOT/'sonnet_files/taper_exploration/ideal_inductor_audit_2026-10-23'


def prepare():
    base=ROOT/'sonnet_files/taper_exploration/compact_centered/w500_wire52_gap6_wire_absent/direct/generated_ncell_3.sonx'
    loaded=ROOT/'sonnet_files/taper_exploration/ideal_wire_diagnostic/wire_terminal_control.sonx'
    meta=json.loads((ROOT/'geometry/taper/output/compact_centered/w500_wire52_gap6_wire_absent/direct/geometry.json').read_text())
    rows=[]
    for source,folder,stem in ((base,'01_T1_no_T2_inductor','T1_reference'),
                               (loaded,'02_T1_plus_ideal_5nH','T1_plus_ideal_5nH')):
        dest=DEST/folder;dest.mkdir(parents=True,exist_ok=True);project=dest/(stem+'.sonx')
        if not project.exists():
            tree=ET.parse(source)
            if source==loaded:
                comp=next(tree.iter(tag('Component')))
                comp.set('Type','IND');comp.set('Units','NH');comp.set('Name','L_T2_ideal_5nH')
                value=ET.Element(tag('Value'));value.text='5';comp.insert(0,value)
            for q in tree.iter(tag('File')):
                if q.get('Type')=='TOUCH':q.set('Name',stem+'.snp')
            write_xml(tree,project)
            cache_source=source.parent/'sondata'/source.stem
            cache_dest=dest/'sondata'/stem
            if cache_source.exists():shutil.copytree(cache_source,cache_dest,dirs_exist_ok=True)
            if source==base:shutil.copyfile(source.with_suffix('.snp'),project.with_suffix('.snp'))
        require_native_locked(project,meta)
        run([BIN/'soncmd.exe','-Read',project],dest/'native_read.log')
        rows.append(dict(source=str(source),source_sha256=digest(source),project=str(project),
                         project_sha256=digest(project),T1_native_XOR_empty=True,
                         ideal_L_nH=5 if source==loaded else None))
    (DEST/'manifest.json').write_text(json.dumps(rows,indent=2)+'\n')
    return rows


def plots():
    """Fit the exact ordinary two-port exports accompanying the review files."""
    from geometry.taper.ideal_topology_diagnostics import fit
    from geometry.taper.minimum_inductance import OUT as SWEEP_OUT
    import matplotlib.pyplot as plt
    results=[]
    for folder,stem,l in (('01_T1_no_T2_inductor','T1_reference',None),
                          ('02_T1_plus_ideal_5nH','T1_plus_ideal_5nH',5.)):
        dest=DEST/folder;project=dest/(stem+'.sonx');snp=project.with_suffix('.snp')
        n=read_snp(snp);data=dict(f=n.f,actual=n.s)
        cache_check=None
        if l is not None:
            expected=np.load(SWEEP_OUT/'ideal/L5nH/fit/response.npz')
            if not np.allclose(n.f,expected['f'],rtol=0,atol=1):raise ValueError('Different frequency grids')
            cache_check=float(abs(n.s-expected['actual']).max())
            if cache_check>1e-8:raise ValueError('Native ideal inductor differs from Python load')
        prior,old_model,_=fit(l,'current_staged',data=data)
        consistent,new_model,_=fit(l,'selected_with_end_LC',data=data)
        report=dict(project=str(project),project_sha256=digest(project),source_snp=str(snp),snp_sha256=digest(snp),
                    native_ideal_vs_Python_max_all_S_error=cache_check,
                    old_mixed_references=prior,consistent_T1_references=consistent,
                    T1_geometry_and_fitted_T1_parameters_unchanged=True,
                    scope='Diagnostic curve comparison; not a final T2 design or accepted component extraction')
        (dest/'comparison.json').write_text(json.dumps(report,indent=2)+'\n')
        np.savez(dest/'comparison.npz',f=n.f,actual=n.s,old_model=old_model,consistent_model=new_model)
        fig,ax=plt.subplots(figsize=(12,7))
        ax.plot(n.f/1e9,abs(n.s[:,1,0]),color='#b22222',lw=3,label='Sonnet simulation')
        ax.plot(n.f/1e9,abs(old_model),'k--',lw=2,label='Fit using staged-extraction constants')
        ax.plot(n.f/1e9,abs(new_model),color='#246cb5',ls=':',lw=2.7,label='Fit using earlier T1-fit constants')
        title='T1 Reference: T2 Inductor Removed' if l is None else 'T1 with an Ideal 5 nH T2 Inductor'
        ax.set_title(title,fontsize=23,pad=15)
        text=(f"Same frame and coupling electrode; 50-ohm ports\n"
              f"Complex S21 RMS:\n"
              f"Staged constants: {prior['S21_complex_RMS']:.4f}\n"
              f"Earlier T1-fit constants: {consistent['S21_complex_RMS']:.4f}\n"
              f"T1 parameters held fixed; diagnostic fit")
        ax.set_ylim(0,float(abs(n.s[:,1,0]).max())*1.6)
        ax.text(.035,.96,text,transform=ax.transAxes,va='top',fontsize=14,
                bbox=dict(boxstyle='round,pad=.5',fc='white',ec='gray'))
        ax.set_xlabel('Frequency (GHz)',fontsize=19);ax.set_ylabel(r'$|S_{21}|$',fontsize=23)
        ax.tick_params(labelsize=16);ax.legend(fontsize=12,loc='upper right');ax.grid(alpha=.25)
        fig.tight_layout();fig.savefig(dest/'s21_comparison.png',dpi=210);plt.close(fig)
        results.append(report)
    (DEST/'comparison_summary.json').write_text(json.dumps(results,indent=2)+'\n')
    return results


if __name__=='__main__':
    import sys
    if len(sys.argv)>1 and sys.argv[1]=='plots':plots()
    else:
        for row in prepare():print(row['project'],flush=True)
