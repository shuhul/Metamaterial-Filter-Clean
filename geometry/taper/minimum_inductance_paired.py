"""Pair every physical width with an equal-L ideal load of the cached EM fixture.

No field solves: terminate the validated four-port wire-absent control with L.
Both responses use identical frequency samples, downstream constants and fit rules.
This diagnoses model/geometry sensitivity; it does not uniquely identify a cause.
"""
import json
import numpy as np
import skrf
from geometry.taper.compact_stages import ROOT
from geometry.taper.minimum_inductance import OUT
from geometry.taper.ideal_wire_control import load_inductor
from geometry.taper.response_cases import W
from geometry.taper.run_bridged_stage import digest


def compare(width,l,physical):
    from geometry.taper.minimum_inductance_review import fit,PROFILES
    source=ROOT/'sonnet_files/taper_exploration/ideal_wire_diagnostic/wire_terminal_control.snp'
    control=ROOT/'output_data/taper_exploration/ideal_wire_diagnostic/open_control.json'
    if not json.loads(control.read_text())['control_passed']:
        raise ValueError('Ideal terminal control failed')
    n=W['read_network'](source,4)
    folder=OUT/f'meander_w{width}'/'paired_ideal';folder.mkdir(parents=True,exist_ok=True)
    # Exact same independently measured L, never a rounded nearby sweep value.
    s=load_inductor(n.s,n.f,l)
    path=folder/'ideal_loaded.s2p'
    skrf.Network(frequency=skrf.Frequency.from_f(n.f,unit='hz'),s=s,z0=50).write_touchstone(str(path))
    rows={}
    import matplotlib.pyplot as plt
    for profile in PROFILES:
        ideal=fit(path,l,folder/profile,profile,data_label='Cached EM fixture + ideal inductor')
        actual=physical[profile]
        a=np.load(OUT/f'meander_w{width}'/'review'/profile/'response.npz')
        b=np.load(folder/profile/'response.npz')
        if not np.allclose(a['f'],b['f'],rtol=0,atol=1):
            raise ValueError('Paired frequency grids differ')
        fig,axes=plt.subplots(1,2,figsize=(15,6),sharex=True,sharey=True)
        for ax,data,r,label in ((axes[0],a,actual,'Physical meander'),(axes[1],b,ideal,'Ideal inductor in EM fixture')):
            ax.plot(data['f']/1e9,abs(data['actual']),color='#b22222',lw=2.5,label=label)
            ax.plot(data['f']/1e9,abs(data['model']),'k--',lw=2,label='Same fixed downstream circuit + fit')
            ax.set_title(f"{label}\nRMS {r['complex_RMS']:.5f}; {'PASS' if r['assessment']['pass_curve'] else 'FAIL'}",fontsize=16)
            ax.set_xlabel('Frequency (GHz)',fontsize=15);ax.tick_params(labelsize=12)
            ax.grid(alpha=.25);ax.legend(fontsize=10)
        axes[0].set_ylabel(r'$|S_{21}|$',fontsize=20)
        fig.suptitle(f'Equal-L Comparison: {l:.5f} nH — {profile}',fontsize=18)
        fig.tight_layout();fig.savefig(folder/profile/'paired_s21.png',dpi=180);plt.close(fig)
        rows[profile]=dict(physical_RMS=actual['complex_RMS'],ideal_RMS=ideal['complex_RMS'],
                           physical_pass=actual['assessment']['pass_curve'],ideal_pass=ideal['assessment']['pass_curve'],
                           ideal_fit=ideal,plot=str(folder/profile/'paired_s21.png'))
    result=dict(L_nH=l,operator_sha256=digest(source),open_control_sha256=digest(control),
                profiles=rows,interpretation='Ideal failure prevents attributing physical failure solely to the inductor geometry; ideal success isolates an additional physical-geometry effect but not its unique mechanism.')
    (folder/'comparison.json').write_text(json.dumps(result,indent=2)+'\n')
    return result
