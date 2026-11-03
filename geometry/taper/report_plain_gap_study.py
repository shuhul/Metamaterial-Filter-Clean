"""Collect the controlled study and maintain an explicit current handoff."""
from datetime import datetime, timezone
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from circuit_common import read_snp
from geometry.taper.plain_gap_study import ROOT, NATIVE, OUT, BASE_NATIVE, BASE_OUT, configurations
from geometry.taper.response_diagnostics import summarize_study
from geometry.taper.target_iteration import target_response


def findings(report):
    """Write interpretations separately from the machine-readable evidence."""
    base=report['baseline'];goal=base['target_peaks'];p=base['peaks']
    lines=['# Findings: approved geometry, parameter-only investigation','',
        '## What was established','',
        'This study retains the approved topology: flat upper ground, uniform input IDC, straight centered T2 inductor, '
        'and plain 4 um T2/T1 coupling gap. T1 and all downstream native metal are locked. '
        'The three variants change one T2/input parameter each. No topology changes or extra EM probes were introduced.',
        '', 'The acceptance test compares Sonnet complex S21 with the specified fixed-target circuit over 4.65–5.40 GHz. '
        'Only matched fixture delay is optimized. No capacitances, inductances or frequency shifts are fitted for this score.',
        '', 'The baseline lower/upper frequencies must move upward by '
        f"{(goal[0]['frequency_GHz']-p[0]['frequency_GHz'])*1000:.1f} / {(goal[1]['frequency_GHz']-p[1]['frequency_GHz'])*1000:.1f} MHz. "
        f"The linewidths must grow from {p[0]['power_FWHM_MHz']:.2f} / {p[1]['power_FWHM_MHz']:.2f} "
        f"to {goal[0]['power_FWHM_MHz']:.2f} / {goal[1]['power_FWHM_MHz']:.2f} MHz.",
        '', '## Measured geometry sensitivities','']
    for row in report['cases']:
        d=row['changes_vs_baseline'];widths=[v['power_FWHM_MHz'] for v in row['peaks']]
        percentages=[(new/old['power_FWHM_MHz']-1)*100 for new,old in zip(widths,p)]
        lines.append(f"- **{row['case']}**: frequencies {d[0]['frequency_shift_MHz']:+.2f} / {d[1]['frequency_shift_MHz']:+.2f} MHz; "
            f"linewidths {percentages[0]:+.1f}% / {percentages[1]:+.1f}%; fixed-target complex RMS {row['target_complex_RMS']:.5f}.")
    lines+=['', 'Increasing input IDC overlap is the cleanest measured linewidth control in this neighborhood: '
        'it broadened both peaks with negligible frequency shift. Increasing T2 facing-leg thickness moved both peaks '
        'downward with only modest broadening, so thickness increase is the wrong frequency direction here.',
        '', 'These are finite changes, not exact infinitesimal derivatives. Do not extrapolate them over many steps '
        'or add their effects without a validating EM run. A large frequency error compared with the linewidth also makes '
        'the overall RMS a poor stand-alone guide to which local parameter is physically useful.',
        '', '## What the additional S-parameters and fields say','',
        'Common-pole fits to complex S11/S21/S22 reproduce local baseline resonance shapes with RMS below 1e-5, '
        'giving loaded linewidths 2.488 / 2.624 MHz. This confirms that the narrow peaks are present in the EM response; '
        'they are not created by a degenerate capacitance fit. These poles do not separately identify internal/external Q or every lumped capacitor.',
        '', 'At the two baseline peaks, approximately 3.0% / 2.1% of incident power is absent from S11/S21. '
        'Reflection accounts for much more power than that deficit. The deficit includes possible radiation/material loss and numerical/calibration effects; '
        'it is not a direct measurement of conductor loss.',
        '', 'Both T1 and the regular-cell meander carry strong current at both resonances. '
        'The same vertical neck-current component has T1-minus-regular phase of about -179.7 degrees at the lower peak '
        'and -1.7 degrees at the upper peak. This is consistent with two coupled modes with different relative current phases. '
        'The maps do not establish unique capacitances or prove a particular parasitic-coupling mechanism.',
        '', 'The three direct current solves reproduce the saved baseline adaptive S-parameters to within 0.0005 '
        'maximum complex error across all S entries. Geometry, stack and ports are exactly unchanged. '
        'Whole-layout current/charge coordinates were checked against actual native metal: zero significant nonzero samples outside metal. '
        'Rect and Line exports failed validation and are explicitly rejected; unknown charge token also silently produced current, '
        'so surface charge was re-exported with token q and verified Coulomb/sq m units.',
        '', '## Circuit-only guide','',
        '`theory_sensitivity.json` changes one electrical capacitance by +10% around the theoretical circuit. '
        'Input coupling broadens both peaks by roughly 12% with <0.3 MHz frequency movement. '
        'T1 C0 and T2/T1 Cg move the modes and redistribute their linewidths. '
        'These are model-level trends; a geometry change may alter several effective elements together.',
        '', '## Next decision','',
        'The small controlled changes have not met RMS <=0.05. They also do not prove the topology cannot work. '
        'Use input-IDC overlap to control overall loading and a separate parameter to move the frequencies. '
        'Frame-width increase moved the peaks by less than 0.4 MHz, so it is a weak frequency control over the tested interval. '
        'A useful next parameter-only test is the plain T2/T1 gap (for example 4 to 6 um), keeping the same facing-leg topology '
        'and every other parameter fixed. The present study deliberately kept that gap at 4 um and has not simulated this next step. '
        'Increasing input overlap further can then address the linewidth separately, with another validating EM run. '
        'If these controls cannot cover the required frequency/linewidth changes, examine other permitted dimensions '
        'one at a time before considering any topology change. T1 changes and added structures are not part of this study.',
        '', 'No variant has been promoted to final_designs. The original baseline standalone T2 inductance was 0.73423 nH '
        '(0.78% below 0.74 nH); candidate environment-dependent inductances have not been independently re-extracted in this sensitivity study.',
        '', 'See README.md for the numerical table, parameter_comparison.png for the curves, '
        'and baseline_currents/current_charge_maps.png for the fields.']
    (OUT/'FINDINGS.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')


def main():
    report=summarize_study();base=read_snp(BASE_NATIVE.with_suffix('.snp'))
    completed=len(report['cases'])==len(configurations())
    if completed:findings(report)
    fig,axs=plt.subplots(1,3,figsize=(19,5),gridspec_kw={'width_ratios':[1.4,1,1]})
    curves=[('Baseline',base,'#777777')]
    colors=['#3a78ad','#bf6836','#5d9055']
    for row,color in zip(report['cases'],colors):
        curves.append((row['case'],read_snp(NATIVE/row['case']/'idc/generated_ncell_3.snp'),color))
    for ax in axs:
        for name,net,color in curves:ax.plot(net.f/1e9,abs(net.s[:,1,0]),color=color,lw=1.7,label=name)
        ax.plot(base.f/1e9,abs(target_response(base.f)),color='#202a35',ls='--',lw=2,label='Fixed target circuit')
        ax.set_xlabel('Frequency (GHz)',fontsize=14);ax.set_ylim(0,1.02);ax.grid(alpha=.2);ax.tick_params(labelsize=12)
    axs[0].set_xlim(4.65,5.4);axs[0].set_ylabel('|S21|',fontsize=18)
    axs[1].set_xlim(4.82,4.91);axs[2].set_xlim(5.055,5.13)
    for ax,title in zip(axs,['Full Scoring Band','Lower Resonance','Upper Resonance']):ax.set_title(title,fontsize=16)
    axs[0].legend(fontsize=10,loc='upper left')
    fig.suptitle('Approved Geometry: One-Parameter Response Changes',fontsize=22)
    fig.tight_layout();fig.savefig(OUT/'parameter_comparison.png',dpi=170);fig.savefig(OUT/'parameter_comparison.pdf');plt.close(fig)
    now=datetime.now(timezone.utc).isoformat()
    path=ROOT/'geometry/taper/taper2_status.json';status=json.loads(path.read_text())
    archive=OUT/'status_before_parameter_investigation.json'
    if not archive.exists():archive.write_text(json.dumps(status,indent=2)+'\n')
    # Archive prior mixed exploration rather than presenting rejected topologies
    # and old fit criteria as the active best design.
    status=dict(schema_version=3,updated_utc=now,
        state='plain_gap_sensitivity_complete' if completed else 'plain_gap_sensitivity_running',
        geometry_approved=True,final_T2_accepted=False,campaign=str(NATIVE.relative_to(ROOT)),
        baseline_project=str(BASE_NATIVE.relative_to(ROOT)),latest_plot=str((OUT/'parameter_comparison.png').relative_to(ROOT)),
        report=str((OUT/'README.md').relative_to(ROOT)),field_plot=str((OUT/'baseline_currents/current_charge_maps.png').relative_to(ROOT)),
        completed_cases=[r['case'] for r in report['cases']],
        geometry_constraints=dict(top_ground_edge='straight',T2_T1_gap_um=4,
            T2_T1_coupler='plain legs, no IDC',facing_lengths_um=[584,584],
            input_IDC='equal buses, uniform top-start fingers',T1_and_downstream='original native metal locked',grid_um=[2,2]),
        acceptance=dict(metric='complex_S21_RMS_against_fixed_target_circuit',maximum=.05,
            band_GHz=[4.65,5.4],free_parameter='matched fixture delay only',extracted_capacitances='heuristics'),
        baseline_complex_RMS=json.loads((BASE_OUT/'target_comparison.json').read_text())['complex_RMS'],
        cases=[dict(case=r['case'],complex_RMS=r['target_complex_RMS']) for r in report['cases']],
        next_step='Use measured parameter trends to choose next geometry parameters; no topology change authorized',
        historical_status_archive=str(archive.relative_to(ROOT)))
    path.write_text(json.dumps(status,indent=2)+'\n')
    NATIVE.mkdir(parents=True,exist_ok=True)
    (NATIVE/'README.md').write_text('# Approved plain-gap investigation\n\n'
        'The user approved this topology and authorized controlled parameter tests. '
        'No upper-ground notch, added T2/T1 IDC, pads, or other topology changes. '
        'T1 and downstream native metal are checked against the verified design.\n\n'
        'The three cases each change one parameter: input IDC span 140→156 um; '
        'T2 facing-leg thickness 20→24 um at fixed outer edge and 4 um gap; '
        'T2 frame width 400→416 um with the 60 um straight inductor centered.\n\n'
        'Run: `.venv\\Scripts\\python.exe -B -m geometry.taper.plain_gap_study`. '
        'Completed, hash-verified runs are reused. One guarded four-thread EM job at a time.\n\n'
        f'Results: `{OUT.relative_to(ROOT)}/README.md`; `parameter_comparison.png`; '
        '`baseline_currents/current_charge_maps.png`; numerical provenance in JSON beside each output.\n\n'
        'Baseline current diagnostics use exactly the approved native Geometry XML and three direct frequencies. '
        'Current and charge export is offline, using installed soncmd -JXYExport. '
        'Only Whole export files passed coordinate/metal-support checks. Older Rect/Line exports are rejected. '
        'The quantity token for charge is q; unknown charge token silently yielded current. '
        'Whole coordinates map as native Y=1500-exported Y. See field_validation.json.\n\n'
        'Sonnet current-density documentation: '
        'https://www.sonnetsoftware.com/support/doc_19.52/Topics/Current_Density_Plots/Current_Density_Plots.htm\n\n'
        'Export XML documentation (19.54; syntax verified locally in 19.52): '
        'https://www.sonnetsoftware.com/support/doc_19.54/Topics/Current_Density_Plots/Export_Data.htm\n',encoding='utf-8')
    old=ROOT/'sonnet_files/taper_exploration/target_response_2026-10-25/README.md'
    addition='## Superseding approved parameter investigation'
    if addition not in old.read_text(encoding='utf-8'):
        with old.open('a',encoding='utf-8') as handle:
            handle.write(f'\n\n{addition}\n\nThe user approved the plain 4 um-gap geometry and said "Ok go". '
                f'The active investigation is now `{NATIVE.relative_to(ROOT)}`; see its README. '
                'Earlier requests to await geometry approval are historical.\n')
    print(json.dumps(status,indent=2))


if __name__=='__main__':main()
