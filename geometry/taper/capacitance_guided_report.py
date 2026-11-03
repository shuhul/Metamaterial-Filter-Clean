"""Refresh the capacitance-guided campaign's readable evidence and handoff."""
from datetime import datetime,timezone
import json
import numpy as np
from circuit_common import read_snp
from geometry.taper.capacitance_guided_iteration import ROOT,OUT,NATIVE,network


def main():
    lines=['# Capacitance-guided iteration','',
        '## Rules','',
        'Same approved topology: flat top ground, uniform input IDC, centered straight T2 wire, '
        'plain T2/T1 facing legs. No extra IDC, ground notch, shield, pads, or pedestal. '
        'T1 and downstream native metal are checked against the verified design. '
        'All native vertices remain on the 2 um grid. One guarded four-thread EM job at a time.',
        '', 'Acceptance remains complex S21 RMS <=0.05 against the fixed target L/C circuit, '
        '4.65–5.40 GHz, with matched fixture delay only. Extraction is a geometry heuristic.',
        '', '## Extraction of existing simulations','',
        'Five taper/junction capacitances and separate left/right matched fixture delays were fitted to '
        'complex S11/S21/S22. All inductances and regular-cell parameters were fixed to the documented '
        'measured reference used in earlier extractions. These reference values differ slightly from the '
        'nominal target values used for acceptance; extraction never changes the target curve.',
        '', 'The narrow band still permits substantially different circuit parameter solutions. '
        'This is why the individual C0_T2/input_Cg values cannot be treated as certified physical measurements. '
        'A more stable heuristic across the recorded solutions is T1 C0 + T2/T1 Cg, approximately '
        '164–167 fF versus the target sum 156.7 fF. That suggests reducing assembled T1 capacitive loading. '
        'The measured IDC-length perturbation independently supports increasing input loading to broaden the peaks.',
        '', '| Saved case | Joint S RMS | C0 T2 (fF) | Input Cg (fF) | C0 T1 (fF) | T2/T1 Cg (fF) | Sum (fF) |',
        '|---|---:|---:|---:|---:|---:|---:|']
    for name in ('baseline','input_span156','facing_thickness24','frame_width416'):
        p=OUT/name/'extraction/extraction.json'
        if not p.exists():continue
        r=json.loads(p.read_text());v=r['values']
        lines.append(f"| {name} | {r['all_channel_complex_RMS']:.5f} | {v['C0_T2_fF']:.1f} | {v['Cg_input_T2_fF']:.1f} | {v['C0_T1_fF']:.1f} | {v['Cg_T2_T1_fF']:.1f} | {v['C0_T1_fF']+v['Cg_T2_T1_fF']:.1f} |")
        net=read_snp(r['provenance']['source']);y=np.column_stack([net.s[:,0,0],net.s[:,1,0],net.s[:,1,1]])
        # Re-score every saved alternative on the SAME full grid. Older raw
        # extraction files compared residual norms from different sample counts;
        # these rows supersede that historical similar-residual label.
        branches=[list(v.values())+r['fixture_delay_ps']]+r['similar_residual_start_solutions']
        scores=[dict(values=dict(zip(list(v)+['left_ps','right_ps'],b)),
            full_grid_joint_complex_RMS=float(np.sqrt(np.mean(abs(network(net.f,b)-y)**2)))) for b in branches]
        (p.parent/'extraction_branches.json').write_text(json.dumps(dict(
            note='Alternatives rescored on identical full grid; supersedes historical similar-residual list label.',
            branches=scores),indent=2)+'\n')
    lines+=['', 'The table lists numerical best solutions, including high-C/negative-delay branches. '
        'Those solutions are not accepted physical capacitor values. See each extraction_branches.json '
        'for other recorded solutions with modestly worse curve errors and very different individual capacitors.',
        '', '## Current maps','',
        'Baseline maps show strong T1 and regular-meander participation at both peaks. Their neck currents '
        'are approximately 180 degrees apart at the lower mode and aligned at the upper mode. '
        'This confirms coupled modes and provides an independent qualitative check. '
        'It does not uniquely identify the capacitances. Whole-layout fields passed native-metal support validation; '
        'rectangular and line export formats did not and are excluded.',
        '', '## New geometry runs','',
        '| Case | Status | Fixed-target complex RMS |', '|---|---|---:|']
    rows=[]
    for path in sorted(OUT.glob('*/idc/run_status.json')):
        r=json.loads(path.read_text());case=path.parents[1].name;s=path.parent/'target_comparison.json'
        score=json.loads(s.read_text()) if s.exists() else None
        rows.append(dict(case=case,status=r['status'],score=score))
        intent_path=path.parent/'design_intent.json'
        metadata=ROOT/'geometry/taper/output'/NATIVE.name/case/'idc/geometry.json'
        if intent_path.exists() and metadata.exists():
            from dataclasses import asdict
            from geometry.taper.second_layout import from_dict
            from geometry.taper.plain_gap_study import baseline,leaves
            intent=json.loads(intent_path.read_text());old=leaves(asdict(baseline()))
            actual=leaves(asdict(from_dict(json.loads(metadata.read_text())['second'])))
            intent['changed_parameters']={k:dict(baseline=old[k],candidate=v) for k,v in actual.items() if old[k]!=v}
            if 'slightly larger plain gap' in intent.get('rationale',''):
                intent['previous_generic_rationale']=intent['rationale']
                intent['rationale']='Parameter-only iteration guided by extracted loading and measured response sensitivities; see changed_parameters and campaign README.'
            intent_path.write_text(json.dumps(intent,indent=2)+'\n')
        text=f"{score['complex_RMS']:.5f}" if score else 'pending'
        lines.append(f"| {case} | {r['status']} | {text} |")
    lines+=['', 'First candidate: span256_gap6 increases input IDC overlap from 112 to 228 um and '
        'increases the existing plain T2/T1 gap from 4 to 6 um. Every other geometry parameter stays at baseline. '
        'This combines the measured input-loading direction with the extracted excess-T1-loading heuristic.',
        '', 'Each native folder contains generated_ncell_3.sonx and geometry_review.png. '
        'Each output folder contains solver logs/hashes, fixed-target S21 comparisons, and response diagnostics.',
        '', 'No extracted individual capacitor value is used as a completion criterion. '
        'A successful candidate will receive an independent final inductance check before promotion to final_designs.']
    lines+=['', '### Response features','',
        '| Case | Peaks (GHz) | Power FWHM (MHz) | Reflection-pole input proxies (MHz) |',
        '|---|---|---|---|']
    for row in rows:
        path=OUT/row['case']/'idc/response_diagnostics.json'
        if not path.exists():continue
        d=json.loads(path.read_text());peaks=', '.join(f"{p['frequency_GHz']:.4f}" for p in d['peaks'])
        widths=', '.join(f"{p['power_FWHM_MHz']:.3f}" if p['power_FWHM_MHz'] else 'undefined' for p in d['peaks'])
        loading=', '.join(f"{p['residue_magnitudes'][0]:.3f}" for p in d['local_poles'])
        lines.append(f"| {row['case']} | {peaks} | {widths} | {loading} |")
    lines+=['', 'The reflection-pole residues are loading proxies near isolated resonances, '
        'not independently certified external decay rates or unique capacitor measurements.']
    lines+=['', '## Narrow-frame follow-up','',
        'The first gap6/span256 trial overshot the peak frequencies and remained narrow. '
        'The next parameter-only experiment reduces T2 frame width from 400 to 200 um, '
        'shortening the top-bar paths between IDC, wire, and T1 coupler. '
        'The unchanged shared top ground remains flat. T1 and downstream shapes are still exactly locked. '
        'This is a test of a distributed-path hypothesis, not a demonstrated cause.',
        '', 'The narrower frame needs its own wire calibration. The unused meander.w metadata is reduced '
        'only to satisfy legacy meander-parameter validation; straight_inductor=True still renders one uniform rectangle.',
        '', '| Isolated case | L (nH) | Error from 0.74 nH | Held-out Y error (S) |',
        '|---|---:|---:|---:|']
    for path in sorted(OUT.glob('*/idc_Lcal/results.json')):
        r=json.loads(path.read_text())
        lines.append(f"| {path.parents[1].name} | {r['L0_nH']:.6f} | {r['target_error_percent']:+.3f}% | {r['held_out_max_Y_error_S']:.2e} |")
    lines+=['', '## What the later tests established','',
        'The 200 um frame with a 48 um wire calibrated to 0.74084 nH, but its full-response linewidths '
        'fell to 2.054/1.965 MHz. Narrowing the frame did not solve the loading problem.',
        '', 'The 400 um frame with a 60 um coupling leg and 2 um gap calibrated to 0.75038 nH. '
        'Its lower linewidth increased to 6.662 MHz, but the modes shifted to 4.7314/5.0453 GHz and the upper '
        'linewidth remained 3.049 MHz. Stronger plain coupling alone therefore does not meet the target.',
        '', 'The next larger test uses a 900 um frame, 40 um facing leg, 6 um gap, and 256 um IDC span. '
        'The frame-width sensitivity motivates this loading test; extrapolation is uncertain and the full EM curve decides. '
        'A 108 um wire calibrated too high (0.75866 nH), so it was replaced with a uniform 128 um wire '
        '(0.73806 nH, -0.262%). No full response was run for the failed 108 um calibration. '
        'The unused meander-width metadata is 400 um in these cases, but no meander is rendered.',
        '', 'The 900 um frame improved target RMS to 0.09700, with modes at 4.8864/5.1045 GHz and '
        'widths 5.947/4.583 MHz. Its joint extraction gave input Cg about 717 fF and effective T2/T1 Cg '
        'about 29 fF, while T2 C0 remained poorly identified near a bound. The next candidate uses '
        'a 1100 um frame, 176 um straight wire, IDC span 440 um, 6 um gap and 50 um facing leg. '
        'The larger IDC targets a modeled input-coupling region near 1.3 pF; this is a heuristic loading adjustment, '
        'not a revision of the 459 fF electrical target. Currents are saved in that full sweep.',
        '', '## Existing junction ground-position parameter','',
        'x_junction_ground is now available in SecondParams: positive moves the existing T1/regular ground strip '
        'toward the regular cell. It defaults to zero and is omitted from legacy metadata, so every archived default '
        'GDS still rebuilds identically. No simulation has used an offset yet at the time of adding this parameter. '
        'The motivation is to obtain a separate environmental-loading control without altering either cell shape. '
        'Native locking allows changes only inside the existing junction gap; GDS/native generation still validates '
        'the intended strip geometry there. Tests verify unchanged geometry outside that gap, flat upper ground, '
        'and grid/clearance bounds.',
        '', '## Conditional mutual-inductance hypothesis','',
        'Positive mutual inductance between the T2 wire and T1 can oppose capacitive cross admittance: '
        'Y01 = -j*w*Cg + j*M/[w*(L2*L1-M²)]. A capacitor-only model can then trade an apparent lower Cg '
        'against higher shunt C0. Assuming physical T1 C0=100 fF, one continuous recorded fit family '
        'illustrates M around 0.19 nH at 5 GHz. This is not an independent M measurement; other fit branches '
        'give different values. mutual_cancellation_hypothesis.json records the assumptions and alternatives. '
        'It motivates testing wire separation while acceptance remains the original fixed-target circuit.',
        '', '## Reference methodology','',
        'Frank Yang et al., Appendix C.2.b, describe isolated ideal-capacitor inductance sweeps, progressive '
        'two-capacitance extraction, then validation of the assembled taper/regular response. '
        'The current full two-port extraction is a heuristic because recorded staged constants are fixture dependent. '
        'Acceptance uses the assembled fixed-target response, not those fitted values. '
        '[Paper](https://arxiv.org/pdf/2607.05760), pp. 23–24.']
    OUT.mkdir(parents=True,exist_ok=True);(OUT/'README.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
    (NATIVE/'README.md').write_text('# Capacitance-guided parameter iterations\n\n'
        f'Read `{OUT.relative_to(ROOT)}/README.md` for methodology, numerical evidence, current status, and limitations.\n',encoding='utf-8')
    path=ROOT/'geometry/taper/taper2_status.json';status=json.loads(path.read_text())
    archive=OUT/'previous_status.json'
    if not archive.exists():archive.write_text(json.dumps(status,indent=2)+'\n')
    status.update(updated_utc=datetime.now(timezone.utc).isoformat(),state='capacitance_guided_iteration',
        campaign=str(NATIVE.relative_to(ROOT)),report=str((OUT/'README.md').relative_to(ROOT)),
        geometry_constraints=dict(top_ground_edge='straight',T2_T1_coupler='plain facing legs, gap varied as a parameter',
            T1_and_downstream='original native metal locked',input_IDC='equal buses, uniform fingers',grid_um=[2,2]),
        current_campaign_cases=rows,final_T2_accepted=False)
    complete=[r for r in rows if r['score'] is not None]
    if complete:
        best=min(complete,key=lambda r:r['score']['complex_RMS']);name=best['case']
        status.update(best_completed_case=name,best_completed_target_RMS=best['score']['complex_RMS'],
            latest_plot=str((OUT/name/'idc/s21_target_comparison.png').relative_to(ROOT)))
    if rows:
        last=next((r for r in rows if r['status']=='running'),rows[-1]);name=last['case']
        status.update(current_case=name,latest_project=str((NATIVE/name/'idc/generated_ncell_3.sonx').relative_to(ROOT)),
            latest_geometry_preview=str((NATIVE/name/'idc/geometry_review.png').relative_to(ROOT)))
    path.write_text(json.dumps(status,indent=2)+'\n')
    print(json.dumps([dict(case=r['case'],status=r['status'],RMS=r['score']['complex_RMS'] if r['score'] else None) for r in rows],indent=2))


if __name__=='__main__':main()
