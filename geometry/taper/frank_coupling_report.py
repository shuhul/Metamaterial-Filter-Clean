"""Save the control-test conclusions and a draft update for Frank; no EM."""
import json
from pathlib import Path
import numpy as np
from geometry.taper.frank_coupling_controls import ROOT,NATIVE,OUT,CASES
from geometry.taper.run_bridged_stage import digest


def run():
    import matplotlib.pyplot as plt
    cases=(*CASES,'bridged_T1_existing')
    data={c:json.loads((OUT/c/'joint_port_fit/fit.json').read_text()) for c in cases}
    for r in data.values():
        candidates=r['multistart_solutions'];cost=min(q['cost'] for q in candidates)
        basin=[q['x'] for q in candidates if q['success'] and q['cost'] <= cost*1.01+1e-12]
        r['best_basin_start_count']=len(basin)
        r['best_basin_parameter_spread']=np.ptp(np.array(basin),axis=0).tolist() if basin else None
    validation=json.loads((OUT/'minimal_feed10_direct_validation/validation.json').read_text())
    primary=data['minimal_feed10'];frame=data['frame_without_inductor'];bridge=data['bridged_T1_existing']
    plain=json.loads((OUT/'minimal_feed10/N4_independent/fit.json').read_text())
    transfer=json.loads((OUT/'joint_reference_transfer/comparison.json').read_text())
    titles=dict(minimal_feed10='Coupling Electrode Only: 10 um Feed',minimal_feed14='Coupling Electrode Only: 14 um Feed',
                frame_without_inductor='T1 Frame Retained, Inductor Removed',bridged_T1_existing='Complete T1, Direct Feed')
    fig,axes=plt.subplots(2,2,figsize=(16,11))
    for ax,case in zip(axes.flat,cases):
        r=data[case];a=np.load(OUT/case/'joint_port_fit/response.npz')
        ax.plot(a['f']/1e9,abs(a['actual'][:,1,0]),color='#b22222',lw=2.8,label='Sonnet')
        ax.plot(a['f']/1e9,abs(a['model'][:,1,0]),'k--',lw=2.2,label='Joint complex S11/S21/S22 fit')
        text=f"Cg = {r['Cg_junction_fF']:.3f} fF\nRegular C0 = {r['C0_regular_fF']:.3f} fF\nS21 RMS = {r['S21_complex_RMS']:.4f}"
        ax.text(.03,.96,text,transform=ax.transAxes,va='top',fontsize=14,bbox=dict(boxstyle='round,pad=.4',fc='white',ec='gray'))
        ax.set_title(titles[case],fontsize=17,pad=12);ax.set_xlabel('Frequency (GHz)',fontsize=15)
        ax.set_ylabel(r'$|S_{21}|$',fontsize=20);ax.set_ylim(0,1.55);ax.tick_params(labelsize=13)
        ax.set_xlim(4.98,5.04)
        ax.legend(fontsize=10,loc='upper right');ax.grid(alpha=.25)
    fig.suptitle('T1-to-Regular Coupling: Controlled Geometry Comparison',fontsize=23)
    fig.tight_layout();fig.savefig(NATIVE/'coupling_controls.png',dpi=190);plt.close(fig)
    a=np.load(OUT/'minimal_feed10/N4_independent/response.npz')
    b=np.load(OUT/'minimal_feed10/joint_port_fit/response.npz')
    fig,ax=plt.subplots(figsize=(12,7))
    ax.plot(a['f']/1e9,abs(a['actual'][:,1,0]),color='#b22222',lw=3,label='Sonnet: coupling electrode only')
    ax.plot(a['f']/1e9,abs(a['model'][:,1,0]),'k--',lw=2.2,label='Cg-only fit; regular C0 fixed at 188.878 fF')
    ax.plot(b['f']/1e9,abs(b['model'][:,1,0]),color='#246cb5',ls=':',lw=2.7,label='Joint port fit; regular C0 and output Cg released')
    ax.set_xlim(4.99,5.04);ax.set_ylim(0,1.6)
    old=plain['fixed_regular']['C0_fF'];percent=100*(primary['C0_regular_fF']/old-1)
    text=(f"Fixed regular C0: Cg = {plain['Cg_junction_fF']:.3f} fF; S21 RMS = {plain['S21_complex_RMS']:.4f}\n"
          f"Joint fit: Cg = {primary['Cg_junction_fF']:.3f} fF; S21 RMS = {primary['S21_complex_RMS']:.4f}\n"
          f"Joint regular C0 = {primary['C0_regular_fF']:.3f} fF ({percent:+.2f}% from N=4 reference)")
    ax.text(.03,.96,text,transform=ax.transAxes,va='top',fontsize=13.5,bbox=dict(boxstyle='round,pad=.5',fc='white',ec='gray'))
    ax.set_title('Fixed Regular-Cell Parameters Bias the Coupling Extraction',fontsize=22,pad=14)
    ax.set_xlabel('Frequency (GHz)',fontsize=18);ax.set_ylabel(r'$|S_{21}|$',fontsize=23);ax.tick_params(labelsize=15)
    ax.legend(fontsize=11,loc='upper right',bbox_to_anchor=(1,.72));ax.grid(alpha=.25)
    fig.tight_layout();fig.savefig(NATIVE/'regular_reference_sensitivity.png',dpi=190);plt.close(fig)
    normal=np.load(OUT/'full_T1_free_junction/N4_independent/response.npz')
    forced=np.load(OUT/'joint_reference_transfer/full_T1_fixed_control/response.npz')
    released=np.load(OUT/'joint_reference_transfer/full_T1_free_junction/response.npz')
    fig,ax=plt.subplots(figsize=(12,7))
    ax.plot(normal['f']/1e9,abs(normal['actual'][:,1,0]),color='#b22222',lw=3,label='Sonnet: full T1 restored')
    ax.plot(normal['f']/1e9,abs(normal['model'][:,1,0]),color='#666666',ls='--',lw=2,label='Original regular references; junction Cg free')
    ax.plot(released['f']/1e9,abs(released['model'][:,1,0]),color='#246cb5',ls=':',lw=2.5,label='Isolated regular references; junction Cg free')
    ax.plot(forced['f']/1e9,abs(forced['model'][:,1,0]),'k--',lw=2.3,label='Isolated regular references and Cg fixed')
    ax.set_ylim(0,.58);ax.set_title('Coupling-Only Values Applied to Full T1',fontsize=23,pad=15)
    ax.set_xlabel('Frequency (GHz)',fontsize=18);ax.set_ylabel(r'$|S_{21}|$',fontsize=23);ax.tick_params(labelsize=15)
    ax.legend(fontsize=11,loc='upper right');ax.grid(alpha=.25)
    fig.tight_layout();fig.savefig(NATIVE/'full_T1_transfer.png',dpi=190);plt.close(fig)
    rows='\n'.join(f"| {titles[c]} | {data[c]['Cg_junction_fF']:.4f} | {data[c]['C0_regular_fF']:.4f} | {data[c]['Cg_output_fF']:.4f} | {data[c]['S21_complex_RMS']:.5f} |" for c in cases)
    draft=(f"I ran the coupling-only control, keeping the T1-to-unit coupling leg, gap/ground strip and unit-cell geometry unchanged, "
           f"but removing T1's meander and opposite capacitor leg. With the old fixed unit-cell C0=188.878 fF, a Cg-only fit gives "
           f"{plain['Cg_junction_fF']:.2f} fF but misses the linewidth/height (complex S21 RMS {plain['S21_complex_RMS']:.4f}). "
           f"A joint complex S11/S21/S22 fit gives Cg={primary['Cg_junction_fF']:.2f} fF, unit C0={primary['C0_regular_fF']:.2f} fF "
           f"and output Cg={primary['Cg_output_fF']:.2f} fF, with S21 RMS {primary['S21_complex_RMS']:.4f}. "
           f"Keeping T1's capacitor frame but removing only its meander instead gives effective Cg={frame['Cg_junction_fF']:.2f} fF. "
           f"With the complete T1 directly fed, the same joint fit gives {bridge['Cg_junction_fF']:.2f} fF. "
           f"So the fixed L/C concern was useful: the extracted Cg depends strongly on both the fixed reference values and which T1 metal remains. "
           f"However, transferring the stripped-layout values unchanged back to full T1 fails (S21 RMS 0.108), "
           f"so I cannot call the isolated 5.9 fF value the correct full-design Cg. The results point to additional loading/coupling "
           f"being absorbed into the scalar parameters. Two extra direct EM points confirm the narrow control resonance "
           f"to better than 3e-5 in all S parameters. Separately, the corrected same-layout ideal-5-nH test does not worsen the fit.")
    (NATIVE/'message_for_Frank.md').write_text(draft+'\n')
    summary=dict(joint_controls=data,primary_fixed_reference=plain,full_T1_transfer=transfer,direct_validation=validation,
                 production_parameters_changed=False,accepted_T2_design=False,
                 full_test_count=103,full_tests_passed=True,
                 selected_T1_sha256=digest(ROOT/'geometry/taper/selected_design.json'),
                 conclusion='Scalar parameters are fixture dependent; no certified transferable single Cg replacement.')
    (OUT/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    text=f'''# Frank's coupling-control tests — 23 October 2026

## Result

The requested stripped-T1 test was run. Its extracted coupling differs from
the complete T1 model, but it is **not valid to transfer that isolated number
unchanged into the full design**. Removing other T1 metal changes the fitted
effective coupling and regular-cell capacitance even though the retained
coupling leg and regular geometry are exactly preserved. This confirms a
parameter-transfer problem, not one uniquely identified wrong capacitor.
No production constants or accepted geometry were changed; no final T2 is claimed.

## Controlled geometries

All projects retain the exact T1 right leg (6 um wide, 310 um long), 56 um
gap with its 4 um grounded strip, regular resonator and feed-connected output
end. Native polygon XOR is empty from the retained coupling leg through the
entire right-hand layout. Same 50-ohm ports, 400 um reference planes, stack,
materials and 2 x 2 um grid. No T2 or input IDC appears in these controls.

- `minimal_feed10/generated_ncell_2.sonx`: remove the meander and opposite T1
  capacitor leg; use a 10 um feed up to the retained right coupling leg.
- `minimal_feed14/generated_ncell_2.sonx`: same removal, retain the original
  14 um top-bar thickness to check that feed-width change separately.
- `frame_without_inductor/generated_ncell_2.sonx`: retain the entire T1 cap
  frame and remove only its meander. This isolates the effect of that frame.
- Existing bridged T1 is reused from
  `../bridged_feed_tests/01_t1_input_bridged/generated_ncell_2.sonx`.
- Full T1 is reused from
  `../final_designs/finger536_t20_right310_input6_j56_final/generated_ncell_2.sonx`.

The generated filename's `ncell_2` is inherited metadata compatibility, not a
claim that the removed T1 still resonates. The stripped tests have one weakly
coupled regular resonance plus the feed-connected output-end LC.

## Extraction and comparison

The initial test fits only the coupling Cg and total matched fixture delay,
holding independently extracted regular L/C and output coupling fixed. For the
10 um stripped case it yields Cg={plain['Cg_junction_fF']:.5f} fF but S21 RMS
{plain['S21_complex_RMS']:.5f}, maximum error {plain['S21_complex_max']:.5f}. Its
peak height/linewidth mismatch makes that scalar value uncertified.

Releasing the regular C by less than 1% dramatically improves the fit. The
table below uses the **same joint complex S11/S21/S22 objective** in every row:
three electrical unknowns (junction Cg, regular C0, output Cg) and two matched
port delays. Both regular inductances stay fixed at {LR_text()} nH and the
output-end C stays at 188.8781163861302 fF. The complete directly fed T1 row
also includes its known L=6.014105932449249 nH and C0=101.4920484340242 fF.
No T1 parameter is released in these joint-control fits.

| Geometry | Effective junction Cg (fF) | Regular C0 (fF) | Output Cg (fF) | Complex S21 RMS |
| --- | ---: | ---: | ---: | ---: |
{rows}

These are effective parameters of the declared models, **not unique static
electrode capacitances**. The feed-width change alone shifts Cg about 2.2%.
The much larger frame/meander changes cannot be ignored when transferring
isolated circuit values. The extra paths have not been uniquely separated
into electric coupling, mutual inductance or distributed-ground effects.
S12 is a reciprocal consistency check, not independent physical validation.
The fits' source hashes, bounds, multistart solutions and scaled Jacobian
condition numbers are saved in their JSON results.
Some initial starts miss the narrow resonance and terminate at bounds or the
iteration limit; they are not accepted. The best low-error basin was reached
from multiple starts and is reported separately in `summary.json`.

## Restoring full T1

The full T1 Sonnet data and geometry remain untouched. As a diagnostic, its
L stays fixed and C0/input coupling/delay are refitted. Setting the isolated
junction Cg={primary['Cg_junction_fF']:.5f} fF, regular C0={primary['C0_regular_fF']:.5f}
fF and output Cg={primary['Cg_output_fF']:.5f} fF gives complex S21 RMS **0.10836**.
Letting junction Cg vary while keeping that same isolated regular-reference
tuple gives RMS **0.02252**, with Cg about 7.8059 fF. Using the original N=4
regular references and freeing junction Cg gives RMS **0.00829**, Cg about
8.1855 fF. These tests are explicitly diagnostics; refitted T1 values were
not promoted into the production selection.

Thus simply replacing the old full-design coupling by the isolated 5.9 fF
value **does not solve the complete T1 model**, and cannot certify T2.
Frank's concern about incorrect fixed inputs is supported; a claim that the
whole linewidth problem is one wrong geometric capacitor is not yet supported.

## Numerical/resource verification

Three new ABS sweeps: 4.75–5.3 GHz, 5501 output samples, four actual field
frequencies each; zero native errors/warnings. One solver at a time, four
threads, with the existing resource guard. No CPU stall/overload was required
to stop a run. Two additional **direct** field points at 5.0123 and 5.0163 GHz
check the primary peak's half-power flanks. Sonnet's direct export also retains
the four cached points, so it has six rows. Maximum discrepancy against ABS:
S21 {validation['S21_max_error']:.3g}; all S {validation['all_S_max_error']:.3g}.
The 4.060 MHz control linewidth is therefore not an interpolation artifact.
Native logs and resource records accompany each case. Absolute certainty or
fabrication accuracy cannot be claimed from one fixed mesh/equivalent model.

## Plots and data

- `coupling_controls.png`: four joint-fit control comparisons.
- `regular_reference_sensitivity.png`: failure of the one-Cg fit and improvement
  when the small regular-C reference discrepancy is allowed.
- `full_T1_transfer.png`: failure of the stripped-layout tuple to describe full T1.
- Each native case has `geometry_review.png`, project and matching `.snp`.
- Fit JSON/NPZ and validation:
  `output_data/taper_exploration/frank_coupling_controls_2026-10-23/`.
- `message_for_Frank.md`: draft text only; **nothing has been sent to Frank**.

Code: `geometry/taper/frank_coupling_controls.py`,
`frank_coupling_analysis.py`, `frank_coupling_report.py`.
New geometry/linewidth checks: `tests/test_frank_coupling_controls.py` and
`tests/test_minimum_inductance.py`. Test output is saved with the results.
The full Python suite passes **103 tests**; the selected T1 manifest hash is
unchanged. All EM jobs have exited.

The next step should validate a circuit description that retains the extra
coupling/loading from the full neighboring geometry, rather than assume one
isolated Cg remains invariant as structures are added. No further EM job is
needed merely to repeat the now-confirmed scalar-transfer failure.
'''
    (NATIVE/'README.md').write_text(text)
    return summary


def LR_text():
    from geometry.taper.staged_fitting import LR
    return str(LR)


if __name__=='__main__':run()
