"""Offline circuit sensitivity from completed EM fits; no solver or acceptance.

The local mapping can fail away from the measured geometries. Save predictions
as guidance only, then compare every chosen geometry against fresh EM.
"""
import json
import numpy as np
from scipy.optimize import least_squares,minimize_scalar
from geometry.taper.taper2_v3 import ROOT, CAMPAIGN, transmission
from geometry.specifications import TAPER_1, TAPER_2, UNIT_CELL


def main():
    base = ROOT / 'output_data/taper_exploration' / CAMPAIGN
    def load(name):
        return json.loads((base / name / 'response/extraction/fit.json').read_text())
    d = load('t2_c264_top10_w152_pad532_idc34_n72_j112')
    e = load('t2_c264_top10_w152_pad532_idc34_n72_j112_rg64tg4')
    h = load('t2_c264_top10_w152_pad532_idc34_n72_j112_rg62tg6')
    name = 'clear344_p300_d60_i34n70_j136_r62'
    ref = load(name)
    c1_ground = (d['C_fF'][1] - e['C_fF'][1]) / 2
    c1_clear = (d['C_fF'][1] - h['C_fF'][1]) / 2
    gu_ground = (d['G_fF'][2] - e['G_fF'][2]) / 2
    gu_clear = (d['G_fF'][2] - h['G_fF'][2]) / 2
    data = np.load(base / name / 'response/comparison.npz')
    f = data['f'][::5]
    cells = [TAPER_2, TAPER_1, UNIT_CELL, UNIT_CELL]
    target = transmission(f, [v.L0_nH for v in cells], [v.C0_fF for v in cells],
                          [v.Cg_input_fF for v in cells], 0)
    rows = []
    for ground in (2, 4, 6):
        for clear in (28, 30, 32, 34, 36, 38, 40):
            c1 = ref['C_fF'][1] + c1_ground*(ground-6) + c1_clear*(clear-28)
            gu = ref['G_fF'][2] + gu_ground*(ground-6) + gu_clear*(clear-28)
            def residual(x):
                delta = transmission(f, ref['L_nH'], [x[0], c1, *ref['C_fF'][2:]],
                                     [x[1], x[2], gu, ref['G_fF'][3]], x[3]) - target
                return np.r_[delta.real, delta.imag]
            result = least_squares(residual, [80,180,8,25],
                       bounds=([0,10,.1,-100], [400,700,40,100]),
                       x_scale=[100,150,10,20], max_nfev=250)
            rows.append(dict(ground_um=ground, clearance_um=clear,
                        estimated_C1=c1, estimated_Gunit=gu,
                        model_RMS=float(np.sqrt(np.mean(result.fun**2)*2)),
                        required_C2_G2_G1_delay=result.x.tolist()))
    report = dict(scope='Local fit derivative estimate, not an EM prediction or acceptance',
                  rows=sorted(rows, key=lambda r:r['model_RMS']))
    (base / 'shared_junction_circuit_guidance.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report['rows'][:8], indent=2))


def ideal_guidance(name):
    """Find circuit directions around any completed case; no geometry claim."""
    base = ROOT / 'output_data/taper_exploration' / CAMPAIGN / name / 'response'
    data = np.load(base / 'comparison.npz')
    ref = json.loads((base / 'extraction/fit.json').read_text())
    f = data['f'][::5]
    cells = [TAPER_2, TAPER_1, UNIT_CELL, UNIT_CELL]
    target = transmission(f, [v.L0_nH for v in cells], [v.C0_fF for v in cells],
                          [v.Cg_input_fF for v in cells], 0)
    rows=[]
    for c1 in (ref['C_fF'][1], 144, 146, 148, 150, 152, 154):
        def residual(x):
            delta = transmission(f, ref['L_nH'], [x[0],c1,*ref['C_fF'][2:]],
                                 [x[1],x[2],x[3],ref['G_fF'][3]],x[4])-target
            return np.r_[delta.real,delta.imag]
        starts=[[ref['C_fF'][0],ref['G_fF'][0],ref['G_fF'][1],8.2,0],
                [50,200,10,8.2,20],[140,150,15,8.2,10]]
        fits=[least_squares(residual,x,bounds=([0,1,.1,4,-100],[400,700,50,12,100]),
                           x_scale=[100,150,10,8,20],max_nfev=400) for x in starts]
        result=min(fits,key=lambda r:r.fun@r.fun)
        rows.append(dict(C1_fF=c1, model_RMS=float(np.sqrt(np.mean(result.fun**2)*2)),
                         C2_G2_G1_Gunit_delay=result.x.tolist()))
    report=dict(scope='Circuit-only directions; not an EM prediction or acceptance', rows=rows)
    (base/'circuit_guidance.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))
    return report


def ledger(active_case=None):
    """Summarize actual completed target comparisons, never predictions."""
    base = ROOT / 'output_data/taper_exploration' / CAMPAIGN
    rows = []
    for path in base.glob('*/response/target_comparison.json'):
        score = json.loads(path.read_text())
        name = path.parent.parent.name
        fitted = json.loads((path.parent/'extraction/fit.json').read_text())
        measured = json.loads((path.parent.parent/'Lcal/results.json').read_text())
        geometry = json.loads((ROOT/'geometry/taper/output'/CAMPAIGN/name/'response/geometry.json').read_text())
        rows.append(dict(case=name, target_RMS=score['complex_RMS'],
                    target_pass=score['pass'], fit_RMS=fitted['complex_RMS'],
                    fit_at_bound=bool(any(fitted['active_bounds'])),
                    L_nH=measured['L0_nH'], L_error_percent=measured['target_error_percent'],
                    unit=geometry['unit'], coupler=geometry['coupler'],
                    snp_sha256=score['source_sha256'], plot=str(path.parent/'s21_fit.png')))
    rows.sort(key=lambda row:row['target_RMS'])
    report=dict(scope='Completed Sonnet results only; heuristic fits are not physical component acceptance',
                active_case=active_case, trials=rows)
    (base/'completed_trials.json').write_text(json.dumps(report,indent=2)+'\n')
    lines=['# Completed T2 trials', '', 'All target scores use the same fixed nominal circuit, 4.65–5.40 GHz and 7501 samples.',
           'The only target-curve adjustment is one matched delay. Predictions and unrun layouts are excluded.', '',
           '| Case | Target RMS | Fit RMS | L (nH) | L error | Fit at bound |',
           '|---|---:|---:|---:|---:|---|']
    for row in rows:
        lines.append(f"| {row['case']} | {row['target_RMS']:.6f} | {row['fit_RMS']:.6f} | {row['L_nH']:.7f} | {row['L_error_percent']:+.3f}% | {row['fit_at_bound']} |")
    (base/'completed_trials.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
    print(json.dumps(dict(active_case=active_case,best=rows[0]['case'],best_target_RMS=rows[0]['target_RMS'],completed=len(rows)),indent=2))
    return report


def junction_prediction(name):
    """Transport measured junction sensitivities; only an EM test can accept."""
    base=ROOT/'output_data/taper_exploration'/CAMPAIGN
    def load(case):return json.loads((base/case/'response/extraction/fit.json').read_text())
    def vector(q):return np.array([q['C_fF'][0],q['G_fF'][0],q['C_fF'][1],q['G_fF'][1],q['G_fF'][2]])
    d=load('t2_c264_top10_w152_pad532_idc34_n72_j112')
    e=load('t2_c264_top10_w152_pad532_idc34_n72_j112_rg64tg4')
    h=load('t2_c264_top10_w152_pad532_idc34_n72_j112_rg62tg6')
    dg=(vector(d)-vector(e))/2; dc=(vector(d)-vector(h))/2
    ref=load(name); data=np.load(base/name/'response/comparison.npz'); f=data['f'][::5]
    obs=data['observed'][::5]
    geometry=json.loads((ROOT/'geometry/taper/output'/CAMPAIGN/name/'response/geometry.json').read_text())
    coupling=geometry['coupler']; old_g=coupling.get('t_regular_ground',6)
    old_c=(coupling.get('d_regular',66)-old_g)/2
    cells=[TAPER_2,TAPER_1,UNIT_CELL,UNIT_CELL]
    target=transmission(f,[v.L0_nH for v in cells],[v.C0_fF for v in cells],[v.Cg_input_fF for v in cells],0)
    original=transmission(f,ref['L_nH'],ref['C_fF'],ref['G_fF'],ref['delay_ps'])
    rows=[]
    for ground in (2,4,6,8):
        for clear in range(20,41,2):
            values=vector(ref)+dg*(ground-old_g)+dc*(clear-old_c)
            if np.any(values<=0):continue
            cp=[values[0],values[2],*ref['C_fF'][2:]]
            gp=[values[1],values[3],values[4],ref['G_fF'][3]]
            curve=obs*transmission(f,ref['L_nH'],cp,gp,ref['delay_ps'])/original
            def residual(x):
                delta=curve*np.exp(-2j*np.pi*f*x[0]*1e-12)-target
                return np.r_[delta.real,delta.imag]
            fits=[least_squares(residual,[tau],bounds=([-100],[100]),max_nfev=100) for tau in (-30,0,30)]
            best=min(fits,key=lambda q:q.fun@q.fun)
            rows.append(dict(ground_um=ground,clearance_um=clear,gap_um=ground+2*clear,
                        estimated_C2_G2_C1_G1_Gunit=values.tolist(),
                        predicted_target_RMS=float(np.sqrt(np.mean(best.fun**2)*2))))
    report=dict(scope='Transferred fit sensitivity; not a simulation or acceptance',source=name,
                rows=sorted(rows,key=lambda q:q['predicted_target_RMS']))
    (base/name/'response/junction_proposal.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report['rows'][:8],indent=2));return report


def shortened_pad_guidance(name,delta_c_values=(-15,-20,-25,-30)):
    """Optimize couplings for hypothetical C0 reduction, then require actual EM."""
    base=ROOT/'output_data/taper_exploration'/CAMPAIGN/name/'response'
    proposal=json.loads((base/'junction_proposal.json').read_text())
    ref=json.loads((base/'extraction/fit.json').read_text());data=np.load(base/'comparison.npz')
    f=data['f'][::5]; cells=[TAPER_2,TAPER_1,UNIT_CELL,UNIT_CELL]
    target=transmission(f,[v.L0_nH for v in cells],[v.C0_fF for v in cells],[v.Cg_input_fF for v in cells],0)
    rows=[]
    for source in proposal['rows']:
        if source['ground_um']!=6 or source['gap_um'] not in (54,58,62,66):continue
        for delta_c in delta_c_values:
            c2,g2,c1,g1,gu=source['estimated_C2_G2_C1_G1_Gunit'];c2+=delta_c
            def residual(x):
                curve=transmission(f,ref['L_nH'],[c2,c1,*ref['C_fF'][2:]],
                              [x[0],x[1],gu,ref['G_fF'][3]],x[2])
                delta=curve-target
                return np.r_[delta.real,delta.imag]
            fits=[least_squares(residual,[g2*s,g1,0],bounds=([50,5,-100],[800,40,100]),
                     x_scale=[250,20,20],max_nfev=250) for s in (.8,1,1.2)]
            best=min(fits,key=lambda q:q.fun@q.fun)
            rows.append(dict(ground_um=source['ground_um'],gap_um=source['gap_um'],
                         assumed_C2_reduction_fF=delta_c,C2=c2,C1=c1,Gunit=gu,
                         needed_G2_G1_delay=best.x.tolist(),baseline_G2_G1=[g2,g1],
                         predicted_target_RMS=float(np.sqrt(np.mean(best.fun**2)*2))))
    report=dict(scope='Hypothetical pad-capacitance reduction; not an EM prediction or acceptance',
                rows=sorted(rows,key=lambda r:r['predicted_target_RMS']))
    (base/'short_pad_proposal.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report['rows'][:10],indent=2));return report


def small_refinement(name,observed_case=None):
    """Local measured directions plus isolated-junction slope; fresh EM needed."""
    base=ROOT/'output_data/taper_exploration'/CAMPAIGN
    def load(case):return json.loads((base/case/'response/extraction/fit.json').read_text())
    def vector(q):return np.array([q['C_fF'][0],q['G_fF'][0],q['C_fF'][1],q['G_fF'][1],q['G_fF'][2]])
    ref=load(name)
    m=load('thick64_p200_w152_i16n70_j68_r54'); n=load('thick64_p200_w152_i28n70_j68_r54')
    d=load('t2_c264_top10_w152_pad532_idc34_n72_j112')
    h=load('t2_c264_top10_w152_pad532_idc34_n72_j112_rg62tg6')
    e=load('t2_c264_top10_w152_pad532_idc34_n72_j112_rg64tg4')
    input_slope=(vector(n)-vector(m))/12
    regular_slope=(vector(h)-vector(d))/-4
    ground_slope=(vector(d)-vector(e))/2
    data=np.load(base/name/'response/comparison.npz');f=data['f'];obs=data['observed']
    old_model=transmission(f,ref['L_nH'],ref['C_fF'],ref['G_fF'],ref['delay_ps'])
    cells=[TAPER_2,TAPER_1,UNIT_CELL,UNIT_CELL]
    target=transmission(f,[v.L0_nH for v in cells],[v.C0_fF for v in cells],[v.Cg_input_fF for v in cells],0)
    def uncorrected(din,junction,regular,ground=6):
        values=vector(ref)+input_slope*(din-26)+regular_slope*(regular-ground-60)+ground_slope*(ground-6)
        # Two independent 64-to-6 um electrode tests supply this local scale.
        values[3]*=(106+37.2)/(junction+37.2)
        values[[0,2]]+=.049*(junction-106)
        if np.min(values)<=0:return np.zeros(len(f),complex),values
        cp=[values[0],values[2],*ref['C_fF'][2:]]
        gp=[values[1],values[3],values[4],ref['G_fF'][3]]
        model=transmission(f,ref['L_nH'],cp,gp,ref['delay_ps'])
        return obs*model/old_model,values
    correction=None
    if observed_case is not None:
        geometry=ROOT/'geometry/taper/output'/CAMPAIGN
        source_meta=json.loads((geometry/name/'response/geometry.json').read_text())
        cases=[observed_case] if isinstance(observed_case,str) else list(observed_case)
        steps=[];errors=[]
        for case in cases:
            measured_meta=json.loads((geometry/case/'response/geometry.json').read_text())
            if source_meta['unit']!=measured_meta['unit']:raise ValueError('Correction requires the same body')
            actual=np.load(base/case/'response/comparison.npz')
            if not np.array_equal(f,actual['f']):raise ValueError('Correction frequency grids differ')
            coupling=measured_meta['coupler'];point=np.array([coupling['d'],coupling['d_t1'],coupling.get('d_regular',66),coupling.get('t_regular_ground',6)])
            step=np.array([point[0]-26,point[1]-106,point[2]-point[3]-60,point[3]-6])/np.array([4,10,4,2])
            if step@step==0:raise ValueError('Correction case must differ from the source')
            predicted,_=uncorrected(*point)
            steps.append(step);errors.append(1/actual['observed']-1/predicted)
        correction=np.linalg.pinv(np.array(steps))@np.array(errors)
    def predict(din,junction,regular,ground=6):
        curve,values=uncorrected(din,junction,regular,ground)
        if np.min(values)<=0:return curve,values
        if correction is not None:
            displacement=np.array([din-26,junction-106,regular-ground-60,ground-6])/np.array([4,10,4,2])
            curve=1/(1/curve+displacement@correction)
        return curve,values
    ix=np.arange(0,len(f),5);rows=[];omega=2*np.pi*f[ix]*1e-12
    for ground in (4,6):
        for clearance in (26,28,30,32,34):
            regular=ground+2*clearance
            for din in range(16,45,2):
                for junction in range(72,133,2):
                    curve,values=predict(din,junction,regular,ground)
                    if np.min(values)<=0 or not np.all(np.isfinite(curve)):continue
                    sample=curve[ix]
                    seed=float(np.clip(np.angle(np.vdot(target[ix],sample))/(2*np.pi*np.mean(f[ix]))*1e12,-99,99))
                    objective=lambda tau:float(np.mean(abs(sample*np.exp(-1j*omega*tau)-target[ix])**2))
                    phase=minimize_scalar(objective,bounds=(max(-100,seed-20),min(100,seed+20)),method='bounded')
                    rows.append(dict(input_gap=din,T2_T1_gap=junction,regular_gap=regular,ground_width=ground,
                        predicted_RMS=float(np.sqrt(phase.fun)),phase_ps=float(phase.x),
                        estimated_C2_G2_C1_G1_Gunit=values.tolist()))
    rows.sort(key=lambda q:q['predicted_RMS'])
    for row in rows[:20]:
        curve,_=predict(row['input_gap'],row['T2_T1_gap'],row['regular_gap'],row['ground_width'])
        def residual(x):
            delta=curve*np.exp(-2j*np.pi*f*x[0]*1e-12)-target
            return np.r_[delta.real,delta.imag]
        phase=least_squares(residual,[row['phase_ps']],bounds=([-100],[100]),max_nfev=100)
        row['predicted_RMS']=float(np.sqrt(np.mean(phase.fun**2)*2))
    report=dict(scope='Local parameter proposal only; transferred slopes require actual EM validation',
                source=name,correction_case=observed_case,rows=sorted(rows,key=lambda q:q['predicted_RMS']))
    label='small_refinement' if observed_case is None else 'corrected_small_refinement'
    (base/name/'response'/f'{label}.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report['rows'][:10],indent=2));return report


if __name__ == '__main__':
    main()
