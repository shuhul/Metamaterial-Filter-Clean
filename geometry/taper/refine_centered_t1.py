"""Local geometry guide using EM fits and isolated capacitor sensitivities.

Predictions never accept a design: every proposed layout needs its own EM run.
"""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
import json
import itertools
import numpy as np
from scipy.optimize import least_squares,minimize_scalar
from scipy.signal import find_peaks,peak_widths
from geometry.taper.three_tapers_v3 import ROOT,CAMPAIGN
from geometry.taper.measure_compact_inductance import AP,W
from geometry.taper.progressive_fitting import transmission
from geometry.specifications import TAPER_1,UNIT_CELL,specification

A='plain_w220_p544_g52_j72'
B='plain_w220_p558_g52_j84'
BRACKET='plain_w220_p560_g40_j76'
OUT=ROOT/'output_data/taper_exploration'/CAMPAIGN


def isolated(name):
    project=ROOT/'sonnet_files/taper_exploration'/CAMPAIGN/name/'Lcal/generated_ncell_1.sonx'
    n=W['read_network'](project.with_suffix('.snp'),4)
    s,z=AP['reduce_pairs'](n.s);_,c,error=AP['extract_kc'](n.f,AP['s_to_y'](s,z))
    return np.array([c[2].sum(),-c[2,0],-c[2,1]])


def peak_diagnostics(name):
    """Record peak positions/heights and half-power widths from the dense data."""
    data=np.load(OUT/name/'response/extraction/response.npz')
    f=data['f']
    def describe(s):
        power=abs(s)**2
        indices,_=find_peaks(power,prominence=.01)
        _,_,left,right=peak_widths(power,indices,rel_height=.5)
        return [{'frequency_GHz':float(f[i]/1e9),'magnitude':float(abs(s[i])),
                 'half_power_prominence_width_MHz':float((b-a)*(f[1]-f[0])/1e6)}
                for i,a,b in zip(indices,left,right)]
    result={'scope':'Dense exported response; widths use half power prominence',
            'Sonnet':describe(data['actual'][:,1,0]),
            'fit':describe(data['model_S21']),'fixed_target':describe(data['target_S21'])}
    (OUT/name/'response/peak_diagnostics.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2),flush=True)
    return result


def isolation_dispersion(names):
    """Inspect terminal coupling dispersion, not separately certified Cg values."""
    rows=[]
    for name in names:
        path=ROOT/'sonnet_files/taper_exploration'/CAMPAIGN/name/'Lcal/generated_ncell_1.snp'
        n=W['read_network'](path,4)
        s,z=AP['reduce_pairs'](n.s);y=AP['s_to_y'](s,z)
        c=-y[:,2,:2].imag/(2*np.pi*n.f[:,None])*1e15
        rows.append({'case':name,'frequency_GHz':(n.f/1e9).tolist(),
                     'input_and_dummy_effective_fF':c.tolist()})
    result={'scope':'Raw -Im(Y_node,terminal)/omega; includes distributed effects',
            'cases':rows}
    (OUT/'isolated_coupling_dispersion.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2),flush=True)
    return result


def reference_diagnostic(name):
    """Compare fixed regular references without changing the acceptance model."""
    data=np.load(OUT/name/'response/extraction/response.npz')
    old=json.loads((OUT/name/'response/extraction/fit.json').read_text())
    f=data['f'][::3];obs=data['actual'][::3,1,0]
    rows=[]
    for c,g in ((188.8781163861302,8.241620220231225),(190.4,8.2865)):
        def response(x):
            return transmission(f,[old['L_T1_nH'],4.922723945068658,4.922723945068658],
                                [x[0],c,c],[x[1],x[2],g],x[3])
        def error(x):
            delta=response(x)-obs
            return np.r_[delta.real,delta.imag]
        seed=np.r_[old['C_fF'][0],old['G_fF'],old['delay_ps']]
        fit=least_squares(error,seed,x_scale=[150,15,8,10],
                          bounds=([20,.1,.1,-100],[400,100,30,100]))
        rows.append({'fixed_regular_C0_fF':c,'fixed_output_Cg_fF':g,
                     'free_C0_Ginput_Gunit_delay':fit.x.tolist(),
                     'sampled_fit_RMS':float(np.sqrt(np.mean(abs(response(fit.x)-obs)**2)))})
    result={'scope':'Diagnostic only; acceptance and production fit are unchanged',
            'case':name,'reference_choices':rows}
    (OUT/name/'response/reference_diagnostic.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2),flush=True)
    return result


def fine_length_guide():
    """Use same-input pad data and same-junction input data before trimming legs."""
    b='plain_cw344_w220_p530_g54_j70_c596_tg6'
    c='plain_cw344_w220_p518_g56_j66_c596_tg6'
    d='plain_cw344_w220_p528_g56_j70_c596_tg6'
    short='plain_cw344_w220_p518_g56_j66_c584_tg6'
    read=lambda name:json.loads((OUT/name/'response/extraction/fit.json').read_text())
    fb,fc,fd=map(read,(b,c,d))
    values=lambda q:np.r_[q['C_fF'][0],q['G_fF']]
    vb,vc,vd=map(values,(fb,fc,fd))
    pad=(isolated(d)-isolated(c))/10 # Same input gap and dummy-output geometry.
    input_gap=(vd-vb+2*pad)/2       # Same real junction, pad change removed.
    junction=(vc-vb+12*pad-2*input_gap)/-4
    cap=(isolated(short)-isolated(c))/-12
    cap[2]*=vc[2]/isolated(c)[2] # Scale dummy coupling to the real junction.
    jac=np.stack([pad,input_gap,junction,cap],axis=1)
    read_l=lambda name:json.loads((OUT/name/'Lcal/results.json').read_text())['L0_nH']
    lp=(read_l(d)-read_l(c))/10
    li=(read_l(d)-read_l(b)+2*lp)/2
    lc=(read_l(short)-read_l(c))/-12
    data=np.load(OUT/c/'response/extraction/response.npz')
    f=data['f'][::3];residue=(data['actual'][:,1,0]-data['model_S21'])[::3]
    target=transmission(f,[TAPER_1.L0_nH,UNIT_CELL.L0_nH,UNIT_CELL.L0_nH],
                        [TAPER_1.C0_fF,UNIT_CELL.C0_fF,UNIT_CELL.C0_fF],
                        [TAPER_1.Cg_input_fF,UNIT_CELL.Cg_input_fF,UNIT_CELL.Cg_input_fF],0)
    tau=json.loads((OUT/c/'response/target_comparison.json').read_text())['fixture_delay_ps']
    trials=[]
    for length,gap,cap_l in itertools.product(range(500,566,2),range(50,64,2),range(584,598,2)):
        delta=np.array([length-518,gap-56,0,cap_l-596])
        caps=vc+jac@delta
        l=fc['L_T1_nH']+np.array([lp,li,0,lc])@delta
        prediction=transmission(f,[l,4.922723945068658,4.922723945068658],
                               [caps[0],188.8781163861302,188.8781163861302],
                               [caps[1],caps[2],8.241620220231225],fc['delay_ps'])+residue
        objective=lambda delay:float(np.mean(abs(prediction-target*np.exp(-2j*np.pi*f*delay*1e-12))**2))
        fitted=minimize_scalar(objective,bounds=(tau-20,tau+20),method='bounded')
        trials.append({'pad_l_um':length,'d_input_um':gap,'d_junction_um':66,'cap_l_um':cap_l,
                       'predicted_C0_Ginput_Gunit_fF':caps.tolist(),'predicted_L_nH':float(l),
                       'predicted_target_RMS':float(np.sqrt(fitted.fun))})
    trials.sort(key=lambda q:q['predicted_target_RMS'])
    result={'scope':'Prediction only, full EM required; measured reference residual retained',
            'reference':c,'geometry_columns':['pad_length','input_gap','regular_gap','cap_length'],
            'Jacobian_fF_per_um':jac.tolist(),'capacitance_rows':['C0','input_Cg','regular_Cg'],
            'sources':{'pad_calibration':[c,d],'input_full_responses':[b,d],
                       'junction_full_responses':[b,c],'cap_leg_calibrations':[c,short]},
            'assumption':'Scaled isolated cap-leg sensitivity approximates the assembled change',
            'ranked_quantized_proposals':trials[:12]}
    (OUT/'fine_length_guide.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2),flush=True)
    return result


def guide():
    ca=isolated(A);cb=isolated(B);cc=isolated(BRACKET)
    fa=json.loads((OUT/A/'response/extraction/fit.json').read_text())
    fb=json.loads((OUT/B/'response/extraction/fit.json').read_text())
    va=np.r_[fa['C_fF'][0],fa['G_fF']];vb=np.r_[fb['C_fF'][0],fb['G_fF']]
    # A and B have identical input geometry. Their isolation differs only
    # in pad length; the full response additionally changes the junction gap.
    pad=(cb-ca)/14
    junction=(vb-va-pad*14)/12
    input_gap=(ca-(cc-pad*16))/12
    input_gap[2]=0  # Distant regular junction; this term is not measured here.
    jac=np.stack([pad,input_gap,junction],axis=1)
    f=np.linspace(4.65e9,5.4e9,2501)
    target=transmission(f,[TAPER_1.L0_nH,UNIT_CELL.L0_nH,UNIT_CELL.L0_nH],
                        [TAPER_1.C0_fF,UNIT_CELL.C0_fF,UNIT_CELL.C0_fF],
                        [TAPER_1.Cg_input_fF,UNIT_CELL.Cg_input_fF,UNIT_CELL.Cg_input_fF],0)
    l=fb['L_T1_nH']
    def model(v,tau=0):
        return transmission(f,[l,4.922723945068658,4.922723945068658],
                            [v[0],188.8781163861302,188.8781163861302],
                            [v[1],v[2],8.241620220231225],tau)
    def residual(x):
        d=model(x[:3],x[3])-target
        return np.r_[d.real,d.imag]
    optimal=least_squares(residual,np.r_[vb,0.],bounds=([80,1,1,-50],[220,50,20,50]),
                          x_scale=[150,15,8,20],max_nfev=400)
    delta=np.linalg.lstsq(jac,optimal.x[:3]-vb,rcond=None)[0]
    current=np.array([558.,52.,84.]);continuous=current+delta
    steps=np.array([2,2,4]);nearest=np.round(continuous/steps)*steps
    trials=[]
    for offset in itertools.product((-1,0,1),repeat=3):
        geom=nearest+np.array(offset)*steps
        if not (20<=geom[0]<596 and 8<=geom[1]<=180 and 8<=geom[2]<=160):continue
        caps=vb+jac@(geom-current)
        if np.any(caps<=0):continue
        prediction=model(caps)
        objective=lambda tau:float(np.mean(abs(prediction*np.exp(-2j*np.pi*f*tau*1e-12)-target)**2))
        fit=minimize_scalar(objective,bounds=(-20,20),method='bounded')
        trials.append({'pad_l_um':int(geom[0]),'d_input_um':int(geom[1]),'d_junction_um':int(geom[2]),
                       'predicted_C0_Ginput_Gunit_fF':caps.tolist(),'predicted_target_RMS':float(np.sqrt(fit.fun))})
    trials.sort(key=lambda q:q['predicted_target_RMS'])
    result={'scope':'Local first-order geometry prediction, not EM acceptance',
            'specification':specification(),'references':[A,B,BRACKET],
            'reference_source_hashes':{name:W['digest'](OUT/name/'Lcal/results.json') for name in (A,B,BRACKET)},
            'geometry_columns':['pad_length','input_gap','regular_gap'],
            'capacitance_rows':['C0','input_Cg','regular_Cg'],
            'Jacobian_fF_per_um':jac.tolist(),'electrical_optimum_C0_Ginput_Gunit':optimal.x[:3].tolist(),
            'continuous_geometry_proposal':continuous.tolist(),'rank':int(np.linalg.matrix_rank(jac)),
            'assumptions':['Isolation C0 and dummy-output-Cg derivatives approximate assembled derivatives',
                           'Input-gap effect on the distant regular coupling is negligible',
                           'L held at the reference measured value; proposed geometry must be recalibrated'],
            'ranked_quantized_proposals':trials[:8]}
    (OUT/'local_geometry_guide.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2),flush=True)
    return result


def residual_guide(assembled=False,reference_case=None,geometry_origin=None,vary_cap_length=True):
    """Update input sensitivity from full EM and retain its measured residual."""
    name='plain_w220_p564_g54_j84'
    def read(case):return json.loads((OUT/case/'response/extraction/fit.json').read_text())
    fa,fb,fc=map(read,(A,B,name))
    values=lambda r:np.r_[r['C_fF'][0],r['G_fF']]
    va,vb,vc=map(values,(fa,fb,fc))
    pad=(isolated(B)-isolated(A))/14
    junction=(vb-va-pad*14)/12
    input_gap=(vc-vb-pad*6)/2
    cap_case='plain_w220_p564_g54_j84_c596'
    cap_length=(isolated(cap_case)-isolated(name))/-4
    jac=np.stack([pad,input_gap,junction,cap_length],axis=1)
    lp=(fb['L_T1_nH']-fa['L_T1_nH'])/14
    li=(fc['L_T1_nH']-fb['L_T1_nH']-6*lp)/2
    lj=0. # Isolation output geometry is fixed, so junction is absent there.
    lc=(json.loads((OUT/cap_case/'Lcal/results.json').read_text())['L0_nH']-fc['L_T1_nH'])/-4
    origin=np.array([564,54,84,600])
    if assembled:
        # The new full trial changes pad +6, input -2 and cap length -4.
        # Subtract the already measured directions to isolate its cap effect.
        next_name='plain_w220_p570_g52_j84_c596'
        fd=read(next_name);vd=values(fd)
        jac[:,3]=(vd-vc-6*pad+2*input_gap)/-4
        lc=(fd['L_T1_nH']-fc['L_T1_nH']-6*lp+2*li)/-4
        name=next_name;fc=fd;vc=vd
        origin=np.array([570,52,84,596])
    if reference_case is not None:
        if geometry_origin is None:raise ValueError('Provide the measured reference geometry')
        name=reference_case;fc=read(name);vc=values(fc)
        origin=np.asarray(geometry_origin)
        if reference_case.endswith('_tg6'):
            ga=read('plain_cw344_w220_p530_g54_j70_c596_tg6')
            gb=read('plain_cw344_w220_p518_g56_j66_c596_tg6')
            # Update the gap direction with the two measured 6 um-strip runs.
            # Their pad/input shifts are -12/+2 um and their gap shift is -4.
            jac[:,2]=(values(gb)-values(ga)+12*pad-2*input_gap)/-4
    data=np.load(OUT/name/'response/extraction/response.npz')
    f=data['f'][::3];residual=(data['actual'][:,1,0]-data['model_S21'])[::3]
    target=transmission(f,[TAPER_1.L0_nH,UNIT_CELL.L0_nH,UNIT_CELL.L0_nH],
                        [TAPER_1.C0_fF,UNIT_CELL.C0_fF,UNIT_CELL.C0_fF],
                        [TAPER_1.Cg_input_fF,UNIT_CELL.Cg_input_fF,UNIT_CELL.Cg_input_fF],0)
    old_score=json.loads((OUT/name/'response/target_comparison.json').read_text())
    tau=old_score['fixture_delay_ps'];trials=[]
    ranges=[range(int(origin[0]-20),int(origin[0]+22),2),
            range(int(origin[1]-6),int(origin[1]+8),2),
            range(int(origin[2]-4),int(origin[2]+8),4),
            range(max(560,int(origin[3]-6)),min(604,int(origin[3]+8)),2)]
    if not vary_cap_length:ranges[3]=[int(origin[3])]
    for length,gap,jgap,cap_l in itertools.product(*ranges):
        if not (20<=length<596 and gap>=8 and jgap>=8):continue
        delta=np.array([length,gap,jgap,cap_l])-origin
        caps=vc+jac@delta;l=fc['L_T1_nH']+np.array([lp,li,lj,lc])@delta
        if np.any(caps<=0):continue
        prediction=transmission(f,[l,4.922723945068658,4.922723945068658],
                               [caps[0],188.8781163861302,188.8781163861302],
                               [caps[1],caps[2],8.241620220231225],fc['delay_ps'])+residual
        objective=lambda delay:float(np.mean(abs(prediction-target*np.exp(-2j*np.pi*f*delay*1e-12))**2))
        fitted=minimize_scalar(objective,bounds=(tau-20,tau+20),method='bounded')
        trials.append({'pad_l_um':length,'d_input_um':gap,'d_junction_um':jgap,'cap_l_um':cap_l,
                       'predicted_C0_Ginput_Gunit_fF':caps.tolist(),'predicted_L_nH':float(l),
                       'predicted_target_RMS':float(np.sqrt(fitted.fun))})
    trials.sort(key=lambda q:q['predicted_target_RMS'])
    result={'scope':'Local prediction, not EM acceptance; measured residual retained',
            'reference':name,'Jacobian_fF_per_um':jac.tolist(),
            'reference_geometry':origin.tolist(),
            'geometry_columns':['pad_length','input_gap','regular_gap','cap_length'],
            'cap_length_sensitivity_source':cap_case,
            'assembled_cap_length_sensitivity':assembled,
            'assumption':('Cap-leg sensitivity inferred from combined full EM perturbation after subtracting other measured directions'
                          if assembled else 'Isolated cap-leg sensitivity approximates the assembled sensitivity; acceptance requires full EM'),
            'capacitance_rows':['C0','input_Cg','regular_Cg'],
            'full_EM_input_cross_sensitivity_included':True,
            'source_fit_sha256':{case:W['digest'](OUT/case/'response/extraction/fit.json')
                                 for case in {A,B,'plain_w220_p564_g54_j84',name,
                                              *(['plain_w220_p570_g52_j84_c596'] if assembled else [])}},
            'ranked_quantized_proposals':trials[:10]}
    (OUT/'local_residual_geometry_guide.json').write_text(json.dumps(result,indent=2)+'\n')
    (OUT/name/'response/local_geometry_guide.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2),flush=True)
    return result


if __name__=='__main__':guide()
