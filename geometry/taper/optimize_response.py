"""T2 response matching: optimize EM against the specified ideal filter.

No external-S21 fit is used to decide whether the geometry meets the goal.
Internal-port calibration supplies a local search direction, not acceptance.
"""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
from pathlib import Path
import runpy,json
import numpy as np
from scipy.optimize import least_squares
from scipy.signal import find_peaks
from geometry.specifications import UNIT_CELL,TAPER_1,TAPER_2,specification
from geometry.taper.fitting import chain
from geometry.taper.locked_taper2 import FIXED

ROOT=Path(__file__).resolve().parents[2]
INV=ROOT/'investigations/taper2_fit_2026-10-22'
OUT=ROOT/'output_data/taper_exploration/response_optimization_2026-10-22'
OUT.mkdir(parents=True,exist_ok=True)
A=runpy.run_path(str(INV/'analyze_series.py'))
THEORY_VERSION='nominal_with_loaded_output_end_v2'

def append_output_shunt(s,admittance):
    """Cascade a shunt admittance at port 2, normalized to 50 ohms."""
    y=50*np.asarray(admittance);a=-y/(2+y);b=2/(2+y)
    den=1-s[:,1,1]*a;v=s.copy()
    v[:,0,0]=s[:,0,0]+s[:,0,1]*a*s[:,1,0]/den
    v[:,1,1]=a+b*s[:,1,1]*b/den
    v[:,1,0]=s[:,1,0]*b/den;v[:,0,1]=s[:,0,1]*b/den
    return v

def theory(f,nominal=True,*,include_output_end=True):
    if nominal:
        ls=[TAPER_2.L0_nH,TAPER_1.L0_nH,UNIT_CELL.L0_nH]
        cs=[TAPER_2.C0_fF,TAPER_1.C0_fF,UNIT_CELL.C0_fF]
        gs=[TAPER_2.Cg_input_fF,TAPER_1.Cg_input_fF,UNIT_CELL.Cg_input_fF,UNIT_CELL.Cg_input_fF]
    else:
        ls=[.74,FIXED['L_taper_nH'],FIXED['L_regular_nH']]
        cs=[100,FIXED['C_taper_fF'],FIXED['C_regular_fF']]
        gs=[459,FIXED['Cg_input_fF'],FIXED['Cg_junction_fF'],FIXED['Cg_output_fF']]
    s=chain(2*np.pi*np.asarray(f),ls,cs,gs)
    if include_output_end:
        # The final Unit_Cell output is its LC node, directly connected to
        # port 2. Its inductor is damped by the port but remains in the circuit.
        ell=UNIT_CELL.L0_nH if nominal else FIXED['L_regular_nH']
        cap=UNIT_CELL.C0_fF if nominal else FIXED['C_regular_fF']
        admittance=2j*np.pi*f*cap*1e-15+1/(2j*np.pi*f*ell*1e-9)
        s=append_output_shunt(s,admittance)
    # Effective fixture delay already established with the unchanged T1 fixture.
    tl=FIXED['fixture_delay_left_ps']*1e-12;tr=FIXED['fixture_delay_right_ps']*1e-12
    s[:,0,0]*=np.exp(-4j*np.pi*f*tl);s[:,1,1]*=np.exp(-4j*np.pi*f*tr)
    s[:,1,0]*=np.exp(-2j*np.pi*f*(tl+tr));s[:,0,1]=s[:,1,0]
    return s

def baseline(f):
    m=np.load(INV/'series_impedance_model.npz');r=f/5e9
    z=((r[:,None]**m['powers'])@m['coefficients']).reshape(len(f),6,6)
    return z,m['static_C_fF']

def local_response(f,z,c,x):
    # x = ground C of whole T2 island, input mutual C, T2-T1 mutual C,
    # added series inductance at the T2 gap. These are local surrogate controls.
    cc=c.copy();d=x[2]+c[1,2]
    cc[1,1]+=d;cc[2,2]+=d;cc[1,2]-=d;cc[2,1]-=d
    if len(x)>4:cc[0,0]+=x[4]-c[0].sum()
    zz=A['change_two_capacitances'](z,f,cc,x[:2])
    # change_two_capacitances uses cc as its baseline: add the junction change
    # relative to the ORIGINAL matrix explicitly as well.
    b=A['conductor_mapping']()
    delta=b.T@(np.linalg.inv(cc*1e-15)-np.linalg.inv(c*1e-15))@b
    zz+=delta[None,:,:]/(2j*np.pi*f[:,None,None])
    zz[:,2,2]+=2j*np.pi*f*x[3]*1e-9
    return A['short_internal_z'](zz)

def optimize():
    f=np.linspace(4.5e9,5.5e9,1001);z,c=baseline(f);target=theory(f)
    ref=np.array([c[1].sum(),-c[0,1],-c[1,2],0.])
    def residual(x):
        d=local_response(f,z,c,x)[:,1,0]-target[:,1,0]
        return np.r_[d.real,d.imag]
    def denominator_residual(x):
        prediction=local_response(f,z,c,x)[:,1,0]
        d=(1/prediction-1/target[:,1,0])*np.maximum(abs(target[:,1,0]),.02)**2
        return np.r_[d.real,d.imag]
    rng=np.random.default_rng(29)
    seeds=[ref,ref*[1,.7,1,1],ref*[1,1.3,1,1]]
    seeds += [ref+[rng.uniform(-10,10),rng.uniform(-120,120),rng.uniform(-15,15),rng.uniform(-.04,.04)] for _ in range(9)]
    smooth=[least_squares(denominator_residual,x,bounds=([10,30,5,-.65],[500,1500,200,5]),
                        x_scale=[100,459,57,.74],max_nfev=250,ftol=1e-10) for x in seeds[:4]]
    fits=[least_squares(residual,q.x,bounds=([10,30,5,-.65],[500,1500,200,5]),
                        x_scale=[100,459,57,.74],max_nfev=250,ftol=1e-10) for q in smooth]
    best=min(fits,key=lambda q:q.fun@q.fun)
    report={'specification':specification(),'objective':'Direct complex S21 agreement with nominal theoretical circuit; fixed T1 metal',
            'local_surrogate_reference':ref.tolist(),'local_surrogate_best':best.x.tolist(),
            'reference_RMSE':float(np.sqrt(np.mean(abs(residual(ref))**2)*2)),
            'best_RMSE':float(np.sqrt(np.mean(abs(residual(best.x))**2)*2)),
            'physical_EM_acceptance':False,
            'starts':[{'x':q.x.tolist(),'cost':float(q.fun@q.fun)} for q in fits]}
    (OUT/'surrogate_direction.json').write_text(json.dumps(report,indent=2)+'\n')
    np.savez(OUT/'surrogate_direction.npz',f=f,target=target,baseline=A['short_internal_z'](z),predicted=local_response(f,z,c,best.x))
    print(json.dumps(report,indent=2),flush=True)

if __name__=='__main__':optimize()
