"""Offline magnetic-coupling sensitivity of the calibrated bare field model.

Defined loop coefficients are design sensitivities, not uniquely extracted
ideal mutual inductances. This diagnostic never releases T1 fit parameters.
"""
import json
import numpy as np
from scipy.optimize import least_squares
from geometry.taper.plain_series import OUT,A,change_caps
from geometry.taper.compact_fit import direct_predictor


def search():
    out=OUT/'w500_top10_wire52_gap6/series_probes';m=np.load(out/'coupled_model.npz')
    f=np.linspace(4.65e9,5.4e9,1501);co=m['coefficients'];powers=m['powers'];cap=m['static_C_fF']
    z=((f[:,None]/5e9)**powers@co).reshape(len(f),6,6)
    obs=direct_predictor(f,.7327695472953403)([100,56.7,5])
    def predict(x):
        zz=change_caps(z,f,cap,x[:2])
        for j,cross in ((3,.10416484567579731),(4,.01911955440805398),(5,.006275668516268254)):
            delta=2j*np.pi*f*(x[2]-1)*cross*1e-9
            zz[:,2,j]+=delta;zz[:,j,2]+=delta
        return A['short_internal_z'](zz)[:,1,0]*np.exp(-2j*np.pi*f*x[3]*1e-12)
    def residual(x,smooth=False):
        p=predict(x);d=(1/p-1/obs)*np.maximum(abs(obs),.02)**2 if smooth else p-obs
        return np.r_[d.real,d.imag]
    lo=[1,20,-1,-20];hi=[2000,160,2,20];scale=[200,70,1,5]
    seeds=[[c,g,q,0] for c in (50,300,600) for g in (60,90) for q in (0,1)]
    smooth=[least_squares(lambda x:residual(x,True),x,bounds=(lo,hi),x_scale=scale,max_nfev=120) for x in seeds]
    fits=[least_squares(residual,q.x,bounds=(lo,hi),x_scale=scale,max_nfev=200)
          for q in sorted(smooth,key=lambda q:np.linalg.norm(residual(q.x)))[:4]]
    q=min(fits,key=lambda q:q.fun@q.fun);p=predict(q.x)
    result=dict(diagnostic_only=True,physical_design_accepted=False,T1_parameters_fitted=False,
        controls=['whole_source_ground_C_fF','whole_side_mutual_C_fF','defined_loop_cross_coefficient_scale','delay_ps'],
        x=q.x.tolist(),complex_S21_RMS=float(np.sqrt(np.mean(abs(p-obs)**2))),max_S21_error=float(abs(p-obs).max()),
        scope='Local field-model sensitivity; coefficient scaling is not an independently validated physical geometry')
    (out/'magnetic_direction.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2),flush=True);return result


if __name__=='__main__':search()
