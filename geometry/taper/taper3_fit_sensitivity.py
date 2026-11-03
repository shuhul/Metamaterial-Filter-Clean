"""Local identifiability diagnostics on saved T3 fits; no Sonnet execution."""
import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'
os.environ['OMP_NUM_THREADS'] = '1'
import json
import numpy as np
from scipy.optimize import least_squares
from geometry.taper.taper3 import OUT
from geometry.taper.progressive_fitting import transmission


def analyze(case, *, stage='response', fit_folder='fit'):
    folder = OUT/case/stage/fit_folder
    fit = json.loads((folder/'fit.json').read_text())
    data = np.load(folder/'response.npz')
    indices = np.unique(np.linspace(0, len(data['f'])-1, 2001).astype(int))
    f, actual = data['f'][indices], data['actual'][indices, 1, 0]
    declared = fit['provenance']['free_capacitances']
    new_free = [(k,i) for k,i in declared if (k,i) in (('C',0),('G',0),('G',1))]
    names = {('C',0):'C0_T3_fF',('G',0):'input_Cg_fF',('G',1):'T3_to_T2_Cg_fF'}
    labels = [names[q] for q in new_free]+['delay_ps']
    x = np.array([fit[k+'_fF'][i] for k,i in new_free]+[fit['delay_ps']])
    count = len(x)
    scales = np.array([max(v,8.) for v in x[:-1]]+[20.])
    lower = np.array([0. if k=='C' else .1 for k,i in new_free]+[-100.])
    upper = np.array([1000. if k=='C' else 3000. for k,i in new_free]+[100.])

    def residual(values):
        cs, gs = fit['C_fF'].copy(), fit['G_fF'].copy()
        for (kind,index),value in zip(new_free,values[:-1]):
            (cs if kind=='C' else gs)[index] = value
        delta = transmission(f,fit['L_nH'],cs,gs,values[-1])-actual
        return np.r_[delta.real,delta.imag]

    step = np.eye(count)*scales*1e-4
    jac = np.column_stack([(residual(x+q)-residual(x-q))/(2e-4) for q in step])
    _, singular, directions = np.linalg.svd(jac, full_matrices=False)
    covariance = np.linalg.pinv(jac.T@jac)
    correlation = covariance/np.sqrt(np.outer(np.diag(covariance), np.diag(covariance)))
    probes = []
    for index in range(count-1):
        other = [j for j in range(count) if j != index]
        for factor in (.8, .95, 1.05, 1.2):
            held = x[index]*factor
            if not lower[index] < held < upper[index]:
                continue

            def expand(values):
                trial = x.copy(); trial[index] = held; trial[other] = values
                return trial

            result = least_squares(lambda v: residual(expand(v)), x[other],
                bounds=(lower[other], upper[other]), x_scale=scales[other],
                max_nfev=400, ftol=1e-11, xtol=1e-11, gtol=1e-11)
            trial = expand(result.x)
            probes.append(dict(held_parameter=index, factor=factor, values=trial.tolist(),
                complex_RMS=float(np.sqrt(np.mean(abs(
                    residual(trial)[:len(f)]+1j*residual(trial)[len(f):])**2))),
                success=bool(result.success)))
    report = dict(case=case, stage=stage, fit_folder=fit_folder,
        parameters=labels,
        values=x.tolist(), scales=scales.tolist(), sample_count=len(f),
        baseline_complex_RMS=float(np.linalg.norm(residual(x))/np.sqrt(len(f))),
        scaled_singular_values=singular.tolist(),
        scaled_condition_number=float(singular[0]/singular[-1]),
        linearized_parameter_correlation=correlation.tolist(),
        weakest_scaled_direction=directions[-1].tolist(), local_profile_probes=probes,
        scope='Local numerical sensitivity with older values and all L fixed. '
              'Not statistical confidence, global uniqueness, or an independent physical extraction. '
              'Profiles are local reoptimizations; residual model error is not random noise.')
    (folder/'local_sensitivity.json').write_text(json.dumps(report, indent=2)+'\n')
    return report


if __name__ == '__main__':
    from geometry.taper.run_taper3 import CASE
    print(json.dumps(analyze(CASE), indent=2))
