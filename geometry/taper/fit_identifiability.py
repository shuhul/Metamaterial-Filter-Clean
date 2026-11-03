"""Local sensitivity and reoptimized profiles for saved complex-S21 fits.

Profiles are deterministic diagnostics, not statistical confidence intervals:
EM/model residuals are not independent measurement noise.
"""
import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'
os.environ['OMP_NUM_THREADS'] = '1'
import json
from pathlib import Path
import numpy as np
from scipy.optimize import least_squares
from geometry.taper.progressive_fitting import transmission


def investigate(folder):
    folder = Path(folder)
    r = json.loads((folder/'fit.json').read_text())
    data = np.load(folder/'response.npz'); f = data['f']; obs = data['actual'][:, 1, 0]
    free = r['free_capacitances']; ls = r['L_nH']; cs = r['C_fF']; gs = r['G_fF']
    values = np.array([(cs if k == 'C' else gs)[i] for k, i in free]+[r['delay_ps']])
    scales = np.r_[np.maximum(values[:-1], 1), 20.]
    lower = np.array([0 if k == 'C' else .1 for k, _ in free]+[-100.])/scales
    upper = np.array([1000 if k == 'C' else 3000 for k, _ in free]+[100.])/scales
    def model(x, ix):
        cp = list(cs); gp = list(gs); v = x*scales
        for (k, i), q in zip(free, v): (cp if k == 'C' else gp)[i] = q
        return transmission(f[ix], ls, cp, gp, v[-1])
    def residual(x, ix):
        delta = model(x, ix)-obs[ix]
        return np.r_[delta.real, delta.imag]
    full = np.arange(len(f)); ix = np.unique(np.linspace(0, len(f)-1, min(2001, len(f))).astype(int))
    x = values/scales
    columns = []
    for j in range(len(x)):
        delta = np.zeros(len(x)); delta[j] = 1e-5
        columns.append((residual(x+delta, full)-residual(x-delta, full))/2e-5)
    jac = np.array(columns).T
    _, sv, directions = np.linalg.svd(jac/np.sqrt(len(f)), full_matrices=False)
    profiles = []
    for j, (kind, index) in enumerate(free):
        others = [q for q in range(len(x)) if q != j]
        rows = [dict(factor=1., fixed_value=float(values[j]), complex_RMS=r['complex_RMS'],
                     other_values=values.tolist(), converged=True)]
        # Track each profile continuously from the optimum. Jumping directly
        # by 10% can lose a 3 MHz resonance and produce a false local minimum.
        for factors in ((.9, .75, .5), (1.1, 1.25, 1.5)):
            previous = x.copy(); previous_factor = 1.
            for factor in factors:
                for step in np.linspace(previous_factor, factor,
                                        int(np.ceil(abs(factor-previous_factor)/.005))+1)[1:]:
                    fixed = float(np.clip(x[j]*step, lower[j]+1e-10, upper[j]-1e-10))
                    def expand(q):
                        v = previous.copy(); v[j] = fixed; v[others] = q
                        return v
                    fitted = least_squares(lambda q: residual(expand(q), ix), previous[others],
                        bounds=(lower[others], upper[others]), max_nfev=100,
                        ftol=1e-9, xtol=1e-9, gtol=1e-9)
                    previous = expand(fitted.x)
                previous_factor = factor
                rows.append(dict(factor=factor, fixed_value=float(previous[j]*scales[j]),
                    complex_RMS=float(np.linalg.norm(residual(previous, full))/np.sqrt(len(f))),
                    other_values=(previous*scales).tolist(), converged=bool(fitted.success)))
        rows.sort(key=lambda q:q['factor'])
        profiles.append(dict(parameter=[kind, index], fitted_value=float(values[j]), rows=rows))
    result = dict(source=str(folder/'fit.json'), fixed_inductances=ls, free_capacitances=free,
        parameter_scales=scales.tolist(), RMS_scaled_singular_values=sv.tolist(),
        scaled_condition=float(sv[0]/sv[-1]), weakest_scaled_direction=directions[-1].tolist(),
        base_complex_RMS=r['complex_RMS'], profiles=profiles,
        scope='Reoptimize every other free value at each fixed-capacitance profile point; no statistical confidence claim.')
    (folder/'identifiability.json').write_text(json.dumps(result, indent=2)+'\n')
    return result


if __name__ == '__main__':
    import sys
    r = investigate(sys.argv[1])
    print('Condition:', r['scaled_condition'])
    for p in r['profiles']:
        print(p['parameter'], p['fitted_value'], [(v['factor'], round(v['complex_RMS'], 6)) for v in p['rows']])
