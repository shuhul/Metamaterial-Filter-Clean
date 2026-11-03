"""Offline sensitivity check; never runs Sonnet or changes the primary fit.

Compare weighting and bandwidth because an accurate-looking transmission
curve alone does not guarantee unique capacitances in this strongly loaded IDC.
"""
from pathlib import Path
import json
import numpy as np
from scipy.optimize import least_squares
from circuit_common import read_snp
from geometry.taper.second_experiment import model_response


def check(output):
    output = Path(output)
    primary = json.loads((output / 'fit_network.json').read_text())
    net = read_snp(output / 'response.snp')
    lower = [20, 60, 150, 100, 10, 2, 2, 0, 0]
    upper = [300, 160, 230, 1600, 140, 25, 25, 20, 20]
    rows = []
    for stop, weight, step in ((9, 1, 2), (9, 3, 2), (9, 10, 2), (6, 3, 2), (9, 3, 1)):
        ix = np.flatnonzero(net.f <= stop * 1e9)[::step]
        f, obs = net.f[ix], net.s[ix]
        def residual(x):
            pred = model_response(f, x, primary['L_t2_nH'], l1=primary['L_t1_nH'])
            a = weight * (pred[:, 1, 0] - obs[:, 1, 0])
            b = pred[:, 0, 0] - obs[:, 0, 0]
            return np.r_[a.real, a.imag, b.real, b.imag]
        fit = least_squares(residual, primary['x'], bounds=(lower, upper),
                            x_scale='jac', max_nfev=180, ftol=1e-10)
        pred = model_response(f, fit.x, primary['L_t2_nH'], l1=primary['L_t1_nH'])
        row = dict(stop_GHz=stop, s21_weight=weight, sample_step=step, x=fit.x.tolist(),
                   success=bool(fit.success), active_bounds=fit.active_mask.tolist(),
                   s21_complex_rmse=float(np.sqrt(np.mean(abs(pred[:, 1, 0]-obs[:, 1, 0])**2))),
                   s11_complex_rmse=float(np.sqrt(np.mean(abs(pred[:, 0, 0]-obs[:, 0, 0])**2))))
        rows.append(row)
        print(json.dumps(row), flush=True)
    result = {'purpose': 'Sensitivity diagnostic, not a replacement target extraction', 'fits': rows}
    (output / 'fit_stability.json').write_text(json.dumps(result, indent=2) + '\n')
    return result
