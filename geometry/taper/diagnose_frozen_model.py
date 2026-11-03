"""Offline sensitivity experiment; never used to accept or promote a design.

Compare the requested two-variable model with alternative hypotheses. Releasing
one saved value here diagnoses a failure; it does NOT change the production fit,
the geometry, or the selected first-taper manifest. No EM is launched.
"""
import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'
os.environ['OMP_NUM_THREADS'] = '1'
import json
from pathlib import Path
import numpy as np
from scipy.optimize import least_squares
from geometry.taper.locked_taper2 import ROOT, FIXED, response
from geometry.taper.fitting import chain
from circuit_common import read_snp


def diagnose(name):
    out = ROOT/'output_data/taper_exploration'/name
    trial = json.loads((out/'trial_locked.json').read_text())
    net = read_snp(out/'response.snp')
    mask = (net.f >= 4.5e9) & (net.f <= 5.5e9)
    f, observed = net.f[mask], net.s[mask, 1, 0]
    l2 = trial['fit']['L_t2_nH']
    records = {}
    for hypothesis in ('junction', 'T1_C0', 'regular_C0', 'extra_delay'):
        def predict(x):
            c1, cr, cg21 = FIXED['C_taper_fF'], FIXED['C_regular_fF'], FIXED['Cg_input_fF']
            delay = FIXED['effective_fixture_delay_ps']
            if hypothesis == 'junction': cg21 = x[2]
            if hypothesis == 'T1_C0': c1 = x[2]
            if hypothesis == 'regular_C0': cr = x[2]
            if hypothesis == 'extra_delay': delay += x[2]
            s = chain(2*np.pi*f, [l2, FIXED['L_taper_nH'], FIXED['L_regular_nH']],
                      [x[0], c1, cr], [x[1], cg21, FIXED['Cg_junction_fF'], FIXED['Cg_output_fF']])
            return s[:, 1, 0]*np.exp(-2j*np.pi*f*delay*1e-12)
        def residual(x):
            d = predict(x)-observed
            return np.r_[d.real, d.imag]
        ranges = {'junction': (10, 120), 'T1_C0': (80, 125),
                  'regular_C0': (175, 205), 'extra_delay': (-100, 100)}
        lo, hi = ranges[hypothesis]
        starts = [np.array([c, g, z]) for c in (20, 100, 220)
                  for g in (100, 459, 800) for z in np.linspace(lo, hi, 9)]
        saved = {'junction': FIXED['Cg_input_fF'], 'T1_C0': FIXED['C_taper_fF'],
                 'regular_C0': FIXED['C_regular_fF'], 'extra_delay': 0}[hypothesis]
        starts.extend([np.r_[trial['fit']['x'], saved], np.array([100, 459, saved])])
        starts.sort(key=lambda x: float(residual(x)@residual(x)))
        fits = [least_squares(residual, x, bounds=([0, 30, lo], [400, 1500, hi]),
                              x_scale='jac', max_nfev=250) for x in starts[:8]]
        best = min(fits, key=lambda r: r.fun@r.fun)
        records[hypothesis] = {'x': best.x.tolist(), 'active_bounds': best.active_mask.tolist(),
                               'complex_rmse': float(np.sqrt(np.mean(abs(predict(best.x)-observed)**2))),
                               'max_complex_error': float(np.max(abs(predict(best.x)-observed)))}
    # A separate counterfactual: hold BOTH T2 capacitances exactly at target
    # and see whether a downstream capacitance shift could explain the data.
    # This explicitly violates the frozen assumption, so it cannot be selected.
    def shifted_downstream(x):
        s = chain(2*np.pi*f, [l2, FIXED['L_taper_nH'], FIXED['L_regular_nH']],
                  [100, x[0], x[1]], [459, FIXED['Cg_input_fF'],
                  FIXED['Cg_junction_fF'], FIXED['Cg_output_fF']])
        return s[:, 1, 0]*np.exp(-2j*np.pi*f*FIXED['effective_fixture_delay_ps']*1e-12)
    def shifted_residual(x):
        d = shifted_downstream(x)-observed
        return np.r_[d.real, d.imag]
    grid = [np.array([c1, cr]) for c1 in np.linspace(80, 125, 46)
            for cr in np.linspace(175, 205, 61)]
    grid.sort(key=lambda x: float(shifted_residual(x)@shifted_residual(x)))
    fits = [least_squares(shifted_residual, x, bounds=([80, 175], [125, 205]),
                         x_scale='jac', max_nfev=250) for x in grid[:12]]
    best = min(fits, key=lambda r: r.fun@r.fun)
    records['downstream_caps_with_T2_at_targets'] = {
        'C0_T1_fF': float(best.x[0]), 'C0_regular_fF': float(best.x[1]),
        'active_bounds': best.active_mask.tolist(),
        'complex_rmse': float(np.sqrt(np.mean(abs(shifted_downstream(best.x)-observed)**2))),
        'max_complex_error': float(np.max(abs(shifted_downstream(best.x)-observed)))}
    target = response(f, [100, 459], l2)[:, 1, 0]
    result = {'purpose': 'Offline diagnostic only. Alternative hypotheses are NOT accepted fits.',
              'production_fit_unchanged': True, 'case': name,
              'target_model_complex_rmse': float(np.sqrt(np.mean(abs(target-observed)**2))),
              'hypotheses': records}
    (out/'frozen_model_diagnostic.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result, indent=2), flush=True)
    return result


if __name__ == '__main__':
    from geometry.taper.run_locked_taper2 import CASE_NAME
    diagnose(CASE_NAME)
