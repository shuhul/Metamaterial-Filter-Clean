"""Correct the local response model using the latest measured EM secant.

Predictions only guide the next geometry; they are never final acceptance.
"""
import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'
os.environ['OMP_NUM_THREADS'] = '1'
import itertools
import json
import numpy as np
from scipy.optimize import differential_evolution, minimize_scalar
from geometry.taper.response_cases import OUT, A, theory
from geometry.taper.refine_physical_response import features
from geometry.taper.tune_wide_response import impedance

CASE = 'response_match_span262_rise20_right558_tab116'


def main():
    previous = np.load(OUT / 'local_tab_response_model.npz')
    f = previous['f']
    z = impedance(CASE, f)
    derivatives = previous['derivatives'].copy()
    delta = np.array([20., -20., 56., 0., 116.])
    scales = np.array([40., 24., 32., 4., 100.])
    direction = delta / scales
    error = z - previous['z'] - np.einsum('p,pfij->fij', delta, derivatives)
    derivatives += (direction / (direction @ direction) / scales)[:, None, None, None] * error
    target = theory(f)[:, 1, 0]
    wanted = features(f, target)

    def predict(x):
        return A['short_internal_z'](z + np.einsum('p,pfij->fij', x, derivatives))[:, 1, 0]

    def score(x):
        measured = features(f, predict(x))
        if measured is None:
            return 1e6
        # Frequency, peak height and power width all contribute. Width and
        # height scales are the actual 2% tolerances, not arbitrary weights.
        residual = (measured - wanted) / np.r_[.001, .001, .02 * wanted[2:]]
        return float(residual @ residual)

    bounds = [(-20, 40), (-16, 16), (-28, 28), (0, 2), (-40, 48)]
    fit = differential_evolution(score, bounds, seed=180, popsize=10,
                                 maxiter=140, workers=1, polish=True, tol=1e-5)
    grid = [2, 2, 2, 1, 4]
    options = [np.unique([step * np.floor(v / step), step * np.ceil(v / step)])
               for v, step in zip(fit.x, grid)]
    rows = []
    for row in itertools.product(*options):
        x = np.array(row)
        s = predict(x)
        feat = features(f, s)
        if feat is None:
            continue
        delay = minimize_scalar(lambda t: np.mean(abs(s - target * np.exp(-2j*np.pi*f*t*1e-12))**2),
                                bounds=(-30, 30), method='bounded')
        rows.append(dict(delta=x.tolist(), score=score(x), features=feat.tolist(),
                         frequency_errors_MHz=((feat[:2]-wanted[:2])*1000).tolist(),
                         height_width_errors_percent=(100*(feat[2:]/wanted[2:]-1)).tolist(),
                         complex_RMSE=float(np.sqrt(delay.fun)), delay_ps=float(delay.x)))
    rows.sort(key=lambda row: row['score'])
    result = dict(base=CASE, measured_features=features(f, A['short_internal_z'](z)[:,1,0]).tolist(),
                  target_features=wanted.tolist(), continuous_delta=fit.x.tolist(),
                  predicted_score=float(fit.fun), candidates=rows[:10], physical_acceptance=False)
    np.savez(OUT / 'refined_tab_response_model.npz', f=f, z=z, derivatives=derivatives)
    (OUT / 'refined_tab_response_direction.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result, indent=2), flush=True)


if __name__ == '__main__':
    main()
