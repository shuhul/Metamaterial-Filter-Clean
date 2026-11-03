r"""Refit cached plain-wire data: free taper C's, fixed L's and regular cells.

Run with: .venv\Scripts\python.exe -B -m geometry.taper.refit_plain_t2
No geometry changes or Sonnet simulations are performed.
"""
import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'
os.environ['OMP_NUM_THREADS'] = '1'
import json
from pathlib import Path
import numpy as np
from scipy.optimize import least_squares
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from circuit_common import read_snp
from geometry.taper.staged_fitting import circuit
from geometry.taper.run_bridged_stage import digest


def run():
    root = Path(__file__).resolve().parents[2]
    case = root / 'output_data/taper_exploration/compact_centered/w500_top10_wire52_gap6/direct'
    old = json.loads((case / 'fit/fit.json').read_text())
    source = Path(old['source'])
    net = read_snp(source)
    f, actual = net.f, net.s[:, 1, 0]
    fixed = old['fixed']
    ls = [old['L_T2_nH'], fixed['L_T1_nH'], fixed['L_regular_nH'], fixed['L_regular_nH']]
    names = ['C0_T1_fF', 'Cg_T2_T1_fF', 'Cg_T1_regular_fF', 'fixture_delay_ps']
    base = np.array([fixed['C0_T1_fF'], old['Cg_T2_T1_fF'], fixed['Cg_T1_regular_fF'], old['fixture_delay_ps']])
    lower, upper = np.array([1, 1, .1, -100]), np.array([1000, 1000, 100, 100])
    scales = [100, 57, 8, 10]

    def predict(x, frequency):
        caps = [0, x[0], fixed['C0_regular_fF'], fixed['C0_regular_fF']]
        gaps = [None, x[1], x[2], fixed['Cg_output_fF']]
        return circuit(frequency, ls, caps, gaps)[:, 1, 0] * np.exp(-2j*np.pi*frequency*x[3]*1e-12)

    # Retain every sample for final optimization; multistarts use a smaller grid.
    ix = np.unique(np.linspace(0, len(f)-1, min(len(f), 1801)).astype(int))
    def residual(x, indices, inverse=False):
        obs = actual[indices]
        model = predict(x, f[indices])
        delta = ((1/model - 1/obs)*np.maximum(abs(obs), .03)**2
                 if inverse else model - obs)
        return np.r_[delta.real, delta.imag]

    rng = np.random.default_rng(20261023)
    seeds = [base]
    for _ in range(11):
        seed = base.copy()
        seed[:3] *= np.exp(rng.normal(0, .4, 3))
        seed[3] += rng.normal(0, 8)
        seeds.append(np.clip(seed, lower+1e-7, upper-1e-7))
    candidates = []
    for seed in seeds:
        pre = least_squares(residual, seed, args=(ix, True), bounds=(lower, upper), x_scale=scales, max_nfev=300)
        for start in (seed, pre.x):
            fit = least_squares(residual, start, args=(ix,), bounds=(lower, upper), x_scale=scales,
                                max_nfev=500, ftol=1e-11, xtol=1e-11, gtol=1e-11)
            candidates.append(fit)
    best = min(candidates, key=lambda q: q.fun @ q.fun)
    best = least_squares(residual, best.x, args=(np.arange(len(f)),), bounds=(lower, upper),
                         x_scale=scales, max_nfev=800, ftol=1e-12, xtol=1e-12, gtol=1e-12)
    model = predict(best.x, f)
    delta = model - actual
    out = case / 'fit_free_taper_capacitances'
    out.mkdir(parents=True, exist_ok=True)
    report = dict(source=str(source), source_sha256=digest(source), code_sha256=digest(__file__),
                  free_parameters=dict(zip(names, map(float, best.x))),
                  fixed_inductances_nH=dict(zip(['T2', 'T1', 'regular', 'output_end'], ls)),
                  fixed_capacitances_fF=dict(C0_T2=0, C0_regular=fixed['C0_regular_fF'],
                                             C0_output_end=fixed['C0_regular_fF'], Cg_output=fixed['Cg_output_fF']),
                  objective='Unweighted complex S21 real and imaginary residuals, all original frequency samples',
                  S21_complex_RMSE=float(np.sqrt(np.mean(abs(delta)**2))),
                  S21_max_complex_error=float(abs(delta).max()),
                  S21_magnitude_RMSE=float(np.sqrt(np.mean((abs(model)-abs(actual))**2))),
                  previous_complex_RMSE=old['S21_complex_RMSE'], active_bounds=best.active_mask.tolist(),
                  optimizer_success=bool(best.success), samples=len(f),
                  parameter_extraction_accepted=False,
                  caveat='Transmission agreement alone does not establish unique physical capacitances.',
                  multistart_costs=sorted(float(q.fun@q.fun) for q in candidates))
    (out/'fit.json').write_text(json.dumps(report, indent=2)+'\n')
    np.savez(out/'response.npz', f=f, actual=actual, model=model)
    fig, ax = plt.subplots(figsize=(13, 7.5))
    ax.plot(f/1e9, abs(actual), color='#b22222', lw=2.8, label='Sonnet: straight T2 wire, no IDC')
    ax.plot(f/1e9, abs(model), 'k--', lw=2.3, label='Complex S21 fit: taper capacitances free')
    text = (f'$C_{{0,T1}}$ = {best.x[0]:.3f} fF\n'
            f'$C_{{g,T2-T1}}$ = {best.x[1]:.3f} fF\n'
            f'$C_{{g,T1-unit}}$ = {best.x[2]:.3f} fF\n'
            f'Complex S21 RMS = {report["S21_complex_RMSE"]:.5f}\n'
            'All L and regular-cell parameters fixed\n'
            f'$C_{{0,T2}}$ = 0; delay = {best.x[3]:.2f} ps')
    ax.text(.025, .97, text, transform=ax.transAxes, va='top', fontsize=14,
            bbox=dict(boxstyle='round,pad=.55', fc='white', ec='#999999', alpha=.95))
    ax.set(xlabel='Frequency (GHz)', ylabel=r'$|S_{21}|$', ylim=(0, 1.35), title='Straight-Wire T2: Free Taper Capacitance Fit')
    ax.xaxis.label.set_size(18); ax.yaxis.label.set_size(21); ax.title.set_size(22)
    ax.tick_params(labelsize=15); ax.grid(alpha=.22); ax.legend(fontsize=13, loc='upper right')
    fig.tight_layout(); fig.savefig(out/'s21_fit_linear.png', dpi=190); plt.close(fig)
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    run()
