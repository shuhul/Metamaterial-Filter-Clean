"""Replay completed full-filter metrics and plots from published data; no EM."""
from pathlib import Path
from datetime import datetime
import hashlib
import json
import re
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / 'final_designs/final_full_filter'


def read(path):
    lines = path.read_text().splitlines()
    assert '# GHZ S MA R 50.00000' in lines
    return np.array([[float(v) for v in line.split()] for line in lines
                     if line.strip() and not line.lstrip().startswith(('!', '#'))])


def main():
    direct = read(BASE / 'response_direct.snp')
    dense = read(BASE / 'response_interpolated.snp')
    f, m = dense[:, 0], dense[:, 3]
    assert direct.shape == (123, 9) and dense.shape == (2001, 9)
    assert np.allclose(f[[0, -1]], [4.5, 5.5]) and np.all(np.diff(f) > 0)
    log = (BASE / 'evidence/native_completion.txt').read_text()
    export = (BASE / 'evidence/interpolation_export.txt').read_text(encoding='utf-16')
    assert 'Em errors detected: 0' in log and 'Em simulation completed' in log
    assert 'Total time for 0 frequencies:' in export and 'Em errors detected: 0' in export
    rows = []
    for line in log.splitlines():
        try:
            values = list(map(float, line.split()))
        except ValueError:
            continue
        if len(values) == 9 and 4.5 <= values[0] <= 5.5:
            rows.append(values)
    rows = np.array(sorted(rows))
    assert len(rows) == 123 and np.allclose(rows[:, 0], direct[:, 0])
    log_error = float(np.max(abs(rows[:, 3] - direct[:, 3])))
    assert log_error < 1e-6
    idx = np.flatnonzero(np.diff((m >= .8).astype(int)))
    edges = [float(f[i]+(.8-m[i])*(f[i+1]-f[i])/(m[i+1]-m[i])) for i in idx]
    intervals = list(zip(edges[::2], edges[1::2]))
    lo, hi = max(intervals, key=lambda pair: pair[1]-pair[0])
    conv = re.findall(r'ABS: SNR = ([\d.]+) dB \(([\d.]+) dB Goal\) for (\S+) at ([\d.]+) GHz', log)
    resources = json.loads((BASE / 'evidence/resources_summary.json').read_text())
    # Native completion timestamp is Pacific standard time on this run date.
    started = datetime.fromisoformat(resources['started_utc'])
    completed = datetime.fromisoformat('2026-11-10T04:08:02-08:00')
    metrics = dict(status='native_EM_complete_and_native_interpolation_export_verified',
                   direct_solve_count=123, interpolated_samples=2001,
                   frequency_GHz=[float(f[0]), float(f[-1])], export_step_MHz=.5,
                   native_completion_Pacific=completed.isoformat(),
                   start_utc=started.isoformat(), elapsed_to_native_completion_seconds=(completed-started).total_seconds(),
                   resource_runner_elapsed_seconds=resources['elapsed_seconds'],
                   native_errors=0, native_warnings=0, final_ABS_SNR_dB=float(conv[-1][0]),
                   ABS_goal_dB=float(conv[-1][1]), interpolation_added_EM_solves=0,
                   thread_count=36, subsections=41299, native_peak_memory_GB=14,
                   amplitude_cutoff=.8, cutoff_dB=float(20*np.log10(.8)),
                   main_passband_GHz=[lo, hi], main_passband_width_MHz=(hi-lo)*1000,
                   target_passband_GHz=[4.8, 5.2],
                   crossing_method='Linear crossing estimates on native Sonnet ABS 0.5 MHz export; not independent dense solves',
                   direct_export_log_magnitude_max_difference=log_error,
                   limitations=['Main contiguous band only; no claim entire 4.8-5.2 GHz meets cutoff',
                                'ABS convergence is interpolation evidence, not mesh or fabrication convergence',
                                'One narrow above-cutoff peak has unresolved width on this export grid'],
                   hashes={n: hashlib.sha256((BASE/n).read_bytes()).hexdigest()
                           for n in ['project.sonx', 'response_direct.snp', 'response_interpolated.snp']})
    (BASE / 'metrics.json').write_text(json.dumps(metrics, indent=2)+'\n', encoding='utf-8')
    figures = BASE / 'figures'
    fig, ax = plt.subplots(figsize=(11, 5.8))
    fig.subplots_adjust(left=.09, right=.97, bottom=.22, top=.88)
    ax.axvspan(lo, hi, color='#319566', alpha=.13, label='Main achieved passband')
    ax.plot(f, m, color='#555555', lw=1.8, label='Simulated |S21|')
    for j, x in enumerate([4.8, 5.2]):
        ax.axvline(x, color='#2865ad', ls='--', lw=1.5, label='Target edges' if j == 0 else None)
        ax.text(x, .035, f'Target {x:.1f} GHz', rotation=90, ha='right', va='bottom', color='#2865ad', fontsize=10)
    for j, x in enumerate([lo, hi]):
        ax.axvline(x, color='#c47816', ls=':', lw=1.5, label='Main passband edges' if j == 0 else None)
    ax.axhline(.8, color='#319566', ls='-.', lw=1.2)
    ax.text(4.705, .815, 'Cutoff: |S21| = 0.8', color='#21744b', fontsize=11)
    ax.annotate('', xy=(lo, 1.075), xytext=(hi, 1.075), arrowprops={'arrowstyle': '<->', 'color': '#21744b'})
    ax.text((lo+hi)/2, 1.095, f'Main achieved band: {lo:.4f} to {hi:.4f} GHz ({(hi-lo)*1000:.1f} MHz)',
            ha='center', va='bottom', fontsize=11, color='#21744b')
    ax.set(xlim=(4.7, 5.3), ylim=(0, 1.19), xlabel='Frequency (GHz)', ylabel='|S21|',
           title='Full-filter simulation: target and main achieved passband')
    ax.set_yticks(np.arange(0, 1.01, .2)); ax.grid(axis='y', alpha=.15)
    ax.legend(loc='upper center', bbox_to_anchor=(.5, -.15), ncol=4, frameon=False, fontsize=10)
    for ext in ['png', 'pdf']:
        fig.savefig(figures/f'annotated_main_passband_linear.{ext}', dpi=180)
    plt.close(fig)
    fig, ax = plt.subplots(figsize=(10, 4))
    # The logs do not specify a fixed solve-index offset for each SNR report.
    ax.plot(np.arange(1, len(conv)+1), [float(v[0]) for v in conv], 'o-', ms=3)
    ax.axhline(40, color='#c73731', ls='--', label='40 dB goal')
    ax.set(xlabel='ABS SNR check number', ylabel='Reported ABS SNR (dB)',
           title='Completed sweep: interpolation convergence (not stopband attenuation)')
    ax.legend(); ax.grid(alpha=.2); fig.tight_layout()
    fig.savefig(figures/'completed_convergence.png', dpi=160); plt.close(fig)
    print(json.dumps(metrics, indent=2))


if __name__ == '__main__':
    main()
