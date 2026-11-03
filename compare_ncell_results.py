"""Fit baseline/generated four-cell results using the existing extraction model."""

from pathlib import Path
import json

import matplotlib.pyplot as plt
import numpy as np

import constants as cfg
from circuit_common import magnitude_db, passband_window, plot_s21, read_snp, save_figure
from circuit_models import s21_n_cells
from extract_circuit_params_ncell import extract_capacitances_ncell, MASK_HALF_WIDTH_HZ


def main():
    output = cfg.OUTPUT_DATA_DIR / "generated_ncell_4_comparison"
    output.mkdir(parents=True, exist_ok=True)
    sources = {
        "baseline": cfg.SONNET_FILES_DIR / "filter_v6_cells_n=4.snp",
        "generated": cfg.SONNET_FILES_DIR / "generated_ncell_4" / "generated_ncell_4.snp",
    }
    results, networks = {}, {}
    inductance = cfg.TARGET_L0_NH * 1e-9
    for label, source in sources.items():
        print(f"\n--- {label}: {source} ---", flush=True)
        net = read_snp(source)
        networks[label] = net
        c0, cg = extract_capacitances_ncell(inductance, 4, source)
        window, peaks = passband_window(net, MASK_HALF_WIDTH_HZ)
        fitted = s21_n_cells(2*np.pi*net.f, c0, cg, inductance, 4)
        mask = (net.f >= window[0]) & (net.f <= window[1])
        rmse = np.sqrt(np.mean((magnitude_db(fitted[mask], clip=True) -
                               magnitude_db(net.s[mask,1,0], clip=True))**2))
        results[label] = {
            "source": str(source), "L0_nH_fixed": cfg.TARGET_L0_NH,
            "C0_fF": float(c0), "Cg_fF": float(cg), "fit_rmse_db": float(rmse),
            "peaks_GHz": (peaks/1e9).tolist(), "samples": len(net.f),
            "frequency_range_GHz": [float(net.f[0]/1e9), float(net.f[-1]/1e9)],
        }
        plot_s21(net, fitted, window, output / label, inductance, c0, cg, 4)
        print(json.dumps(results[label], indent=2), flush=True)

    changes = {key: 100*(results['generated'][key]/results['baseline'][key]-1)
               for key in ('C0_fF', 'Cg_fF')}
    a, b = networks['baseline'], networks['generated']
    low, high = max(a.f[0], b.f[0]), min(a.f[-1], b.f[-1])
    common = a.f[(a.f >= low) & (a.f <= high)]
    delta_db = np.interp(common, b.f, magnitude_db(b.s[:,1,0])) - np.interp(
        common, a.f, magnitude_db(a.s[:,1,0]))
    results['comparison'] = {
        "parameter_change_percent": changes,
        "common_range_GHz": [float(low/1e9), float(high/1e9)],
        "s21_magnitude_rmse_db": float(np.sqrt(np.mean(delta_db**2))),
        "note": "L0 is fixed, not extracted. Geometry and reference planes differ; compare magnitude first.",
    }
    if len(results['baseline']['peaks_GHz']) == len(results['generated']['peaks_GHz']):
        results['comparison']['sorted_peak_shifts_MHz'] = (
            1000*(np.array(results['generated']['peaks_GHz'])-
                  np.array(results['baseline']['peaks_GHz']))).tolist()
    fig, axes = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
    for label, net in networks.items():
        axes[0].plot(net.f/1e9, magnitude_db(net.s[:,1,0]), label=label)
        axes[1].plot(net.f/1e9, magnitude_db(net.s[:,0,0]), label=label)
    axes[0].set(ylabel='|S21| (dB)', ylim=(-60, 2), title='Four-cell Sonnet comparison')
    axes[1].set(xlabel='Frequency (GHz)', ylabel='|S11| (dB)', ylim=(-50, 2), xlim=(low/1e9,high/1e9))
    for ax in axes:
        ax.grid(True, alpha=.3)
        ax.legend()
    save_figure(fig, output / 'comparison.png')
    (output / 'results.json').write_text(json.dumps(results, indent=2)+'\n', encoding='utf-8')
    print(json.dumps(results['comparison'], indent=2))


if __name__ == '__main__':
    main()
