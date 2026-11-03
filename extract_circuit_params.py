"""Single-cell inductance sweep and complex-S21 capacitance extraction."""

import argparse
from pathlib import Path
import matplotlib.pyplot as plt
from scipy.stats import linregress
import numpy as np

import constants as cfg
from circuit_common import (
    CG_SEEDS_FF, find_resonance, fit_capacitance_seeds, plot_s21,
    read_snp, report_parameters, run_simulation, save_figure, write_sonx, single_sweep_resonance,
)
from circuit_models import resonance_c0_seed, s21_single


# Single-cell settings. Only this workflow edits/sweeps C_add.
FILE_NAME = "filter_v6"
PROJECT_FILE = cfg.SONNET_FILES_DIR / f"{FILE_NAME}.sonx"
SNP_FILE = cfg.SONNET_FILES_DIR / f"{FILE_NAME}.snp"
OUTPUT_DIR = cfg.OUTPUT_DATA_DIR / FILE_NAME
RUN_SIMULATION = True
RUN_INDUCTANCE = True
# Fixed C0 and weak coupling mainly affect the intercept; dominance is not required.
C_ADD_SWEEP_PF = (0.02, 0.03, 0.04, 0.05)
CAPACITANCE_C_ADD_PF = 0.0001
MASK_HALF_WIDTH_HZ = 1.50e9
MAX_NFEV = 5000


def extract_inductance(project_file, snp_file, em_exe, output_dir):
    """Sweep C_add (pF), fit 1/omega^2 = L0 * C_add + intercept, and return H."""
    frequencies = []
    for capacitance in C_ADD_SWEEP_PF:
        write_sonx({"C_add": f"{capacitance:.3f}"}, project_file)
        fresh = run_simulation(project_file, em_exe)
        net = read_snp(snp_file)
        frequency = single_sweep_resonance(net)
        frequencies.append(frequency)
        print(f"  C_add={capacitance * 1e3:.1f} fF; f_res={frequency / 1e9:.4f} GHz; "
              f"{'fresh' if fresh else 'cached'}")
    x = np.asarray(C_ADD_SWEEP_PF) * 1e-12
    if not np.all(np.diff(frequencies) < 0):
        raise ValueError("Resonance must decrease with increasing C_add")
    y = 1 / (2 * np.pi * np.asarray(frequencies)) ** 2
    slope, intercept, r_value, *_ = linregress(x, y)
    if slope <= 0:
        raise ValueError("Inductance slope must be positive")
    print(f"L0 = {slope * 1e9:.4f} nH (R^2 = {r_value**2:.6f})")
    fig, ax = plt.subplots(figsize=(7, 5))
    ax.plot(x * 1e15, y * 1e21, "ko", label="Sonnet data")
    x_fit = np.linspace(x[0], x[-1], 200)
    ax.plot(x_fit * 1e15, (intercept + slope * x_fit) * 1e21, "r--",
            label=f"Linear fit (R^2={r_value**2:.5f})")
    ax.set(xlabel="C_add (fF)", ylabel="1/omega^2 (10^-21 s^2)",
           title=f"Inductance extraction: L0 = {slope * 1e9:.3f} nH")
    ax.legend()
    ax.grid(True, linestyle=":", alpha=0.7)
    save_figure(fig, Path(output_dir) / "inductance_fit.png")
    return slope


def extract_capacitances(L_0, project_file, snp_file, em_exe, verbose=True,
                         *, simulate=True):
    """Return (C0, Cg) in fF for a fixed L_0 in H by fitting complex S21.

    Simulation sets C_add to the configured small nonzero value. The analytical
    model neglects this capacitor. simulate=False only reads existing data.
    """
    if simulate:
        write_sonx({"C_add": str(CAPACITANCE_C_ADD_PF)}, project_file)
        fresh = run_simulation(project_file, em_exe)
        if verbose:
            print(f"Simulation: {'fresh' if fresh else 'cached'}")
    net = read_snp(snp_file)
    peak = find_resonance(net)
    mask = ((net.f >= peak - MASK_HALF_WIDTH_HZ)
            & (net.f <= peak + MASK_HALF_WIDTH_HZ))
    omega = 2 * np.pi * net.f[mask]
    s21 = net.s[:, 1, 0][mask]
    data = np.hstack([s21.real, s21.imag])

    def model(omega, c0, cg):
        response = s21_single(omega, c0 * 1e-15, cg * 1e-15, L_0)
        return np.hstack([response.real, response.imag])

    seeds = [(max(1.0, resonance_c0_seed(peak, cg * 1e-15, L_0) * 1e15), cg)
             for cg in CG_SEEDS_FF]
    if verbose:
        print(f"Resonance: {peak / 1e9:.4f} GHz; "
              f"fit window +/-{MASK_HALF_WIDTH_HZ / 1e6:.0f} MHz")
    return fit_capacitance_seeds(model, omega, data, seeds, MAX_NFEV, verbose)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--no-simulate", action="store_true",
                        help="Fit existing .snp with fixed L0; skip Sonnet and project edits")
    args = parser.parse_args(argv)
    simulate = RUN_SIMULATION and not args.no_simulate
    if RUN_INDUCTANCE and simulate:
        inductance = extract_inductance(PROJECT_FILE, SNP_FILE,
                                       cfg.EM_EXE, OUTPUT_DIR)
    else:
        inductance = cfg.TARGET_L0_NH * 1e-9
        print(f"Using fixed L0 = {cfg.TARGET_L0_NH} nH")
    c0, cg = extract_capacitances(inductance, PROJECT_FILE, SNP_FILE,
                                 cfg.EM_EXE, simulate=simulate)
    report_parameters(inductance, c0, cg)
    net = read_snp(SNP_FILE)
    peak = find_resonance(net)
    window = (peak - MASK_HALF_WIDTH_HZ, peak + MASK_HALF_WIDTH_HZ)
    fitted = s21_single(2 * np.pi * net.f, c0 * 1e-15, cg * 1e-15, inductance)
    plot_s21(net, fitted, window, OUTPUT_DIR, inductance, c0, cg)


if __name__ == "__main__":
    main()
