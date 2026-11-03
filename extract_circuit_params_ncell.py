"""Multi-cell Purcell filter: fixed-L0 capacitance extraction using an ABCD cascade."""

import argparse
import json
from pathlib import Path
import numpy as np

import constants as cfg
from circuit_common import (
    CG_SEEDS_FF, find_resonance, fit_capacitance_seeds, magnitude_db, passband_window,
    plot_s21, read_snp, report_parameters, run_simulation,
)
from circuit_models import s21_n_cells


# Multi-cell settings. L0 is always taken directly from the shared target.
N_CELLS = 4
FILE_NAME = f"filter_v6_cells_n={N_CELLS}"
PROJECT_FILE = cfg.SONNET_FILES_DIR / f"{FILE_NAME}.sonx"
SNP_FILE = cfg.SONNET_FILES_DIR / f"{FILE_NAME}.snp"
OUTPUT_DIR = cfg.OUTPUT_DATA_DIR / FILE_NAME
RUN_SIMULATION = False
MASK_HALF_WIDTH_HZ = 0.5e9
MAX_NFEV = 8000


def extract_capacitances_ncell(L_0, N, snp_file, verbose=True):
    """Fit clipped S21 magnitude in dB; return (C0, Cg) in fF for L_0 in H."""
    net = read_snp(snp_file)
    (low, high), peaks = passband_window(net, MASK_HALF_WIDTH_HZ)
    if len(peaks) >= 2:
        c0_phys = max(1.0, 1e15 / (L_0 * (2 * np.pi * peaks[-1]) ** 2))
        total_c = 1e15 / (L_0 * (2 * np.pi * peaks[0]) ** 2)
        cg_phys = max(0.1, (total_c - c0_phys) / 4)
    else:
        center = find_resonance(net)
        c0_phys = max(1.0, 1e15 / (L_0 * (2 * np.pi * center) ** 2))
        cg_phys = 5.0
    mask = (net.f >= low) & (net.f <= high)
    omega = 2 * np.pi * net.f[mask]
    data = magnitude_db(net.s[:, 1, 0][mask], clip=True)

    def model(omega, c0, cg):
        return magnitude_db(s21_n_cells(omega, c0, cg, L_0, N), clip=True)

    seeds = [(cfg.TARGET_C0_FF, cfg.TARGET_CG_FF)]
    seeds.append((c0_phys, cg_phys))
    seeds.extend((c0_phys, cg) for cg in CG_SEEDS_FF)
    if verbose:
        print(f"N={N}; found {len(peaks)} peaks; "
              f"fit window {low / 1e9:.4f}-{high / 1e9:.4f} GHz ({mask.sum()} points)")
    return fit_capacitance_seeds(model, omega, data, seeds, MAX_NFEV, verbose)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--no-simulate", action="store_true",
                        help="Fit the existing .snp without running Sonnet")
    parser.add_argument("--snp", type=Path, default=SNP_FILE)
    parser.add_argument("--output-dir", type=Path, default=OUTPUT_DIR)
    parser.add_argument("--n-cells", type=int, default=N_CELLS)
    parser.add_argument("--l0-nh", type=float, default=cfg.TARGET_L0_NH,
                        help="Fixed inductance in nH; defaults to target, never implicitly reads old results")
    args = parser.parse_args(argv)
    if args.n_cells < 1 or not np.isfinite(args.l0_nh) or args.l0_nh <= 0:
        parser.error("n-cells and l0-nh must be positive")
    if args.snp != SNP_FILE and RUN_SIMULATION and not args.no_simulate:
        parser.error("Use --no-simulate when fitting an explicit --snp file")
    if RUN_SIMULATION and not args.no_simulate:
        fresh = run_simulation(PROJECT_FILE, cfg.EM_EXE)
        print(f"Simulation: {'fresh' if fresh else 'cached'}")
    inductance = args.l0_nh * 1e-9
    c0, cg = extract_capacitances_ncell(inductance, args.n_cells, args.snp)
    report_parameters(inductance, c0, cg, args.n_cells)
    net = read_snp(args.snp)
    window, _ = passband_window(net, MASK_HALF_WIDTH_HZ)
    fitted = s21_n_cells(2 * np.pi * net.f, c0, cg, inductance, args.n_cells)
    plot_s21(net, fitted, window, args.output_dir, inductance, c0, cg, args.n_cells)
    (args.output_dir / 'fit_parameters.json').write_text(json.dumps({
        'source': str(args.snp.resolve()), 'n_cells': args.n_cells,
        'L0_nH_fixed': args.l0_nh, 'C0_fF': float(c0), 'Cg_fF': float(cg)
    }, indent=2)+'\n')


if __name__ == "__main__":
    main()
