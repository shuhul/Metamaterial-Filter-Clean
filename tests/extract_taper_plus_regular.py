"""
Taper cell + regular unit cell S21 fitting.

Topology:  Port1 → Cg3 → [L03‖C03 shunt] → Cg1 → [L01‖C01 shunt] → Cg1 → Port2

Known (fixed):
  - Regular unit cell: L01, C01, Cg1   (from single-cell extraction)
  - Taper inductance:  L03              (from taper isolation extraction)

Fitted:
  - C03, Cg3  (taper shunt cap and coupling cap)

Usage:
  python fixed_scripts/extract_taper_plus_regular.py
"""

import sys, math
import numpy as np
import matplotlib.pyplot as plt
import subprocess, shutil, tempfile, os
import skrf as rf
from scipy.optimize import curve_fit
from scipy.signal import find_peaks

sys.stdout.reconfigure(encoding='utf-8')

# ── Configuration ─────────────────────────────────────────────────────────────
EM_EXE       = r"C:\Program Files\Sonnet Software\19.52.2025.261\bin\em.exe"
PROJECT_FILE = r"S:\shuhul\sonnet\20260409_MM_resonator_example_FY\filter_new_v6_taper1_unit.sonx"
SNP_FILE     = r"S:\shuhul\sonnet\20260409_MM_resonator_example_FY\filter_new_v6_taper1_unit.snp"

RUN_SIMULATION = False

# Regular unit cell — all known
L01_nH = 4.9194
C01_FF  = 191.192
Cg1_FF  = 8.286

# Taper inductance — known from isolation extraction
L03_nH = 6.06

# Reference from taper isolation fit — used as seed and for deviation reporting
C03_REF_FF = 100.0
Cg3_REF_FF = 56.7

MASK_HALF_WIDTH = 1.5e9   # wide window — taper resonator has very wide passband
Z_0 = 50.0

# Diagnostic: pin C03 at isolation reference, fit only Cg3 + Cg1.
# If Cg1 drops to ~8.3 fF, the +10% C03 drift is the root cause.
# If Cg1 stays high, the junction geometry genuinely requires more coupling.
FIX_C03_AT_ISOLATION = True
# ─────────────────────────────────────────────────────────────────────────────

L01 = L01_nH * 1e-9
L03 = L03_nH * 1e-9


# ── I/O helpers ───────────────────────────────────────────────────────────────

def read_snp(snp_file):
    tmp = tempfile.mktemp(suffix=".s2p")
    shutil.copy(snp_file, tmp)
    try:
        return rf.Network(tmp)
    finally:
        os.unlink(tmp)


def run_simulation(project_file):
    result = subprocess.run([EM_EXE, "-v", project_file],
                            capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"Sonnet em failed:\n{result.stderr}")
    return "Begin electromagnetic simulation" in result.stdout


# ── ABCD model ────────────────────────────────────────────────────────────────

def _s21_taper_plus_regular(omega_arr, C03_fF, Cg3_fF, Cg1_fF):
    """
    Port1 → Cg3 → [L03‖C03 shunt] → Cg1 → [L01‖C01 shunt] → Cg1 → Port2

    Cg1 is the regular unit cell coupling cap at the junction — fitted separately
    since the geometry at the taper/regular interface may differ slightly.
    """
    C03 = C03_fF * 1e-15
    Cg3 = Cg3_fF * 1e-15
    C01 = C01_FF  * 1e-15
    Cg1 = Cg1_fF  * 1e-15

    out = np.empty(len(omega_arr), dtype=complex)
    for i, w in enumerate(omega_arr):
        Z_g3  = 1.0 / (1j * w * Cg3)
        Y_lc3 = 1j * w * C03 + 1.0 / (1j * w * L03)
        Z_g1  = 1.0 / (1j * w * Cg1)
        Y_lc1 = 1j * w * C01 + 1.0 / (1j * w * L01)

        M_Cg3 = np.array([[1.0, Z_g3 ], [0.0,   1.0]], dtype=complex)
        M_LC3 = np.array([[1.0, 0.0  ], [Y_lc3, 1.0]], dtype=complex)
        M_Cg1 = np.array([[1.0, Z_g1 ], [0.0,   1.0]], dtype=complex)
        M_LC1 = np.array([[1.0, 0.0  ], [Y_lc1, 1.0]], dtype=complex)

        M = M_Cg3 @ M_LC3 @ M_Cg1 @ M_LC1 @ M_Cg1
        A, B, C, D = M[0,0], M[0,1], M[1,0], M[1,1]
        out[i] = 2.0 / (A + B / Z_0 + C * Z_0 + D)
    return out


# ── Fitting ───────────────────────────────────────────────────────────────────

def _s21_taper_plus_regular_2p(omega_arr, Cg3_fF, Cg1_fF, C03_fixed_fF):
    """Same topology but C03 pinned — only Cg3 and Cg1 are free."""
    return _s21_taper_plus_regular(omega_arr, C03_fixed_fF, Cg3_fF, Cg1_fF)


def extract_taper_capacitances(snp_file, verbose=True):
    mode = "C03 fixed at isolation ref" if FIX_C03_AT_ISOLATION else "C03, Cg3, Cg1 free"
    if verbose:
        print(f"── Taper + regular cell fit ({mode}) ────────")

    net     = read_snp(snp_file)
    s21_mag = np.abs(net.s[:, 1, 0])

    # Find peaks and build fit window tightly around the two-peak region
    peaks_idx, _ = find_peaks(s21_mag, prominence=0.005)
    if len(peaks_idx) >= 2:
        peak_freqs = net.f[peaks_idx]
        f_lower = peak_freqs[0]
        f_upper = peak_freqs[-1]
        if verbose:
            print(f"  Found {len(peaks_idx)} peaks: "
                  + ", ".join(f"{f/1e9:.4f}" for f in peak_freqs) + " GHz")
        # Tight window — 400 MHz below lowest peak, 400 MHz above highest
        f_lo_edge = f_lower - 0.4e9
        f_hi_edge = f_upper + 0.4e9
    else:
        f_ctr = net.f[peaks_idx[0]] if len(peaks_idx) == 1 else net.f[np.argmax(s21_mag)]
        f_lo_edge = f_ctr - MASK_HALF_WIDTH
        f_hi_edge = f_ctr + MASK_HALF_WIDTH
        if verbose:
            print(f"  {len(peaks_idx)} peak(s) found, centre at {f_ctr/1e9:.4f} GHz")

    mask = (net.f >= f_lo_edge) & (net.f <= f_hi_edge)
    if verbose:
        print(f"  Fit window: {f_lo_edge/1e9:.3f}–{f_hi_edge/1e9:.3f} GHz  "
              f"({mask.sum()} points)")

    omega = 2 * np.pi * net.f[mask]
    y     = s21_mag[mask]   # linear amplitude — peaks weighted naturally

    if FIX_C03_AT_ISOLATION:
        # ── 2-param fit: C03 pinned at isolation reference ────────────────────
        C03_pin = C03_REF_FF
        if verbose:
            print(f"  C03 pinned at {C03_pin:.1f} fF (isolation ref)")

        def fit_wrapper(omega, Cg3_fF, Cg1_fF):
            return np.abs(_s21_taper_plus_regular_2p(omega, Cg3_fF, Cg1_fF, C03_pin))

        def residual_2p(cg3, cg1):
            return float(np.sum((fit_wrapper(omega, cg3, cg1) - y) ** 2))

        Cg3_grid = np.linspace(Cg3_REF_FF * 0.50, Cg3_REF_FF * 1.80, 12)
        Cg1_grid = np.linspace(4.0, 18.0, 10)

        grid_results = []
        for cg3 in Cg3_grid:
            for cg1 in Cg1_grid:
                r = residual_2p(cg3, cg1)
                grid_results.append((r, cg3, cg1))
        grid_results.sort(key=lambda x: x[0])

        if verbose:
            print(f"  Grid search ({len(grid_results)} pts), "
                  f"best: Cg3={grid_results[0][1]:.1f}, Cg1={grid_results[0][2]:.1f} fF  "
                  f"res={grid_results[0][0]:.3e}")

        bounds_lo = [Cg3_REF_FF * 0.40, 3.0]
        bounds_hi = [Cg3_REF_FF * 2.00, 20.0]

        best_popt, best_res = None, np.inf
        for _, cg3, cg1 in grid_results[:8]:
            try:
                popt, _ = curve_fit(
                    fit_wrapper, omega, y,
                    p0=[cg3, cg1],
                    bounds=(bounds_lo, bounds_hi),
                    max_nfev=15000,
                )
                res = float(np.sum((fit_wrapper(omega, *popt) - y) ** 2))
                if res < best_res:
                    best_res  = res
                    best_popt = popt
                    if verbose:
                        print(f"  Seed (Cg3={cg3:.1f}, Cg1={cg1:.1f}) → "
                              f"Cg3={popt[0]:.3f}, Cg1={popt[1]:.3f}  res={res:.3e}  ✓")
            except Exception:
                pass

        if best_popt is None:
            raise RuntimeError("All seeds failed (2-param mode)")

        C03_fF = C03_pin
        Cg3_fF, Cg1_fF = best_popt

    else:
        # ── 3-param fit: C03, Cg3, Cg1 all free ─────────────────────────────
        def fit_wrapper(omega, C03_fF, Cg3_fF, Cg1_fF):
            s21_c = _s21_taper_plus_regular(omega, C03_fF, Cg3_fF, Cg1_fF)
            return np.abs(s21_c)

        def residual(C03_fF, Cg3_fF, Cg1_fF):
            return float(np.sum((fit_wrapper(omega, C03_fF, Cg3_fF, Cg1_fF) - y) ** 2))

        C03_grid = np.linspace(C03_REF_FF * 0.75, C03_REF_FF * 1.25, 6)
        Cg3_grid = np.linspace(Cg3_REF_FF * 0.50, Cg3_REF_FF * 1.80, 8)
        Cg1_grid = np.linspace(4.0, 18.0, 6)

        grid_results = []
        for c03 in C03_grid:
            for cg3 in Cg3_grid:
                for cg1 in Cg1_grid:
                    r = residual(c03, cg3, cg1)
                    grid_results.append((r, c03, cg3, cg1))
        grid_results.sort(key=lambda x: x[0])

        if verbose:
            print(f"  Grid search ({len(grid_results)} pts), "
                  f"best grid residual: {grid_results[0][0]:.3e} "
                  f"at C03={grid_results[0][1]:.1f}, Cg3={grid_results[0][2]:.1f}, "
                  f"Cg1={grid_results[0][3]:.1f} fF")

        bounds_lo = [C03_REF_FF * 0.60, Cg3_REF_FF * 0.40, 3.0]
        bounds_hi = [C03_REF_FF * 1.40, Cg3_REF_FF * 2.00, 20.0]

        best_popt, best_res = None, np.inf
        for _, c03, cg3, cg1 in grid_results[:8]:
            try:
                popt, _ = curve_fit(
                    fit_wrapper, omega, y,
                    p0=[c03, cg3, cg1],
                    bounds=(bounds_lo, bounds_hi),
                    max_nfev=15000,
                )
                res = float(np.sum((fit_wrapper(omega, *popt) - y) ** 2))
                if res < best_res:
                    best_res  = res
                    best_popt = popt
                    if verbose:
                        print(f"  Seed ({c03:.0f}, {cg3:.1f}, {cg1:.1f}) → "
                              f"C03={popt[0]:.2f}, Cg3={popt[1]:.3f}, Cg1={popt[2]:.3f}  "
                              f"res={res:.3e}  ✓")
            except Exception:
                pass

        if best_popt is None:
            raise RuntimeError("All seeds failed")

        C03_fF, Cg3_fF, Cg1_fF = best_popt
    if verbose:
        print(f"  ➜  C03={C03_fF:.3f} fF  Cg3={Cg3_fF:.3f} fF  Cg1={Cg1_fF:.3f} fF")
        print("──────────────────────────────────────────────────────────")

    return C03_fF, Cg3_fF, Cg1_fF


# ── Main ──────────────────────────────────────────────────────────────────────

if RUN_SIMULATION:
    fresh = run_simulation(PROJECT_FILE)
    print(f"Simulation: {'fresh' if fresh else 'cached'}")

C03, Cg3, Cg1_fit = extract_taper_capacitances(SNP_FILE)

dC03 = (C03  - C03_REF_FF) / C03_REF_FF * 100
dCg3 = (Cg3  - Cg3_REF_FF) / Cg3_REF_FF * 100
dCg1 = (Cg1_fit - Cg1_FF)  / Cg1_FF     * 100

print(f"\nFitted (taper cell):")
print(f"  L03   = {L03_nH:.4f} nH  (fixed)")
print(f"  C03   = {C03:.3f} fF  (ref {C03_REF_FF:.1f}, {dC03:+.2f}%)")
print(f"  Cg3   = {Cg3:.3f} fF  (ref {Cg3_REF_FF:.2f}, {dCg3:+.2f}%)")
print(f"  Cg1   = {Cg1_fit:.3f} fF  (ref {Cg1_FF:.3f}, {dCg1:+.2f}%)")

# ── Plot ──────────────────────────────────────────────────────────────────────

net        = read_snp(SNP_FILE)
freqs      = net.f
s21_raw    = net.s[:, 1, 0]
s21_mag    = np.abs(s21_raw)

peaks_idx, _ = find_peaks(s21_mag, prominence=0.005)
if len(peaks_idx) >= 2:
    f_plot_lo = net.f[peaks_idx[0]]  - MASK_HALF_WIDTH
    f_plot_hi = net.f[peaks_idx[-1]] + MASK_HALF_WIDTH
else:
    f_ctr = net.f[peaks_idx[0]] if len(peaks_idx) else net.f[np.argmax(s21_mag)]
    f_plot_lo = f_ctr - MASK_HALF_WIDTH
    f_plot_hi = f_ctr + MASK_HALF_WIDTH

omega_full  = 2 * np.pi * freqs
s21_fit     = _s21_taper_plus_regular(omega_full, C03, Cg3, Cg1_fit)
s21_ref     = _s21_taper_plus_regular(omega_full, C03_REF_FF, Cg3_REF_FF, Cg1_FF)

s21_data_db = 20 * np.log10(s21_mag + 1e-12)
s21_fit_db  = 20 * np.log10(np.abs(s21_fit)  + 1e-12)
s21_ref_db  = 20 * np.log10(np.abs(s21_ref)  + 1e-12)

plt.rcParams.update({'font.size': 13})
fig, ax = plt.subplots(figsize=(9, 5))

ax.plot(freqs / 1e9, s21_data_db, 'r-',  lw=2,   label='Sonnet simulation')
ax.plot(freqs / 1e9, s21_fit_db,  'k--', lw=2,   label='ABCD fit')
ax.plot(freqs / 1e9, s21_ref_db,  'b:',  lw=1.5, label='Reference (isolation vals)')

ax.set_xlim([f_plot_lo / 1e9 - 0.2, f_plot_hi / 1e9 + 0.2])
ax.set_ylim([-60, 5])
ax.set_xlabel('Frequency (GHz)')
ax.set_ylabel('|S21| (dB)')
fit_label = "C03 pinned at isolation" if FIX_C03_AT_ISOLATION else "C03, Cg3, Cg1 free"
ax.set_title(f'Taper + Regular cell  —  ABCD fit  [{fit_label}]')
ax.grid(True, linestyle=':', alpha=0.7)
ax.legend(loc='lower right')

text = (f"$C_{{03}} = {C03:.1f}$ fF  (ref {C03_REF_FF:.1f}, {dC03:+.2f}%)\n"
        f"$C_{{g3}} = {Cg3:.2f}$ fF  (ref {Cg3_REF_FF:.2f}, {dCg3:+.2f}%)\n"
        f"$C_{{g1}} = {Cg1_fit:.2f}$ fF  (ref {Cg1_FF:.2f}, {dCg1:+.2f}%)\n"
        f"$L_{{03}} = {L03_nH:.3f}$ nH  (fixed)")
ax.text(0.05, 0.95, text, transform=ax.transAxes, fontsize=11,
        verticalalignment='top',
        bbox=dict(boxstyle='round', facecolor='white', alpha=0.9, edgecolor='gray'))

plt.tight_layout()
suffix = "_c03fixed" if FIX_C03_AT_ISOLATION else ""
out_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), f'taper_plus_regular_fit{suffix}.png')
plt.savefig(out_path, dpi=300)
print(f"Saved {out_path}")

# ── Linear scale plot ─────────────────────────────────────────────────────────
fig2, ax2 = plt.subplots(figsize=(9, 5))

ax2.plot(freqs / 1e9, s21_mag,            'r-',  lw=2,   label='Sonnet simulation')
ax2.plot(freqs / 1e9, np.abs(s21_fit),    'k--', lw=2,   label='ABCD fit')
ax2.plot(freqs / 1e9, np.abs(s21_ref),    'b:',  lw=1.5, label='Reference (isolation vals)')

ax2.set_xlim([f_plot_lo / 1e9 - 0.2, f_plot_hi / 1e9 + 0.2])
ax2.set_ylim([0, 1.05])
ax2.set_xlabel('Frequency (GHz)')
ax2.set_ylabel('|S21| (linear)')
ax2.set_title(f'Taper + Regular cell  —  ABCD fit (linear)  [{fit_label}]')
ax2.grid(True, linestyle=':', alpha=0.7)
ax2.legend(loc='lower right')
ax2.text(0.05, 0.95, text, transform=ax2.transAxes, fontsize=11,
         verticalalignment='top',
         bbox=dict(boxstyle='round', facecolor='white', alpha=0.9, edgecolor='gray'))

plt.tight_layout()
out_path2 = os.path.join(os.path.dirname(os.path.abspath(__file__)), f'taper_plus_regular_fit_linear{suffix}.png')
plt.savefig(out_path2, dpi=300)
print(f"Saved {out_path2}")
