"""Shared Sonnet I/O, fitting, reporting, and plotting utilities."""

from pathlib import Path
import re
import shutil
import subprocess
import tempfile

import matplotlib

# These scripts save PNG files; they do not need a GUI/Tk installation.
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import skrf as rf
from scipy.optimize import curve_fit
from scipy.signal import find_peaks

import constants as cfg

# Shared numerical and plotting controls, owned by these helpers.
CG_SEEDS_FF = (0.5, 1.0, 2.0, 5.0, 10.0, 20.0)
CAPACITANCE_BOUNDS_FF = ((0.1, 0.1), (2000.0, 2000.0))
RESONANCE_PROMINENCE = 0.02
PASSBAND_PROMINENCE = 0.005
DB_EPSILON = 1e-12
FIT_DB_MIN = -80.0
FIT_DB_MAX = 0.0
PLOT_DPI = 300


def set_param(content, name, value):
    """Replace a named Sonnet Variable value, failing if it is absent."""
    pattern = rf'(<Variable\b[^>]*\bName="{re.escape(name)}"[^>]*\bValue=")[^"]*(")'
    content, count = re.subn(pattern, lambda m: f"{m[1]}{value}{m[2]}", content)
    if not count:
        raise ValueError(f"Parameter '{name}' not found in .sonx")
    return content


def write_sonx(params, project_file):
    """Update project variables in place, as required by a Sonnet sweep."""
    path = Path(project_file)
    content = path.read_text(encoding="utf-8")
    for name, value in params.items():
        content = set_param(content, name, value)
    path.write_text(content, encoding="utf-8")


def run_simulation(project_file, em_exe=None):
    """Run Sonnet; return True for a fresh run, False for cached output."""
    project_file = Path(project_file).resolve()
    from geometry.taper.solver_guard import run as guarded_solver
    log_path = project_file.with_suffix('.em.log')
    with log_path.open('w') as log:
        result = guarded_solver([str(em_exe or cfg.EM_EXE), "-v", str(project_file)],
                                cwd=project_file.parent, stdout=log)
    text = log_path.read_text(errors='replace')
    if result.returncode:
        raise RuntimeError(f"Sonnet em failed:\n{text}")
    return "Begin electromagnetic simulation" in text


def read_snp(snp_file):
    """Load a two-port Sonnet .snp using a safely managed temporary .s2p copy."""
    with tempfile.TemporaryDirectory() as directory:
        temp_path = Path(directory) / "network.s2p"
        shutil.copyfile(snp_file, temp_path)
        net = rf.Network(str(temp_path))
    if net.nports != 2 or len(net.f) < 2:
        raise ValueError("Expected a two-port network with at least two frequency samples")
    if (not np.all(np.isfinite(net.f)) or not np.all(np.isfinite(net.s))
            or np.any(net.f <= 0) or np.any(np.diff(net.f) <= 0)):
        raise ValueError("S-parameters must be finite with positive, increasing frequencies")
    if not np.allclose(net.z0, cfg.Z_0):
        raise ValueError(f"Circuit extraction expects {cfg.Z_0:g} ohm port normalization")
    return net


def find_resonance(net):
    magnitude = np.abs(net.s[:, 1, 0])
    peaks, props = find_peaks(magnitude, prominence=RESONANCE_PROMINENCE)
    index = peaks[np.argmax(props["prominences"])] if len(peaks) else np.argmax(magnitude)
    return net.f[index]


def single_sweep_resonance(net):
    """Require one resolved interior resonance for the C_add slope method."""
    magnitude = np.abs(net.s[:, 1, 0])
    index = int(np.argmax(magnitude))
    peaks, _ = find_peaks(magnitude, prominence=RESONANCE_PROMINENCE)
    if index in (0, len(net.f)-1):
        raise ValueError("Resonance is at the sweep boundary; expand the frequency sweep")
    if len(peaks) != 1 or peaks[0] != index:
        raise ValueError("Inductance extraction requires one prominent resonance to track")
    return float(net.f[index])


def passband_window(net, half_width_hz):
    """Return window edges and peaks, using the dominant peak as fallback."""
    indices, _ = find_peaks(np.abs(net.s[:, 1, 0]), prominence=PASSBAND_PROMINENCE)
    peaks = net.f[indices]
    if len(peaks) >= 2:
        return (peaks[0] - half_width_hz, peaks[-1] + half_width_hz), peaks
    center = find_resonance(net)
    return (center - half_width_hz, center + half_width_hz), peaks


def magnitude_db(s21, clip=False):
    values = 20 * np.log10(np.abs(s21) + DB_EPSILON)
    return np.clip(values, FIT_DB_MIN, FIT_DB_MAX) if clip else values


def fit_capacitance_seeds(model, omega, data, seeds, max_nfev, verbose=True):
    """Fit each (C0, Cg) seed in fF and return the least-squares winner."""
    if len(omega) < 2:
        raise ValueError("Fit window must contain at least two frequency samples")
    best, best_residual = None, np.inf
    for c0, cg in seeds:
        try:
            fitted, _ = curve_fit(model, omega, data, p0=[c0, cg],
                                  bounds=CAPACITANCE_BOUNDS_FF, max_nfev=max_nfev)
            residual = float(np.sum((model(omega, *fitted) - data) ** 2))
            if residual < best_residual:
                best, best_residual = fitted, residual
                if verbose:
                    print(f"  Seed ({c0:.2f}, {cg:.2f}) fF -> "
                          f"C0={fitted[0]:.3f}, Cg={fitted[1]:.3f}; residual={residual:.3e}")
        except (RuntimeError, ValueError, FloatingPointError) as error:
            if verbose:
                print(f"  Seed ({c0:.2f}, {cg:.2f}) fF failed: {error}")
    if best is None:
        raise RuntimeError("All capacitance fit seeds failed; check data, bounds, and fit window")
    return tuple(best)


def report_parameters(L0, C0_fF, Cg_fF, n=1):
    """Print independent component references and the model's implied band."""
    high = 1 / (2 * np.pi * np.sqrt(L0 * C0_fF * 1e-15))
    low = 1 / (2 * np.pi * np.sqrt(L0 * (C0_fF + 4 * Cg_fF) * 1e-15))
    ratio = Cg_fF / C0_fF
    print(f"\nN = {n}")
    print(f"L0 = {L0 * 1e9:.4f} nH (component target {cfg.TARGET_L0_NH})")
    print(f"C0 = {C0_fF:.3f} fF (component target {cfg.TARGET_C0_FF})")
    print(f"Cg = {Cg_fF:.3f} fF (component target {cfg.TARGET_CG_FF})")
    print(f"Cg/C0 = {ratio * 100:.2f}% (band target {cfg.TARGET_CG_C0_RATIO * 100:.2f}%)")
    print(f"  {'Increase' if ratio < cfg.TARGET_CG_C0_RATIO else 'Decrease'} Cg or "
          f"{'decrease' if ratio < cfg.TARGET_CG_C0_RATIO else 'increase'} C0")
    print(f"Implied passband: {low / 1e9:.4f}-{high / 1e9:.4f} GHz "
          f"(BW = {(high - low) / 1e6:.0f} MHz)")
    print(f"Target band: {cfg.TARGET_FL_HZ / 1e9:.4f}-{cfg.TARGET_FH_HZ / 1e9:.4f} GHz "
          f"(BW = {(cfg.TARGET_FH_HZ - cfg.TARGET_FL_HZ) / 1e6:.0f} MHz)")


def save_figure(fig, path):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(path, dpi=PLOT_DPI)
    plt.close(fig)
    print(f"Saved {path}")


def plot_s21(net, fitted, window, output_dir, L0, C0_fF, Cg_fF, n=1):
    """Use the same plot style and target comparisons for every cell count."""
    low, high = window
    plt.rcParams.update({"font.size": 13})
    fig, ax = plt.subplots(figsize=(9, 5))
    ax.plot(net.f / 1e9, magnitude_db(net.s[:, 1, 0]), "r-", lw=2, label="Sonnet simulation")
    ax.plot(net.f / 1e9, magnitude_db(fitted), "k--", lw=2, label="Circuit fit")
    ax.axvspan(low / 1e9, high / 1e9, color="blue", alpha=0.1, label="Fit window")
    ax.axvspan(cfg.TARGET_FL_HZ / 1e9, cfg.TARGET_FH_HZ / 1e9,
               color="green", alpha=0.06, label="Target band")
    ax.set(xlim=(low / 1e9 - 0.5, high / 1e9 + 0.5), ylim=(-60, 5),
           xlabel="Frequency (GHz)", ylabel="|S21| (dB)",
           title=f"N = {n} cell Purcell filter; L0 = {L0 * 1e9:.3f} nH")
    ax.grid(True, linestyle=":", alpha=0.7)
    ax.legend(loc="lower right")
    dc0 = (C0_fF / cfg.TARGET_C0_FF - 1) * 100
    dcg = (Cg_fF / cfg.TARGET_CG_FF - 1) * 100
    annotation = (f"C0 = {C0_fF:.1f} fF (target {cfg.TARGET_C0_FF:.1f}, {dc0:+.2f}%)\n"
                  f"Cg = {Cg_fF:.2f} fF (target {cfg.TARGET_CG_FF:.4f}, {dcg:+.2f}%)")
    ax.text(0.05, 0.95, annotation, transform=ax.transAxes, fontsize=11, va="top",
            bbox=dict(boxstyle="round", facecolor="white", alpha=0.9, edgecolor="gray"))
    save_figure(fig, Path(output_dir) / "s21_fit.png")
