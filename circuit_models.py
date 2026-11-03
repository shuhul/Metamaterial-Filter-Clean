"""Circuit models: SI units for the single cell, fF inputs for the cascade."""


import numpy as np
import constants as cfg


def s21_single(omega, C_0, C_g, L_0):
    """Complex S21 for one ideal resonator; C_0/C_g in F and L_0 in H."""
    term_A    = 1 + (C_0 / C_g) - (1 / (omega**2 * L_0 * C_g))
    imag_part = (
          (omega * C_0 * cfg.Z_0)
        - (cfg.Z_0 / (omega * L_0))
        - (1 / (omega * C_g * cfg.Z_0)) * (2 + (C_0 / C_g) - (1 / (omega**2 * L_0 * C_g)))
    )
    return 2 / (2 * term_A + 1j * imag_part)


def resonance_c0_seed(f_res_hz, C_g, L_0):
    """Solve Im(S21)=0 numerically for C_0."""
    omega = 2 * np.pi * f_res_hz

    def imag_part(C_0):
        return (
            omega * C_0 * cfg.Z_0
            - cfg.Z_0 / (omega * L_0)
            - (1 / (omega * C_g * cfg.Z_0)) * (2 + (C_0 / C_g) - (1 / (omega**2 * L_0 * C_g)))
        )

    C_lo, C_hi = 1e-15, 1000e-15
    f_lo, f_hi = imag_part(C_lo), imag_part(C_hi)
    if f_lo * f_hi > 0:
        return max(1e-15, 1.0 / (L_0 * omega**2))

    for _ in range(60):
        C_mid = (C_lo + C_hi) / 2
        if imag_part(C_mid) * f_lo < 0:
            C_hi = C_mid
        else:
            C_lo = C_mid

    return max(1e-15, (C_lo + C_hi) / 2)


def s21_n_cells(omega_arr, C0_fF, Cg_fF, L0, N, Z0=None):
    """
    S21 for an N-cell Purcell filter via ABCD matrix cascade.

    Topology: Cg | (L0‖C0 | Cg)^N
      M_total = M_Cg · (M_LC · M_Cg)^N

    Verified: N=1 reproduces the closed-form single-cell S21 exactly.
    """
    Z0 = cfg.Z_0 if Z0 is None else Z0
    if not isinstance(N, (int, np.integer)) or N < 1:
        raise ValueError("N must be a positive integer")
    C0 = C0_fF * 1e-15
    Cg = Cg_fF * 1e-15
    out = np.empty(len(omega_arr), dtype=complex)
    for i, w in enumerate(omega_arr):
        Z_g  = 1.0 / (1j * w * Cg)
        Y_lc = 1j * w * C0 + 1.0 / (1j * w * L0)
        M_Cg = np.array([[1.0, Z_g ], [0.0, 1.0 ]], dtype=complex)
        M_LC = np.array([[1.0, 0.0 ], [Y_lc, 1.0]], dtype=complex)
        # Start with the input coupling cap, then append N (LC · Cg) blocks
        M = M_Cg.copy()
        for _ in range(N):
            M = M @ M_LC @ M_Cg
        A, B, C, D = M[0, 0], M[0, 1], M[1, 0], M[1, 1]
        out[i] = 2.0 / (A + B / Z0 + C * Z0 + D)
    return out
