"""Shared targets, paths, and reference impedance; workflow settings live in each script."""

from pathlib import Path
from geometry.specifications import UNIT_CELL

BASE_DIR = Path(__file__).resolve().parent
SONNET_FILES_DIR = BASE_DIR / "sonnet_files"
OUTPUT_DATA_DIR = BASE_DIR / "output_data"
EM_EXE = Path(r"C:\Program Files\Sonnet Software\19.52.2025.261\bin\em.exe")
Z_0 = 50.0
TARGET_L0_NH = UNIT_CELL.L0_nH
TARGET_C0_FF = UNIT_CELL.C0_fF
TARGET_CG_FF = UNIT_CELL.Cg_input_fF
TARGET_FL_HZ = 4.8e9
TARGET_FH_HZ = 5.2e9
TARGET_CG_C0_RATIO = ((TARGET_FH_HZ / TARGET_FL_HZ) ** 2 - 1) / 4
