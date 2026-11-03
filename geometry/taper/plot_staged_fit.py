"""Plot a saved joint fit with target comparisons; never run EM or refit."""
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from geometry.specifications import TAPER_1, TAPER_2, UNIT_CELL


def plot(folder):
    folder = Path(folder)
    result = json.loads((folder / 'results.json').read_text())
    values = result['fitted']
    with np.load(folder / 'idc_response.npz') as data:
        frequency = data['f'] / 1e9
        actual, model = data['actual'], data['model']

    rows = [
        ('L T2', result['fixed_inductances_nH'][1][0], TAPER_2.L0_nH, 'nH'),
        ('C0 T2', values['C0_T2_fF'], TAPER_2.C0_fF, 'fF'),
        ('Cg T2 (input)', values['Cg_input_T2_fF'], TAPER_2.Cg_input_fF, 'fF'),
        ('L T1', result['fixed_inductances_nH'][1][1], TAPER_1.L0_nH, 'nH'),
        ('C0 T1', values['C0_T1_fF'], TAPER_1.C0_fF, 'fF'),
        ('Cg T1 (T2 to T1)', values['Cg_T2_T1_fF'], TAPER_1.Cg_input_fF, 'fF'),
        ('Cg (T1 to unit)', values['Cg_T1_regular_fF'], UNIT_CELL.Cg_input_fF, 'fF'),
    ]
    with plt.rc_context({'font.size': 15, 'axes.labelsize': 20,
                         'xtick.labelsize': 15, 'ytick.labelsize': 15}):
        fig, (ax, panel) = plt.subplots(
            1, 2, figsize=(19, 8), gridspec_kw={'width_ratios': [1.35, 1]})
        ax.plot(frequency, abs(actual), color='#2768aa', lw=3, label='Sonnet EM')
        ax.plot(frequency, abs(model), color='#db6d28', ls='--', lw=2.5,
                label='Equivalent circuit fit')
        ax.set(xlabel='Frequency (GHz)', ylabel=r'$|S_{21}|$',
               xlim=(frequency.min(), frequency.max()), ylim=(0, 1.06))
        ax.grid(alpha=.2)
        ax.spines[['top', 'right']].set_visible(False)
        ax.legend(loc='upper right', fontsize=15, framealpha=.95)
        panel.axis('off')
        panel.text(0, .96, 'Extracted Values and Targets', fontsize=21, weight='bold')
        panel.text(0, .90, 'Value (target; signed difference)', fontsize=15, color='#555555')
        for index, (label, value, target, unit) in enumerate(rows):
            y = .81 - index * .092
            panel.text(0, y, label, fontsize=16, weight='bold')
            panel.text(0, y - .040,
                       f'{value:.3f} {unit}  (target {target:g}; {(value / target - 1) * 100:+.2f}%)',
                       fontsize=15)
        panel.text(0, .105, f"Complex S21 RMS = {result['complex_RMS'][1]:.4f}  (limit 0.05)",
                   fontsize=16, weight='bold')
        panel.text(0, .035,
                   'Capacitances: joint complex S11, S21, S22 fit.\n'
                   'Inductances: independently extracted and held fixed.\n'
                   'Coarse equivalent values; geometry is not target-optimized.',
                   fontsize=12, color='#555555', va='top', linespacing=1.5)
        fig.suptitle('Taper 2 + Taper 1: Transmission Fit', fontsize=27, weight='bold')
        fig.tight_layout(rect=(0, .025, 1, .94), w_pad=3)
        for suffix in ('png', 'pdf'):
            fig.savefig(folder / f's21_fit_linear.{suffix}', dpi=200,
                        facecolor='white', bbox_inches='tight')
        plt.close(fig)
    print((folder / 's21_fit_linear.png').resolve())


if __name__ == '__main__':
    root = Path(__file__).resolve().parents[2]
    plot(root / 'output_data/taper_exploration/clearance_extraction_2026-10-25'
         '/frame700_wire88_left214/staged_joint_with_reflections')
