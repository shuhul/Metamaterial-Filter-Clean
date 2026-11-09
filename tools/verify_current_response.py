"""Replay the published taper evidence without launching Sonnet or changing files."""
from pathlib import Path
import json
import shutil
import tempfile
import numpy as np
import skrf as rf

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / 'final_designs/current_2026_11_08/taper'


def ladder(f, c, l, delay_ps=0):
    w = 2 * np.pi * f
    m = np.tile(np.eye(2, dtype=complex), (len(f), 1, 1))
    for k in range(5):
        z = 1 / (1j * w * c[5 + k] * 1e-15)
        y = 1j * w * c[k] * 1e-15 + 1 / (1j * w * l[k] * 1e-9)
        a = np.tile(np.eye(2, dtype=complex), (len(f), 1, 1))
        a[:, 0, 0] += z * y
        a[:, 0, 1] = z
        a[:, 1, 0] = y
        m = m @ a
    return 2 * np.exp(-1j * w * delay_ps * 1e-12) / (
        m[:, 0, 0] + m[:, 0, 1] / 50 + 50 * m[:, 1, 0] + m[:, 1, 1])


def main():
    result = json.loads((BASE / 'result.json').read_text())
    data = json.loads((BASE / 'local_model/network.json').read_text())
    with tempfile.TemporaryDirectory() as tmp:
        p = Path(tmp) / 'response.s2p'
        shutil.copyfile(BASE / 'project.snp', p)
        native = rf.Network(str(p))
    co = np.array(data['coefficients']['real']) + 1j * np.array(data['coefficients']['imag'])
    z = ((native.f / 5e9)[:, None] ** np.array(data['powers']) @ co).reshape(-1, 7, 7)
    ze = z[:, :2, :2] - z[:, :2, 2:] @ np.linalg.solve(z[:, 2:, 2:], z[:, 2:, :2])
    model = ((ze - 50 * np.eye(2)) @ np.linalg.inv(ze + 50 * np.eye(2)))[:, 1, 0]
    target = ladder(native.f, result['normal_target_C_fF'], result['target_L_nH'], result['total_reference_delay_ps'])
    rms = lambda a, b: float(np.sqrt(np.mean(abs(a - b) ** 2)))
    measured = dict(native_to_comparator_RMS=rms(native.s[:, 1, 0], target),
                    native_to_terminal_model_RMS=rms(native.s[:, 1, 0], model),
                    terminal_model_to_same_comparator_RMS=rms(model, target),
                    samples=len(native.f), EM_started=False)
    assert abs(measured['native_to_comparator_RMS'] - result['actual_target_RMS']) < 1e-9
    assert abs(measured['native_to_terminal_model_RMS'] - result['model_S21_RMS']) < 1e-9
    assert measured['native_to_comparator_RMS'] < .0025
    for value, nominal in [('normal_target_C_fF', 'nominal_C_fF'), ('target_L_nH', 'nominal_L_nH')]:
        assert np.max(abs(np.array(result[value]) / np.array(result[nominal]) - 1)) <= .0025 + 1e-12
    print(json.dumps(measured, indent=2))


if __name__ == '__main__':
    main()
