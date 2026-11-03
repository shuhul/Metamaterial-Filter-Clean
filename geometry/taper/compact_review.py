"""Summarize compact-stage measurements without launching Sonnet.

This is an investigation index, not a design acceptance decision. In particular,
an isolated inductance result cannot certify the assembled transmission response.
"""
import json
from datetime import datetime, timezone
from geometry.taper.compact_stages import OUT, NATIVE, GDS


def main():
    rows = []
    for folder in sorted(OUT.iterdir()):
        if not folder.is_dir() or not (GDS / folder.name).exists():
            continue
        row = {"case": folder.name, "accepted": False}
        meta = next((GDS / folder.name / stage / 'geometry.json' for stage in ('idc','idc_Lcal','Lcal','direct')
                     if (GDS / folder.name / stage / 'geometry.json').exists()), GDS / folder.name / 'missing.json')
        if meta.exists():
            row["parameters"] = json.loads(meta.read_text())["second"]
        row['calibrations'] = {}
        for stage in ('Lcal', 'idc_Lcal'):
            path = folder / stage / 'results.json'
            if path.exists():row['calibrations'][stage] = json.loads(path.read_text())
        calibration = folder / ('idc_Lcal' if 'idc_Lcal' in row['calibrations'] else 'Lcal') / "results.json"
        if calibration.exists():
            result = json.loads(calibration.read_text())
            row.update(L_nH=result["L0_nH"], L_error_percent=result["target_error_percent"],
                       calibration=str(calibration.resolve()))
        direct = folder / "direct"
        fits = sorted(direct.glob("fit*/fit.json")) if direct.exists() else []
        row["fits"] = [dict(path=str(path.resolve()), **json.loads(path.read_text())) for path in fits]
        project = NATIVE / folder.name / "direct" / "generated_ncell_3.sonx"
        row["direct_project"] = str(project.resolve()) if project.exists() else None
        idc_project = NATIVE / folder.name / 'idc/generated_ncell_3.sonx'
        row['idc_project'] = str(idc_project.resolve()) if idc_project.exists() else None
        row['user_rejected'] = (NATIVE / folder.name / 'REJECTED.md').exists()
        idc = folder / 'idc'
        row['idc_fits'] = [dict(path=str(path.resolve()), **json.loads(path.read_text()))
                           for path in sorted(idc.glob('fit*/fit.json'))] if idc.exists() else []
        idc_marker = idc / 'run_status.json'
        row['complete_IDC_response'] = idc_marker.exists() and json.loads(idc_marker.read_text()).get('status') == 'complete'
        idc_screen = idc / 'screen_completed.json'
        row['IDC_screen_record'] = json.loads(idc_screen.read_text()) if idc_screen.exists() else None
        component = folder / 'idc_check/results.json'
        row['independent_IDC_check'] = json.loads(component.read_text()) if component.exists() else None
        row["complete_full_response"] = (direct / "run_status.json").exists() and json.loads(
            (direct / "run_status.json").read_text()).get("status") == "complete"
        screen = direct / "screen_completed.json"
        row["screen_record"] = json.loads(screen.read_text()) if screen.exists() else None
        for key in ("two_point_diagnostic", "paired_two_point_diagnostic"):
            path = direct / (key + ".json")
            if path.exists():
                row[key] = json.loads(path.read_text())
        rows.append(row)
    report = dict(updated_utc=datetime.now(timezone.utc).isoformat(),
                  accepted_T2_design=None, cases=rows)
    path = OUT / "investigation_index.json"
    path.write_text(json.dumps(report, indent=2) + "\n")
    print(path.resolve())
    for row in rows:
        fits = row["fits"] + row['idc_fits']
        print(row["case"], "L=", row.get("L_nH"), "direct full=", row["complete_full_response"],
              "IDC full=", row['complete_IDC_response'], "rejected=", row['user_rejected'],
              "fit RMS=", [(q["S21_complex_RMSE"], q.get("Cg_T2_T1_fF",q.get('Cg_input_T2_fF'))) for q in fits])


if __name__ == "__main__":
    main()
