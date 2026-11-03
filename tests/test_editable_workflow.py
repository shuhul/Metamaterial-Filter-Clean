"""Offline checks for the editable entry point; no native EM execution."""
from dataclasses import replace
from pathlib import Path
import contextlib
import io
import json
import shutil
import tempfile
import unittest
from unittest.mock import patch
import xml.etree.ElementTree as ET

from geometry import run_filter as workflow
from geometry import prepare_inductance
from geometry.prepare_sonnet import tag

BASE=Path(__file__).resolve().parents[1]
SELECTED=BASE/'sonnet_files/generated_ncell_4_w186_d56/generated_ncell_4.sonx'
SNP=SELECTED.with_suffix('.snp')


class EditableWorkflowTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory(prefix='filter_workflow_test_')
        self.root=Path(self.temp.name)
        self.patch=patch.multiple(workflow,RUN_NAME='case',GEOMETRY_ROOT=self.root/'geometry',
            PROJECT_ROOT=self.root/'projects',RESULT_ROOT=self.root/'results',
            GENERATE_GEOMETRY=False,PREPARE_SONNET=False,SIMULATE_INDUCTANCE=False,
            SIMULATE_ARRAY=False,FIT_S21=True,INDUCTANCE_RESULT_FILE=None)
        self.patch.start()
        self.output=contextlib.redirect_stdout(io.StringIO());self.output.__enter__()
    def tearDown(self):
        self.output.__exit__(None,None,None)
        self.patch.stop();self.temp.cleanup()
    def fixtures(self):
        stem=workflow.GEOMETRY_ROOT/'case';stem.parent.mkdir(parents=True)
        with patch.object(workflow,'save_preview'):
            workflow.generate(workflow.PARAMS,stem)
        folder=workflow.PROJECT_ROOT/'case'/'array';folder.mkdir(parents=True)
        project=folder/'generated_ncell_4.sonx'
        shutil.copyfile(SELECTED,project);shutil.copyfile(SNP,project.with_suffix('.snp'))
        workflow.save_json(folder/'simulation_record.json',{
            **workflow.input_identity(project,stem),'snp_sha256':workflow.digest(project.with_suffix('.snp')),
            'origin':'offline fixture copied from the previously completed selected run'})
        return project,stem
    def test_geometry_only_never_calls_native_tools_or_fit(self):
        with patch.multiple(workflow,GENERATE_GEOMETRY=True,FIT_S21=False), \
             patch.object(workflow,'save_preview'), patch.object(workflow,'prepare_sonnet') as native, \
             patch.object(workflow,'simulate_array') as sim,patch.object(workflow,'extract_inductance') as sweep:
            workflow.main()
        native.assert_not_called();sim.assert_not_called();sweep.assert_not_called()
        meta=json.loads((workflow.GEOMETRY_ROOT/'case.json').read_text())
        self.assertEqual(meta['parameters']['unit']['meander']['w'],186)
        self.assertEqual(meta['parameters']['d_cell'],56)
        self.assertEqual(json.loads((workflow.RESULT_ROOT/'case'/'run_status.json').read_text())['state'],'complete')
    def test_refit_saved_data_reproduces_selected_capacitances_without_simulation(self):
        self.fixtures()
        with patch.object(workflow,'simulate_array') as sim,patch.object(workflow,'extract_inductance') as sweep:
            workflow.main()
        sim.assert_not_called();sweep.assert_not_called()
        output=workflow.RESULT_ROOT/'case'
        fit=json.loads((output/'fit_parameters.json').read_text())
        self.assertAlmostEqual(fit['Cg_fF'],8.241620220231225,places=5)
        self.assertTrue((output/'s21_fit.png').exists())
    def test_stale_geometry_is_rejected_before_fitting(self):
        self.fixtures()
        with patch.object(workflow,'PARAMS',replace(workflow.PARAMS,d_cell=54)),patch.object(workflow,'fit_array') as fit:
            with self.assertRaisesRegex(ValueError,'differs from PARAMS'):workflow.main()
        fit.assert_not_called()
    def test_changed_project_or_export_is_rejected_before_fitting(self):
        project,_=self.fixtures()
        with project.with_suffix('.snp').open('a') as stream:stream.write('\n! changed\n')
        with patch.object(workflow,'fit_array') as fit:
            with self.assertRaisesRegex(ValueError,'S-parameter file changed'):workflow.main()
        fit.assert_not_called()
    def test_contradictory_saved_fit_flags_rejected(self):
        with patch.object(workflow,'GENERATE_GEOMETRY',True):
            with self.assertRaisesRegex(ValueError,'saved-data fitting'):workflow.check_settings()
    def test_existing_geometry_is_not_overwritten(self):
        self.fixtures()
        with patch.multiple(workflow,GENERATE_GEOMETRY=True,PREPARE_SONNET=True,SIMULATE_ARRAY=True):
            with self.assertRaisesRegex(ValueError,'already has geometry'):workflow.main()
    def test_full_orchestration_uses_new_extracted_L0(self):
        def native(n_cells,*,gds_path,output_dir):
            output_dir.mkdir(parents=True)
            shutil.copyfile(SELECTED,output_dir/f'generated_ncell_{n_cells}.sonx')
        def capacitor(**kwargs):
            shutil.copyfile(kwargs['project_path'],kwargs['output_path'])
        def extract(**kwargs):
            kwargs['output_dir'].mkdir(parents=True)
            self.assertFalse(kwargs['fit_capacitance'])
            self.assertEqual(tuple(kwargs['c_add_sweep_pf']),(.02,.03,.04,.05))
            self.assertIn('inductance',str(kwargs['project_file']))
            return {'L0_nH':4.93,'R_squared':.99999}
        def simulate(project,stem,analysis):
            shutil.copyfile(SNP,project.with_suffix('.snp'))
            workflow.save_json(project.parent/'simulation_record.json',{
                **workflow.input_identity(project,stem),'snp_sha256':workflow.digest(project.with_suffix('.snp'))})
        def fit(args):
            self.assertEqual(float(args[args.index('--l0-nh')+1]),4.93)
            self.assertIn('--no-simulate',args)
            folder=Path(args[args.index('--output-dir')+1])
            workflow.save_json(folder/'fit_parameters.json',{'L0_nH_fixed':4.93,'C0_fF':188.8,'Cg_fF':8.2})
        with patch.multiple(workflow,GENERATE_GEOMETRY=True,PREPARE_SONNET=True,
                            SIMULATE_INDUCTANCE=True,SIMULATE_ARRAY=True), \
             patch.object(workflow,'save_preview'),patch.object(workflow,'prepare_sonnet',side_effect=native), \
             patch.object(workflow,'prepare_inductance',side_effect=capacitor), \
             patch.object(workflow,'check_project_geometry'), \
             patch.object(workflow,'extract_inductance',side_effect=extract) as sweep, \
             patch.object(workflow,'simulate_array',side_effect=simulate) as sim, \
             patch.object(workflow,'fit_array',side_effect=fit):
            workflow.main()
        self.assertEqual(sweep.call_count,1);self.assertEqual(sim.call_count,1)
        output=workflow.RESULT_ROOT/'case'
        self.assertEqual(json.loads((output/'fit_parameters.json').read_text())['inductance_provenance']['type'],'new_extraction')
        measured,_=workflow.saved_inductance();self.assertEqual(measured,4.93)
    def test_default_inductance_extractor_preserves_native_tuner_override(self):
        import extract_generated_inductance as extraction
        gold=BASE/'output_data/meander_tuning_ref400/best_refined'
        out=self.root/'sweep'
        def cached(value,label):
            shutil.copyfile(gold/f'{label}.sonx',out/f'{label}.sonx')
            shutil.copyfile(gold/f'{label}.snp',out/f'{label}.snp')
            net=extraction.read_snp(out/f'{label}.snp')
            return net,extraction.single_sweep_resonance(net)
        with patch.multiple(extraction,OUTPUT=out,PROJECT=self.root/'unused.sonx'), \
             patch.object(extraction,'simulate',side_effect=cached):
            result=extraction.main(fit_capacitance=False)
        self.assertAlmostEqual(result['L0_nH'],4.922723945068658,places=7)
    def test_forced_fixed_inductance_ignores_saved_source(self):
        with patch.multiple(workflow,USE_SAVED_INDUCTANCE=False,FIXED_L0_NH=4.9,
                            INDUCTANCE_RESULT_FILE=Path('does_not_exist.json')):
            value,record=workflow.saved_inductance()
        self.assertEqual(value,4.9);self.assertEqual(record['type'],'user_fixed')
    def test_generalized_component_preparation_uses_supplied_paths_and_sweep(self):
        folder=self.root/'single';folder.mkdir()
        source=BASE/'sonnet_files/generated_ncell_1/generated_ncell_1.sonx'
        project=folder/'base.sonx';shutil.copyfile(source,project)
        target=folder/'custom_inductance.sonx'
        with patch.object(prepare_inductance,'run') as native:
            prepare_inductance.main(gds_path=BASE/'geometry/output/ncell_1.gds',project_path=project,
                                    output_path=target,sweep=(4.35,4.95,1500))
        self.assertEqual(native.call_count,1)
        tree=ET.parse(target)
        sweep=next(s for ss in tree.findall(f'{tag("Sweeps")}/{tag("Set")}') if ss.get('On')=='TRUE' for s in ss.findall('.//'+tag('Sweep')))
        self.assertEqual(sweep.get('Target'),'1500')
        self.assertEqual(sweep.get('Start'),'4.35')
        export=tree.find('.//'+tag('OutputFiles')+'//'+tag('File'))
        self.assertEqual(export.get('Name'),'custom_inductance.snp')
        self.assertEqual(len(tree.findall('.//'+tag('Component'))),1)


if __name__=='__main__':unittest.main()
