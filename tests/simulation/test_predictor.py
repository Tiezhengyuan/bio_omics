'''
unit tests for Predictor: subprocess handling and command builders
- pytest -s tests/simulation/test_predictor.py
'''
import asyncio
from pathlib import Path
from unittest.mock import patch

from ..helper import *
from src.bioomics.simulation.predictor import Predictor


@ddt
class TestPredictor(TestCase):

    def test_run_subprocess_echo(self):
        res = Predictor.run_subprocess(['echo', 'hello world'],
                                       tool_name='echo')
        self.assertEqual(res.returncode, 0)
        self.assertIn('hello world', res.stdout)

    def test_run_subprocess_failure_does_not_raise(self):
        res = Predictor.run_subprocess(['ls', '/definitely/not/a/path'],
                                       tool_name='ls')
        self.assertNotEqual(res.returncode, 0)

    def test_haddock3_builds_run_command(self):
        # regression: haddock3() used to call the (missing) _subprocess
        # helper instead of Predictor.run_subprocess()
        with patch.object(Predictor, 'run_subprocess') as mocked:
            out = Predictor.haddock3({'cgi_file': 'restraints.cgi'})
        mocked.assert_called_once_with(
            ['haddock3', 'restraints.cgi'], tool_name='HadDock3')
        self.assertEqual(out, mocked.return_value)

    def test_haddock3_missing_cgi_key(self):
        with self.assertRaises(KeyError):
            Predictor.haddock3({})

    def test_colabfold_batch_command(self):
        params = {
            'use_templates': True, 'save_all': True, 'relax': True,
            'max_iterations': 200, 'model_type': 'alphafold2_ptm',
            'num_recycle': 5, 'colabfold_bin_dir': '/tmp/bin',
            'tools_dir': Path('/tmp'),
        }
        row = {'input_file': 'query.fasta', 'outdir': '/tmp/out'}
        with patch.object(Predictor, 'colabfold_subprocess') as mocked:
            Predictor.colabfold_batch(row, params)
        # colabfold_subprocess(params, cmd): the command is arg #2
        cmd = mocked.call_args[0][1]
        self.assertEqual(cmd, ['colabfold_batch',
            '--num-models', '5',
            '--model-order', '1,2,3,4,5',
            '--random-seed', '42',
            '--templates',
            '--save-all',
            '--amber', '--use-gpu-relax',
            '--relax-max-iterations', '200',
            '--model-type', 'alphafold2_ptm',
            '--num-recycle', '5',
            'query.fasta', '/tmp/out'])

    def test_colabfold_search_command(self):
        params = {'gpu': True, 'colabfold_db': '/databases/uniclust30'}
        row = {'query_faa': 'query.fasta', 'outdir': '/tmp/msa'}
        with patch.object(Predictor, 'colabfold_subprocess') as mocked:
            Predictor.colabfold_search(row, params)
        cmd = mocked.call_args[0][1]
        self.assertEqual(cmd, ['colabfold_search',
            '--threads', '8',
            '--use-env', '0',
            '--use-templates', '1',
            '--db1', 'uniref30_2302_db',
            '--db2', 'pdb100_230517',
            '--db3', 'colabfold_envdb_202108_db',
            '--gpu', '1',
            'query.fasta', '/databases/uniclust30', '/tmp/msa'])

    def test_run_command_async_success(self):
        ok = asyncio.run(Predictor.run_command(['echo', 'async ok'], 'echo'))
        self.assertTrue(ok)

    def test_run_command_async_failure(self):
        ok = asyncio.run(
            Predictor.run_command(['ls', '/definitely/not/a/path'], 'ls'))
        self.assertFalse(ok)
