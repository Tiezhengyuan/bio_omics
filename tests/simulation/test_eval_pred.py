'''
unit tests for EvalPred
- parse_prodigy / TMalign: pure logic + mocked subprocess, run everywhere
- rmsd (pyrosetta) and prodigy (prodigy binary): live checks
- pytest -s tests/simulation/test_eval_pred.py
'''
import os
from pathlib import Path
from unittest import skipUnless
from unittest.mock import patch

from ..helper import *
from src.bioomics.simulation.eval_pred import EvalPred

PDB_1RZI = os.path.join(DATA_DIR, 'PDB', '1RZI.pdb')
PDB_6BSX = os.path.join(DATA_DIR, 'PDB', '6BSX_0_A_0.pdb')
MMALIGN_DIR = Path(src_dir, 'bioomics', 'tools', 'external')
VENV_BIN = os.path.join(os.path.dirname(src_dir), '.venv', 'bin')
PRODIGY_BIN = os.path.join(VENV_BIN, 'prodigy')

# verbatim output of `prodigy tests/data/PDB/1RZI.pdb`
PRODIGY_OUTPUT = (
    '[+] No. of intermolecular contacts: 1681\n'
    '[+] No. of charged-charged contacts: 75.0\n'
    '[+] No. of charged-polar contacts: 230.0\n'
    '[+] No. of charged-apolar contacts: 200.0\n'
    '[+] No. of polar-polar contacts: 278.0\n'
    '[+] No. of apolar-polar contacts: 527.0\n'
    '[+] No. of apolar-apolar contacts: 371.0\n'
    '[+] Percentage of apolar NIS residues: 37.34\n'
    '[+] Percentage of charged NIS residues: 19.19\n'
    '[++] Predicted binding affinity (kcal.mol-1):    -98.5\n'
    '[++] Predicted dissociation constant (M) at 25.0C:  5.8e-73'
)


@ddt
class TestEvalPred(TestCase):

    def test_parse_prodigy_full(self):
        res = EvalPred().parse_prodigy(PRODIGY_OUTPUT)
        self.assertTrue(res['contacts'])
        self.assertEqual(res['intermolecular_contacts'], 1681.0)
        self.assertEqual(res['charged-charged_contacts'], 75.0)
        self.assertEqual(res['charged-polar_contacts'], 230.0)
        self.assertEqual(res['charged-apolar_contacts'], 200.0)
        self.assertEqual(res['polar-polar_contacts'], 278.0)
        self.assertEqual(res['apolar-polar_contacts'], 527.0)
        self.assertEqual(res['apolar-apolar_contacts'], 371.0)
        self.assertEqual(res['apolar_NIS_residues'], 37.34)
        self.assertEqual(res['charged_NIS_residues'], 19.19)
        self.assertEqual(res['binding_affinity'], -98.5)
        self.assertEqual(res['dissociation_constant'], 5.8e-73)

    def test_parse_prodigy_no_contacts(self):
        res = EvalPred().parse_prodigy('*** No contacts found ***')
        self.assertEqual(res, {'contacts': False})

    def test_parse_prodigy_empty(self):
        # no markers at all -> contacts still reported as True
        self.assertEqual(EvalPred().parse_prodigy(''), {'contacts': True})

    @skipUnless(os.path.exists(PRODIGY_BIN), 'prodigy binary not available')
    def test_prodigy_live(self):
        # prodigy() invokes the console script by name: make sure the
        # virtualenv's bin dir is on PATH
        if VENV_BIN not in os.environ['PATH'].split(os.pathsep):
            os.environ['PATH'] = VENV_BIN + os.pathsep + os.environ['PATH']
        res = EvalPred().prodigy(Path(PDB_1RZI))
        self.assertIsNotNone(res)
        self.assertEqual(res['intermolecular_contacts'], 1681.0)
        self.assertEqual(res['binding_affinity'], -98.5)
        self.assertEqual(res['dissociation_constant'], 5.8e-73)

    def test_prodigy_failing_call_returns_none(self):
        self.assertIsNone(EvalPred().prodigy(Path('/not/a/pdb.pdb')))

    def test_tmalign_without_binary(self):
        # bin_dir=None -> 'None/MMalign' cannot run -> caught -> {}
        self.assertEqual(EvalPred().TMalign('a.pdb', 'b.pdb'), {})

    @skipUnless((MMALIGN_DIR / 'MMalign').exists(), 'MMalign not available')
    def test_tmalign_live(self):
        res = EvalPred(MMALIGN_DIR).TMalign(PDB_1RZI, PDB_6BSX)
        self.assertNotEqual(res, {})
        self.assertIn('TM-score', res)
        self.assertIn('RMSD', res)
        self.assertGreater(float(res['TM-score']), 0)

    def test_rmsd_self(self):
        res = EvalPred().rmsd(PDB_6BSX, PDB_6BSX)
        self.assertIsNotNone(res)
        self.assertAlmostEqual(res, 0.0)

    def test_rmsd_missing_file(self):
        self.assertIsNone(EvalPred().rmsd('/no/such.pdb', '/no/such.pdb'))

    def test_batch_rmsd(self):
        pool = [{'pdb1': PDB_6BSX, 'pdb2': PDB_6BSX}]
        rows = list(EvalPred().batch_rmsd(pool))
        self.assertEqual(len(rows), 1)
        self.assertAlmostEqual(rows[0]['rmsd'], 0.0)
