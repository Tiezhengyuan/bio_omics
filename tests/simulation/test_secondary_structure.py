'''
unit tests for SecondaryStructure
- DSSP-free parts (clash score, percentages) run everywhere
- DSSP-dependent parts skip when no mkdssp / dssp binary is installed
- fixture: tests/data/PDB/1I8I_0_A-B.pdb (chains A: 106 res, B: 119 res)
- pytest -s tests/simulation/test_secondary_structure.py
'''
import os
import shutil
from unittest import skipUnless

from ..helper import *
from src.bioomics.simulation.secondary_structure import SecondaryStructure

PDB_COMPLEX = os.path.join(DATA_DIR, 'PDB', '1I8I_0_A-B.pdb')
PDB_SINGLE = os.path.join(DATA_DIR, 'PDB', '6BSX_0_A_0.pdb')

DSSP_BIN = shutil.which('mkdssp') or shutil.which('dssp')


@ddt
class TestSecondaryStructure(TestCase):

    def setUp(self):
        self.ss = SecondaryStructure(PDB_COMPLEX)

    def test_initial_state(self):
        self.assertEqual(dict(self.ss.ss_counts), {})
        self.assertEqual(self.ss.i_plddt(), 0)
        self.assertEqual(self.ss.ss_plddt(), 0)
        self.assertEqual(self.ss.total_percentages(), (0.0, 0.0, 0.0))

    def test_calc_clash_score_ca_only(self):
        # no CA atoms of the two chains come within 2.4 A
        self.assertEqual(
            self.ss.calc_clash_score(threshold=2.4, only_ca=True), 0)

    def test_calc_clash_score_state_reset(self):
        # repeated calls must not accumulate neighbours from previous runs
        first = self.ss.calc_clash_score(threshold=2.4, only_ca=True)
        second = self.ss.calc_clash_score(threshold=2.4, only_ca=True)
        self.assertEqual(first, second)

    def test_calc_clash_score_all_atom(self):
        # the two-chain complex has one heavy-atom contact below 2.4 A
        self.assertEqual(
            self.ss.calc_clash_score(threshold=2.4, only_ca=False), 1)

    def test_calc_clash_score_single_chain(self):
        # a single-chain structure has no inter-chain contacts
        single = SecondaryStructure(PDB_SINGLE)
        self.assertEqual(single.calc_clash_score(only_ca=False), 0)

    @data(
        [[100, 25, 50], (25.0, 50.0, 25.0)],
        [[200, 30, 70], (15.0, 35.0, 50.0)],
        [[0, 0, 0], (0.0, 0.0, 0.0)],
    )
    @unpack
    def test_calculate_percentages(self, counts, expect):
        self.assertEqual(SecondaryStructure.calculate_percentages(*counts),
                         expect)

    @skipUnless(DSSP_BIN, 'mkdssp / dssp binary not available')
    def test_count_residues(self):
        self.ss.count_residues(DSSP_BIN, chain_id='B', atom_distance_cutoff=4.0)
        total = sum(self.ss.ss_counts.values())
        self.assertEqual(total, 119)    # all residues of chain B
        helix, sheet, loop = self.ss.total_percentages()
        self.assertAlmostEqual(helix + sheet + loop, 100.0)
        # interface residues must be a subset of the chain residues
        interface_total = sum(self.ss.interface_ss_counts.values())
        self.assertLessEqual(interface_total, total)
        self.assertGreaterEqual(self.ss.i_plddt(), 0)
        self.assertLessEqual(self.ss.ss_plddt(), 1)
