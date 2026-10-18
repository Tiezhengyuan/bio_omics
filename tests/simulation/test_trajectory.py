'''
unit tests for the Trajectory class (interface trajectory analysis)
- fixture: tests/data/PDB/1I8I_0_A-B.pdb
    chain A: 106 residues (trajectory chain), chain B: 119 residues
- pytest -s tests/simulation/test_trajectory.py
'''
import os

from ..helper import *
from src.bioomics.simulation.trajectory import Trajectory

TRAJ_PDB = os.path.join(DATA_DIR, 'PDB', '1I8I_0_A-B.pdb')


@ddt
class TestTrajectory(TestCase):

    def setUp(self):
        self.traj = Trajectory(TRAJ_PDB)

    def test_retrieve_traj(self):
        '''chain A of the trajectory file holds the trajectory residues'''
        residues = self.traj.retrieve_traj()
        self.assertEqual(len(residues), 106)
        # explicit chain id agrees with the default
        self.assertEqual(len(self.traj.retrieve_traj('A')), 106)

    @data(['A', 106], ['B', 119])
    @unpack
    def test_retriev_target(self, chain_ids, expect_residues):
        residues = self.traj.retriev_target(TRAJ_PDB, chain_ids)
        self.assertEqual(len(residues), expect_residues)

    def test_retriev_target_multi_chain(self):
        residues = self.traj.retriev_target(TRAJ_PDB, 'A,B')
        self.assertEqual(len(residues), 225)

    def test_ca_residues(self):
        residues = self.traj.retrieve_traj()
        cas = Trajectory.ca_residues(residues)
        self.assertEqual(len(cas), len(residues))
        self.assertEqual({a.get_name() for a in cas}, {'CA'})

    def test_target_pdb_rmsd_self_zero(self):
        '''the trajectory against itself has zero rmsd'''
        rmsd = self.traj.target_pdb_rmsd(TRAJ_PDB, 'A,B')
        self.assertEqual(rmsd, 0.0)

    def test_target_pdb_rmsd_different_conformation(self):
        '''chain A vs chain B: different conformations -> positive rmsd'''
        rmsd = self.traj.target_pdb_rmsd(TRAJ_PDB, 'B')
        self.assertGreater(rmsd, 1.0)
        self.assertLess(rmsd, 100.0)

    def test_hotspot_residues(self):
        '''binder chain B contacts 16 residues of trajectory chain A'''
        hs = self.traj.hotspot_residues(
            binder_chain='B', atom_distance_cutoff=4.0, traj_chain='A')
        self.assertEqual(len(hs), 16)
        self.assertEqual(sorted(hs.keys()),
            [337, 339, 343, 344, 345, 347, 359, 395, 400, 404,
             405, 406, 407, 408, 410, 411])
        # values are one-letter amino-acid codes
        for aa in hs.values():
            self.assertEqual(len(aa), 1)
            self.assertTrue(aa.isupper())

    def test_hotspot_residues_swapped_chains(self):
        '''regression: a traj_chain other than 'A' must not raise
        (unfold_entities() used to receive the chain id as entity level)'''
        hs = self.traj.hotspot_residues(
            binder_chain='A', atom_distance_cutoff=4.0, traj_chain='B')
        self.assertEqual(len(hs), 21)
