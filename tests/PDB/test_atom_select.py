'''
Unit tests for src/bioomics/PDB/atom_select.py

1I8I (tests/data/PDB/1I8I.pdb) contains ATOM records for chains A, B, C
and HETATM water residues (HOH, het flag 'W') on chain A.
'''
import os
import shutil
import tempfile

from Bio.PDB import PDBIO, PDBParser

from ..helper import *

from src.bioomics import (
    ChainAtomSelect, LigandAtomSelect, AllAtomSelect, NonHetSelect)

PDB_FILE = os.path.join(DATA_DIR, 'PDB', '1I8I.pdb')


def chain_ids(structure):
    return [chain.id for model in structure for chain in model]


def het_residues(structure):
    return [
        (chain.id, res)
        for model in structure for chain in model for res in chain
        if res.id[0] != ' '
    ]


class TestAtomSelect(TestCase):

    def setUp(self):
        self.tmp_dir = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.tmp_dir, ignore_errors=True)
        self.parser = PDBParser(QUIET=True)
        self.structure = self.parser.get_structure('1I8I', PDB_FILE)
        self.io = PDBIO()
        self.io.set_structure(self.structure)

    def _save_and_reload(self, select):
        outfile = os.path.join(self.tmp_dir, 'sel.pdb')
        self.io.save(outfile, select)
        return self.parser.get_structure('sel', outfile)

    def test_source_structure_layout(self):
        # sanity: chains A/B/C; all HETATM records are HOH waters
        self.assertEqual(chain_ids(self.structure), ['A', 'B', 'C'])
        water_chains = {chain for chain, _ in het_residues(self.structure)}
        self.assertEqual(water_chains, {'A', 'B', 'C'})
        res_names = {res.get_resname() for _, res in het_residues(self.structure)}
        self.assertEqual(res_names, {'HOH'})

    def test_chain_atom_select_filters_chain_and_hetero(self):
        select = ChainAtomSelect(['A'])
        selected = self._save_and_reload(select)
        self.assertEqual(chain_ids(selected), ['A'])
        # hetero water residues (het flag 'W') are dropped
        self.assertEqual(het_residues(selected), [])
        n_atoms = sum(1 for a in selected.get_atoms())
        self.assertGreater(n_atoms, 0)

    def test_chain_atom_select_accept_chain(self):
        select = ChainAtomSelect(['A'])
        for model in self.structure:
            chains = {chain.id: chain for chain in model}
            break
        self.assertTrue(select.accept_chain(chains['A']))
        self.assertFalse(select.accept_chain(chains['B']))
        self.assertFalse(select.accept_chain(chains['C']))

    def test_ligand_atom_select_keeps_hetero_only(self):
        select = LigandAtomSelect(['A'])
        selected = self._save_and_reload(select)
        self.assertEqual(chain_ids(selected), ['A'])
        # only water residues survive, all flagged as hetero
        residues = [res.id for res in list(selected)[0]['A']]
        self.assertGreater(len(residues), 0)
        self.assertTrue(all(res[0] != ' ' for res in residues))

    def test_all_atom_select_keeps_residues(self):
        # accept_residue keeps any residue with a (truthy) id
        select = AllAtomSelect(('A',))
        for model in self.structure:
            chain_a = list(model)[0]
            break
        residues = list(chain_a)
        self.assertTrue(all(select.accept_residue(res) for res in residues))

    def test_non_het_select_drops_waters(self):
        select = NonHetSelect()
        n_standard, n_waters = 0, 0
        for model in self.structure:
            for chain in model:
                for res in chain:
                    if select.accept_residue(res):
                        n_standard += 1
                    if res.id[0] != ' ':
                        n_waters += 1
        self.assertGreater(n_standard, 0)
        self.assertGreater(n_waters, 0)
        # every water is rejected
        for model in self.structure:
            for chain in model:
                for res in chain:
                    if res.id[0] != ' ':
                        self.assertFalse(select.accept_residue(res))
