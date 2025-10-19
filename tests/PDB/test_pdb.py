import os
from unittest import TestCase
from ddt import ddt, data, unpack
from pprint import pprint

from src.bioomics import Pdb

TESTS_DIR = os.path.dirname(os.path.dirname(__file__))
DATA_DIR = os.path.join(TESTS_DIR, 'data', 'PDB')
RESULTS_DIR = os.path.join(TESTS_DIR, 'results')

@ddt
class TestParsePdb(TestCase):

    def setUp(self):
        self.c = {}
        for file_name in ['7S8L.pdb',]:
            pdb_id = os.path.splitext(file_name)[0]
            pdb_file = os.path.join(DATA_DIR, file_name)
            self.c[pdb_id] = Pdb(pdb_file, RESULTS_DIR, True)
    
    @data(
        ['7S8L', 'succeed'],
    )
    @unpack
    def test_export_fasta(self, pdb_id, expect):
        fa_file = os.path.join(RESULTS_DIR, pdb_id + '.fa')
        result = self.c[pdb_id].export_fasta(fa_file)
        assert result == expect

    @data(
        ['7S8L',],
    )
    @unpack
    def test_export_chains_fasta(self, pdb_id):
        result = self.c[pdb_id].export_chains_fasta()

    @data(
        ['7S8L', ['R',], os.path.join(RESULTS_DIR, '7S8L_R.pdb')],
        ['7S8L', ['E', 'C'], os.path.join(RESULTS_DIR, '7S8L_E_C.pdb')],
    )
    @unpack
    def test_export_chain_pdb(self, pdb_id, chain_ids, expect):
        result = self.c[pdb_id].export_chain_pdb(chain_ids, overwrite=True)
        assert result == expect

    @data(
        ['7S8L', '',],
    )
    @unpack
    def test_split_by_chain(self, pdb_id, expect):
        result = self.c[pdb_id].split_by_chain(overwrite=True)
        # print(result)


    '''
    @data(
        ['7S8L', '',],
    )
    @unpack
    def test_(self, pdb_id, expect):
        result = self.c[pdb_id].export_chain_pdb(chain_ids, overwrite=True)
        assert result == expect
    '''