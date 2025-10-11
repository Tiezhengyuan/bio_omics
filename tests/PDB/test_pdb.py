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
            self.c[pdb_id] = Pdb(pdb_file)
    
    @data(
        ['7S8L', 6],
    )
    @unpack
    def test_export_fasta(self, pdb_id, expect):
        file_name = pdb_id + '.fa'
        result = self.c[pdb_id].export_fasta(file_name, RESULTS_DIR)
        outfile = os.path.join(RESULTS_DIR, file_name)
        assert result[0] == outfile
        assert len(result[1]) == expect

    @data(
        ['7S8L', '/home/yuan/bio/bio_omics/tests/results/7S8L_pdb.pickle', '7S8L'],
    )
    @unpack
    def test_pdb_obj(self, pdb_id, expect_outfile, structure_id):
        pdb_file = os.path.join(DATA_DIR, f'{pdb_id}.pdb')
        pfile = Pdb.build_pdb_obj(pdb_file, RESULTS_DIR)
        assert pfile == expect_outfile
        result = Pdb.load_pdb_obj(pfile)
        assert result.structure_id == structure_id


    '''
    @data(
        ['7S8L', '',],
    )
    @unpack
    def test_(self, input, expect):

        assert result == expect
    '''