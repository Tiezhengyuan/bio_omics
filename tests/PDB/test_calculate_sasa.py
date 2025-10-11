
import os
import pickle
from ddt import ddt, data, unpack
from unittest import TestCase, skip

from src.bioomics import CalculateSasa

TESTS_DIR = os.path.dirname(os.path.dirname(__file__))
DATA_DIR = os.path.join(TESTS_DIR, 'data', 'PDB')
RESULTS_DIR = os.path.join(TESTS_DIR, 'results')

@ddt
class TestCalculateSasa(TestCase):

    
    @data(
        ['7S8L', (1096, 20)],
    )
    @unpack
    def test_freesasa(self, pdb_id, expect):
        pdb_file = os.path.join(DATA_DIR, f'{pdb_id}.pdb')
        c = CalculateSasa(pdb_file)
        result = c.cal_freesasa(pdb_file)
        assert result.shape == expect

    @data(
        ['7S8L', 'A', (1096, 20)],
    )
    @unpack
    def test_freesasa(self, pdb_id, chain_id, expect):
        pdb_file = os.path.join(DATA_DIR, f'{pdb_id}.pdb')
        c = CalculateSasa(pdb_file)
        result = c.cal_sasa_df(0, chain_id)
        print(result)
        # assert result.shape == expect