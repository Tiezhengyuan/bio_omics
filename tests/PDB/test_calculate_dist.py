
import os
from ddt import ddt, data, unpack
from unittest import TestCase, skip

from src import CalculateDist

TESTS_DIR = os.path.dirname(os.path.dirname(__file__))
DATA_DIR = os.path.join(TESTS_DIR, 'data', 'PDB')
RESULTS_DIR = os.path.join(TESTS_DIR, 'results')

@ddt
class TestCalculateDist(TestCase):

    def setUp(self):
        self.c = CalculateDist()
    
    def test_(self):
        pass