'''
Test class IntegrateUniProtKBProtein
'''
from ..helper import *
from src.bioomics import IntegrateUniProtKBProtein

@ddt
class TestIntegrateUniProtKBProtein(TestCase):

    @skip('heavy: integrates the full local UniProtKB/NCBI data')
    def test_uniprotkb_protein(self):
        IntegrateUniProtKBProtein(DIR_DATA, False)()
    
