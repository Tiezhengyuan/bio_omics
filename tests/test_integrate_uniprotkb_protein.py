'''
Test class ConnectNCBI
'''
from tests.helper import *
from src.bioomics import IntegrateUniProtKBProtein

@ddt
class TestIntegrateUniProtKBProtein(TestCase):

    def test_uniprotkb_protein(self):
        IntegrateUniProtKBProtein(DIR_DATA, False)()
    
