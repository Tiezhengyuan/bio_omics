'''
Test class ConnectNCBI
'''
from tests.helper import *
from src.bioomics import IntegrateID

@ddt
class TestIntegrateID(TestCase):

    def test_ncbi_protein(self):
        IntegrateID(DIR_DATA).ncbi_protein()

    def test_uniprotkb_protein(self):
        IntegrateID(DIR_DATA).uniprotkb_protein()
    
