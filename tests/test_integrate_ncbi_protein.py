'''
Test class ConnectNCBI
'''
from tests.helper import *
from src.bioomics import IntegrateNCBIProtein

@ddt
class TestIntegrateNCBIProtein(TestCase):

    def test_ncbi_protein(self):
        IntegrateNCBIProtein(DIR_DATA, False)()

    
