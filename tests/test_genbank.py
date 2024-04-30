'''
Test class ConnectNCBI
'''
from tests.helper import *
from src.bioomics import GenBank

@ddt
class TestGenBank(TestCase):


    def test_process_epitope(self):
        entity_path = os.path.join(DIR_DATA, 'epitope')
        GenBank(DIR_DATA).process_epitope(entity_path)

