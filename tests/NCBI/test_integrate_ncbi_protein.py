'''
Test class ConnectNCBI
'''
from tests.helper import *
from src.bioomics import IntegrateNCBIProtein

@ddt
class TestIntegrateNCBIProtein(TestCase):

    # downloads gene_refseq_uniprotkb_collab.gz and integrates millions of
    # accession rows - network-heavy, skipped
    @skip
    def test_ncbi_protein(self):
        IntegrateNCBIProtein(DIR_DATA, False)()

    
