'''
Test class ConnectNCBI
'''
from tests.helper import *
from src.bioomics import ParseID

@ddt
class TestParseID(TestCase):

    def setUp(self) -> None:
        self.c = ParseID(nrows=10)

    @data(
        ['#NCBI_protein_accession', "AP_000001.1"],
        ["UniProtKB_protein_accession", "A0A9W3HR54"],
    )
    @unpack
    def test_gene_mapp(self, index_key, expect):
        infile = os.path.join(DIR_DATA, 'NCBI', 'gene', 'DATA', 'gene_refseq_uniprotkb_collab.gz')
        res = self.c.gene_map(infile, index_key)
        assert expect in res
    
