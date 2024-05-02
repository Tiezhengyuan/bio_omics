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
        [False, "AP_000001.1"],
        [True, "A0A9W3HR54"],
    )
    @unpack
    def test_parse_uniprotkb(self, by_uniprotkb, expect):
        infile = os.path.join(DIR_DATA, 'NCBI', 'gene', 'DATA', 'gene_refseq_uniprotkb_collab.gz')
        res = self.c.parse_uniprotkb(infile, by_uniprotkb)
        assert expect in res
    
