'''
Test class ConnectNCBI
'''
from tests.helper import *
from src.bioomics import ParseID

from unittest import skipUnless

# the collab mapping dump must be downloaded from NCBI FTP (multi-GB set) -
# only run when the file is already available locally
COLLAB_FILE = os.path.join(DIR_DATA, 'NCBI', 'gene', 'DATA',
    'gene_refseq_uniprotkb_collab.gz')

@ddt
class TestParseID(TestCase):

    def setUp(self) -> None:
        self.c = ParseID(nrows=10)

    @skipUnless(os.path.isfile(COLLAB_FILE),
        'gene_refseq_uniprotkb_collab.gz not available locally')
    @data(
        ['#NCBI_protein_accession', "AP_000001.1"],
        ["UniProtKB_protein_accession", "A0A9W3HR54"],
    )
    @unpack
    def test_gene_mapp(self, index_key, expect):
        res = self.c.gene_map(COLLAB_FILE, index_key)
        assert expect in res

