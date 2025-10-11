'''
Test class 
'''
from ..helper import *
from src.bioomics import UniProt

@ddt
class TestUniProtEpitopes(TestCase):

    def setUp(self):
        self.c = UniProt(RESULTS_DIR)

    @skip
    def test_parse_uniprot_sprot_epitope(self):
        dat_gz = os.path.join(DIR_DATA, 'UniProt', 'uniprot_sprot.dat.gz')
        parser = self.c.parse_dat(dat_gz)
        res = self.c.parse_epitope(parser)

