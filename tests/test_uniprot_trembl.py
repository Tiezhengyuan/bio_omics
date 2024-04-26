'''
Test class 
'''
from .helper import *
from src.bioomics import UniProtTrembl


@ddt
class TestUniProtTrembl(TestCase):

    def setUp(self):
        self.c = UniProtTrembl(DIR_DATA, False)

    def test_parse_epitope(self):
        dat_gz = self.c.download_dat()
        assert dat_gz == os.path.join(DIR_DATA, 'UniProt', 'uniprot_trembl.dat.gz')

        parser = self.c.parse_dat(dat_gz)
        assert parser

        entity_data = self.c.parse_epitope(parser)
        assert len(entity_data) > 0
        

