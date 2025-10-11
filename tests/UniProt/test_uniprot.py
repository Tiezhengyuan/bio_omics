'''
Test class 
'''
from ..helper import *
from src.bioomics import Uniprot

@ddt
class TestUniprot(TestCase):

    @data(
        ['uniprot_sprot.dat', 'Q6GZX4'],
        ['uniprot_sprot.dat.gz', 'Q6GZX4'],
        ['uniprot_trembl.dat', 'A0A7C4TNZ0'],
        ['uniprot_trembl.dat.gz', 'A0A7C4TNZ0'],
    )
    @unpack
    def test_dat_iter(self, file_name, expect):
        dat_file = os.path.join(DATA_DIR, 'UniProt', file_name)
        res = Uniprot(dat_file).dat_iter()
        rec = next(res)
        assert rec.id == expect

