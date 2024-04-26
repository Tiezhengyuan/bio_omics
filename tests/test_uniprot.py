'''
Test class 
'''
from .helper import *
from src.bioomics import UniProt


@ddt
class TestUniProt(TestCase):

    def setUp(self):
        self.c = UniProt(DIR_DATA, False)

    '''
    def test_download_uniprot_sprot_dat(self):
        res = self.c.download_uniprot_sprot_dat()
        assert res == os.path.join(DIR_DATA, 'UniProt', 'uniprot_sprot.dat.gz')


    def test_parse_dat(self):
        dat_gz = os.path.join(DIR_DATA, 'UniProt', 'uniprot_sprot.dat.gz')
        res = self.c.parse_dat(dat_gz)
        rec = next(res)
        assert res
        # print(rec, '\n\n')
        # print(dir(rec))
        # print(rec.dbxrefs, type(rec.dbxrefs))

    def test_parse_uniprot_sprot_epitope(self):
        dat_gz = os.path.join(DIR_DATA, 'UniProt', 'uniprot_sprot.dat.gz')
        parser = self.c.parse_dat(dat_gz)
        res = self.c.parse_epitope(parser)
    '''

    def test_download_uniprot_trembl_dat(self):
        res = self.c.download_uniprot_trembl_dat()
        assert res[0] == os.path.join(DIR_DATA, 'UniProt')
        assert res[1] == os.path.join(DIR_DATA, 'UniProt', 'uniprot_trembl.dat.gz')

        # print(rec, '\n\n')
        # print(dir(rec))
        # ft = rec.features[-1]
        # print(ft.location.start, ft.location.end,ft.id)
        # print(ft.qualifiers, type(ft.qualifiers))
