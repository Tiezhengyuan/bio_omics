'''
Test class 
'''
from ..helper import *
from src.bioomics import DownloadUniprot

@ddt
class TestDownloadUniprot(TestCase):

    def setUp(self):
        DIR_DATA = '/home/yuan/data/'
        self.c = DownloadUniprot(DIR_DATA, False)

    @skip
    def test_download_uniprot_sprot(self):
        res = self.c.uniprot_sprot()
        assert os.path.basename(res[0]) == 'uniprot_sprot.dat.gz'
        assert os.path.basename(res[0]) == 'uniprot_sprot.fasta.gz'
        assert os.path.basename(res[0]) == 'uniprot_sprot.xml.gz'

    @skip
    def test_download_uniprot_trembl(self):
        res = self.c.uniprot_trembl()
        assert res[0] == os.path.join(DIR_DATA, 'UniProt', 'uniprot_trembl.dat.gz')


