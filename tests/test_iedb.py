'''
Test class 
'''
from tests.helper import *

from src.bioomics import IEDB

@ddt
class TestIEDB(TestCase):
    def setUp(self):
        self.conn = IEDB(DIR_DATA, None, True)

    @skip
    @data(
        ['antigen', 'antigen_full_v3.zip'],
    )
    @unpack
    def test_csv(self, itype, file_name):
        res = self.conn.download_csv(itype)
        assert res == os.path.join(DIR_DATA, file_name)


    def test_pull_csv(self):
        types = ['antigen', 'epitope', 'tcell', 'bcell', 'mhc', 
            'reference', 'receptor', 'tcr', 'bcr', 'iedb']
        for i in types:
            self.conn.download_csv(i)
    @skip
    def test_pull(self):
        res = IEDB(DIR_DATA, None, False).pull('mhc')
        assert 'records' in res

# https://www.iedb.org/downloader.php?file_name=doc/epitope_full_v3_json.zip