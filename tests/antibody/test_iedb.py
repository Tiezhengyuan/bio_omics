'''
Test class 
'''
from tests.helper import *

from src.bioomics import IEDB

@ddt
class TestIEDB(TestCase):
    def setUp(self):
        self.conn = IEDB(DIR_DATA, None, True)

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

    def test_pull(self):
        res = IEDB(DIR_DATA, None, False).pull('mhc')
        assert 'records' in res

    def test_antigen_csv(self):
        res = self.conn.antigen_csv()
        assert len(res) == 79546

    def test_antigen_json(self):
        res = self.conn.antigen_json()
        assert len(res) == 79546

    def test_epitope_csv(self):
        res = self.conn.epitope_csv()
        assert len(res) == 2226136
    
    def test_epitope_json(self):
        res = self.conn.epitope_json()
        assert len(res) == 2226136
        
    def test_tcr_csv(self):
        res = self.conn.tcr_csv()
        assert len(res) == 204754

    def test_tcr_json(self):
        res = self.conn.tcr_json()
        assert len(res) == 204754
    
    def test_bcr_csv(self):
        res = self.conn.bcr_csv()
        assert len(res) == 5127
    
    def test_bcr_json(self):
        res = self.conn.bcr_json()
        assert len(res) == 5127
    
    def test_reference_csv(self):
        res = self.conn.reference_csv()
        assert len(res) == 24687
    
    def test_reference_json(self):
        res = self.conn.reference_json()
        assert len(res) == 24687

    def test_tcell_csv(self):
        res = self.conn.tcell_csv()
        assert len(res) == 512545

    def test_tcell_json(self):
        res = self.conn.tcell_json()
        assert len(res) == 512545

    def test_bcell_csv(self):
        res = self.conn.bcell_csv()
        assert len(res) == 1403929

    def test_bcell_json(self):
        res = self.conn.bcell_json()
        assert len(res) == 1403929

    def test_mhc_csv(self):
        res = self.conn.mhc_csv()
        assert len(res) == 4803074
        
    def test_mhc_json(self):
        res = self.conn.mhc_json()
        assert len(res) == 4803074

