'''
Test class IEDB
'''
from tests.helper import *

from src.bioomics import IEDB

@ddt
class TestIEDB(TestCase):
    def setUp(self):
        # overwrite=False: reuse already downloaded IEDB dumps
        self.conn = IEDB(DIR_DATA, False)

    @data(
        ['antigen', 'antigen_full_v3.zip'],
    )
    @unpack
    def test_csv(self, itype, file_name):
        res = self.conn.download_csv(itype)
        # dumps are stored in DIR_DATA/IEDB
        assert res == os.path.join(DIR_DATA, 'IEDB', file_name)

    # NOTE: the full IEDB dumps (epitope/tcell/bcell/mhc/...) are
    # multi-GB downloads - the tests below are network-heavy and skipped.
    @skip
    def test_pull_csv(self):
        types = ['antigen', 'epitope', 'tcell', 'bcell', 'mhc',
            'reference', 'receptor', 'tcr', 'bcr', 'iedb']
        for i in types:
            self.conn.download_csv(i)

    @skip
    def test_pull(self):
        res = IEDB(DIR_DATA, False).pull('mhc')
        assert 'records' in res

    def test_antigen_csv(self):
        res = self.conn.antigen_csv()
        # current IEDB antigen_full_v3 dump size
        assert len(res) == 125129

    def test_antigen_json(self):
        res = self.conn.antigen_json()
        assert len(res) == 125129

    @skip
    def test_epitope_csv(self):
        res = self.conn.epitope_csv()
        assert len(res) == 2226136

    @skip
    def test_epitope_json(self):
        res = self.conn.epitope_json()
        assert len(res) == 2226136

    @skip
    def test_tcr_csv(self):
        res = self.conn.tcr_csv()
        assert len(res) == 204754

    @skip
    def test_tcr_json(self):
        res = self.conn.tcr_json()
        assert len(res) == 204754

    @skip
    def test_bcr_csv(self):
        res = self.conn.bcr_csv()
        assert len(res) == 5127

    @skip
    def test_bcr_json(self):
        res = self.conn.bcr_json()
        assert len(res) == 5127

    @skip
    def test_reference_csv(self):
        res = self.conn.reference_csv()
        assert len(res) == 24687

    @skip
    def test_reference_json(self):
        res = self.conn.reference_json()
        assert len(res) == 24687

    @skip
    def test_tcell_csv(self):
        res = self.conn.tcell_csv()
        assert len(res) == 512545

    @skip
    def test_tcell_json(self):
        res = self.conn.tcell_json()
        assert len(res) == 512545

    @skip
    def test_bcell_csv(self):
        res = self.conn.bcell_csv()
        assert len(res) == 1403929

    @skip
    def test_bcell_json(self):
        res = self.conn.bcell_json()
        assert len(res) == 1403929

    @skip
    def test_mhc_csv(self):
        res = self.conn.mhc_csv()
        assert len(res) == 4803074

    @skip
    def test_mhc_json(self):
        res = self.conn.mhc_json()
        assert len(res) == 4803074
