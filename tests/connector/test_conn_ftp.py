'''
Test class 
'''
from tests.helper import *

from src.bioomics import ConnFTP

@ddt
class TestConnFTP(TestCase):

    def setUp(self):
        self.c = ConnFTP('ftp.ncbi.nlm.nih.gov')

    @data(
        ['blast', True],
        ['robts.txt', False]
    )
    @unpack
    def test_is_dir(self, name, expect):
        ftp = self.c.connect()
        res = self.c._is_dir(ftp, name)

    @data(
        [
            'geo/series/GSE3nnn/GSE3341/matrix',
            [('geo/series/GSE3nnn/GSE3341/matrix', 'GSE3341_series_matrix.txt.gz')],
        ],
    )
    @unpack
    def test_list_files(self, path, expect):
        res = self.c.list_files(path)
        assert res == expect

    @data(
        ['geo/series/GSE3nnn/GSE3341/matrix', '.gz', 1],
        ['geo/series/GSE3nnn/GSE3341/', '.gz', 0],
        ['geo/series/GSE3nnn/GSE3341/', None, 0],
    )
    @unpack
    def test_download_files(self, path, pattern, expect):
        res = self.c.download_files(path, pattern, DIR_DATA)
        assert len(res) == expect


    @data(
        ['gene/tools/', 'README', None, 'README'],
        ['gene/DATA/', 'gene_orthologs.gz', False, 'gene_orthologs.gz'],
        ['gene/DATA/', 'gene_orthologs.gz', True, 'gene_orthologs'],
    )
    @unpack
    def test_download_file(self, ftp_path, file_name, run_gunzip, expect):
        res = self.c.download_file(ftp_path, file_name, DIR_DATA, run_gunzip)
        expect = os.path.join(DIR_DATA, expect)
        assert res == expect



    # @skip
    # @data(
    #     ['gene/tools/', None, 3],
    #     # ['pubmed/updatefiles/', 'gz', 0],
    # )
    # @unpack
    # def test_download_tree(self, path, pattern, expect):
    #     res = self.c.download_tree(DIR_DATA, path, pattern)
    #     assert len(res) >= expect