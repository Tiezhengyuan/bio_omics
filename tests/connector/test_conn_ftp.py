'''
Test class 
'''
from tests.helper import *

from src.bioomics import ConnFTP

URL = 'ftp.ncbi.nlm.nih.gov'

@ddt
class TestConnFTP(TestCase):

    @data(
        ['blast', True],
        ['robts.txt', False]
    )
    @unpack
    def test_is_dir(self, name, expect):
        c = ConnFTP(URL, False, False)
        ftp = c.connect()
        res = c._is_dir(ftp, name)
        assert res == expect

    @data(
        [
            'geo/series/GSE3nnn/GSE3341/matrix',
            [('geo/series/GSE3nnn/GSE3341/matrix', 'GSE3341_series_matrix.txt.gz')],
        ],
    )
    @unpack
    def test_list_files(self, path, expect):
        c = ConnFTP(URL, False, False)
        res = c.list_files(path)
        assert res == expect

    @data(
        ['geo/series/GSE3nnn/GSE3341/matrix', '.gz', 1],
        ['geo/series/GSE3nnn/GSE3341/', '.gz', 0],
        ['geo/series/GSE3nnn/GSE3341/', None, 0],
    )
    @unpack
    def test_download_files(self, path, pattern, expect):
        c = ConnFTP(URL, True, False)
        res = c.download_files(path, pattern, DIR_DATA)
        assert len(res) == expect


    @data(
        ['gene/tools/', 'README', False, 'README'],
        ['gene/DATA/', 'gene_orthologs.gz', False, 'gene_orthologs.gz'],
        ['gene/DATA/', 'gene_orthologs.gz', True, 'gene_orthologs'],
    )
    @unpack
    def test_download_file(self, ftp_path, file_name, run_gunzip, expect):
        '''
        NB: run_gunzip moved from download_file args to the constructor
        in the current API.
        '''
        c = ConnFTP(URL, True, run_gunzip)
        res = c.download_file(ftp_path, file_name, DIR_DATA)
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