'''
Test class 
'''
from tests.helper import *

from connector.connect_ftp import ConnectFTP

@ddt
class Test_(TestCase):

    def setUp(self):
        self.endpoint = 'ftp.ncbi.nlm.nih.gov'
        self.c = ConnectFTP(self.endpoint)

    @skip
    @data(
        ['geo/series/GSE3nnn/GSE3341/matrix', True],
        ['geo/series/GSE3nnn/GSE3341/matrix/GSE3341_series_matrix.txt.gz', False],
        ['', True],
    )
    @unpack
    @mock.patch.dict(os.environ, env)
    def test_is_dir(self, path, expect):
        res = self.c.is_dir(path)
        assert res == expect


    @skip
    @data(
        ['geo/series/GSE3nnn/GSE3341/matrix', '.gz', 1],
        ['geo/series/GSE3nnn/GSE3341/', '.gz', 0],
        ['geo/series/GSE3nnn/GSE3341/', None, 0],
    )
    @unpack
    @mock.patch.dict(os.environ, env)
    def test_download_files(self, path, pattern, expect):
        res = self.c.download_files(path, pattern)
        assert len(res) == expect


    @skip
    @data(
        # ['GSE3341', 'geo/series/GSE3nnn/GSE3341/', None, 3],
        ['baseline', 'pubmed/updatefiles/', 'gz', 3],
    )
    @unpack
    @mock.patch.dict(os.environ, env)
    def test_download_tree(self, local_name, path, pattern, expect):
        res = self.c.download_tree(local_name, path, pattern)
        assert len(res) >= expect

    
    @skip
    @data(
        # ['pub/taxonomy/new_taxdump/', 'new_taxdump.zip', None, True],
        # ['pub/taxonomy/new_taxdump/', 'wrong_name.zip', None, False],
        # ['pubmed/baseline/', 'pubmed23n0001.xml.gz', None, True],
        ['pubmed/updatefiles/', 'pubmed23n1237.xml.gz', None, True],
    )
    @unpack
    @mock.patch.dict(os.environ, env)
    def test_download_file(self, ftp_path, file_name, local_path, expect):
        res = self.c.download_file(ftp_path, file_name, local_path)
        assert res == expect

