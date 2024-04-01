'''
Test class 
'''
from .helper import *
from src.bioomics import ConnFTPlib as cf

@ddt
class TestConnFTPlib(TestCase):

    @skip
    @data(
        ['ftp.ncbi.nlm.nih.gov', '/pubmed', None],
        ['ftp.ncbi.nlm.nih.gov', 'pubmed', None],
        # wrong ftp path
        ['ftp.ncbi.nlm.nih.gov', 'wrong_path', 550],
        # wrong ftp url
        ['wrong_ftp', 'pubmed', 11001],
    )
    @unpack
    def test_connect_ftp(self, endpoint, ftp_path, expect):
        res = cf.connect_ftp(endpoint, ftp_path)
        assert res.get('errno') == expect

    @skip
    @data(
        ['ftp.ncbi.nlm.nih.gov', '/pubmed'],
    )
    @unpack
    def test_list_contents(self, endpoint, ftp_path):
        cf.list_contents(endpoint, ftp_path)

    @skip
    @data(
        ['ftp.ncbi.nlm.nih.gov', '/pubmed', 'gz'],
    )
    @unpack
    def test_scan_tree(self, endpoint, ftp_path, pattern):
        res = cf.scan_tree(endpoint, ftp_path, pattern)

    @data(
        ['ftp.ncbi.nlm.nih.gov', '/pubmed', 'gz'],
    )
    @unpack
    def test_download_tree(self, endpoint, ftp_path, pattern):
        '''
        download all files
        '''
        local_path = os.path.join(DIR_DATA, 'NCBI')
        cf.download_tree(endpoint, ftp_path, pattern, local_path)


    @skip
    @data(
        ['ftp.ncbi.nlm.nih.gov', '/pubmed',
            {'/pubmed/baseline', '/pubmed/pubmedcommons', '/pubmed/updatefiles'},
            {'deleted.pmids.gz', 'id_list.txt.gz', 'J_Entrez.gz', 'J_Medline.gz'},
        ],
    )
    @unpack
    def test_parse_contents(self, endpoint, ftp_path, expect_dirs, expect_files):
        res = cf.parse_contents(endpoint, ftp_path, 'gz')
        assert set(res['dirs']) == expect_dirs
        assert res['current_path'] == ftp_path
        assert set(res['files']) == expect_files


    @skip
    @data(
        ['/pubmed/updatefiles', 'pubmed23n1237.xml.gz', True],
        ['/pubmed/updatefiles', 'wrong_file_name', False],
    )
    @unpack
    def test_download_file(self, ftp_path, file_name, expect):
        endpoint = 'ftp.ncbi.nlm.nih.gov'
        res = cf.download_file(endpoint, ftp_path, file_name, DIR_DATA)
        assert res == expect

    @skip
    @data(
        ['/pubmed/baseline', 'gz', False],
        ['/pubmed', 'gz', False],
        ['/wrong_path', 'gz', True],
        ['/pubmed', 'unmatched', True],
    )
    @unpack
    def test_retrieve_file_names(self, path, pattern, expect):
        endpoint = 'ftp.ncbi.nlm.nih.gov'
        res = cf.retrieve_file_names(endpoint, path, pattern)
        assert (len(res)==0) == expect
