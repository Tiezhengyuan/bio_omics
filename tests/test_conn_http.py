'''
Test class 
'''
from tests.helper import *

from src.bioomics import ConnHTTP

@ddt
class TestConnHTTP(TestCase):
    
    @data(
        ["http://httpbin.org", False, '/robots.txt', 200, \
            os.path.join(DIR_DATA, 'robots.txt')],
        # wrong url
        ["http://wrong.org", False, '/xx.txt', 501, None],
        # wrong file name
        ["http://httpbin.org", False, '/xx.txt', 404, None],
    )
    @unpack
    def test_text(self, url, overwrite, file_name, \
            expect_status, expect_local_file):
        conn = ConnHTTP(url, DIR_DATA, overwrite)
        res, local_file = conn.download_file(file_name)
        assert res.status == expect_status
        assert local_file == expect_local_file

    @data(
        [
            "https://www.iedb.org/downloader.php",
            '?file_name=doc/antigen_full_v3.zip', 
            'antigen_full_v3.zip',
            os.path.join(DIR_DATA, 'antigen_full_v3.zip')
        ],
    )
    @unpack
    def test_zip(self, url, end_point, file_name, expect_local_file):
        conn = ConnHTTP(url, DIR_DATA, True)
        _, local_file = conn.download_file(end_point, file_name)
        assert local_file == expect_local_file



