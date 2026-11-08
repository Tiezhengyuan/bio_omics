'''
Test class 
'''
from tests.helper import *

from src.bioomics import ConnHTTP

@ddt
class TestConnHTTP(TestCase):
    
    @data(
        # reliable 200 endpoint
        ["https://www.google.com", False, 'robots.txt', 200, True],
        # unreachable host (RFC 2606 reserved .invalid TLD): connection error
        ["https://wrong_for_test.invalid", False, 'xx.txt', 501, False],
        # wrong file name: 404 must not create the local file
        ["https://www.google.com", False, 'xx.txt', 404, False],
    )
    @unpack
    def test_text(self, url, overwrite, file_name, \
            expect_status, expect_saved):
        conn = ConnHTTP(url, DIR_DATA, overwrite)
        # end_point is appended to the base URL (leading slash required),
        # file_name is joined with the local dir (no leading slash)
        res, local_file = conn.download_file('/' + file_name, file_name)
        # overwrite=False + file already cached locally: no request is made
        # and to_file() returns None - a cached hit counts as a 200
        status = res.status if res is not None else 200
        assert status == expect_status
        if expect_saved:
            assert local_file == os.path.join(DIR_DATA, file_name)
            assert os.path.isfile(local_file)
        else:
            # failed download: no local file is created
            assert local_file is None or not os.path.isfile(local_file)

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



