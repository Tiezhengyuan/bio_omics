'''
Test class DownloadExpasy: FTP of ExPASy
'''
from ..helper import *
from src.bioomics import DownloadExpasy
from src.bioomics.connector.conn_ftp import ConnFTP


@ddt
class TestExpasy(TestCase):

    def setUp(self):
        self.c = DownloadExpasy(DIR_DATA, False)

    def test_init(self):
        # local dir is <local_dir>/ExPASy
        assert self.c.local_dir == os.path.join(DIR_DATA, 'ExPASy')
        assert self.c.url == 'ftp.expasy.org'
        assert self.c.overwrite is False

    def test_download_data(self):
        # mock the FTP download tree: the real ExPASy databases directory
        # is far too large to download in a unit test
        with mock.patch.object(ConnFTP, 'download_tree',
                return_value=['databases/somestruct.dat']) as m:
            res = self.c.download_data()
        assert res == ['databases/somestruct.dat']
        m.assert_called_once_with(
            local_path=self.c.local_dir,
            endpoint='databases',
        )

