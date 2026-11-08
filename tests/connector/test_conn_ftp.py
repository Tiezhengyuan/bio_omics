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


import tempfile
from ftplib import error_perm


@ddt
class TestConnFTPRetry(TestCase):
    '''
    Behavior merged from the former ConnFTPlib (conn_ftplib.py):
    connection/download retry with backoff, FTP error triage and local
    file handling. Verified with mocks, no network needed.
    '''

    def test_connect_retries_transient_errors(self):
        '''transient connection errors are retried try_times times'''
        ftp = mock.MagicMock()
        with mock.patch('src.bioomics.connector.conn_ftp.FTP',
                side_effect=[OSError('t1'), OSError('t2'), ftp]), \
            mock.patch('src.bioomics.connector.conn_ftp.time.sleep') as sleep:
            c = ConnFTP(URL, False, False, try_times=3)
            res = c.connect('pubmed')
            assert res is ftp
            assert sleep.call_count == 2

    def test_connect_fails_after_retries(self):
        '''connect raises IOError after exhausting the retries'''
        with mock.patch('src.bioomics.connector.conn_ftp.FTP',
                side_effect=OSError('connection down')), \
            mock.patch('src.bioomics.connector.conn_ftp.time.sleep'):
            c = ConnFTP(URL, False, False, try_times=3)
            with self.assertRaises(IOError):
                c.connect('pubmed')

    def test_connect_wrong_path_no_retry(self):
        '''a wrong FTP path (error_perm 550) raises without retrying'''
        ftp_cls = mock.MagicMock()
        ftp = mock.MagicMock()
        ftp.cwd.side_effect = error_perm('550 No such directory')
        ftp_cls.return_value = ftp
        with mock.patch('src.bioomics.connector.conn_ftp.FTP', ftp_cls), \
            mock.patch('src.bioomics.connector.conn_ftp.time.sleep') as sleep:
            c = ConnFTP(URL, False, False, try_times=3)
            with self.assertRaises(IOError):
                c.connect('wrong_path')
            assert ftp_cls.call_count == 1
            assert sleep.call_count == 0

    def test_download_file_retries_and_cleans_partial(self):
        '''a failed download is retried and leaves no partial file'''
        with tempfile.TemporaryDirectory() as local_path:
            with mock.patch('src.bioomics.connector.conn_ftp.FTP',
                    side_effect=OSError('net down')), \
                mock.patch('src.bioomics.connector.conn_ftp.time.sleep'):
                c = ConnFTP(URL, True, False, try_times=2)
                res = c.download_file(
                    'pubmed/updatefiles', 'README.gz', local_path)
                assert res is None
                assert not os.path.isfile(
                    os.path.join(local_path, 'README.gz'))

    def test_download_file_skips_existing(self):
        '''existing local files are skipped when overwrite is False'''
        with tempfile.TemporaryDirectory() as local_path:
            local_file = os.path.join(local_path, 'README')
            with open(local_file, 'w') as f:
                f.write('data')
            with mock.patch('src.bioomics.connector.conn_ftp.FTP') as ftp_cls:
                c = ConnFTP(URL, False, False)
                res = c.download_file('gene/tools', 'README.gz', local_path)
                # 'README.gz' maps to local 'README' which exists: no FTP call
                assert res == local_file
                ftp_cls.assert_not_called()

    def test_download_file_gunzip(self):
        '''run_gunzip decompresses the downloaded .gz file'''
        with tempfile.TemporaryDirectory() as local_path:
            with mock.patch('src.bioomics.connector.conn_ftp.FTP'), \
                mock.patch('src.bioomics.connector.conn_ftp.gunzip') as gz:
                c = ConnFTP(URL, True, True)
                res = c.download_file('gene/tools', 'README.gz', local_path)
                expect = os.path.join(local_path, 'README')
                assert res == expect
                gz.assert_called_once_with(
                    os.path.join(local_path, 'README.gz'), '-f')

    def test_download_tree_broken_path_is_skipped(self):
        '''an unreachable directory doesn't break the whole tree download'''
        with mock.patch.object(ConnFTP, 'connect',
                side_effect=IOError("Failure: can't connect")):
            c = ConnFTP(URL, False, False)
            with tempfile.TemporaryDirectory() as local_path:
                res = c.download_tree(local_path, endpoint='wrong_path')
                assert res == []