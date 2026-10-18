'''
RNAcentral
download the sequence files of the integrated databases from the FTP server

note: the test downloads pirbase.fasta.gz (~9 MB) from the FTP server;
lncbook.fasta.gz (~180 MB) is too large for a unit test
'''
import os
from unittest import mock

from Bio import SeqIO

from ..helper import *
from src.bioomics import RNACentral


@ddt
class TestRNACentral(TestCase):

    @data(
        ['pirbase.fasta', os.path.join(DATA_DIR, 'RNACentral', 'pirbase.fasta')],
        ['wrong', None],
    )
    @unpack
    def test_download_sequence(self, file_name, expect):
        # pirbase.fasta is downloaded from <file_name>.gz and gunzipped,
        # an unknown database has no file on the server
        local_path, local_file = RNACentral(DATA_DIR, True).download_sequence(file_name)
        assert local_path == os.path.join(DATA_DIR, 'RNACentral')
        assert local_file == expect

    def test_download_sequence_content(self):
        # the gunzipped file is valid fasta
        local_file = os.path.join(DATA_DIR, 'RNACentral', 'pirbase.fasta')
        with open(local_file) as f:
            records = list(SeqIO.parse(f, 'fasta'))
        assert len(records) > 0

    def test_download_sequence_cached(self):
        # an existing file is returned without any FTP access
        local_file = os.path.join(DATA_DIR, 'RNACentral', 'pirbase.fasta')
        with mock.patch('src.bioomics.connector.conn_ftp.FTP') as m:
            m.side_effect = AssertionError('FTP must not be touched for cached files')
            local_path, got = RNACentral(DATA_DIR, False).download_sequence('pirbase.fasta')
        assert local_path == os.path.join(DATA_DIR, 'RNACentral')
        assert got == local_file

    def test_download_sequence_gz_cached(self):
        # a locally cached .gz file is returned as-is
        gz_file = os.path.join(DATA_DIR, 'RNACentral', 'cached.fasta.gz')
        with open(gz_file, 'wb') as f:
            f.write(b'test')
        local_path, got = RNACentral(DATA_DIR, False).download_sequence('cached.fasta.gz')
        assert got == gz_file
        os.remove(gz_file)
