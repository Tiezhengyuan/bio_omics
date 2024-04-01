'''
'''

from .helper import *
from src.bioomics import RNACentral

@ddt
class TestRNACentral(TestCase):

    @data(
        ['lncbook.fasta', os.path.join(DIR_DATA, 'RNACentral', 'lncbook.fasta')],
        ['pirbase.fasta', os.path.join(DIR_DATA, 'RNACentral', 'pirbase.fasta')],
        ['wrong', None],
    )
    @unpack
    def test_download_sequence(self, file_name, expect):
        local_path, local_file = RNACentral(DIR_DATA, True).download_sequence(file_name)
        assert local_path == os.path.join(DIR_DATA, 'RNACentral')
        assert local_file == expect