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
        res = RNACentral(DIR_DATA, False).download_sequence(file_name)
        assert res == expect