'''
MiRBase
'''
from .helper import *
from src.bioomics import Mirbase

class TestMirbase(TestCase):
    def setUp(self):
        self.c = Mirbase(DIR_DATA, False)

    def test_download_hairpin(self):
        local_path, local_file = self.c.download_hairpin()
        assert local_path == os.path.join(DIR_DATA, 'miRBase')
        assert local_file == os.path.join(DIR_DATA, 'miRBase', 'hairpin.fa')

    def test_download_mature(self):
        local_path, local_file = self.c.download_mature()
        assert local_path == os.path.join(DIR_DATA, 'miRBase')
        assert local_file == os.path.join(DIR_DATA, 'miRBase', 'mature.fa')

