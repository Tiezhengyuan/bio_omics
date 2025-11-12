from ..helper import *
from src.bioomics import CalculateDist

@ddt
class TestCalculateDist(TestCase):

    def setUp(self):
        self.c = {}
        for pdb_id, file_name in [('7S8L', '7S8L.pdb'),]:
            pdb_file = os.path.join(DATA_DIR, 'PDB', file_name)
            self.c[pdb_id] = CalculateDist(pdb_file)
    
    @data(
        ['7S8L', 'A', 'B', 6, 221],
        ['7S8L', 'C', 'D', 337, 51],
        ['7S8L', 'A', 'Wrong', None, None],
    )
    @unpack
    def test_cal_dist(self, input, chain1_id, chain2_id, expect1, expect2):
        g1, g2 = self.c[input].cal_dist(chain1_id, chain2_id)
        if g1 is None or g2 is None:
            assert g1 == expect1
            assert g2 == expect2
        else:
            assert len(g1) == expect1
            assert len(g2) == expect2
        