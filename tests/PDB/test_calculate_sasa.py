
from ..helper import *
from src.bioomics import SasaCalculate


@ddt
class TestSasaCalculate(TestCase):

    
    @data(
        ['7S8L', None, (1096, 12)],
        ['7S8L', True, (1096, 12)],
    )
    @unpack
    def test_freesasa(self, pdb_id, hetatm, expect):
        pdb_file = os.path.join(DATA_DIR, 'PDB', f'{pdb_id}.pdb')
        c = SasaCalculate(pdb_file)
        result = c.cal_freesasa(None, hetatm)
        assert result.shape == expect

