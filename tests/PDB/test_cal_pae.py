from ..helper import *
from src.bioomics import CalculatePae

print(DATA_DIR)

@ddt
class TestCalculatePae(TestCase):

    def test_(self):
        pdb_file1 = os.path.join(DATA_DIR, 'PDB', '2OX8_A.pdb')
        pdb_file2 = os.path.join(DATA_DIR, 'PDB', 'AF-Q8IZS7-F1-model_v6.pdb')
        p = CalculatePae(pdb_file1, pdb_file2)
        pae_matrix = p.calculate_pae()
        stats = p.calculate_pae_statistics(pae_matrix)
        print(stats)


