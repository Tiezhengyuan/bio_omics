from ..helper import *
from src.bioomics import CalculatePlddt

print(DATA_DIR)

@ddt
class TestCalculatePlddt(TestCase):

    @data(
        ['2OX8_A.pdb', (606, 36.46), (734, 24.96)]
    )
    @unpack
    def test_monomer_ca_plddt(self, input, expect_first, expect_last):
        pdb_file = os.path.join(DATA_DIR, 'PDB', input)
        result = CalculatePlddt(pdb_file).monomer_ca_plddt()
        assert tuple(result.iloc[0]) == expect_first
        assert tuple(result.iloc[-1]) == expect_last

    @data(
        ['2OX8_A.pdb', 28.79],
        ['7S8L.pdb', 59.14],
        ['AF-Q181Q7-F1-model_v6.pdb', 95.05],
    )
    @unpack
    def test_plddt_statistics(self, input, expect_avg_plddt):
        pdb_file = os.path.join(DATA_DIR, 'PDB', input)
        cp = CalculatePlddt(pdb_file)
        plddt = cp.monomer_ca_plddt()
        result = cp.plddt_statistics(plddt)
        assert round(result['mean_plddt'], 2) == expect_avg_plddt
