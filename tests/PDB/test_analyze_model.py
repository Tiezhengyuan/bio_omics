from ..helper import *
from src.bioomics import AnalyzeModel

@ddt
class TestAnalyzeModel(TestCase):

    @data(
        # identical chain seq but different experimental data
        ['103M_0_A_0.pdb', '8F9I_0_A_0.pdb', True],
    )
    @unpack
    def test_tmalign(self, pdb1, pdb2):
        bin_dir = '/home/yuan/bio/TMalign-20180426/'
        pdb1 = os.path.join(DATA_DIR, 'PDB',  pdb1)
        pdb2 = os.path.join(DATA_DIR, 'PDB',  pdb2)
        res = AnalyzeModel.tmalign(bin_dir, pdb1, pdb2)
        print(res)
        assert res['TM1'] == res['TM2']

