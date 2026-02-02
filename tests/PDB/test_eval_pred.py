'''
- pytest -s tests/PDB/test_eval_pred.py
'''
from ..helper import *
from src.bioomics import EvalPred

@ddt
class TestEvalPred(TestCase):

    @skip
    @data(
        # identical chain seq but different experimental data
        ['103M_0_A_0.pdb', '8F9I_0_A_0.pdb', True],
    )
    @unpack
    def test_tmalign(self, pdb1, pdb2):
        bin_dir = '/home/yuan/bio/TMalign-20180426/'
        pdb1 = os.path.join(DATA_DIR, 'PDB',  pdb1)
        pdb2 = os.path.join(DATA_DIR, 'PDB',  pdb2)
        res = EvalPred().tmalign(bin_dir, pdb1, pdb2)
        assert res['TM1'] == res['TM2']

    @data(
        # contacts exist
        ['1RZI.pdb', -98.5, 5.8e-73, True,],
        # no contact due to single chain
        ['6BSX_0_A_0.pdb', None, None, False],
    )
    @unpack
    def test_prodigy(self, pdb, expect_affinity, expect_dc, expect_contacts):
        pdb_file = os.path.join(DATA_DIR, 'PDB',  pdb)
        res = EvalPred().prodigy(pdb_file)
        assert res.get('binding_affinity') == expect_affinity
        assert res.get('dissociation_constant') == expect_dc
        assert res['contacts'] == expect_contacts

