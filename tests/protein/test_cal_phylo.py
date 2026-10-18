'''
example: pytest -s tests/protein/test_cal_phylo.py
'''
from ..helper import *

from src.bioomics import CalPhylo

@ddt
class TestCalPhylo(TestCase):

    def setUp(self):
        self.indir = os.path.join(DATA_DIR, 'MSA')
        self.c = CalPhylo(self.indir)

    def test_msa_distance(self):
        msa = os.path.join(self.indir, '1A6U_0_L.afa')
        scores = self.c.msa_distance(msa)
        assert set(scores) == {'identity', 'blosum62'}
        for dist in scores.values():
            # distance of sequences 1..9 to the first one
            assert len(dist) == 9
            assert dist.name == '1A6U_0_L'
        identity = scores['identity']
        assert all(0 <= v <= 1 for v in identity)
        assert round(identity.iloc[0], 6) == 0.633028

    def test_calc_scan(self):
        res = self.c.calc_scan()
        assert set(res) == {'identity', 'blosum62'}
        for df in res.values():
            # one .afa file in the directory
            assert len(df) == 1
            assert list(df.index) == ['1A6U_0_L']
        # identity distances vs the first sequence
        identity = res['identity']
        assert round(identity.iloc[0]['IGLV1-40*02'], 6) == 0.633028
        assert round(identity.iloc[0]['IGLV1-36*01'], 6) == 0.669725

    def test_calc_scan_skips_other_extensions(self):
        res = self.c.calc_scan()
        # only .afa files are scanned
        assert all(len(df) == 1 for df in res.values())