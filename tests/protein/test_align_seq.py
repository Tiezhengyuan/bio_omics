'''
example: pytest -s tests/protein/test_align_seq.py
'''
from contextlib import redirect_stdout
import io

from ..helper import *

from src.bioomics import AlignSeq

@ddt
class TestAlignSeq(TestCase):

    @classmethod
    def setUpClass(cls):
        cls.seq_dict = {
            'AAAA': ['s1'],
            'AAAC': ['s2'],
            'GGGG': ['s3'],
        }
        cls.aln = AlignSeq.pairwise_align(cls.seq_dict)

    def test_pairwise_align(self):
        assert len(self.aln) == 3
        def get_rec(id_pair):
            return [r for r in self.aln if {r['id1'], r['id2']} == id_pair][0]
        # one mismatch: 3 matches of 4 residues
        rec = get_rec({'s1', 's2'})
        assert rec['seq1'] == 'AAAA' and rec['seq2'] == 'AAAC'
        assert rec['score'] == 3
        assert rec['ratio'] == 0.75
        # no common residue
        rec = get_rec({'s1', 's3'})
        assert rec['score'] == 0
        assert rec['ratio'] == 0.0

    def test_view_alignment(self):
        buf = io.StringIO()
        with redirect_stdout(buf):
            AlignSeq.view_alignment(self.aln)
        out = buf.getvalue()
        assert 'seq1=4, seq2=4' in out
        assert 'Score=3' in out

    def test_muscle_fail(self):
        a = AlignSeq(verbose=True)
        outfile = os.path.join(RESULTS_DIR, 'muscle.faa')
        # input file does not exist: muscle fails, no exception raised
        assert a.muscle('no_exist.faa', outfile) == 'fail'
        assert not os.path.isfile(outfile)