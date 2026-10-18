'''
example: pytest -s tests/protein/test_analyze_seq.py
'''
from contextlib import redirect_stdout
import io

from ..helper import *

from src.bioomics import AnalyzeSeq

@ddt
class TestAnalyzeSeq(TestCase):

    def setUp(self):
        # alignment free of ambiguity for the cases below
        self.a = AnalyzeSeq('AAAA')

    def test_align_pairwise(self):
        assert self.a.align_pairwise('AAAA').score == 4
        # one mismatch and one indel are scored the same (both 0)
        assert self.a.align_pairwise('AAHA').score == 3
        assert self.a.align_pairwise('AA').score == 2

    def test_iter_align(self):
        items = [
            {'name': 'exact', 'query_seq': 'AAAA', 'query_val': [1, 2, 3, 4]},
            {'name': 'one_mismatch', 'query_seq': 'AAHA', 'query_val': [1, 2, 3, 4]},
            {'name': 'short', 'query_seq': 'AA', 'query_val': [7, 9]},
            {'name': 'no_val', 'query_seq': 'AAAA'},
        ]
        res = self.a.iter_align(items)
        assert len(res) == 4
        rec = res[0]
        assert rec['name'] == 'exact'
        assert rec['score'] == 4
        assert rec['ratio'] == 1.0
        # exact match: values are unchanged
        assert list(rec['adj_val']) == [1, 2, 3, 4]
        rec = res[1]
        assert rec['score'] == 3
        assert rec['ratio'] == 0.75
        # H is skipped: value of the last A (4) fills the gap position
        assert list(rec['adj_val']) == [1, 2, 4, 1]
        rec = res[2]
        assert rec['score'] == 2
        assert rec['ratio'] == round(2*2/6, 4)
        # target positions outside the alignment hold the minimum value
        assert list(rec['adj_val']) == [7, 9, 7, 7]
        # no query_val: no adj_val
        assert 'adj_val' not in res[3]

    def test_parse_values(self):
        b = AnalyzeSeq('AAAAAA')
        # aligned columns take query values, the rest the minimum
        assert list(b.parse_values('AAAA', [1, 2, 3, 4])) == [1, 2, 3, 4, 1, 1]
        assert list(b.parse_values('AAAAAA', [5, 6, 7, 8, 9, 10])) \
            == [5, 6, 7, 8, 9, 10]

    def test_print_alignments(self):
        self.a.iter_align([{'name': 'exact', 'query_seq': 'AAAA'}])
        buf = io.StringIO()
        with redirect_stdout(buf):
            self.a.print_alignments()
        out = buf.getvalue()
        assert 'AAAA' in out
        assert 'target' in out
        assert 'query' in out