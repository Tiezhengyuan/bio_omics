'''
Unit tests for src/bioomics/PDB/analyze_complex.py (AnalyzeComplex)
All functions are pure; no files or network access needed.
'''
import pandas as pd

from ..helper import *

from src.bioomics.PDB.analyze_complex import AnalyzeComplex


@ddt
class TestCalSpan(TestCase):

    @data(
        ['1', '10', 9],
        # numbers are also allowed
        [1, 10, 9],
        # first digit run of each string is used
        ['res5x', 'res25x', 20],
        ['10', '10', 0],
    )
    @unpack
    def test_cal_span(self, a, b, expect):
        self.assertEqual(AnalyzeComplex.cal_span(a, b), expect)


def sig_df(res_nos):
    return pd.DataFrame({'res_no': res_nos})


class TestGroupSignificantPoints(TestCase):

    def test_groups(self):
        res = AnalyzeComplex.group_significant_points(
            sig_df([1, 2, 5, 6, 100]), 1)
        self.assertEqual(res, [
            [(0, 1), (1, 2)],
            [(2, 5), (3, 6)],
            [(4, 100)],
        ])

    def test_larger_max_span_merges_groups(self):
        res = AnalyzeComplex.group_significant_points(sig_df([1, 2, 6]), 4)
        self.assertEqual(res, [[(0, 1), (1, 2), (2, 6)]])

    def test_empty_df(self):
        res = AnalyzeComplex.group_significant_points(sig_df([]), 1)
        self.assertEqual(res, [])


class TestGetSpans(TestCase):

    def test_spans_filtered(self):
        # spans between consecutive residues, gaps <= 1 removed
        self.assertEqual(AnalyzeComplex.get_spans(sig_df([1, 2, 3, 10])), [7])
        self.assertEqual(
            AnalyzeComplex.get_spans(sig_df([1, 2, 10, 20])), [8, 10])

    def test_consecutive_only_returns_empty(self):
        self.assertEqual(AnalyzeComplex.get_spans(sig_df([1, 2, 3])), [])


class TestGetFragments(TestCase):

    def test_fragments(self):
        # [1,2,3] is one fragment of length 3, [10,11] of length 2
        self.assertEqual(
            AnalyzeComplex.get_fragments(sig_df([1, 2, 3, 10, 11])), [3, 2])

    def test_single_fragment_of_two(self):
        self.assertEqual(AnalyzeComplex.get_fragments(sig_df([4, 5])), [2])

    def test_single_residue_returns_empty(self):
        self.assertEqual(AnalyzeComplex.get_fragments(sig_df([1])), [])

    def test_consecutive_only_returns_empty(self):
        self.assertEqual(AnalyzeComplex.get_fragments(sig_df([1])), [])


class TestPrintMeta(TestCase):

    def test_print_meta(self):
        data = [
            {
                'type': 'antibody-antigen',
                'chains': [('L', 'light chain'), ('A', 'antigen chain')],
            },
        ]
        # smoke test: prints pairwise chain type and id
        AnalyzeComplex.print_meta(data)


class TestBuildDataset(TestCase):

    def setUp(self):
        self.ab1 = [{'seq': 'AAA'}, {'seq': 'C'}]
        self.rec1 = [{'seq': 'RR'}, {'seq': 'DE'}]
        self.data = [('ab1', self.ab1, 'rec1', self.rec1)]

    def test_build_dataset(self):
        res = AnalyzeComplex.build_dataset(self.data, 1)
        self.assertEqual(res, [{
            'input': 'RR | DE',
            'input_label': 'rec1',
            'output': 'AAA | C',
            'output_label': 'ab1',
        }])

    def test_min_fragment_filter(self):
        res = AnalyzeComplex.build_dataset(self.data, 2)
        self.assertEqual(res, [{
            'input': 'RR | DE',
            'input_label': 'rec1',
            'output': 'AAA',
            'output_label': 'ab1',
        }])

    def test_empty_antibody_skipped(self):
        res = AnalyzeComplex.build_dataset(
            [('ab2', [], 'rec1', self.rec1)], 1)
        self.assertEqual(res, [])

    def test_empty_receptor_skipped(self):
        res = AnalyzeComplex.build_dataset(
            [('ab1', self.ab1, 'rec2', [])], 1)
        self.assertEqual(res, [])
