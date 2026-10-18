'''
Test src/bioomics/common/parse_text.py
'''
from ..helper import *

from src.bioomics import ParseText


class TestParseText(TestCase):

    # ----- split_block -----
    def test_split_block(self):
        block = [
            'Alignments',
            'seq1 line A',
            '',
            'seq2 line A',
            'seq2 line B',
            '',
            '',
            'after the double blank',
        ]
        res = ParseText.split_block(block, 'Alignments')
        self.assertEqual(res, [['seq1 line A'], ['seq2 line A', 'seq2 line B']])

    def test_split_block_trailing_single_blank(self):
        # regression: block[i+1] raised IndexError at the end of the block
        block = ['Alignments', 'seq1 line A', '']
        res = ParseText.split_block(block, 'Alignments')
        self.assertEqual(res, [['seq1 line A']])

    def test_split_block_header_missing(self):
        self.assertEqual(ParseText.split_block(['a', '', ''], 'Nothing'), [])

    def test_split_block_igblast_shape(self):
        # usage seen in src/bioomics/protein/parse_igblast.py
        block = [
            'Alignments',
            'Query_1  50 nucleotides',
            '',
            'Query_2  60 nucleotides',
            '',
            '',
        ]
        res = ParseText.split_block(block, 'Alignments')
        self.assertEqual(res, [
            ['Query_1  50 nucleotides'],
            ['Query_2  60 nucleotides'],
        ])

    # ----- parse_lines_items -----
    def test_parse_lines_items(self):
        lines = [
            'AB   CD',
            'EF   GH',
        ]
        res = ParseText.parse_lines_items(lines)
        self.assertEqual(res, [['AB', 'CD'], ['EF', 'GH']])

    def test_parse_lines_items_three_columns(self):
        lines = [
            'A    BB    CCC',
            'D    EE    FFF',
        ]
        res = ParseText.parse_lines_items(lines)
        self.assertEqual(res, [['A', 'BB', 'CCC'], ['D', 'EE', 'FFF']])

    def test_parse_lines_items_empty(self):
        self.assertEqual(ParseText.parse_lines_items([]), [])
