'''
Unit tests for src/bioomics/antibody/process_seq.py (ProcessSeq)

The antibody ProcessSeq mirrors the protein one, except that insert_chr
applies insertions in idx_chr insertion order (no sorting).
'''
import os

from ..helper import *

from src.bioomics.antibody.process_seq import ProcessSeq


@ddt
class TestInsertChr(TestCase):

    @data(
        ['ABCD', {1: '-'}, 'A-BCD'],
        # insertions applied in idx_chr insertion order, with offset shifting
        ['ABCD', {1: '-', 2: '+'}, 'A-B+CD'],
        # unsorted keys are NOT sorted here; the '-' insertion shifts '+'
        ['ABCD', {2: '+', 1: '-'}, 'AB-+CD'],
        # insertion at both ends
        ['AB', {0: '[', 2: ']'}, '[AB]'],
        # empty sequence
        ['', {0: 'X'}, 'X'],
    )
    @unpack
    def test_insert_chr(self, seq, idx_chr, expect):
        self.assertEqual(ProcessSeq.insert_chr(seq, idx_chr), expect)


@ddt
class TestEncompass(TestCase):

    @data(
        # markers shift the insert positions of later markers
        ['AABBCC', [(1, 4)], 'A[ABB]CC'],
        # spans (0,2) and (4,6) enclose 'AA' and 'CC'
        ['AABBCCDD', [(0, 2), (4, 6)], '[AA]BB[CC]DD'],
        # no spans: unchanged
        ['AABB', [], 'AABB'],
    )
    @unpack
    def test_encompass(self, seq, idx_pool, expect):
        self.assertEqual(ProcessSeq.encompass(seq, idx_pool, '[', ']'), expect)

    @data(
        [
            'QVQLVQSGAEVKKPGASVKVSCKASGYTFT',
            ['GYTFT'],
            'QVQLVQSGAEVKKPGASVKVSCKAS[GYTFT]',
        ],
        # regexp patterns are allowed
        [
            'QVQLVQSGAEVKKPGASVKVSCKASGYTFT',
            ['G[A-Z]{3}T'],
            'QVQLVQSGAEVKKPGASVKVSCKAS[GYTFT]',
        ],
        # no match: unchanged
        ['AAAA', ['ZZZ'], 'AAAA'],
    )
    @unpack
    def test_encompass_substr(self, seq, substr_pool, expect):
        self.assertEqual(
            ProcessSeq.encompass_substr(seq, substr_pool, '[', ']'), expect)


@ddt
class TestMaskSeq(TestCase):

    @data(
        # default labels are <CDR1> ... <CDR5>
        ['AA[BB]CC[DD]EE', None, None, 'AA<CDR1>CC<CDR2>EE'],
        # custom labels are stringified
        ['AA[BB]CC[DD]EE', None, [1, 2], 'AA1CC2EE'],
        # custom regexp
        ['AA(BB)CC(DD)EE', r'\(.*?\)', None, 'AA<CDR1>CC<CDR2>EE'],
    )
    @unpack
    def test_mask_seq(self, seq, regexp, labels, expect):
        self.assertEqual(ProcessSeq.mask_seq(seq, regexp, labels), expect)

    def test_no_mask_returns_none(self):
        self.assertIsNone(ProcessSeq.mask_seq('AABBCC'))


@ddt
class TestDigitCdrMask(TestCase):

    @data(
        # trailing <CDRn> marker is stripped
        ['QQ<CDR1>', 'QQ'],
        # internal markers become digits
        ['QQ<CDR1>EE<CDR2>', 'QQ1EE'],
        ['<CDR3>AAAA', '3AAAA'],
        # empty/falsy passthrough
        ['', ''],
        ['ACDEFG', 'ACDEFG'],
    )
    @unpack
    def test_digit_cdr_mask(self, seq, expect):
        self.assertEqual(ProcessSeq.digit_cdr_mask(seq), expect)


@ddt
class TestTrimCdr3(TestCase):

    @data(
        # trims at the 3rd [CDR] block
        ['AA[BB]CC[DD]EE[FF]GG', 'AA[BB]CC[DD]EE'],
        # fewer than 3 blocks: unchanged
        ['AA[BB]CC', 'AA[BB]CC'],
        ['AABBCC', 'AABBCC'],
        # empty/falsy passthrough
        ['', ''],
    )
    @unpack
    def test_trim_cdr3(self, seq, expect):
        self.assertEqual(ProcessSeq.trim_cdr3(seq), expect)


class TestToFasta(TestCase):

    def setUp(self):
        import tempfile
        import shutil
        self.tmp_dir = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.tmp_dir, ignore_errors=True)
        self.records = [
            {'id': 'seq1', 'seq': 'AAA', 'description': 'first'},
            {'id': 'seq2', 'seq': 'CCC'},
        ]

    def test_write_fasta_file(self):
        from Bio import SeqIO
        outfile = os.path.join(self.tmp_dir, 'seqs.fa')
        self.assertTrue(ProcessSeq.to_fasta(self.records, outfile))
        with open(outfile) as f:
            text = f.read()
        self.assertIn('>seq1', text)
        self.assertIn('>seq2', text)
        with open(outfile) as f:
            records = list(SeqIO.parse(f, 'fasta'))
        self.assertEqual([r.id for r in records], ['seq1', 'seq2'])
        self.assertEqual(str(records[0].seq), 'AAA')
        self.assertEqual(str(records[1].seq), 'CCC')

    def test_no_outfile_returns_false(self):
        self.assertFalse(ProcessSeq.to_fasta(self.records))
