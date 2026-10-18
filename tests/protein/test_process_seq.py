'''
example: pytest -s tests/protein/test_process_seq.py
'''
from ..helper import *

from src.bioomics import ProcessSeq

@ddt
class TestProcessSeq(TestCase):

    def test_insert_chr(self):
        assert ProcessSeq.insert_chr('ACGTACGT', {0: '['}) == '[ACGTACGT'
        assert ProcessSeq.insert_chr('ACGTACGT', {4: ']'}) == 'ACGT]ACGT'
        assert ProcessSeq.insert_chr('ACGTACGT', {4: '[', 8: ']'}) == 'ACGT[ACGT]'
        # unsorted indices are applied in ascending order
        assert ProcessSeq.insert_chr('ACGTACGT', {8: ']', 4: '['}) == 'ACGT[ACGT]'

    def test_encompass(self):
        assert ProcessSeq.encompass('ACGTACGT', [(4, 8)], '[', ']') == 'ACGT[ACGT]'
        assert ProcessSeq.encompass('ACGTACGT', [(1, 3), (5, 7)], '<', '>') \
            == 'A<CG>TA<CG>T'

    def test_encompass_substr(self):
        # only the first occurrence of each substring is enclosed
        assert ProcessSeq.encompass_substr('ACGTACGT', ['CGT'], '[', ']') \
            == 'A[CGT]ACGT'
        assert ProcessSeq.encompass_substr('ACGTACGT', ['ACG'], '[', ']') \
            == '[ACG]TACGT'
        # substring not found: sequence unchanged
        assert ProcessSeq.encompass_substr('ACGT', ['TTT'], '[', ']') == 'ACGT'

    def test_mask_seq(self):
        assert ProcessSeq.mask_seq('A[XY]B[Z]C') == 'A<CDR1>B<CDR2>C'
        # custom labels
        assert ProcessSeq.mask_seq('A[XY]B[Z]C', None, [7, 9]) == 'A7B9C'
        # custom regexp
        assert ProcessSeq.mask_seq('A..B..C', r'\.\.', ['x', 'y']) == 'AxByC'
        # nothing to mask
        assert ProcessSeq.mask_seq('ABC') is None

    def test_digit_cdr_mask(self):
        assert ProcessSeq.digit_cdr_mask('A<CDR1>B<CDR2>C') == 'A1B2C'
        # the trailing mask digit is removed
        assert ProcessSeq.digit_cdr_mask('A<CDR1>B<CDR3>') == 'A1B'
        assert ProcessSeq.digit_cdr_mask('') == ''
        assert ProcessSeq.digit_cdr_mask(None) is None

    def test_trim_cdr3(self):
        seq = 'QVQL[GYTFTGYY]IHW[INPNRGGT]NYC[AR]GQG'
        assert ProcessSeq.trim_cdr3(seq) == 'QVQL[GYTFTGYY]IHW[INPNRGGT]NYC'
        # less than three CDR: unchanged
        assert ProcessSeq.trim_cdr3('A[XY]B') == 'A[XY]B'
        assert ProcessSeq.trim_cdr3('') == ''
        assert ProcessSeq.trim_cdr3(None) is None

    def test_to_fasta(self):
        outfile = os.path.join(RESULTS_DIR, 'process_seq.faa')
        data = [{'id': 'a', 'seq': 'ACD'}, {'id': 'b', 'seq': 'EEE'}]
        assert ProcessSeq.to_fasta(data, outfile) is True
        with open(outfile) as f:
            assert f.read() == '>a\nACD\n>b\nEEE\n'
        # no outfile: only records are built
        assert ProcessSeq.to_fasta(data) is False

    def test_to_colabfold_csv(self):
        outfile = os.path.join(RESULTS_DIR, 'colabfold.csv')
        data = {'complex1': ['AAA', 'DDD']}
        assert ProcessSeq.to_colabfold_csv(data, outfile) == 'succeed'
        with open(outfile) as f:
            assert f.read() == 'id,sequence\ncomplex1,AAA:DDD\n'

    def test_to_colabfold_fasta(self):
        outfile = os.path.join(RESULTS_DIR, 'colabfold.faa')
        data = {'complex1': ['AAA', 'DDD']}
        assert ProcessSeq.to_colabfold_fasta(data, outfile) == 'succeed'
        with open(outfile) as f:
            assert f.read() == '>complex1\nAAA:\nDDD'