'''
6PUH: Structure of human MAIT A-F7 TCR in complex with human MR1-Ribityl-less
6P8M: Crystal Structure of Antibody P-p3b3 A60C Heavy Chain 
    in Complex with 426c HIV-1 gp120 core G459C
'''
from ..helper import *

from src.bioomics import ParseIgblast

@ddt
class TestParseIgblast(TestCase):

    def setUp(self):
        self.c = {}
        _dir = os.path.join(DATA_DIR, 'igblastp')
        for prefix in ('6PUH', '6P8M', '6EDZ'):
            fa_file = os.path.join(_dir, f'{prefix}_chain_sequences.faa')
            aln_file = os.path.join(_dir, f'{prefix}_chain_sequences.vregion.igblastp3')
            self.c[prefix] = ParseIgblast(aln_file, fa_file)
    @data(
        ['6P8M', 3],
        ['6PUH', 8],
        ['6EDZ', 4],
    )
    @unpack
    def test_iter_query(self, input, expect):
        result = self.c[input].iter_query()
        result = len(list(result))
        assert result == expect

    @data(
        ['6PUH', 8, 'Query= 6PUH_0_A'],
        ['6P8M', 3, 'Query= 6P8M_0_C'],
        ['6EDZ', 4, 'Query= 6EDZ_0_A_0'],
    )
    @unpack
    def test_iter_fmt3(self, input, expect, expect_first):
        p = self.c[input]
        result = [i for i in p.iter_fmt3()]
        assert len(result) == expect
        assert result[0][0] == expect_first

    @data(
        ['6PUH', 0, '6PUH_0_A'],
        ['6PUH', -1, '6PUH_0_H'],
        ['6P8M', 0, '6P8M_0_C'],
        ['6P8M', -1, '6P8M_0_L'],
    )
    @unpack
    def test_query_id(self, input, idx, expect):
        p = self.c[input]
        blocks = [i for i in p.iter_fmt3()]
        result = p.query_id(blocks[idx])
        assert result == expect

    @data(
        ['6PUH', 0, 265],
        ['6PUH', -1, 98],
        ['6P8M', 0, 333],
        ['6P8M', -1, 204],
    )
    @unpack
    def test_query_length(self, input, idx, expect):
        p = self.c[input]
        blocks = [i for i in p.iter_fmt3()]
        result = p.query_length(blocks[idx])
        assert result == expect

    @data(
        ['6PUH', {'specie': None, 'fasta': 8, 'aln': 8, 'pairwise': 8, 'json': 8}],
        ['6P8M', {'specie': None, 'fasta': 3, 'aln': 3, 'pairwise': 3, 'json': 3}],
    )
    @unpack
    def test_analyze_aln(self, input, expect):
        result = self.c[input].analyze_aln()
        assert result == expect

    @data(
        # not hit found
        ['6P8M', 0, []],
        # signficant
        ['6P8M', 1, [
            {'allele_name': 'IGHV1-2*02|Homo_sapiens', 'product': 'unnamed protein product', 
                'bits-score': '202', 'e-value': '1e-68', 'gene_name': 'IGHV1-2', 
                'gene_family': 'IGHV1', 'isotype': 'IG', 'chain_type': 'H'}, 
            {'allele_name': 'IGHV1-2*07|Homo_sapiens', 'product': 'unnamed protein product', 
                'bits-score': '200', 'e-value': '7e-68', 'gene_name': 'IGHV1-2', 
                'gene_family': 'IGHV1', 'isotype': 'IG', 'chain_type': 'H'}, 
            {'allele_name': 'IGHV1-2*04|Homo_sapiens', 'product': 'unnamed protein product', 
                'bits-score': '199', 'e-value': '2e-67', 'gene_name': 'IGHV1-2', 
                'gene_family': 'IGHV1', 'isotype': 'IG', 'chain_type': 'H'}]
        ],
        ['6P8M', 1, [
            {'allele_name': 'IGHV1-2*02|Homo_sapiens', 'product': 'unnamed protein product', 
                'bits-score': '202', 'e-value': '1e-68', 'gene_name': 'IGHV1-2', 
                'gene_family': 'IGHV1', 'isotype': 'IG', 'chain_type': 'H'}, 
            {'allele_name': 'IGHV1-2*07|Homo_sapiens', 'product': 'unnamed protein product', 
                'bits-score': '200', 'e-value': '7e-68', 'gene_name': 'IGHV1-2', 
                'gene_family': 'IGHV1', 'isotype': 'IG', 'chain_type': 'H'}, 
            {'allele_name': 'IGHV1-2*04|Homo_sapiens', 'product': 'unnamed protein product', 
                'bits-score': '199', 'e-value': '2e-67', 'gene_name': 'IGHV1-2', 
                'gene_family': 'IGHV1', 'isotype': 'IG', 'chain_type': 'H'}]
        ],
    )
    @unpack
    def test_sig_alignment(self, input, idx, expect):
        p = self.c[input]
        blocks = [i for i in p.iter_fmt3()]
        result = p.sig_alignment(blocks[idx])
        assert result == expect

    @data(
        ['6P8M', 0, [], [], []],
        ['6P8M', 1, 
            ['<-------FR1-IMGT--------><CDR1-I><---FR2-IMGT----><CDR2-I><--------------FR3-IMGT---------', '----->  '], 
            ['QVQLVQSGAEVKKPGASVKVSCKASGYTFTGYYIHWVRQAPGQGLEWMGWINPNRGGTNYCQKFQGRVTMTRDTSISTAYMELSRLRSDD', 'TAVYYCAR'], 
            [1, 90, 91, 98]
        ],
        ['6P8M', 2, 
            ['<------FR1-IMGT-------><-CDR1-IMGT><---FR2-IMGT----><C><-------------FR3-IMGT-------------', '>   '],
            ['MSQSPSSLAVSVGEKVTMSCKSSQ-------QKNYLAWYQQKPGQSPKLLIYWASTRESGVPDRFTGSGSGTDFTLTISSVKAEDLAVYY', 'CQQY'],
            [1, 83, 84, 87],
        ],
    )
    @unpack
    def test_alignments(self, input, idx, expect_region, expect_aln_seq, expect_pos):
        p = self.c[input]
        blocks = [i for i in p.iter_fmt3()]
        region, aln_seq, pos = p.alignments(blocks[idx])
        assert region == expect_region
        assert aln_seq == expect_aln_seq
        assert pos == expect_pos

    @data(
        ['6P8M', 0, {}],
        ['6P8M', 1, {'start': 0, 'end': 98, 
            'gap_seq': 'QVQLVQSGAEVKKPGASVKVSCKASGYTFTGYYIHWVRQAPGQGLEWMGWINPNRGGTNYCQKFQGRVTMTRDTSISTAYMELSRLRSDDTAVYYCAR', 
            'seq': 'QVQLVQSGAEVKKPGASVKVSCKASGYTFTGYYIHWVRQAPGQGLEWMGWINPNRGGTNYCQKFQGRVTMTRDTSISTAYMELSRLRSDDTAVYYCAR', 
            'gap_seq_cdr': 'QVQLVQSGAEVKKPGASVKVSCKAS[GYTFTGYY]IHWVRQAPGQGLEWMGW[INPNRGGT]NYCQKFQGRVTMTRDTSISTAYMELSRLRSDDTAVYYC[AR]', 
            'seq_cdr': 'QVQLVQSGAEVKKPGASVKVSCKAS[GYTFTGYY]IHWVRQAPGQGLEWMGW[INPNRGGT]NYCQKFQGRVTMTRDTSISTAYMELSRLRSDDTAVYYC[AR]', 
            'gap_seq_cdr_mask': 'QVQLVQSGAEVKKPGASVKVSCKAS1IHWVRQAPGQGLEWMGW2NYCQKFQGRVTMTRDTSISTAYMELSRLRSDDTAVYYC', 
            'gap_seq_cdr_trim': 'QVQLVQSGAEVKKPGASVKVSCKAS[GYTFTGYY]IHWVRQAPGQGLEWMGW[INPNRGGT]NYCQKFQGRVTMTRDTSISTAYMELSRLRSDDTAVYYC'}
        ],
        ['6P8M', 2, {'start': 0, 'end': 87, 
            'gap_seq': 'MSQSPSSLAVSVGEKVTMSCKSSQ-------QKNYLAWYQQKPGQSPKLLIYWASTRESGVPDRFTGSGSGTDFTLTISSVKAEDLAVYYCQQY', 
            'seq': 'MSQSPSSLAVSVGEKVTMSCKSSQQKNYLAWYQQKPGQSPKLLIYWASTRESGVPDRFTGSGSGTDFTLTISSVKAEDLAVYYCQQY', 
            'gap_seq_cdr': 'MSQSPSSLAVSVGEKVTMSCKSSQ-------QKNYLAWYQQKPGQSPKLLIY[WAS]TRESGVPDRFTGSGSGTDFTLTISSVKAEDLAVYYC[QQY]', 
            'seq_cdr': 'MSQSPSSLAVSVGEKVTMSCKSSQQKNYLAWYQQKPGQSPKLLIY[WAS]TRESGVPDRFTGSGSGTDFTLTISSVKAEDLAVYYC[QQY]', 
            'gap_seq_cdr_mask': 'MSQSPSSLAVSVGEKVTMSCKSSQ-------QKNYLAWYQQKPGQSPKLLIY1TRESGVPDRFTGSGSGTDFTLTISSVKAEDLAVYYC', 
            'gap_seq_cdr_trim': 'MSQSPSSLAVSVGEKVTMSCKSSQ-------QKNYLAWYQQKPGQSPKLLIY[WAS]TRESGVPDRFTGSGSGTDFTLTISSVKAEDLAVYYC[QQY]'}
        ],
    )
    @unpack
    def test_gap_seq(self, input, idx, expect):
        p = self.c[input]
        blocks = [i for i in p.iter_fmt3()]
        region, aln_seq, pos = p.alignments(blocks[idx])
        result = p.gap_seq(region, aln_seq, pos)
        assert result == expect

    @data(
        ['6P8M', 0, {}],
        ['6P8M', 1, {'CDR1-IMGT': 'GYTFTGYY', 'CDR2-IMGT': 'INPNRGGT', 'CDR3-IMGT': 'AR',
            'FR1-IMGT': 'QVQLVQSGAEVKKPGASVKVSCKAS', 'FR2-IMGT': 'IHWVRQAPGQGLEWMGW',
            'FR3-IMGT': 'NYCQKFQGRVTMTRDTSISTAYMELSRLRSDDTAVYYC',
            'Total': 'QVQLVQSGAEVKKPGASVKVSCKASGYTFTGYYIHWVRQAPGQGLEWMGWINPNRGGTNYCQKFQGRVTMTRDTSISTAYMELSRLRSDDTAVYYCARGINSDYTWDFQHWGQGTLVTVSSASTKGPSVFPLAPSGGTAALGCLVKDYFPEPVTVSWNSGALTSGVHTFPAVLQSSGLYSLSSVVTVPSSSLGTQTYICNVNHKPSNTKVDKRVEP'}
        ],
        ['6P8M', 2, {'CDR1-IMGT': 'QQKNY', 'CDR2-IMGT': 'WAS', 'CDR3-IMGT': 'QQY',
            'FR1-IMGT': 'MSQSPSSLAVSVGEKVTMSCKSS', 'FR2-IMGT': 'LAWYQQKPGQSPKLLIY',
            'FR3-IMGT': 'TRESGVPDRFTGSGSGTDFTLTISSVKAEDLAVYYC',
            'Total': 'MSQSPSSLAVSVGEKVTMSCKSSQQKNYLAWYQQKPGQSPKLLIYWASTRESGVPDRFTGSGSGTDFTLTISSVKAEDLAVYYCQQYETLGSGTKLEIKRTVAAPSVFIFPPSDEQLKSGTASVVCLLNNFYPREAKVQWKVDNALQSGNSQESVTEQDSKDSTYSLSSTLTLSKADYEKHKVYACEVTHQGLSSPVTKSFNRG'}
        ],
        ['6EDZ', 0, {}],
    )
    @unpack
    def test_alignment_summary(self, input, idx, expect):
        p = self.c[input]
        records = [i for i in p.iter_query()]
        query_seq = str(records[idx].seq)
        blocks = [i for i in p.iter_fmt3()]
        alignment = p.alignment_summary(blocks[idx], query_seq)
        result = {k:v['seq'] for k,v in alignment.items()}
        assert result == expect

    @data(
        ['6P8M', {'IGHV1', None, 'IGKV8'}],
        ['6PUH', {'IGKV14S22', 'TRBV6', 'IGHV1', 'TRAV1'}],
    )
    @unpack
    def test_retrieve_hitss(self, input, expect):
        p = self.c[input]
        result = p.retrieve_hits()
        gene_family = set([i.get('gene_family') for i in result])
        assert gene_family == expect
    

    @data(
        ['6P8M',
            {'6P8M_0_H', '6P8M_0_L', '6P8M_0_C'}, 
            {'Total', 'CDR2-IMGT', 'CDR1-IMGT', 'FR2-IMGT', \
                'FR1-IMGT', 'FR3-IMGT', 'CDR3-IMGT', None}
        ],
        ['6PUH',
            {'6PUH_0_D', '6PUH_0_E', '6PUH_0_H', '6PUH_0_B', '6PUH_0_F', \
                '6PUH_0_G', '6PUH_0_A', '6PUH_0_C'},
            {'Total', 'CDR2-IMGT', 'CDR1-IMGT', 'FR2-IMGT', 'FR1-IMGT', \
                'FR3-IMGT', 'CDR3-IMGT'}
        ],
        ['6EDZ', set(), set()]
    )
    @unpack
    def test_retrieve_regions(self, input, expect_chain_ids, expect_regions):
        p = self.c[input]
        result = p.retrieve_regions()
        chain_ids = set([i['query_id'] for i in result])
        assert chain_ids == expect_chain_ids
        regions = set([i.get('region') for i in result])
        assert regions == expect_regions
