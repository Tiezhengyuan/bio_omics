import os
from unittest import TestCase, skip
from ddt import ddt, data, unpack
from Bio.Seq import Seq
# from Bio.SeqRecord import SeqRecord

from src.bioomics import ImgtLigmdb

TESTS_DIR = os.path.dirname(os.path.dirname(__file__))
DATA_DIR = os.path.join(TESTS_DIR, 'data', 'IMGT')

@ddt
class TestImgtLigmdb(TestCase):

    def setUp(self):
        data_dir = '/home/yuan/data/IMGT'
        self.obj = ImgtLigmdb(data_dir)
    
    @skip
    def test_A00673(self):
        infile = os.path.join(DATA_DIR, 'LIGM-DB', 'A00673.dat')
        iter_block = self.obj.iter_imgt_dat(infile)
        
        block = next(iter_block)
        assert block[0] == [
            'ID   A00673; SV 1; linear; unassigned DNA; PAT; SYN; 45 BP.'
        ]
        # accession number
        assert block[1] == ['AC   A00673;']
        assert block[2] == [
            'DT   24-MAY-2002 (Rel. 200221-5, arrived in LIGM-DB)',
            'DT   14-JAN-2013 (Rel. 201303-1, Last updated, Version 4)'
        ]
        assert block[3] == [
            'DE   Artificial sequence for plasmid pSV-V-NP gamma-SNase'
        ]
        assert block[4] == [
            'KW   antigen receptor; Immunoglobulin superfamily (IgSF); immunoglobulin (IG);',
            'KW   IG-Heavy; patent; regular; unassigned DNA; synthetic construct.'
        ]
        # organism species
        assert block[5] == [
            'OS   synthetic construct (synthetic construct)',
            'OC   other sequences; artificial sequences; synthetic construct.'
        ]
        # reference number, reference authors, reference title
        assert block[6] == [
            'RN   [1]', 'RP   1-45', 'RA   Neuberger M.S.;',
            'RT   "PRODUCTION OF CHIMERIC ANTIBODIES";',
            'RL   Patent number WO8601533-A/2, 13-MAR-1986.',
            'RL   CellTech Ltd..'
        ]
        assert block[7] == ['DR   ENA; A00673.']
        assert block[8] == [
            'CC   IMGT/LIGM-DB annotation level: keyword level'
        ]
        # feature header, feature table
        assert block[9] == [
            'FH   Key             Location/Qualifiers (from EMBL)',
            'FH',
            'FT   source          1..45',
            'FT                   /organism="synthetic construct"',
            'FT                   /mol_type="unassigned DNA"',
            'FT                   /db_xref="taxon:32630"',
            'FT   exon            8..>45',
            'FT   intron          <1..7'
        ]
        assert block[10] == [
            'SQ   Sequence 45 BP; 14 A; 13 C; 7 G; 11 T; 0 other;',
            '     tcatcagctc ctaacctcga ggatccaaca gtatatagtg caact                        45'
        ]

        record = self.obj.parse_record(block)
        assert record['id'] == ['A00673', 'SV 1', 'linear', \
            'unassigned DNA', 'PAT', 'SYN', '45 BP.']
        assert record['date'] == [
            '24-MAY-2002 (Rel. 200221-5, arrived in LIGM-DB)',
            '14-JAN-2013 (Rel. 201303-1, Last updated, Version 4)'
        ]
        assert record['description'] == \
            'Artificial sequence for plasmid pSV-V-NP gamma-SNase'
        assert record['keywords'] == [
            'antigen receptor',
            'Immunoglobulin superfamily (IgSF)',
            'immunoglobulin (IG)',
            'IG-Heavy',
            'patent',
            'regular',
            'unassigned DNA',
            'synthetic construct'
        ]
        assert record['features'] == [
            {'name': 'source', 'start': '1', 'end': '45', \
                'organism': 'synthetic construct', \
                'mol_type': 'unassigned DNA', \
                'db_xref': 'taxon:32630'},
            {'name': 'exon', 'start': '8', 'end': '>45'},
            {'name': 'intron', 'start': '<1', 'end': '7'}  
        ]
        assert record['sequence']['header'] == ['Sequence 45 BP', \
            '14 A', '13 C', '7 G', '11 T', '0 other;']
        assert len(record['sequence']['seq']) == 45

    @skip
    def test_Z99606(self):
        infile = os.path.join(DATA_DIR, 'LIGM-DB', 'Z99606.dat')
        iter_block = self.obj.iter_imgt_dat(infile)
        block = next(iter_block)
        assert block[0] == ['ID   Z99606; SV 1; linear; gDNA; STD; HUM; 319 BP.']
        assert block[-1] == [
            'SQ   Sequence 319 BP; 72 A; 88 C; 92 G; 67 T; 0 other;',
            '     agctacagca gtggggcgca ggactgttga agccttcgga gaccctgtcc ctcacctgcg        60',
            '     ctgtctatgg tgggtccttc agtggttact actggagctg gatccgccag cccccaggga       120',
            '     aggggctgga gtggattggg gaaatcaatc atagtggaag caccaactac aacccgtccc       180', '     tctagagtcg agtcaccata tcagtagaca cgtccaagaa ccagttctcc ctgaagctga       240', '     gctctgtgac cgccacggac acggctgtgt attactgtgc gagaggccgt ctaacaatct       300',
            '     ttgactactg gggccaagg                                                    319'
        ]

        record = self.obj.parse_record(block)
        assert record['id'] == ['Z99606', 'SV 1', 'linear', \
            'gDNA', 'STD', 'HUM', '319 BP.']
        assert record['description'] == 'Homo sapiens rearranged gene' + \
            ' for Ig heavy chain, VH4.21-D-JH4, clone F4J'
        assert record['features'][0] =={
            'name': 'source',
            'start': '1',
            'end': '319',
            'organism': 'Homo sapiens',
            'mol_type': 'genomic DNA',
            'dev_stage': 'fetal',
            'clone': 'F4J',
            'cell_type': 'B cell',
            'tissue_type': 'gut',
            'rearranged': True,
            'note': 'PCR amplified genomic DNA',
            'db_xref': 'taxon:9606'
        }
        assert record['features'][1] == {
            'name': 'V_region',
            'start': '1',
            'end': '319',
            'product': 'immunoglobulin rearranged heavy chain, VH4.21-D-JH4'
        }
        assert record['features'][2] == {
            'name': 'V_segment',
            'start': '1',
            'end': '284',
            'note': 'VH region'
        }
        assert record['features'][3] == {
            'name': 'D_segment',
            'start': '285',
            'end': '298',
            'note': 'DH/N nucleotide region'
        }
        assert record['features'][4] == {
            'name': 'J_segment',
            'start': '299',
            'end': '319',
            'note': 'JH region'
        }
        assert record['sequence']['header'] == ['Sequence 319 BP', \
            '72 A', '88 C', '92 G', '67 T', '0 other;']
        assert len(record['sequence']['seq']) == 319

    @skip
    def test_A02856(self):
        infile = os.path.join(DATA_DIR, 'LIGM-DB', 'A02856.dat')
        iter_block = self.obj.iter_imgt_dat(infile)
        block = next(iter_block)
        record = self.obj.parse_record(block)

        assert record['features'][0] == {
            'name': 'source',
            'start': '1',
            'end': '894',
            'organism': 'synthetic construct',
            'mol_type': 'unassigned DNA',
            'db_xref': 'taxon:32630'
        }
        assert record['features'][1] == {
            'name': 'CDS',
            'start': '184',
            'end': '888',
            'transl_table': '11',
            'product': 'ZZ-IGF-1',
            'protein_id': 'CAA00282.1',
            'translation': 'MKKKNIYSIRKLGVGIASVTLGTLLISGXVTPAANAAQHDEAVDN' +\
                'KFNKEQQNAFYEILHLPNLNEEQRNAFIQSLKDDPSQSANLLAEAKKLNDAQAPKVDNK' +\
                'FNKEQQNAFYEILHLPNLNEEQRNAFIQSLKDDPSQSANLLAEAKKLNDAQAPKVDANS' +\
                'NGPETLCGAELVDPLQFVCGDRGFYFNKPTAYGSSSRRAPQTGIVDECCFRSCDLRRLEMYCAPLKPAKSA'
        }
        assert record['sequence']['header'] == ['Sequence 894 BP', \
            '295 A', '185 C', '157 G', '256 T', '1 other;']
        assert len(record['sequence']['seq']) == 894

    @skip
    def test_A16547(self):
        infile = os.path.join(DATA_DIR, 'LIGM-DB', 'A16547.dat')
        iter_block = self.obj.iter_imgt_dat(infile)
        block = next(iter_block)
        record = self.obj.parse_record(block)

        result = [i for i in record['features'] if i['name']=='V-REGION'][0]
        assert result == {
            'name': 'V-REGION',
            'start': '95',
            'end': '379',
            'translation': 'DAGVIQSPRHEVTEMGQEVTLRCKPISGHNSLFWYRQTMMR' +\
                'GLELLIYFNNNVPIDDSGMPEDRFSAKMPNASFSTLKIQPSEPRDSAVYFCASS',
            'CDR_length': '[5.6.14]',
            'gene': 'TRBV12-3',
            'putative_limit': "3' side"
        }
        start, end = int(result['start']) - 1, int(result['end'])
        nt_seq = record['sequence']['seq'][start:end]
        aa_seq = Seq(nt_seq).translate()
        assert aa_seq == result['translation']

        result = [i for i in record['features'] if i['name']=='CDR1-IMGT'][0]
        assert result == {
            'name': 'CDR1-IMGT',
            'start': '173',
            'end': '187',
            'translation': 'SGHNS',
            'AA_IMGT': 'AA 27 to 38, AA 30, 31, 32, 33, 34, 35 and 36 are missing'
        }
        start, end = int(result['start']) - 1, int(result['end'])
        nt_seq = record['sequence']['seq'][start:end]
        assert nt_seq == 'TCAGGCCACAACTCC'
        aa_seq = Seq(nt_seq).translate()
        assert aa_seq == result['translation']

    def test_AF173903(self):
        infile = os.path.join(DATA_DIR, 'LIGM-DB', 'AF173903.dat')
        iter_block = self.obj.iter_imgt_dat(infile)
        block = next(iter_block)
        record = self.obj.parse_record(block)

        result = [i for i in record['features'] if i['name']=='V-GENE'][0]
        assert result['IMGT_allele'] == 'IGHV3-13*01'
        assert result['IMGT_gene'] == 'IGHV3-13'

        result = [i for i in record['features'] if i['name']=='V-REGION'][0]
        assert result == {
            'name': 'V-REGION',
            'start': '157',
            'end': '456',
            'translation': 'EVQLVESGGGLVQPGGSLRLSCAASGFTFSNYYMHWVRQAQGK' +\
                'GLEWVGLIRNKANSYTTEYAAAVKGRFTISRDDSKNTLYLQMSSLKTEDTALYYCTK',
            'CDR_length': '[8.10.2]',
            'FR_length': '[25.17.38.G]',
            'functional': True,
            'IMGT_allele': 'IGHV3-13*01',
            'IMGT_gene': 'IGHV3-13'
        }
        start, end = int(result['start']) - 1, int(result['end'])
        nt_seq = record['sequence']['seq'][start:end]
        aa_seq = Seq(nt_seq).translate()
        assert aa_seq == result['translation']
