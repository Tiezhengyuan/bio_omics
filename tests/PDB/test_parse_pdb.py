import os
from unittest import TestCase
from ddt import ddt, data, unpack
from pprint import pprint

from ..helper import *
from src.bioomics import Pdb


@ddt
class TestParsePdb(TestCase):

    def setUp(self):
        self.c = {}
        for file_name in ['7S8L.pdb',]:
            pdb_file = os.path.join(DATA_DIR, 'PDB', file_name)
            pdb_id = os.path.splitext(file_name)[0]
            self.c[pdb_id] = Pdb(pdb_file, RESULTS_DIR, True)
    
    @data(
        ['7S8L', {
            '1': {'misc': '', 'molecule': 'mas-related g-protein coupled receptor member x2', 'chain': 'r', 'engineered': 'yes'}, 
            '2': {'misc': '', 'molecule': 'cortistatin 14', 'chain': 'a', 'engineered': 'yes'}, 
            '3': {'misc': '', 'molecule': 'gs-mini-gq chimera', 'chain': 'b', 'engineered': 'yes'}, 
            '4': {'misc': '', 'molecule': 'guanine nucleotide-binding protein g(i)/g(s)/g(t) subunitbeta-1 ', 
                'chain': 'c', 'synonym': 'transducin beta chain 1', 'engineered': 'yes'}, 
            '5': {'misc': '', 'molecule': 'guanine nucleotide-binding protein g(i)/g(s)/g(o) subunitgamma-2 ', 
                'chain': 'd', 'synonym': 'g gamma-i', 'engineered': 'yes'}, 
            '6': {'misc': '', 'molecule': 'scfv16', 'chain': 'e', 'engineered': 'yes'}
        }],
    )
    @unpack
    def test_get_header(self, input, expect_compound):
        result = self.c[input].get_header()
        assert result['compound'] == expect_compound

    @data(
        ['7S8L', {0: {'R': 250, 'A': 6, 'B': 221, 'C': 337, 'D': 51, 'E': 231}}],
    )
    @unpack
    def test_iter_residues(self, input, expect):
        result = {}
        for model, chain, res in self.c[input].iter_res():
            if model.id not in result:
                result[model.id] = {}
            if chain.id not in result[model.id]:
                result[model.id][chain.id] = 0
            result[model.id][chain.id] +=1
        assert result == expect

    @data(
        ['7S8L', [{'model_no': 0, 'chains': [
            {'chain_no': 'R', 'AA residues': 250, 'seq': 'KETLIPVFLILFIALVGLVGNGFVLWLLGFRMRRNAFSVYVLSLAGADF'+\
                'LFLCFQIINCLVYLSNFFCSISINFPSFFTTVMTCAYLAGLSMLSTVSTERCLSVLWPIWYRCRRPRHLSAVVCVLLWALSLLLSILEGK'+\
                'FCGFLFSDGDSGWCQTFDFITAAWLIFLFMVLCGSSLALLVRILCGLTRLYLTILLTVLVFLLCGLPFGIQWFLILWIWDVLFCHIHPVS'+\
                'VVLSSLNSSANPIIYFFVGSF'}, 
            {'chain_no': 'A', 'AA residues': 6, 'seq': 'KNFFWK'}, {'chain_no': 'B', 'AA residues': 221, 
                'seq': 'VSAEDKAAAERSKMIDKNLREDGEKARRTLRLLLLGADNSGKSTIVKGIFETKFQVDKVNFHMFDVGRRKWIQCFNDVTAIIF'+\
                'VVDSSDYNRLQEALNDFKSIWNNRWLRTISVILFLNKQDLLAEKVLAGKSKIEDYFPEFARYTTPEDATPEPGEDPRVTRAKYFIRKEFV'+\
                'DISTASGDGRHICYPHFTCAVDTENARRIFNDCKDIILQMNLREYNLV'}, 
            {'chain_no': 'C', 'AA residues': 337, 'seq': 'LDQLRQEAEQLKNQIRDARKACADATLSQITNNIDPVGRIQMRTRRTLR'+\
                'GHLAKIYAMHWGTDSRLLVSASQDGKLIIWDSYTTNKVHAIPLRSSWVMTCAYAPSGNYVACGGLDNICSIYNLKTREGNVRVSRELAGHT'+\
                'GYLSCCRFLDDNQIVTSSGDTTCALWDIETGQQTTTFTGHTGDVMSLSLAPDTRLFVSGACDASAKLWDVREGMCRQTFTGHESDINAICF'+\
                'FPNGNAFATGSDDATCRLFDLRADQELMTYSHDNIICGITSVSFSKSGRLLLAGYDDFNCNVWDALKADRAGVLAGHDNRVSCLGVTDDGM'+\
                'AVATGSWDSFLKIWN'}, 
            {'chain_no': 'D', 'AA residues': 51, 'seq': 'QARKLVEQLKMEANIDRIKVSKAAADLMAYCEAHAKEDPLLTPVPASENPF'}, 
            {'chain_no': 'E', 'AA residues': 231, 'seq': 'VQLVESGGGLVQPGGSRKLSCSASGFAFSSFGMHWVRQAPEKGLEWVAYI'+\
                'SSGSGTIYYADTVKGRFTISRDDPKNTLFLQMTSLRSEDTAMYYCVRSIYYYGSSPFDFWGQGTTLTVSSDIVMTQATSSVPVTPGESVSI'+\
                'SCRSSKSLLHSNGNTYLYWFLQRPGQSPQLLIYRMSNLASGVPDRFSGSGSGTAFTLTISRLEAEDVGVYYCMQHLEYPLTFGAGTKLEL'}]}]
        ],
    )
    @unpack
    def test_get_chains(self, input, expect):
        result = self.c[input].get_chains()
        assert result == expect

    @data(
        ['7S8L', ['R', 'A', 'B', 'C', 'D', 'E'],],
    )
    @unpack
    def test_chain_no(self, input, expect):
        result = self.c[input].get_chain_no()
        print(result)
        assert result == expect

    '''
    @data(
        ['7S8L', '',],
    )
    @unpack
    def test_(self, input, expect):

        assert result == expect
    '''