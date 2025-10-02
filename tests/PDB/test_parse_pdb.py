import os
from unittest import TestCase
from ddt import ddt, data, unpack
from pprint import pprint

from src.bioomics import Pdb

TESTS_DIR = os.path.dirname(os.path.dirname(__file__))
DATA_DIR = os.path.join(TESTS_DIR, 'data', 'PDB')

@ddt
class TestParsePdb(TestCase):

    def setUp(self):
        self.c = {}
        for file_name in ['7S8L.pdb',]:
            pdb_id = os.path.splitext(file_name)[0]
            args = {
                'pdb_file': os.path.join(DATA_DIR, file_name)
            }
            self.c[pdb_id] = Pdb(args)
    
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
        ['7S8L', [{'model_id': 0, 'chains': [{'chain_id': 'R', 'AA residues': 250, 'seq': 'KETLIPVFLILFIALVGLVGNGFVLWLLGFRMRR' + \
            'NAFSVYVLSLAGADFLFLCFQIINCLVYLSNFFCSISINFPSFFTTVMTCAYLAGLSMLSTVSTERCLSVLWPIWYRCRRPRHLSAVVCVLLWALSLLLSILEGKFCGFLFS' + \
            'DGDSGWCQTFDFITAAWLIFLFMVLCGSSLALLVRILCGLTRLYLTILLTVLVFLLCGLPFGIQWFLILWIWDVLFCHIHPVSVVLSSLNSSANPIIYFFVGSF'}, 
            {'chain_id': 'A', 'AA residues': 6, 'seq': 'KNFFWK'}, {'chain_id': 'B', 'AA residues': 221, 'seq': 'VSAEDKAAAERSK' + \
                'MIDKNLREDGEKARRTLRLLLLGADNSGKSTIVKGIFETKFQVDKVNFHMFDVGRRKWIQCFNDVTAIIFVVDSSDYNRLQEALNDFKSIWNNRWLRTISVILFLNKQ' + \
                'DLLAEKVLAGKSKIEDYFPEFARYTTPEDATPEPGEDPRVTRAKYFIRKEFVDISTASGDGRHICYPHFTCAVDTENARRIFNDCKDIILQMNLREYNLV'}, 
            {'chain_id': 'C', 'AA residues': 337, 'seq': 'LDQLRQEAEQLKNQIRDARKACADATLSQITNNIDPVGRIQMRTRRTLRGHLAKIYAMHWGTDSRLL' + \
                'VSASQDGKLIIWDSYTTNKVHAIPLRSSWVMTCAYAPSGNYVACGGLDNICSIYNLKTREGNVRVSRELAGHTGYLSCCRFLDDNQIVTSSGDTTCALWDIETGQQTT' + \
                'TFTGHTGDVMSLSLAPDTRLFVSGACDASAKLWDVREGMCRQTFTGHESDINAICFFPNGNAFATGSDDATCRLFDLRADQELMTYSHDNIICGITSVSFSKSGRLLL' + \
                'AGYDDFNCNVWDALKADRAGVLAGHDNRVSCLGVTDDGMAVATGSWDSFLKIWN'}, 
            {'chain_id': 'D', 'AA residues': 51, 'seq': 'QARKLVEQLKMEANIDRIKVSKAAADLMAYCEAHAKEDPLLTPVPASENPF'}, 
            {'chain_id': 'E', 'AA residues': 231, 'seq': 'VQLVESGGGLVQPGGSRKLSCSASGFAFSSFGMHWVRQAPEKGLEWVAYISSGSGTIYYADTVKGRF' + \
                'TISRDDPKNTLFLQMTSLRSEDTAMYYCVRSIYYYGSSPFDFWGQGTTLTVSSDIVMTQATSSVPVTPGESVSISCRSSKSLLHSNGNTYLYWFLQRPGQSPQLLIYR' + \
                'MSNLASGVPDRFSGSGSGTAFTLTISRLEAEDVGVYYCMQHLEYPLTFGAGTKLEL'}]}]],
    )
    @unpack
    def test_get_chains(self, input, expect):
        result = self.c[input].get_chains()
        assert result == expect


    '''
    @data(
        ['7S8L', '',],
    )
    @unpack
    def test_(self, input, expect):

        assert result == expect
    '''