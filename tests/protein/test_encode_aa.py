'''
example: pytest -s tests/protein/test_encode_aa.py
'''
import numpy as np

from ..helper import *

from src.bioomics import EncodeAA

@ddt
class TestEncodeAA(TestCase):

    @classmethod
    def setUpClass(cls):
        cls.e = EncodeAA()

    def test_physical_chemical_properties(self):
        assert self.e.hydrophobicity_ph7('AW') == [41, 97]
        assert self.e.hydrophobicity('AW') == [0.61, 2.65]
        assert self.e.polarity('AW') == [8.1, 5.4]
        assert self.e.polarizability('AW') == [0.046, 0.409]
        assert self.e.van_der_Waals_volume('AW') == [1, 8.08]

    def test_freq_1aa(self):
        assert self.e.freq_1aa('AAAG') == {'A': 0.75, 'G': 0.25}

    def test_freq_2aa(self):
        # two pairs in a 3-residue sequence, normalized by sequence length
        assert self.e.freq_2aa('AAG') == {'freq_AA': 0.333, 'freq_AG': 0.333}

    def test_frequency_aa(self):
        res = self.e.frequency_aa('AAG')
        # [seq, label] + 20 single aa + 20x20 double aa
        assert len(res) == 2 + 20 + 400
        assert res[0] == 'AAG'
        # single aa frequencies, rounded to 3 decimals
        assert res[2 + self.e.aa.index('A')] == '0.667'
        assert res[2 + self.e.aa.index('G')] == '0.333'
        assert res[2 + self.e.aa.index('C')] == '0'
        # double aa frequencies: freq_2aa keys carry the freq_ prefix
        assert res[2 + 20 + self.e.aa2.index('AA')] == '0.333'
        assert res[2 + 20 + self.e.aa2.index('AG')] == '0.333'
        assert res[2 + 20 + self.e.aa2.index('CC')] == '0'

    def test_physical_chemical(self):
        res = self.e.physical_chemical('AAG', 1)
        assert res['seq'] == 'AAG'
        assert res['label'] == 1
        assert res['mean_hydro_ph7'] == round(np.mean([41, 41, 0]), 3)
        assert res['median_VDWV'] == round(np.median([1, 1, 0]), 3)
        assert res['variance_hydro'] == round(np.var([0.61, 0.61, 0.07]), 3)

    def test_physical_chemical_text(self):
        res = self.e.physical_chemical_text('AA')
        # seq | mean van der Waals volume
        assert res == 'A A | 1'

    def test_existing(self):
        res = self.e.existing('AAAG')
        assert res[0] == 'AAAG'
        assert res[1] is None
        # single aa: A exists, C does not
        assert res[2 + self.e.aa.index('A')] == '1'
        assert res[2 + self.e.aa.index('C')] == '0'
        # double aa: AA exists, CC does not
        assert res[2 + 20 + self.e.aa2.index('AA')] == '1'
        assert res[2 + 20 + self.e.aa2.index('CC')] == '0'

    def test_vector_1d(self):
        res = self.e.vector_1d('AAG')
        # 5 properties per residue
        assert len(res) == 15
        assert list(res[:3]) == [41, 41, 0]
        # label is appended when given
        assert len(self.e.vector_1d('AAG', 1)) == 16

    def test_vector_2d(self):
        res = self.e.vector_2d('AAG')
        assert res.shape == (3, 5)
        assert list(res[:, 0]) == [41, 41, 0]

    def test_aa_smiles(self):
        # the SMILES codes of constants.SMILES_1 are single atom chars
        # 'GAG': glycine is dropped (no side chain), the first atom char is
        # stripped by the leading-char rule
        assert self.e.aa_smiles('GAG') == ''
        assert self.e.aa_smiles('AA') == ' C'
        # unknown residues are not supported
        assert self.e.aa_smiles('AGX') is None