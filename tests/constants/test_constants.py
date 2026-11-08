'''
Sanity tests for the constants in src/bioomics/constants
'''
from ..helper import *

import importlib

# `from src.bioomics.constants import amino_acids` would bind the amino_acids
# LIST (re-exported by the package __init__ and shadowing the submodule name),
# so resolve the module through sys.modules instead.
aa_mod = importlib.import_module('src.bioomics.constants.amino_acids')
from src.bioomics.constants import iedb as iedb_const
from src.bioomics.constants.constants import AA, PROPERTY, SMILES, SMILES_1

STANDARD_AA = set('ACDEFGHIKLMNPQRSTVWY')


@ddt
class TestAminoAcids(TestCase):

    def test_20_standard_amino_acids(self):
        self.assertEqual(set(aa_mod.amino_acids), STANDARD_AA)

    def test_longer_names_covers_standard(self):
        self.assertEqual(
            set(aa_mod.longer_names.values()), STANDARD_AA)
        for three in aa_mod.longer_names.keys():
            self.assertEqual(len(three), 3)

    def test_one_letter_names_covers_standard(self):
        self.assertEqual(set(aa_mod.one_letter_names.keys()), STANDARD_AA)

    def test_longer_and_one_letter_names_are_inverse(self):
        for three, one in aa_mod.longer_names.items():
            self.assertEqual(aa_mod.one_letter_names[one], three)

    def test_short_to_long_roundtrip(self):
        # short_to_long mirrors one_letter_names (one-letter -> three-letter)
        for one, three in aa_mod.one_letter_names.items():
            self.assertEqual(aa_mod.short_to_long[one], three)

    def test_hpcg_covers_standard(self):
        self.assertEqual(set(aa_mod.hpcg.keys()), STANDARD_AA)
        self.assertLessEqual(set(aa_mod.hpcg.values()), {'H', 'C', 'P', 'G'})


class TestConstants(TestCase):

    def test_aa_covers_standard(self):
        self.assertEqual(set(AA.keys()), STANDARD_AA)

    def test_aa_fields_and_abbreviations(self):
        for one, rec in AA.items():
            self.assertEqual(set(rec.keys()), {'abb', 'hydrophobicity'})
            self.assertEqual(rec['abb'], aa_mod.one_letter_names[one].lower())

    def test_property_covers_standard(self):
        self.assertEqual(set(PROPERTY.keys()), STANDARD_AA)

    def test_property_matches_hydrophobicity_ph7(self):
        for one in STANDARD_AA:
            self.assertEqual(
                PROPERTY[one]['hydrophobicity_ph7'], AA[one]['hydrophobicity'])

    def test_smiles_covers_standard(self):
        self.assertEqual(set(SMILES.keys()), STANDARD_AA)
        self.assertEqual(set(SMILES_1.keys()), STANDARD_AA)
        for one, smi in SMILES.items():
            self.assertIsInstance(smi, str)
            self.assertGreater(len(smi), 0)


@ddt
class TestIedbConstants(TestCase):

    ENTITY_COLUMNS = (
        'antigen', 'epitope', 'tcr', 'bcr', 'reference',
        'tcell', 'bcell', 'mhc',
    )

    @data(*ENTITY_COLUMNS)
    def test_entity_columns(self, entity):
        cols = getattr(iedb_const, entity)
        self.assertIsInstance(cols, dict)
        self.assertGreater(len(cols), 0)
        for k, v in cols.items():
            self.assertIsInstance(k, str)
            self.assertTrue(k)
            self.assertTrue(v is None or isinstance(v, (str, int, float)))
