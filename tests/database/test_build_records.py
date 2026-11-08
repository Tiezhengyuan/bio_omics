'''
Unit tests for src/bioomics/database/build_records.py (BuildRecords)

Only the pure record-building helper is tested; table insertion requires
a live MySQL database (see tests/database/test_insert.py).
'''
from types import SimpleNamespace

from ..helper import *

from src.bioomics.database.build_records import BuildRecords


class TestBuildRecords(TestCase):

    def setUp(self):
        self.br = BuildRecords(False)

    def test_init(self):
        self.assertFalse(self.br.verbose)
        self.assertTrue(BuildRecords(True).verbose)

    def test_uniprot_alphafold_plddt(self):
        plddt_data = {
            'uniprot_acc': 'P02649',
            'pfile': '/tmp/AF-P02649-F1-model_v6.pdb',
            'pLDDT': [
                {
                    'model': SimpleNamespace(id=0),
                    'chain': SimpleNamespace(id='A'),
                    'stats': {'mean_plddt': 78.5},
                },
                {
                    'model': SimpleNamespace(id=0),
                    'chain': SimpleNamespace(id='B'),
                    'stats': {},
                },
            ],
        }
        res = self.br.uniprot_alphafold_plddt(plddt_data)
        self.assertEqual(res, [
            ('P02649', 0, 'A', 78.5, '/tmp/AF-P02649-F1-model_v6.pdb'),
            # missing mean_plddt becomes None
            ('P02649', 0, 'B', None, '/tmp/AF-P02649-F1-model_v6.pdb'),
        ])

    def test_uniprot_alphafold_plddt_empty(self):
        res = self.br.uniprot_alphafold_plddt({
            'uniprot_acc': 'P12345',
            'pfile': '/tmp/x.pdb',
            'pLDDT': [],
        })
        self.assertEqual(res, [])
