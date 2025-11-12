from ..helper import *
from src.bioomics import Pdb


@ddt
class TestParsePdb(TestCase):

    def setUp(self):
        self.c = {}
        for pdb_id, file_name in [('7S8L', '7S8L.pdb'),]:
            pdb_file = os.path.join(DATA_DIR, 'PDB', file_name)
            obj = Pdb(pdb_file, RESULTS_DIR, True)
            obj.load_structure()
            self.c[pdb_id] = obj
    
    @data(
        ['7S8L', os.path.join(RESULTS_DIR, '7S', '7S8L', 'chain_sequences.faa')],
    )
    @unpack
    def test_export_fasta(self, pdb_id, expect):
        result = self.c[pdb_id].export_fasta()
        assert result[1] == expect

    @data(
        ['7S8L', 6],
    )
    @unpack
    def test_export_chains_fasta(self, pdb_id, expect):
        result = self.c[pdb_id].export_chain_faa()
        assert len(result) == expect

    @data(
        ['7S8L', 6],
    )
    @unpack
    def test_export_chain_pdb(self, pdb_id, expect):
        result = self.c[pdb_id].export_chain_pdb(overwrite=True)
        assert len(result) == expect
