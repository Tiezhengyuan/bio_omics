'''
example: pytest -s tests/protein/test_process_fasta.py
'''
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from ..helper import *

from src.bioomics import ProcessFasta

@ddt
class TestProcessFasta(TestCase):

    def setUp(self):
        self.outdir = os.path.join(RESULTS_DIR, 'fasta')
        os.makedirs(self.outdir, exist_ok=True)
        self.fa_file = os.path.join(self.outdir, 'test.faa')
        if os.path.isfile(self.fa_file):
            os.remove(self.fa_file)

    def test_one_seq_roundtrip(self):
        p = ProcessFasta(self.fa_file)
        assert p.one_seq('ACDEFG', 'seq1') == 'succeed'
        records = p.load_fasta()
        assert len(records) == 1
        assert records[0].id == 'seq1'
        assert str(records[0].seq) == 'ACDEFG'
        assert p.get_first_seq() == 'ACDEFG'

    def test_row_seq(self):
        rows = [
            {'id': 1, 'seq': 'AAA'},
            {'id': 'b', 'seq': 'CCC', 'description': 'some description'},
        ]
        p = ProcessFasta(self.fa_file)
        assert p.row_seq(rows) == 'succeed'
        records = p.load_fasta()
        assert [r.id for r in records] == ['1', 'b']
        # the fasta header is '<id> <description>'
        assert records[1].description == 'b some description'

    def test_wrong_fasta(self):
        p = ProcessFasta(self.fa_file)
        # file does not exist
        assert p.wrong_fasta() is True
        p.one_seq('AAA', 'x')
        assert p.wrong_fasta() is False
        # an empty sequence record makes the fasta wrong
        p.to_fasta([SeqRecord(Seq(''), id='empty'), SeqRecord(Seq('MMP'), id='ok')])
        assert p.wrong_fasta() is True

    def test_clean_fasta(self):
        p = ProcessFasta(self.fa_file)
        p.to_fasta([SeqRecord(Seq(''), id='empty'), SeqRecord(Seq('MMP'), id='ok')])
        # one empty record was removed
        assert p.clean_fasta() is True
        assert p.wrong_fasta() is False
        records = p.load_fasta()
        assert len(records) == 1
        assert records[0].id == 'ok'
        # nothing to remove
        assert p.clean_fasta() is False

    def test_unique_sequences(self):
        df = pd.DataFrame({'id': ['a', 'b', 'c'], 'seq': ['AA', 'BB', 'AA']})
        p = ProcessFasta(self.fa_file)
        p.unique_sequences(df)
        records = p.load_fasta()
        assert len(records) == 2
        # ids of duplicated sequences are joined by comma; groups sorted by seq
        by_id = {r.id: str(r.seq) for r in records}
        assert by_id == {'a,c': 'AA', 'b': 'BB'}