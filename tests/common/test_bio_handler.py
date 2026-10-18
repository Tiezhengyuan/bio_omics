'''
Test src/bioomics/common/bio_handler.py
'''
from ..helper import *
import gzip
import tempfile
from pathlib import Path

from src.bioomics import BioHandler

FASTA_TEXT = '>r1 desc one\nMKVLWA\n>r2 desc two\nAAAAAA\n'

GBK_TEXT = '''LOCUS       SCU49845     5028 bp    DNA     circular BCT 30-JUL-2019
DEFINITION  Saccharomyces cerevisiae TPM1 gene.
ACCESSION   SCU49845
VERSION     SCU49845.1
KEYWORDS    .
SOURCE      Saccharomyces cerevisiae
  ORGANISM  Saccharomyces cerevisiae
COMMENT     test comment.
FEATURES             Location/Qualifiers
     source          1..5028
ORIGIN
        1 atgaaattta aaaaatgcaa
//'''


class TestBioHandler(TestCase):

    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.tmp, ignore_errors=True)

    def _write(self, name, text, gzip_it=False):
        path = os.path.join(self.tmp, name)
        if gzip_it:
            with gzip.open(path, 'wt') as f:
                f.write(text)
        else:
            with open(path, 'w') as f:
                f.write(text)
        return path

    # ----- parse_fasta -----
    def test_parse_fasta_plain(self):
        path = self._write('data.fasta', FASTA_TEXT)
        records = list(BioHandler.parse_fasta(path))
        self.assertEqual([r.id for r in records], ['r1', 'r2'])
        self.assertEqual(str(records[0].seq), 'MKVLWA')

    def test_parse_fasta_gz(self):
        path = self._write('data.fasta.gz', FASTA_TEXT, gzip_it=True)
        records = list(BioHandler.parse_fasta(path))
        self.assertEqual([r.id for r in records], ['r1', 'r2'])

    def test_parse_fasta_accepts_pathlib(self):
        # NB: Path has no .endswith; inputs are coerced to str
        path = Path(self._write('data.fasta', FASTA_TEXT))
        records = list(BioHandler.parse_fasta(path))
        self.assertEqual(len(records), 2)

    # ----- parse_gbk -----
    def test_parse_gbk_plain(self):
        path = self._write('data.gbk', GBK_TEXT)
        record = next(BioHandler.parse_gbk(path))
        self.assertEqual(record.locus, 'SCU49845')
        self.assertEqual(record.accession, ['SCU49845'])
        self.assertTrue(record.sequence.startswith('ATG'))

    def test_parse_gbk_gz(self):
        path = self._write('data.gbk.gz', GBK_TEXT, gzip_it=True)
        record = next(BioHandler.parse_gbk(path))
        self.assertEqual(record.locus, 'SCU49845')

    def test_parse_gbk_accepts_pathlib(self):
        path = Path(self._write('data.gbk', GBK_TEXT))
        record = next(BioHandler.parse_gbk(path))
        self.assertEqual(record.locus, 'SCU49845')
