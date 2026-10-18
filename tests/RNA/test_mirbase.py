'''
MiRBase
download hairpin and mature sequences, the download pages are HTML
'''
import io
from unittest import mock

from Bio import SeqIO

from ..helper import *
from src.bioomics import Mirbase


class TestMirbase(TestCase):
    def setUp(self):
        self.c = Mirbase(DATA_DIR, False)

    def test_local_dir(self):
        # the local dir is <local_dir>/miRBase
        assert self.c.local_dir == os.path.join(DATA_DIR, 'miRBase')

    def test_download_hairpin(self):
        # hairpin.fa exists in the test data: no download with overwrite=False
        local_path, local_file = self.c.download_hairpin()
        assert local_path == os.path.join(DATA_DIR, 'miRBase')
        assert local_file == os.path.join(DATA_DIR, 'miRBase', 'hairpin.fa')
        assert os.path.isfile(local_file)

    def test_download_mature(self):
        # mature.fa exists in the test data: no download with overwrite=False
        local_path, local_file = self.c.download_mature()
        assert local_path == os.path.join(DATA_DIR, 'miRBase')
        assert local_file == os.path.join(DATA_DIR, 'miRBase', 'mature.fa')
        assert os.path.isfile(local_file)

    def test_download_overwrite(self):
        # overwrite=True forces parse_fasta even if the file exists
        c = Mirbase(DATA_DIR, True)
        with mock.patch.object(Mirbase, 'parse_fasta') as m:
            local_path, local_file = c.download_hairpin()
        assert m.call_count == 1
        assert m.call_args[0] == (
            'https://mirbase.org/download/CURRENT/hairpin.fa',
            os.path.join(DATA_DIR, 'miRBase', 'hairpin.fa'))

    def test_parse_fasta(self):
        # the mirbase download pages are HTML: <br> tags are turned into
        # line breaks before the fasta content is parsed; the page starts
        # directly with the first record (no text before it)
        page = (
            '<html><body>'
            '&gt;hsa-let-7a-1 MI0000001 hsa-let-7a-1 stem-loop'
            '<br>UGAGGUAGUAGGUUGUAUAGUU'
            '<br>&gt;hsa-let-7a-2 MI0000002 hsa-let-7a-2 stem-loop'
            '<br>UGAGGUAGUAGGUUGUAUAGUUU'
            '<br></body></html>'
        )
        os.makedirs(RESULTS_DIR, exist_ok=True)
        local_file = os.path.join(RESULTS_DIR, 'mirbase_parse.fa')
        with mock.patch('src.bioomics.RNA.mirbase.urlopen',
                return_value=io.BytesIO(page.encode())):
            self.c.parse_fasta('http://example.com/hairpin.fa', local_file)
        with open(local_file) as f:
            records = list(SeqIO.parse(f, 'fasta'))
        assert [(r.id, str(r.seq)) for r in records] == [
            ('hsa-let-7a-1', 'UGAGGUAGUAGGUUGUAUAGUU'),
            ('hsa-let-7a-2', 'UGAGGUAGUAGGUUGUAUAGUUU'),
        ]
        os.remove(local_file)

