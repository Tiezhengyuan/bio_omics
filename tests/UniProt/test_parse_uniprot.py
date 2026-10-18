'''
Test class ParseUniprot
'''
import gzip
import json
import tempfile
from pathlib import Path

from ..helper import *
from src.bioomics import ParseUniprot

DIR_UNIPROT = os.path.join(DATA_DIR, 'UniProt')

with open(os.path.join(DIR_UNIPROT, 'uniprot_sprot.dat')) as f:
    DAT_TEXT = f.read()

UNIREF_TEXT = \
    '>UniRef50_Q96RD9 Signal sequence\n' \
    'MVLSPADKTN\n' \
    '>UniRef50_Q9ZZZ9 hypothetical protein\n' \
    'MSTNETQW\n'


class TestParseUniprot(TestCase):

    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self._tmp.cleanup)
        base = Path(self._tmp.name)
        self.indir = base / 'in'
        self.outdir = base / 'out'
        self.indir.mkdir()
        self.outdir.mkdir()

        # uniprot_sprot.dat.gz: 3 records, Q96RD9 is the IMGT one
        with gzip.open(self.indir / 'uniprot_sprot.dat.gz', 'wt') as f:
            f.write(DAT_TEXT)

        # uniref50.fasta.gz for filter_uniref50
        with gzip.open(self.indir / 'uniref50.fasta.gz', 'wt') as f:
            f.write(UNIREF_TEXT)

        # sprot copy with one DrugBank dbxref injected into Q6GZX4
        lines = DAT_TEXT.splitlines(keepends=True)
        for i, line in enumerate(lines):
            if line.startswith('DR   EMBL; AY548484;'):
                lines.insert(i + 1, 'DR   DrugBank; DB11111; -.\n')
                break
        with gzip.open(self.indir / 'uniprot_sprot_drug.dat.gz', 'wt') as f:
            f.writelines(lines)

    def test_parse_imgt(self):
        c = ParseUniprot(self.indir)
        outfile = c.parse_imgt(self.outdir)
        assert outfile == os.path.join(self.outdir, 'uniprot_sprot_imgt.json')
        with open(outfile) as f:
            data = json.load(f)
        # Q96RD9 carries 'Immunoglobulin domain' in keywords
        assert [r['id'] for r in data] == ['Q96RD9']
        assert set(data[0].keys()) == {'id', 'name', 'pro_seq', 'molecule_type',
            'organism', 'dbxrefs', 'features', 'keywords'}
        assert 'Immunoglobulin domain' in data[0]['keywords']

    def test_parse_drug(self):
        c = ParseUniprot(self.indir)
        outfile = c.parse_drug(self.outdir, name='sprot_drug')
        assert outfile == os.path.join(self.outdir, 'uniprot_sprot_drug_drug.json')
        with open(outfile) as f:
            data = json.load(f)
        # only Q6GZX4 carries the injected DrugBank dbxref
        assert [r['id'] for r in data] == ['Q6GZX4']
        assert 'DrugBank' in data[0]['dbxrefs']
        assert data[0]['dbxrefs']['DrugBank'] == ['DB11111']

    def test_filter_uniref50(self):
        c = ParseUniprot(self.indir)
        index_meta = {
            'Q96RD9': {'key': 'Q96RD9'},
            'Q9XXXX': {'key': 'Q9XXXX'},  # not in uniref50.fasta.gz
        }
        res = c.filter_uniref50(index_meta)
        assert res == {'Q96RD9': {'protein': 'MVLSPADKTN'}}
