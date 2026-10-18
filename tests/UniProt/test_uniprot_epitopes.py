'''
Test class UniprotEpitope
'''
import gzip
import tempfile

from ..helper import *
from src.bioomics import UniprotEpitope

DAT_FILE = os.path.join(DATA_DIR, 'UniProt', 'uniprot_sprot.dat')

# epitope feature injected into Q6GZX4 (REGION 10..20)
EPITOPE_LINES = [
    'FT   REGION          10..20\n',
    'FT                   /note="epitope: synthetic peptide test"\n',
]


@ddt
class TestUniprotEpitope(TestCase):

    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self._tmp.cleanup)
        self.tmp_dir = self._tmp.name
        self.c = UniprotEpitope(
            os.path.join(self.tmp_dir, 'data'),
            os.path.join(self.tmp_dir, 'results'),
        )
        with open(DAT_FILE) as f:
            self.dat_lines = f.readlines()

    def _write_dat_gz(self, file_name, lines):
        dat_gz = os.path.join(self.tmp_dir, file_name)
        with gzip.open(dat_gz, 'wt') as f:
            f.writelines(lines)
        return dat_gz

    def _dat_lines_with_epitope(self):
        # insert the epitope feature right after Q6GZX4's CHAIN /id line
        lines = list(self.dat_lines)
        for i, line in enumerate(lines):
            if line.startswith('FT                   /id="PRO_0000410512"'):
                lines[i + 1:i + 1] = EPITOPE_LINES
                break
        return lines

    def test_init(self):
        assert self.c.source == 'UniProt'
        assert self.c.meta == {}
        assert self.c.overwrite is False

    def test_parse_epitope_iedb_only(self):
        # no record carries an epitope feature: IEDB accessions present in
        # uniprot get an empty epitopes entry, absent ones are kept
        dat_gz = self._write_dat_gz('uniprot_sprot.dat.gz', self.dat_lines)
        iebd = {'Q96RD9': {'key': 'Q96RD9'}, 'Q9XXXX': {'key': 'Q9XXXX'}}
        entity_data = self.c.parse_epitope(dat_gz, iebd)
        assert list(entity_data.keys()) == ['Q96RD9']
        assert entity_data['Q96RD9']['accession'] == 'Q96RD9'
        assert entity_data['Q96RD9']['epitopes'] == []
        # the found accession is consumed, the absent one is kept
        assert list(iebd.keys()) == ['Q9XXXX']
        assert self.c.meta['uniprot_sprot'] == {
            'proteins_with_epitopes': 0,
            'epitopes': 0,
            'proteins_in_IEDB': 1,
        }

    def test_parse_epitope_with_epitope(self):
        lines = self._dat_lines_with_epitope()
        dat_gz = self._write_dat_gz('uniprot_sprot.dat.gz', lines)
        entity_data = self.c.parse_epitope(dat_gz, {})
        assert list(entity_data.keys()) == ['Q6GZX4']
        assert entity_data['Q6GZX4']['accession'] == 'Q6GZX4'
        epitopes = entity_data['Q6GZX4']['epitopes']
        assert len(epitopes) == 1
        # qualifier values are str (Biopython >= 1.85) or list (before)
        note = epitopes[0]['qualifiers'].get('note')
        if isinstance(note, (list, tuple)):
            note = ' '.join(note)
        assert note == 'epitope: synthetic peptide test'
        assert epitopes[0]['seq_start'] == 9
        assert epitopes[0]['seq_end'] == 20
        assert epitopes[0]['seq'] == entity_data['Q6GZX4']['source']['seq'][9:20]
        assert self.c.meta['uniprot_sprot'] == {
            'proteins_with_epitopes': 1,
            'epitopes': 1,
            'proteins_in_IEDB': 0,
        }

    def test_integrate_epitope(self):
        class FakeIntegrate:
            def __init__(self, rows):
                self.rows = rows
                self.saved = []
                self.added = []
            def scan(self):
                return self.rows
            def save_data(self, json_data, data_source):
                self.saved.append(json_data)
            def add_data(self, data, acc, data_source):
                self.added.append((acc, data))

        self.c.integrate = FakeIntegrate([{'key': 'P_OLD'}])
        entity_data = {
            'P_OLD': {
                'accession': 'P_OLD',
                'source': {'accession': 'P_OLD'},
                'epitopes': [{'note': 'old'}],
            },
            'P_NEW': {
                'accession': 'P_NEW',
                'source': {'accession': 'P_NEW'},
                'epitopes': [{'note': 'new1'}, {'note': 'new2'}],
            },
        }
        data_source = ('UniProt', 'uniprot_sprot', 'epitope')
        count = self.c.integrate_epitope(entity_data, data_source)
        assert count == {
            'epitopes': 3,
            'epitope_proteins': 2,
            'updated_proteins': 1,
            'updated_epitopes': 1,
            'new_epitopes': 2,
            'new_proteins': 1,
        }
        # existing protein: epitopes are merged into the stored json
        assert len(self.c.integrate.saved) == 1
        saved = self.c.integrate.saved[0]
        assert saved['key'] == 'P_OLD'
        assert saved['epitopes']['UniProt'] == [{'note': 'old'}]
        assert saved['UniProt'] == {'accession': 'P_OLD'}
        # new protein: added with source and epitopes
        assert self.c.integrate.added == [(
            'P_NEW',
            {
                'UniProt': {'accession': 'P_NEW'},
                'epitopes': {'UniProt': [{'note': 'new1'}, {'note': 'new2'}]},
            },
        )]

