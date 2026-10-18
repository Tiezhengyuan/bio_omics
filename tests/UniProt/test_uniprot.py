'''
Test class Uniprot: *.dat, *.fasta iterators
'''
import gzip
import tempfile
from pathlib import Path

from ..helper import *
from src.bioomics import Uniprot

DIR_UNIPROT = os.path.join(DATA_DIR, 'UniProt')

FASTA_TEXT = \
    '>UniRef50_A0A7C4TNZ0 Uncharacterized protein 002L n=2 Taxon=Frog virus 3\n' \
    'MKVLWA\n' \
    '>UniRef50_Q96RD9 Signal sequence\n' \
    'MVLSP\n'


@ddt
class TestUniprot(TestCase):

    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self._tmp.cleanup)
        self.tmp_dir = self._tmp.name

    @data(
        ['uniprot_sprot.dat', 'Q6GZX4'],
        ['uniprot_sprot.dat.gz', 'Q6GZX4'],
        ['uniprot_trembl.dat', 'A0A7C4TNZ0'],
        ['uniprot_trembl.dat.gz', 'A0A7C4TNZ0'],
    )
    @unpack
    def test_dat_iter(self, file_name, expect):
        dat_file = os.path.join(DIR_UNIPROT, file_name)
        res = Uniprot(dat_file).dat_iter()
        rec = next(res)
        assert rec.id == expect

    def test_dat_iter_default(self):
        # no infile: uniprot_dir may be a directory (default file
        # uniprot_sprot.dat.gz) or the input file itself
        res = Uniprot(DIR_UNIPROT).dat_iter()
        assert next(res).id == 'Q6GZX4'
        res = Uniprot(os.path.join(DIR_UNIPROT, 'uniprot_trembl.dat')).dat_iter()
        assert next(res).id == 'A0A7C4TNZ0'

    def test_filter_uniprot_sprot_dat(self):
        # records whose accession is NOT in ids are yielded
        res = Uniprot(DIR_UNIPROT).filter_uniprot_sprot_dat(['Q6GZX4'])
        assert [r.id for r in res] == ['Q6GZX3']

    def test_filter_uniprot_trembl_dat(self):
        res = Uniprot(DIR_UNIPROT).filter_uniprot_trembl_dat(['E9LU99'])
        assert [r.id for r in res] == ['A0A7C4TNZ0']

    @data(
        ['.fasta'],
        ['.fasta.gz'],
    )
    @unpack
    def test_parse_fasta(self, suffix):
        # NB: parse_fasta used to crash on infile.surfix (AttributeError)
        infile = Path(self.tmp_dir) / f'uniprot_sprot{suffix}'
        if suffix.endswith('.gz'):
            with gzip.open(infile, 'wt') as f:
                f.write(FASTA_TEXT)
        else:
            infile.write_text(FASTA_TEXT)
        res = list(Uniprot(infile).parse_fasta('UniRef50_'))
        assert [r['id'] for r in res] == ['A0A7C4TNZ0', 'Q96RD9']
        assert [r['sequence'] for r in res] == ['MKVLWA', 'MVLSP']
        assert res[0]['description'].startswith('UniRef50_A0A7C4TNZ0 Uncharacterized protein 002L')

    @data(
        ['.fasta'],
        ['.fasta.gz'],
    )
    @unpack
    def test_fasta_iter(self, suffix):
        infile = Path(self.tmp_dir) / f'uniprot_sprot{suffix}'
        if suffix.endswith('.gz'):
            with gzip.open(infile, 'wt') as f:
                f.write(FASTA_TEXT)
        else:
            infile.write_text(FASTA_TEXT)
        res = list(Uniprot(infile).fasta_iter())
        assert [r.id for r in res] == ['UniRef50_A0A7C4TNZ0', 'UniRef50_Q96RD9']

    def test_scan_uniref50(self):
        # default input: <uniprot_dir>/uniref50.fasta.gz
        with gzip.open(os.path.join(self.tmp_dir, 'uniref50.fasta.gz'), 'wt') as f:
            f.write(FASTA_TEXT)
        res = list(Uniprot(self.tmp_dir).scan_uniref50())
        assert [r['id'] for r in res] == ['A0A7C4TNZ0', 'Q96RD9']
        # uniprot_dir may also be the uniref file itself
        res = list(Uniprot(os.path.join(self.tmp_dir, 'uniref50.fasta.gz')).scan_uniref50())
        assert [r['id'] for r in res] == ['A0A7C4TNZ0', 'Q96RD9']

    def test_uniref50_iter(self):
        infile = os.path.join(self.tmp_dir, 'uniref50.fasta.gz')
        with gzip.open(infile, 'wt') as f:
            f.write(FASTA_TEXT)
        assert [r.id for r in Uniprot(self.tmp_dir).uniref50_iter()] == \
            ['UniRef50_A0A7C4TNZ0', 'UniRef50_Q96RD9']
        # plain (not gzipped) input is supported as well
        plain = os.path.join(self.tmp_dir, 'uniref50.fasta')
        with open(plain, 'w') as f:
            f.write(FASTA_TEXT)
        assert [r.id for r in Uniprot(self.tmp_dir).uniref50_iter(infile=plain)] == \
            ['UniRef50_A0A7C4TNZ0', 'UniRef50_Q96RD9']


