'''
test ProcessBlast with the real ncbi-blast+ tools
example: pytest -s tests/protein/test_process_blast.py
'''
import shutil
from unittest import skipIf

from ..helper import *

from src.bioomics import ProcessBlast

BLAST_DB_BIN = shutil.which('makeblastdb')
BLAST_BIN = shutil.which('blastp')
FAA_FILE = os.path.join(DATA_DIR, 'blast', 'vregion_db.faa')
QUERY_FILE = os.path.join(DATA_DIR, 'blast', 'vregion_query.faa')

@ddt
@skipIf(not (BLAST_DB_BIN and BLAST_BIN), 'ncbi-blast+ is not installed')
class TestProcessBlast(TestCase):

    def setUp(self):
        self.outdir = os.path.join(RESULTS_DIR, 'blast')
        shutil.rmtree(self.outdir, ignore_errors=True)
        os.makedirs(self.outdir, exist_ok=True)
        self.p = ProcessBlast(self.outdir, True)

    def test_build_blastp_db(self):
        db_path = self.p.build_blastp_db(FAA_FILE)
        assert str(db_path) == os.path.join(self.outdir, 'vregion_db')
        assert os.path.isfile(f'{db_path}.pin')

    def test_build_blastp_db_fail(self):
        assert self.p.build_blastp_db('/no_exist.fa') is None

    def test_blastp(self):
        db_path = self.p.build_blastp_db(FAA_FILE)
        outfile = os.path.join(self.outdir, 'query.tsv')
        # default outfmt=7
        assert self.p.blastp(QUERY_FILE, db_path, outfile) == 'succeed'
        assert os.path.isfile(outfile)
        with open(outfile) as f:
            assert f.readline().startswith('# BLASTP')

    def test_blastp_fail(self):
        outfile = os.path.join(self.outdir, 'fail.tsv')
        assert self.p.blastp('/no_exist.fa', 'no_db', outfile) == 'fail'
        assert not os.path.isfile(outfile)

    def test_outdir_fallback(self):
        # outdir=None: database is built next to the fasta file
        workdir = os.path.join(RESULTS_DIR, 'blast_fallback')
        shutil.rmtree(workdir, ignore_errors=True)
        os.makedirs(workdir, exist_ok=True)
        faa = shutil.copy(FAA_FILE, workdir)
        p = ProcessBlast(None)
        db_path = p.build_blastp_db(faa)
        assert str(db_path) == os.path.join(workdir, 'vregion_db')
        assert os.path.isfile(f'{db_path}.pin')

    def test_aln_iter(self):
        db_path = self.p.build_blastp_db(FAA_FILE)
        self.p.blastp(QUERY_FILE, db_path,
            os.path.join(self.outdir, 'query.tsv'))
        files = list(self.p.aln_iter('.tsv'))
        assert len(files) == 1
        assert files[0].endswith('query.tsv')
