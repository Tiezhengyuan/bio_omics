from ..helper import *
from src.bioomics import ProcessBlast

@ddt
class TestProcessBlast(TestCase):

    @skip
    def test_build_blastp_db(self):
        outdir = '/home/yuan/data/UniProt/db'
        p = ProcessBlast(outdir, True)
        fa_file = '/home/yuan/data/UniProt/uniprot_sprot.fasta'
        p.build_blastp_db(fa_file)
