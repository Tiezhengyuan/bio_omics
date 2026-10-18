'''
igblastp executables are not installed: only the failure paths are tested
example: pytest -s tests/protein/test_process_igblast.py
'''
from ..helper import *

from src.bioomics import ProcessIgblast

@ddt
class TestProcessIgblast(TestCase):

    def setUp(self):
        self.p = ProcessIgblast(RESULTS_DIR, True)

    def test_build_pro_db_fail(self):
        # input fasta does not exist: returns None, no exception
        assert self.p.build_pro_db('/no_exist.fa') is None

    def test_igblastp_fail(self):
        outfile = os.path.join(RESULTS_DIR, 'igblastp_fail.out')
        assert self.p.igblastp('/no_exist.fa', 'no_db', outfile) == 'fail'
        assert not os.path.isfile(outfile)