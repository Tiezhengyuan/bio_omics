'''
example: pytest -s tests/protein/test_parse_blast.py
'''
from ..helper import *

from src.bioomics import ParseBlastXml

@ddt
class TestParseBlast(TestCase):

    def test_blastp_xml(self):
        aln_file = '/home/yuan/output/pdb/meta_vdj/4I/4IJ4/4IJ4_0.dregion.blastp5'
        aln_file = '/home/yuan/output/pdb/meta_vdj/4R/4RFL/4RFL_0.dregion.blastp5'
        p = ParseBlastXml(aln_file)
        for rec in p.best_hits():
            print(rec)

