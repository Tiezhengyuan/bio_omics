'''
pytest -s tests/database/test_filter.py
'''
import os
import shutil
import tempfile

from ..helper import *
from src.bioomics import FilterDb, QueryComplex


@ddt
class TestFilter(TestCase):

    def setUp(self) -> None:
        self.tmp_dir = tempfile.mkdtemp()

    def tearDown(self) -> None:
        shutil.rmtree(self.tmp_dir)

    def test_uniprot_align_iter(self):
        # only *.uniprot.blastp5 files are collected
        pdb_dir = os.path.join(self.tmp_dir, 'ab', '1ABC')
        os.makedirs(pdb_dir)
        aln = os.path.join(pdb_dir, 'chain_sequences.uniprot.blastp5')
        with open(aln, 'w'):
            pass
        other = os.path.join(pdb_dir, 'chain_sequences.faa')
        with open(other, 'w'):
            pass

        fd = FilterDb(self.tmp_dir)
        res = list(fd.uniprot_align_iter())
        assert res == [aln]

    def test_uniprot_align_check_iter(self):
        # collect the alignment if its pdb_id is not in the table
        pdb_dir = os.path.join(self.tmp_dir, 'ab', '1ABC')
        os.makedirs(pdb_dir)
        aln = os.path.join(pdb_dir, 'chain_sequences.uniprot.blastp5')
        with open(aln, 'w'):
            pass

        fd = FilterDb(self.tmp_dir, False)
        q = QueryComplex(True)
        in_db = '1ABC' in q.unique_values('align_uniprot', 'pdb_id')

        res = list(fd.uniprot_align_check_iter('align_uniprot'))
        if in_db:
            assert aln not in res
        else:
            assert res == [aln]

    def test_igblastp_align_iter(self):
        # the alignment is skipped when the fasta file is missing
        pdb_dir = os.path.join(self.tmp_dir, 'ab', '1ABC')
        os.makedirs(pdb_dir)
        aln = os.path.join(pdb_dir, 'chain_sequences.vregion.igblastp3')
        with open(aln, 'w'):
            pass

        fd = FilterDb(self.tmp_dir)
        res = list(fd.igblastp_align_iter())
        assert res == []

        # with the fasta file both files are yielded
        faa = os.path.join(pdb_dir, 'chain_sequences.faa')
        with open(faa, 'w'):
            pass
        res = list(fd.igblastp_align_iter())
        assert res == [(aln, faa)]
