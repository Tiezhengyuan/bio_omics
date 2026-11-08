'''
Unit tests for src/bioomics/PDB/pull_pdb.py (PullPdb)

Only the offline scan/metadata helpers are tested here; the download
helpers (rsync/wget/SAbDab/PDBe API) are network-bound and skipped.
'''
import os
import shutil
import tempfile
import time
from datetime import date
from pathlib import Path

from ..helper import *

from src.bioomics import PullPdb


class TestPullPdbScan(TestCase):

    def setUp(self):
        self.tmp_dir = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.tmp_dir, ignore_errors=True)
        self.pdb_dir = Path(self.tmp_dir, 'pdb')
        self.pp = PullPdb(self.pdb_dir)

    def _make_files(self):
        # legacy .ent.gz tree and plain .pdb files
        os.makedirs(os.path.join(self.pdb_dir, 'ab'))
        os.makedirs(os.path.join(self.pdb_dir, 'cd'))
        self.file_abc = os.path.join(self.pdb_dir, 'ab', 'pdb1abc.ent.gz')
        self.file_def = os.path.join(self.pdb_dir, 'cd', '1def.pdb')
        self.file_other = os.path.join(self.pdb_dir, 'cd', 'readme.txt')
        for path in (self.file_abc, self.file_def, self.file_other):
            with open(path, 'w') as f:
                f.write('x')

    def test_init(self):
        self.assertEqual(self.pp.pdb_dir, self.pdb_dir)
        self.assertFalse(self.pp.verbose)

    def test_pdb_id(self):
        cases = [
            ('pdb1abc.ent.gz', '1ABC'),
            ('1abc.cif.gz', '1ABC'),
            ('1abc.pdb', '1ABC'),
            ('xyz.pdb', 'XYZ'),
            ('readme.txt', None),
            ('7s8l.cif', None),
        ]
        for path, expect in cases:
            self.assertEqual(self.pp.pdb_id(path), expect, path)

    def test_iter_raw_pdb_and_get_raw_pdb(self):
        self._make_files()
        res = self.pp.get_raw_pdb()
        self.assertEqual(res, {
            '1ABC': self.file_abc,
            '1DEF': self.file_def,
        })
        self.assertEqual(
            sorted(self.pp.iter_raw_pdb()),
            [('1ABC', self.file_abc), ('1DEF', self.file_def)],
        )

    def test_filter_out_raw_pdb(self):
        self._make_files()
        res = self.pp.filter_out_raw_pdb({'1ABC'})
        # unknown files only
        self.assertEqual(res, {'1DEF': self.file_def})

    def test_filter_in_raw_pdb(self):
        self._make_files()
        res = dict(self.pp.filter_in_raw_pdb({'1ABC'}))
        self.assertEqual(res, {'1ABC': self.file_abc})
        # unknown set yields nothing
        self.assertEqual(dict(self.pp.filter_in_raw_pdb({'9ZZZ'})), {})

    def test_get_new_raw_pdb(self):
        self._make_files()
        today = date.today()
        # freshly created files are reported
        self.assertEqual(
            set(self.pp.get_new_raw_pdb(today)), {'1ABC', '1DEF'})
        # backdate 1DEF: it disappears for today's stamp
        old = time.mktime((2000, 1, 1, 0, 0, 0, 0, 0, -1))
        os.utime(self.file_def, (old, old))
        self.assertEqual(
            set(self.pp.get_new_raw_pdb(today)), {'1ABC'})
        self.assertEqual(
            set(self.pp.get_new_raw_pdb(date(1999, 12, 31))),
            {'1ABC', '1DEF'})

    def test_iter_meta(self):
        os.makedirs(os.path.join(self.pdb_dir, '1abc'))
        meta_file = os.path.join(self.pdb_dir, '1abc', 'meta.json')
        with open(meta_file, 'w') as f:
            json.dump({'x': 1}, f)
        # only files named meta.json are picked up
        with open(os.path.join(self.pdb_dir, '1abc', 'other.json'), 'w') as f:
            json.dump({'y': 2}, f)
        res = list(self.pp.iter_meta())
        self.assertEqual(res, [('1abc', {'x': 1})])


class TestPullPdbNetwork(TestCase):

    def setUp(self):
        self.tmp_dir = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.tmp_dir, ignore_errors=True)
        self.pp = PullPdb(Path(self.tmp_dir), verbose=False)

    @skip('downloads the entire RCSB PDB tree via rsync')
    def test_synchronize_pdb(self):
        self.assertEqual(self.pp.synchronize_pdb(), 'success')

    @skip('network download from files.rcsb.org')
    def test_download_pdb(self):
        self.pp.download_pdb('7S8L')

    @skip('network download from rcsb.org')
    def test_download_fasta(self):
        fasta_file, endpoint = self.pp.download_fasta('7S8L')
        self.assertIsNotNone(fasta_file)

    @skip('network download from opig.stats.ox.ac.uk')
    def test_getpdb(self):
        self.pp.pdb_dir.mkdir(parents=True, exist_ok=True)
        self.assertTrue(self.pp.getpdb('7S8L'))

    @skip('network access to the PDBe API')
    def test_pdb_to_uniprot(self):
        res = self.pp.pdb_to_uniprot('7S8L')
        self.assertIn('UniProt', res)
