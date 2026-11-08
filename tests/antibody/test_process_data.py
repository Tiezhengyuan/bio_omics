'''
Unit tests for src/bioomics/antibody/process_data.py (ProcessData)
All filesystem work happens in temporary directories.
'''
import json
import os
import shutil
import tempfile

from ..helper import *

from src.bioomics.antibody.process_data import ProcessData


class TestProcessData(TestCase):

    def setUp(self):
        self.tmp_dir = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.tmp_dir, ignore_errors=True)
        self.pd = ProcessData(self.tmp_dir)

    def test_init(self):
        # data starts empty; data_dir is stored as-is
        self.assertIsNone(self.pd.data)
        self.assertEqual(self.pd.data_dir, self.tmp_dir)
        self.assertFalse(self.pd.verbose)
        self.assertTrue(ProcessData(self.tmp_dir, True).verbose)

    def test_init_dir_creates_nested(self):
        new_dir = os.path.join(self.tmp_dir, 'a', 'b')
        self.assertIs(ProcessData.init_dir(new_dir), True)
        self.assertTrue(os.path.isdir(new_dir))

    def test_init_dir_existing(self):
        self.assertIs(ProcessData.init_dir(self.tmp_dir), True)

    def test_init_dir_nested_file_blocks_creation(self):
        blocked = os.path.join(self.tmp_dir, 'f.json', 'sub')
        with open(os.path.join(self.tmp_dir, 'f.json'), 'w') as f:
            f.write('{}')
        self.assertIs(ProcessData.init_dir(blocked), False)

    def test_recursive_files_skips_hidden(self):
        visible = os.path.join(self.tmp_dir, 'sub', 'a.txt')
        os.makedirs(os.path.dirname(visible))
        with open(visible, 'w') as f:
            f.write('x')
        with open(os.path.join(self.tmp_dir, '.hidden'), 'w') as f:
            f.write('x')
        with open(os.path.join(self.tmp_dir, 'a.txt'), 'w') as f:
            f.write('x')
        res = sorted(ProcessData.recursive_files(self.tmp_dir))
        self.assertEqual(res, [
            os.path.join(self.tmp_dir, 'a.txt'),
            visible,
        ])

    def test_scan_text(self):
        file_a = os.path.join(self.tmp_dir, 'a.fa')
        file_b = os.path.join(self.tmp_dir, 'b.fa')
        for path, text in [(file_a, 'AB'), (file_b, 'CD')]:
            with open(path, 'w') as f:
                f.write(text)
        with open(os.path.join(self.tmp_dir, '.hidden'), 'w') as f:
            f.write('EF')
        res = sorted(ProcessData.scan_text(self.tmp_dir))
        self.assertEqual(res, [
            (os.path.basename(file_a), 'AB'),
            (os.path.basename(file_b), 'CD'),
        ])

    def test_json_roundtrip_and_sorted_keys(self):
        data = {'b': 1, 'a': [1, 2]}
        outfile = os.path.join(self.tmp_dir, 'data.json')
        self.assertIs(ProcessData.save_json(data, outfile), True)
        with open(outfile) as f:
            text = f.read()
        # save_json stores keys in sorted order
        self.assertLess(text.index('"a"'), text.index('"b"'))
        self.assertEqual(ProcessData.load_json(outfile), data)

    def test_load_json_missing_returns_empty_dict(self):
        missing = os.path.join(self.tmp_dir, 'missing.json')
        self.assertEqual(ProcessData.load_json(missing), {})

    def test_save_json_invalid_path_fails_gracefully(self):
        res = ProcessData.save_json(
            {'a': 1}, os.path.join(self.tmp_dir, 'no_dir', 'x.json'))
        self.assertIs(res, False)

    def test_load_data(self):
        existing = os.path.join(self.tmp_dir, 'exists.json')
        ProcessData.save_json({'k': 'v'}, existing)
        # existing file is loaded and merged into self.data
        self.assertEqual(self.pd.load_data(existing, {'d': 0}), {'k': 'v'})
        self.assertEqual(self.pd.data, {'k': 'v'})
        # missing file falls back to the default
        missing = os.path.join(self.tmp_dir, 'missing.json')
        self.assertEqual(self.pd.load_data(missing, {'d': 0}), {'d': 0})

    def test_retrieve_data(self):
        # retrieve_data updates self.data in place; it must start as a dict
        self.pd.data = {}
        text_iter = [('A1', 'MKVL'), ('A2', 'no protein here'), ('A3', 'ACDE')]
        self.pd.retrieve_data(
            text_iter, lambda text: {'seq': text} if text.startswith(('M', 'A')) else None)
        # every accession gets an entry; failed func results leave it empty
        self.assertEqual(self.pd.data, {
            'A1': {'seq': 'MKVL'},
            'A2': {},
            'A3': {'seq': 'ACDE'},
        })

    def test_to_label_csv(self):
        df = pd.DataFrame({'a': [1, 2], 'b': [3, 4]})
        outfile = os.path.join(self.tmp_dir, 'labels.csv')
        # shuffled with a fixed seed: both rows survive, no header/index
        ProcessData.to_label_csv(df, outfile)
        out = pd.read_csv(outfile, header=None)
        self.assertEqual(out.shape, (2, 2))
        self.assertEqual(
            sorted(map(tuple, out.values.tolist())), [(1, 3), (2, 4)])

    def test_scan_text(self):
        file_a = os.path.join(self.tmp_dir, 'a.fa')
        file_b = os.path.join(self.tmp_dir, 'b.fa')
        for path, text in [(file_a, 'AB'), (file_b, 'CD')]:
            with open(path, 'w') as f:
                f.write(text)
        with open(os.path.join(self.tmp_dir, '.hidden'), 'w') as f:
            f.write('EF')
        # defaults to self.data_dir; yields (basename, text)
        res = sorted(self.pd.scan_text())
        self.assertEqual(res, [('a.fa', 'AB'), ('b.fa', 'CD')])

    @skip('requires live UniProt/NCBI network access')
    def test_pull_list(self):
        pass


if __name__ == '__main__':
    import unittest
    unittest.main()
