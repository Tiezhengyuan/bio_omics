'''
Test src/bioomics/common/process_json.py
'''
from ..helper import *
import tempfile

from src.bioomics import ProcessJson


class TestProcessJson(TestCase):

    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.tmp, ignore_errors=True)

    # ----- load_json / save_json -----
    def test_save_load_roundtrip(self):
        path = os.path.join(self.tmp, 'data.json')
        self.assertTrue(ProcessJson.save_json({'a': 1, 'b': [1, 2]}, path))
        self.assertEqual(ProcessJson.load_json(path), {'a': 1, 'b': [1, 2]})

    def test_load_json_missing_default(self):
        path = os.path.join(self.tmp, 'missing.json')
        self.assertEqual(ProcessJson.load_json(path), {})
        # default=None means "no default given": {} is returned
        self.assertEqual(ProcessJson.load_json(path, default=None), {})
        self.assertEqual(ProcessJson.load_json(path, default={'x': 1}), {'x': 1})

    def test_save_json_non_serializable(self):
        # default=str keeps save_json from raising
        path = os.path.join(self.tmp, 'odd.json')
        self.assertTrue(ProcessJson.save_json({'obj': object()}, path))

    # ----- get_value -----
    def test_get_value_nested(self):
        # regression: Path.isfile() (AttributeError) - correct is is_file()
        path = os.path.join(self.tmp, 'data.json')
        ProcessJson.save_json({'a': {'b': {'c': 42}}}, path)
        self.assertEqual(ProcessJson.get_value(path, ['a', 'b', 'c']), 42)
        self.assertEqual(ProcessJson.get_value(path, ['a', 'b']), {'c': 42})

    def test_get_value_missing(self):
        path = os.path.join(self.tmp, 'data.json')
        ProcessJson.save_json({'a': {'b': 1}}, path)
        self.assertIsNone(ProcessJson.get_value(path, ['a', 'x']))
        self.assertIsNone(ProcessJson.get_value(path, ['a', 'x', 'y']))
        self.assertIsNone(ProcessJson.get_value(path, []))

    def test_get_value_accepts_pathlib(self):
        from pathlib import Path
        path = Path(os.path.join(self.tmp, 'data.json'))
        ProcessJson.save_json({'a': {'b': 1}}, path)
        self.assertEqual(ProcessJson.get_value(path, ['a', 'b']), 1)

    def test_get_value_missing_file(self):
        self.assertIsNone(
            ProcessJson.get_value(os.path.join(self.tmp, 'nope.json'), ['a']))

    # ----- iter_to_json -----
    def test_iter_to_json_chunks(self):
        # NB: accepts a str outdir and creates it (old code required Path)
        outdir = os.path.join(self.tmp, 'chunks', 'nested')
        file_list = ProcessJson.iter_to_json(iter(range(5)), outdir, chunk_size=2)
        self.assertEqual(len(file_list), 3)
        with open(file_list[0]) as f:
            self.assertEqual(json.load(f), [0, 1])
        with open(file_list[-1]) as f:
            self.assertEqual(json.load(f), [4])

    def test_iter_to_json_empty(self):
        outdir = os.path.join(self.tmp, 'empty_out')
        file_list = ProcessJson.iter_to_json(iter([]), outdir, chunk_size=2)
        self.assertEqual(file_list, [])
