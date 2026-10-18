'''
Test src/bioomics/common/process_pickle.py
'''
from ..helper import *
import tempfile

from src.bioomics import ProcessPickle


class TestProcessPickle(TestCase):

    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.tmp, ignore_errors=True)
        self.path = os.path.join(self.tmp, 'data.pkl')

    def test_save_load_roundtrip(self):
        pp = ProcessPickle(self.path)
        self.assertEqual(pp.save_pickle({'a': 1}), 'succeed')
        self.assertEqual(pp.load_pickle(), {'a': 1})

    def test_save_accepts_str_path(self):
        # NB: str input is coerced to Path (is_file/pickle_file)
        pp = ProcessPickle(str(self.path))
        self.assertEqual(pp.save_pickle([1, 2]), 'succeed')
        self.assertEqual(pp.load_pickle(), [1, 2])

    def test_save_data_merges(self):
        pp = ProcessPickle(self.path)
        pp.save_data({'a': 1})
        pp.save_data({'b': 2})
        self.assertEqual(pp.load_pickle(), {'a': 1, 'b': 2})

    def test_load_missing_default(self):
        pp = ProcessPickle(self.path)
        self.assertIsNone(pp.load_pickle())
        self.assertEqual(pp.load_pickle(default={}), {})

    def test_load_corrupted(self):
        pp = ProcessPickle(self.path)
        with open(self.path, 'wb') as f:
            f.write(b'not a pickle')
        self.assertIsNone(pp.load_pickle())
        self.assertEqual(pp.load_pickle(default='fallback'), 'fallback')

    def test_is_good(self):
        pp = ProcessPickle(self.path)
        self.assertFalse(pp.is_good())
        pp.save_pickle({'a': 1})
        self.assertTrue(pp.is_good())
