'''
Test src/bioomics/common/write_tmp.py

WriteTmp writes IO-intensive results into a temporary directory and
moves them to the destination when the context manager exits (rsync).
'''
from ..helper import *
import tempfile
from pathlib import Path
from unittest import skipIf

from src.bioomics import ProcessJson
from src.bioomics.common.write_tmp import WriteTmp

RSYNC = shutil.which('rsync') is not None


@ddt
class TestWriteTmp(TestCase):

    def setUp(self):
        self.base = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.base, ignore_errors=True)
        self.tmp_base = os.path.join(self.base, 'tmp')
        self.dest = os.path.join(self.base, 'dest')
        os.makedirs(self.tmp_base)

    def test_temp_path_before_enter(self):
        # NB: temp_dir_path is initialized to None in __init__ so that a
        # clean RuntimeError is raised instead of AttributeError
        writer = WriteTmp(self.tmp_base, self.dest, 'io')
        with self.assertRaises(RuntimeError):
            writer.temp_path

    def test_destination_dir_created(self):
        writer = WriteTmp(self.tmp_base, self.dest, 'io')
        self.assertEqual(writer.destination_dir, Path(self.dest))
        self.assertTrue(os.path.isdir(self.dest))

    @skipIf(not RSYNC, 'rsync is not installed')
    def test_files_moved_to_destination(self):
        with WriteTmp(self.tmp_base, self.dest, 'unittest') as writer:
            self.assertTrue(writer.temp_dir_path.is_dir())
            self.assertTrue(
                writer.temp_dir_path.name.startswith('unittest_'))
            ProcessJson.save_json({'a': 1}, writer.temp_path / 'a.json')
            (writer.temp_path / 'sub').mkdir()
            ProcessJson.save_json({'b': 2}, writer.temp_path / 'sub' / 'b.json')
        # after the with block: results moved, temporary directory removed
        self.assertTrue(os.path.isfile(os.path.join(self.dest, 'a.json')))
        self.assertTrue(
            os.path.isfile(os.path.join(self.dest, 'sub', 'b.json')))
        self.assertEqual(
            ProcessJson.load_json(os.path.join(self.dest, 'a.json')), {'a': 1})
        self.assertFalse(writer.temp_dir_path.exists())

    @skipIf(not RSYNC, 'rsync is not installed')
    def test_context_exit_without_files(self):
        writer = WriteTmp(self.tmp_base, self.dest, 'io')
        with writer:
            pass
        self.assertFalse(writer.temp_dir_path.exists())


