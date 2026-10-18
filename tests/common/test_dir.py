'''
Test src/bioomics/common/dir.py
'''
from ..helper import *
import tempfile

from src.bioomics import Dir


class TestDir(TestCase):

    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.tmp, ignore_errors=True)

    # ----- init_dir -----
    def test_init_dir_creates_nested(self):
        target = os.path.join(self.tmp, 'a', 'b', 'c')
        self.assertTrue(Dir(target).init_dir())
        self.assertTrue(os.path.isdir(target))

    def test_init_dir_existing(self):
        self.assertTrue(Dir(self.tmp).init_dir())

    def test_init_dir_failure(self):
        # not creatable: returns False instead of raising
        self.assertFalse(Dir('/proc/self/nonexistent/xyz').init_dir())

    # ----- recursive files -----
    def _make_tree(self):
        os.makedirs(os.path.join(self.tmp, 'a'))
        os.makedirs(os.path.join(self.tmp, 'b', 'c'))
        open(os.path.join(self.tmp, 'a', 'f1.txt'), 'w').write('1')
        open(os.path.join(self.tmp, 'b', 'c', 'f2.txt'), 'w').write('2')
        open(os.path.join(self.tmp, 'b', 'f3.log'), 'w').write('3')
        open(os.path.join(self.tmp, '.dotfile'), 'w').write('.')

    def test_recursive_files(self):
        self._make_tree()
        files = sorted(Dir(self.tmp).recursive_files())
        # dotfiles are skipped
        self.assertEqual(files, [
            os.path.join(self.tmp, 'a', 'f1.txt'),
            os.path.join(self.tmp, 'b', 'c', 'f2.txt'),
            os.path.join(self.tmp, 'b', 'f3.log'),
        ])

    def test_filter_recursive_files(self):
        self._make_tree()
        files = list(Dir(self.tmp).filter_recursive_files('.txt'))
        self.assertEqual(len(files), 2)
        self.assertTrue(all(f.endswith('.txt') for f in files))

    def test_recursive_dirs(self):
        self._make_tree()
        # leaf directories only (no sub directories)
        dirs = sorted(Dir(self.tmp).recursive_dirs())
        self.assertEqual(dirs, [
            os.path.join(self.tmp, 'a'),
            os.path.join(self.tmp, 'b', 'c'),
        ])

    def test_clear_dir(self):
        self._make_tree()
        d = Dir(self.tmp)
        d.clear_dir()
        self.assertEqual(list(d.recursive_files()), [])
        # directories are kept
        self.assertTrue(os.path.isdir(os.path.join(self.tmp, 'b', 'c')))

    # ----- compress_files -----
    def test_compress_files(self):
        # regression: info[path].append() raised KeyError on the first match
        txt = os.path.join(self.tmp, 'x.txt')
        log = os.path.join(self.tmp, 'y.log')
        open(txt, 'w').write('compress me')
        open(log, 'w').write('keep me')
        info = Dir(self.tmp).compress_files('.txt')
        self.assertEqual(list(info), [txt])
        self.assertTrue(os.path.isfile(txt + '.gz'))
        self.assertFalse(os.path.isfile(txt))
        self.assertTrue(os.path.isfile(log))

    # ----- format_dir -----
    def test_format_dir_relative(self):
        cwd = os.getcwd()
        os.chdir(self.tmp)
        self.addCleanup(os.chdir, cwd)
        formatted = Dir('rel_sub').format_dir()
        self.assertEqual(formatted, os.path.join(self.tmp, 'rel_sub') + '/')
        self.assertTrue(os.path.isdir(os.path.join(self.tmp, 'rel_sub')))

    def test_format_dir_none(self):
        # regression: indir was dereferenced (.find) before the None check
        self.assertEqual(Dir(None).format_dir(), os.getcwd() + '/')

    def test_format_dir_home(self):
        home = os.getenv('HOME')
        target = os.path.join(home, '.bioomics_test_fmt')
        self.addCleanup(shutil.rmtree, target, ignore_errors=True)
        self.assertEqual(Dir('~/.bioomics_test_fmt').format_dir(), target + '/')
        self.assertTrue(os.path.isdir(target))

    # ----- cascade_dir -----
    def test_cascade_dir(self):
        self.assertEqual(
            Dir.cascade_dir('/tmp/x', 'abcdefgh', 2),
            '/tmp/x/ab/cd/ef/gh')

    def test_cascade_dir_short_id(self):
        self.assertEqual(Dir.cascade_dir('/tmp/x', 'ab', 2), '/tmp/x/ab')
