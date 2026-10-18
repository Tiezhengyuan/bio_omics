'''
Test src/bioomics/common/process_file.py
'''
from ..helper import *
import tempfile

from src.bioomics import ProcessFile


class TestProcessFile(TestCase):

    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.tmp, ignore_errors=True)
        # tree: <tmp>/3GKI/meta.json, <tmp>/1I8I/meta.json, decoys
        for pdb_id in ['3GKI', '1I8I']:
            pdb_dir = os.path.join(self.tmp, pdb_id)
            os.makedirs(pdb_dir)
            with open(os.path.join(pdb_dir, 'meta.json'), 'w') as f:
                json.dump({'pdb_id': pdb_id}, f)
        os.makedirs(os.path.join(self.tmp, 'plain'))
        with open(os.path.join(self.tmp, 'plain', 'other.json'), 'w') as f:
            json.dump({}, f)

    def test_get_meta(self):
        # regression: Dir.recrusive_files() (typo) raised AttributeError
        ready = ProcessFile.get_meta(self.tmp)
        self.assertEqual(sorted(ready), ['1I8I', '3GKI'])
        self.assertTrue(ready['3GKI'].endswith('meta.json'))

    def test_iter_files(self):
        pairs = sorted(ProcessFile.iter_files(self.tmp, 'meta.json'))
        self.assertEqual([pdb_id for pdb_id, _ in pairs], ['1I8I', '3GKI'])
        self.assertTrue(all(path.endswith('meta.json') for _, path in pairs))

    def test_iter_files_no_match(self):
        self.assertEqual(list(ProcessFile.iter_files(self.tmp, 'nope.txt')), [])

    def test_to_label_csv(self):
        df = pd.DataFrame({'x': [1, 2, 3], 'label': [0, 1, 0]})
        outfile = os.path.join(self.tmp, 'labels.csv')
        ProcessFile.to_label_csv(df, outfile)
        with open(outfile) as f:
            lines = [line.strip() for line in f if line.strip()]
        # header=False: 3 data lines only (shuffled deterministically)
        self.assertEqual(len(lines), 3)
        self.assertNotIn('x', lines[0])
        self.assertEqual(sorted(lines), ['1,0', '2,1', '3,0'])

    def test_iter_to_txt(self):
        paths = [os.path.join(self.tmp, str(i)) for i in range(5)]
        # contract: the iterator yields (id, path) pairs
        pairs = [(str(i), p) for i, p in enumerate(paths)]
        file_list = ProcessFile.iter_to_txt(iter(pairs), os.path.join(self.tmp, 'chunk'), chunk_size=2)
        self.assertEqual(len(file_list), 3)
        self.assertTrue(os.path.isfile(os.path.join(self.tmp, 'chunk_0.txt')))
        with open(file_list[0]) as f:
            self.assertEqual(f.read().split('\n'), paths[:2])
        with open(file_list[-1]) as f:
            self.assertEqual(f.read().split('\n'), paths[4:])
