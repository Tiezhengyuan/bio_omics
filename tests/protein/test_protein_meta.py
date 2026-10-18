'''
example: pytest -s tests/protein/test_protein_meta.py
'''
import re
import shutil
from datetime import datetime
from pathlib import Path

from ..helper import *

from src.bioomics import ProteinMeta

@ddt
class TestProteinMeta(TestCase):

    def setUp(self):
        self.result_dir = Path(RESULTS_DIR) / 'protein_meta'
        shutil.rmtree(self.result_dir, ignore_errors=True)
        self.m = ProteinMeta(self.result_dir, 'epitope', 'IEDB')

    def test_paths(self):
        assert self.m.entity_path() == self.result_dir / 'epitope' / 'IEDB'
        assert os.path.isdir(self.m.entity_path())
        assert self.m.meta_path() == self.result_dir / 'epitope' / 'meta'
        assert self.m.meta_file() \
            == self.result_dir / 'epitope' / 'meta' / 'IEDB_meta.json'
        assert self.m.index_file() \
            == self.result_dir / 'epitope' / 'meta' / 'IEDB_index.json'

    def test_get_meta(self):
        meta = self.m.get_meta({'custom': 'v1'})
        assert meta['entity'] == 'epitope'
        assert meta['custom'] == 'v1'
        assert isinstance(meta['start_time'], datetime)
        assert meta['meta_file'] == self.m.meta_file()

    def test_save_meta(self):
        meta = self.m.get_meta({'custom': 'v1'})
        meta_file = self.m.save_meta(meta)
        assert meta_file == self.m.meta_file()
        assert meta_file.is_file()
        data = json.loads(meta_file.read_text())
        assert data['custom'] == 'v1'
        assert isinstance(data['duration'], int)
        # year is formatted with %Y: 4 digits (regression: '%Y' typo)
        assert re.match(r'^\d{2}/\d{2}/\d{4}, \d{2}:\d{2}:\d{2}$',
            data['start_time'])
        assert re.match(r'^\d{2}/\d{2}/\d{4}, \d{2}:\d{2}:\d{2}$',
            data['end_time'])
        # Path objects are serialized
        assert os.path.isdir(data['entity_path'])
        assert data['meta_file'] == str(self.m.meta_file())

    def test_get_meta_loads_saved(self):
        self.m.save_meta(self.m.get_meta({'custom': 'v1'}))
        # saved meta overrides the default one
        meta = self.m.get_meta({'custom': 'v2'})
        assert meta['custom'] == 'v1'
        # start_time is refreshed for the next run
        assert isinstance(meta['start_time'], datetime)
        # meta_file loaded from json is a string
        assert self.m.save_meta(meta) == str(self.m.meta_file())

    def test_get_index_meta(self):
        assert self.m.get_index_meta() == {}
        index_file = self.m.index_file()
        index_file.write_text(json.dumps({'chain': 'H'}))
        assert self.m.get_index_meta() == {'chain': 'H'}