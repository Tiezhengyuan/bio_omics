'''
Test src/bioomics/common/integrate_data.py

layout of the json storage:
    <entity_path>/<sub dirs>/<key>.json  + index file (json)
'''
from ..helper import *
import tempfile

from src.bioomics import IntegrateData


class TestIntegrateData(TestCase):

    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.tmp, ignore_errors=True)
        self.entity = os.path.join(self.tmp, 'entity')
        self.index_file = os.path.join(self.tmp, 'index.json')

    def test_init_from_path(self):
        # NB: callers pass a plain path (genBank, refseq, iedb_antigen...) -
        # the dict-only version crashed with TypeError
        integ = IntegrateData(self.entity)
        self.assertEqual(integ.entity_path, self.entity)
        self.assertEqual(integ.index_meta_file,
                         os.path.join(self.entity, 'index.json'))
        self.assertEqual(integ.index_meta, {})

    def test_init_from_meta_dict(self):
        integ = IntegrateData({
            'entity_path': self.entity,
            'index_file': self.index_file,
        })
        self.assertEqual(integ.entity_path, self.entity)
        self.assertEqual(integ.index_meta_file, self.index_file)

    def test_init_from_meta_dict_default_index(self):
        integ = IntegrateData({'entity_path': self.entity})
        self.assertEqual(integ.index_meta_file,
                         os.path.join(self.entity, 'index.json'))

    def test_get_meta(self):
        integ = IntegrateData(self.entity)
        meta = integ.get_meta({'source': 'UniProt'})
        self.assertEqual(meta, {'source': 'UniProt'})
        self.assertIs(meta, integ.index_meta)

    # ----- add_data / next_id -----
    def test_add_data_auto_numeric_key(self):
        integ = IntegrateData(self.entity)
        json_file = integ.add_data({'v': 1}, source='test')
        self.assertTrue(os.path.isfile(json_file))
        self.assertEqual(integ.index_meta['1']['key'], '1')
        self.assertEqual(integ.index_meta['1']['source'], ['test'])
        with open(json_file) as f:
            self.assertEqual(json.load(f), {'key': '1', 'v': 1})

    def test_add_data_accession_key(self):
        integ = IntegrateData(self.entity)
        json_file = integ.add_data({'v': 1}, key_value='Q6GZX4', source='UniProt')
        # key minus its last 2 chars ('Q6GZ') sliced per 3 chars: Q6G/Z
        self.assertEqual(json_file,
                         os.path.join(self.entity, 'Q6G', 'Z', 'Q6GZX4.json'))
        self.assertTrue(os.path.isfile(json_file))

    def test_add_data_default_source(self):
        integ = IntegrateData(self.entity)
        integ.add_data({'v': 1})
        self.assertEqual(integ.index_meta['1']['source'], ['UNKNOWN'])

    def test_next_id_skips_non_numeric_keys(self):
        # regression: int('Q6GZX4') raised ValueError
        integ = IntegrateData(self.entity)
        integ.add_data({'v': 1}, key_value='Q6GZX4')
        self.assertEqual(integ.next_id(), '1')

    def test_next_id_increments(self):
        integ = IntegrateData(self.entity)
        integ.add_data({'v': 1})
        integ.add_data({'v': 2})
        self.assertEqual(integ.next_id(), '3')

    # ----- save_data / get_data / scan -----
    def test_save_data_updates(self):
        integ = IntegrateData(self.entity)
        json_file = integ.add_data({'v': 1}, key_value='Q6GZX4', source='s1')
        data = integ.get_data('Q6GZX4')
        data['v'] = 2
        out = integ.save_data(data, source='s2')
        self.assertEqual(out, json_file)
        self.assertEqual(integ.get_data('Q6GZX4')['v'], 2)
        # sources are accumulated
        self.assertEqual(integ.index_meta['Q6GZX4']['source'], ['s1', 's2'])

    def test_save_data_recreates_missing_file(self):
        integ = IntegrateData(self.entity)
        json_file = integ.add_data({'v': 1}, key_value='Q6GZX4')
        os.remove(json_file)
        out = integ.save_data({'key': 'Q6GZX4', 'v': 9})
        # fallback re-adds under the same key
        self.assertEqual(out, json_file)
        self.assertEqual(integ.get_data('Q6GZX4')['v'], 9)

    def test_save_data_new_key(self):
        integ = IntegrateData(self.entity)
        out = integ.save_data({'key': 'P12345', 'v': 1})
        self.assertTrue(os.path.isfile(out))
        self.assertEqual(integ.get_data('P12345')['v'], 1)

    def test_get_data_missing(self):
        integ = IntegrateData(self.entity)
        self.assertEqual(integ.get_data('NOPE'), {})

    def test_scan(self):
        integ = IntegrateData(self.entity)
        integ.add_data({'v': 1})
        integ.add_data({'v': 2})
        data = sorted(d['v'] for d in integ.scan())
        self.assertEqual(data, [1, 2])

    # ----- index persistence -----
    def test_save_and_reload_index_meta(self):
        integ = IntegrateData({
            'entity_path': self.entity,
            'index_file': self.index_file,
        })
        integ.add_data({'v': 1}, key_value='Q6GZX4')
        self.assertTrue(integ.save_index_meta())
        # a fresh instance reads the index back
        integ2 = IntegrateData({
            'entity_path': self.entity,
            'index_file': self.index_file,
        })
        self.assertIn('Q6GZX4', integ2.index_meta)
        self.assertEqual(integ2.get_data('Q6GZX4')['v'], 1)

    def test_save_index_meta_empty(self):
        integ = IntegrateData(self.entity)
        self.assertFalse(integ.save_index_meta())

    # ----- path helpers -----
    def test_new_json_path(self):
        integ = IntegrateData(self.entity)
        path = integ.new_json_path('1234567')
        # floor(1234567/1000)='1234' -> 12/34 sub dirs
        self.assertEqual(path,
                         os.path.join(self.entity, '12', '34', '1234567.json'))
        self.assertTrue(os.path.isdir(os.path.join(self.entity, '12', '34')))
