'''
Test src/bioomics/common/utils.py
'''
from ..helper import *
import subprocess

from src.bioomics import Utils


class TestUtils(TestCase):

    def test_count_iterator(self):
        self.assertEqual(Utils.count_iterator(x for x in [1, 2, 3]), 3)
        self.assertEqual(Utils.count_iterator([]), 0)

    def test_update_dict_new_key(self):
        self.assertEqual(Utils.update_dict({}, 'k', 1), {'k': 1})

    def test_update_dict_same_value(self):
        # same value: input unchanged
        self.assertEqual(Utils.update_dict({'k': 1}, 'k', 1), {'k': 1})

    def test_update_dict_scalar_to_list(self):
        # different scalar: old value and new value are merged into a list
        self.assertEqual(Utils.update_dict({'k': 1}, 'k', 2), {'k': [1, 2]})

    def test_update_dict_list_merge(self):
        self.assertEqual(Utils.update_dict({'k': [1, 2]}, 'k', 3), {'k': [1, 2, 3]})

    def test_update_dict_list_dedup(self):
        # value already in list: unchanged
        self.assertEqual(Utils.update_dict({'k': [1, 2]}, 'k', 2), {'k': [1, 2]})

    def test_update_nested_dict(self):
        # regression: every key used to be written flat into input and
        # the function returned None
        inp = {'a': {}}
        res = Utils.update_nested_dict(inp, ['a', 'b'], 1)
        self.assertEqual(inp, {'a': {'b': 1}})
        self.assertIs(res, inp)

    def test_update_nested_dict_creates_intermediates(self):
        inp = {}
        Utils.update_nested_dict(inp, ['x', 'y', 'z'], 'v')
        self.assertEqual(inp, {'x': {'y': {'z': 'v'}}})

    def test_update_nested_dict_replaces_non_dict(self):
        # non-dict on the path is replaced by a dict
        inp = {'x': 'scalar'}
        Utils.update_nested_dict(inp, ['x', 'y'], 1)
        self.assertEqual(inp, {'x': {'y': 1}})

    def test_update_nested_dict_leaf_list_merge(self):
        inp = {}
        Utils.update_nested_dict(inp, ['x', 'y'], 1)
        Utils.update_nested_dict(inp, ['x', 'y'], 2)
        self.assertEqual(inp, {'x': {'y': [1, 2]}})

    def test_update_nested_dict_empty_keys(self):
        inp = {'a': 1}
        self.assertIs(Utils.update_nested_dict(inp, [], 9), inp)
        self.assertEqual(inp, {'a': 1})

    def test_copy_data_no_hostname(self):
        # no DATA_HOSTNAME: skip synchronization
        with mock.patch.dict(os.environ, {'DATA_HOSTNAME': ''}):
            self.assertFalse(Utils.copy_data('/tmp/nowhere'))

    @mock.patch('src.bioomics.common.utils.socket.gethostname',
                return_value='local-host')
    @mock.patch('src.bioomics.common.utils.subprocess.run')
    @mock.patch.dict(os.environ, {
        'DATA_HOSTNAME': 'remote-host',
        'DATA_HOST_USER': 'datauser',
        'DATA_HOST_IP': '10.0.0.1',
    })
    def test_copy_data_rsync(self, m_run, m_hostname):
        m_run.return_value = mock.Mock(stdout='ok', stderr='', returncode=0)
        self.assertTrue(Utils.copy_data('/data/out'))
        cmd = m_run.call_args.args[0]
        self.assertIn('rsync', cmd)
        self.assertIn('datauser@remote-host:/data/out/', cmd)

    @mock.patch('src.bioomics.common.utils.socket.gethostname',
                return_value='remote-host')
    @mock.patch.dict(os.environ, {'DATA_HOSTNAME': 'remote-host'})
    def test_copy_data_identical_hostname(self, m_hostname):
        self.assertFalse(Utils.copy_data('/data/out'))

    @mock.patch('src.bioomics.common.utils.socket.gethostname',
                return_value='local-host')
    @mock.patch('src.bioomics.common.utils.subprocess.run',
                side_effect=subprocess.CalledProcessError(1, 'rsync'))
    @mock.patch.dict(os.environ, {
        'DATA_HOSTNAME': 'remote-host',
        'DATA_HOST_USER': 'datauser',
    })
    def test_copy_data_failure(self, m_run, m_hostname):
        self.assertFalse(Utils.copy_data('/data/out'))
