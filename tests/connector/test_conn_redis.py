'''
Test class
'''
import socket

from tests.helper import *
from unittest import skipIf
# ConnRedis is not re-exported in src.bioomics/__init__: import from its module
from src.bioomics.connector.conn_redis import ConnRedis

# redis is an optional local service: skip redis-dependent tests when offline
def _redis_available() -> bool:
    try:
        with socket.create_connection(('localhost', 6379), timeout=2):
            return True
    except OSError:
        return False

REDIS_AVAILABLE = _redis_available()

@ddt
@skipIf(not REDIS_AVAILABLE, 'redis is not running on localhost:6379')
class TestConnRedis(TestCase):

    def setUp(self):
        self.c = ConnRedis()
        self.c.conn.set('a', 1)
        self.c.conn.set('b', json.dumps({'c':[1,2], 'd':{'e':'f'}}))

    def test_connection(self):
        self.c.is_connected()

    @data(
        ['a', 1],
        ['wrong', None],
        ['b', {'c':[1,2], 'd':{'e':'f'}}],
    )
    @unpack
    def test_get_data(self, key, expect):
        res = self.c.get_data(key)
        assert res == expect

    @data(
        ['a', 1],
        ['b', {'c':[1,2], 'd':{'e':'f'}}],
        ['c', None],
        ['d', [1,2,'a']],
    )
    @unpack
    def test_set_data(self, key, val):
        self.c.set_data(key, val)
        res = self.c.get_data(key)
        assert res == val
    
    @data(
        ['856646', []],
        ['856646', None],
        ['856646', {'symbol': 'CAN1', 'pmid':[34,66]}],
        ['GeneID:856646', {'symbol': 'CAN1'}],
        ['GeneID: 856646', {'symbol': 'CAN1'}],
        ["geneid: 856646 \t ", {'symbol': 'CAN1'}],
    )
    @unpack
    def test_geneid(self, geneid, data):
        c = ConnRedis('geneid')
        c.set_geneid(geneid, data)
        res = c.get_geneid(geneid)
        assert res == data

    @data(
        ['A2BC19', 'abc'],
        ['A2BC19', []],
        ['A2BC19', {'id': 'A2BC19_HELPX'}],
        ['A2BC19', None],
        ['wrong', None],
    )
    @unpack
    def test_uniprotkb_acc(self, acc, data):
        c = ConnRedis('uniprot_acc')
        c.set_uniprotkb_acc(acc, data)
        res = c.get_uniprotkb_acc(acc)
        assert res == data
       
    @data(
        ['NC_000001.11', 'abc'],
        ['NC_000001.11', {'gi': 568815597}],
    )
    @unpack
    def test_ncbi_acc(self, acc, data):
        c = ConnRedis('ncbi_acc')
        c.set_ncbi_acc(acc, data)
        res = c.get_ncbi_acc(acc)
        assert res == data