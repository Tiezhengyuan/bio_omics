'''
Test class 
'''
from .helper import *
from src.bioomics import ConnRedis

@ddt
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
    
    def test_put_data(self):
        key, val = '1', {'id':'2'}
        self.c.set_data(key, val)
        res = self.c.get_data(key)
        assert res == val

        #same key but different values
        val = {'id':'3'}
        self.c.put_data(key, val)
        res = self.c.get_data(key)
        assert res == {'id':['2','3']}

        # duplicaete value
        val = {'id':'3'}
        self.c.put_data(key, val)
        res = self.c.get_data(key)
        assert res == {'id':['2','3']}

        # new key
        key, val = '2', {'id':'3'}
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
        c = ConnectRedis('geneid')
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
        c = ConnectRedis('uniprot_acc')
        c.set_uniprotkb_acc(acc, data)
        res = c.get_uniprotkb_acc(acc)
        assert res == data
       
    @data(
        ['NC_000001.11', 'abc'],
        ['NC_000001.11', {'gi': 568815597}],
    )
    @unpack
    def test_ncbi_acc(self, acc, data):
        c = ConnectRedis('ncbi_acc')
        c.set_ncbi_acc(acc, data)
        res = c.get_ncbi_acc(acc)
        assert res == data