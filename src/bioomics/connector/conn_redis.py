"""
cache data
"""
import re
import redis
import json
from biosequtils import KeyValue

# default database is 0
redis_db = {
    # NCBI
    'geneid':1,
    'gi':2,
    'ncbi_acc':3,

    # UniProt
    'uniprot_acc': 6,

}
class ConnRedis(object):
    def __init__(self, db:str=None, host:str=None, port:int=None):
        if db is None: db = ''
        self.db = redis_db.get(db, 0)
        self.host = host if host else 'localhost'
        self.port = port if port else 6379
    
        # default redis database of ontology is 0
        self.conn = redis.StrictRedis(
            host=self.host,
            port=self.port,
            db=self.db,
            decode_responses=True
         )
    
    def is_connected(self):
        """
        Test redis connection
        """
        try:
            self.conn.set('a', 1)
            self.conn.get('a')
            return True
        except Exception as e:
            print(e)
        return False

    def get_data(self, key:str):
        data = self.conn.get(key)
        if data is None:
            return None
        try:
            return json.loads(data)
        except ValueError as e:
            print(e)
            return {'data_str': data}

    def set_data(self, key:str, val):
        try:
            val = json.dumps(val)
            self.conn.set(key, val)
            return True
        except ValueError as e:
            self.conn.set(key, str(val))
            return True
        except Exception as e:
            return False
        
    def put_data(self, key:str, val):
        old = self.get_data(key)
        insert = KeyValue.merge_dict(old, val) if old else val
        self.set_data(key, insert)
    
    def get_geneid(self, geneid:str):
        '''
        Entrez GeneID: integers
        '''
        integers = re.findall(r'\d+', str(geneid))
        if integers:
            return self.get_data(integers[0])
        return None
    
    def set_geneid(self, geneid:str, data:dict):
        integers = re.findall(r'\d+', str(geneid))
        if integers:
            return self.set_data(integers[0], data)
        return False
        
    def put_geneid(self, geneid:str, data:dict):
        integers = re.findall(r'\d+', str(geneid))
        if integers:
            return self.put_data(integers[0], data)
        return False

    def get_gi(self, gi:str):
        '''
        NCBI GI(GenInfo Identifier): 8 integers
        '''
        integers = re.findall(r'\d+', str(gi))
        if integers:
            return self.get_data(integers[0])
        return None

    def set_gi(self, gi:str, data:dict):
        integers = re.findall(r'\d+', str(gi))
        if integers:
            return self.set_data(integers[0], data)
        return False
        
    def put_gi(self, gi:str, data:dict):
        integers = re.findall(r'\d+', str(gi))
        if integers:
            return self.put_data(integers[0], data)
        return False
      
    def get_uniprotkb_acc(self, acc:str):
        '''
        UniProtKB accession: 6 or 10 alphanumerical characters
        '''
        if len(acc) in (6,10):
            return self.get_data(acc)
        return None

    def set_uniprotkb_acc(self, acc:str, data:dict):
        if len(acc) in (6,10):
            return self.set_data(acc, data)
        return False
        
    def put_uniprotkb_acc(self, acc:str, data:dict):
        if len(acc) in (6,10):
            return self.put_data(acc, data)
        return False

    def get_ncbi_acc(self, acc:str):
        '''
        NCBI sequence accession: xx_xxxxxxx.x
        <two leters>_<6 digits>.<version>
        '''
        return self.get_data(acc)

    def set_ncbi_acc(self, acc:str, data:dict):
        return self.set_data(acc, data)
        
    def put_ncbi_acc(self, acc:str, data:dict):
        return self.put_data(acc, data)
