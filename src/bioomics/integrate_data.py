'''
'''
from biosequtils import Dir
import os
import json
import math
from typing import Iterable

class IntegrateData:
    def __init__(self, entity_path:str):
        '''
        args: entity_path: store integrated data.
        '''
        self.entity_path = entity_path
        Dir(self.entity_path).init_dir()
        self.index_meta = self.get_index_meta()

    def get_index_meta(self) -> dict:
        '''
        self.index_meta
        key: a certain accession
        value: dictionary
        '''
        self.index_meta_file = os.path.join(self.entity_path, 'index_meta.json')
        if os.path.isfile(self.index_meta_file):
            with open(self.index_meta_file, 'r') as f:
                index_meta = json.load(f)
                return index_meta
        return {}

    def save_index_meta(self, input:dict=None) -> bool:
        if isinstance(input, dict):
            self.index_meta.update(input)
        if self.index_meta:
            with open(self.index_meta_file, 'w') as f:
                json.dump(self.index_meta, f, indent=4)
            return True
        return False

    def get_meta(self, updated_meta:dict):
        '''
        meta varies by database
        '''
        meta = {}
        source = meta.get('source', '')
        meta_file = os.path.join(self.entity_path, f"{source}_meta.json")
        if os.path.isfile(meta_file):
            print('get meta.')
            with open(meta_file, 'r') as f:
                meta = json.load(f)
        else:
            meta = {
                'entity_path': self.entity_path,
                'index_meta_file': self.index_meta_file,
                'meta_file': meta_file,
            }
        # update meta
        meta.update(updated_meta)
        return meta

    def save_meta(self, meta:dict):
        with open(meta['meta_file'], 'w') as f:
            json.dump(meta, f, indent=4)
        return meta['meta_file']
        
    def scan(self) -> Iterable:
        '''
        Note: self.index_meta could be updated
        ?? memroy leak
        '''
        files = [i['json_file'] for i in self.index_meta.values() if 'json_file' in i]
        for file in files:
            if os.path.isfile(file):
                with open(file, 'r') as f:
                    data = json.load(f)
                    yield data

    def next_id(self) -> str:
        if self.index_meta:
            ids = [int(i['ID']) for i in self.index_meta.values()]
            return str(max(ids) + 1)
        return '1'
    
    def new_json_path(self, new_id:str) -> str:
        '''
        id = '1234'
        path: ./12/34/1234.json
        '''
        id_prefix = str(math.floor(int(new_id)/1000))
        sub_dirs = [id_prefix[i:i+2] for i in range(0, len(id_prefix), 2)]
        path = os.path.join(self.entity_path, *sub_dirs)
        Dir(path).init_dir()
        json_file = os.path.join(path, f'{new_id}.json')
        return json_file
    
    def add_data(self, data:dict, key_value:str=None):
        '''
        'key' and 'ID' are added into new data
        key is unique id for identification of data
        key could be new_id or accession
        '''
        new_id = self.next_id()
        json_file = self.new_json_path(new_id)
        if key_value is None:
            key_value = new_id
        new_data = {
            'ID': new_id,
            'key': key_value
        }
        new_data.update(data)
        # update index_meta
        self.index_meta[key_value] = {
            'ID': new_id,
            'key': key_value,
            'json_file': json_file,
        }
        with open(json_file, 'w') as f:
            json.dump(new_data, f, indent=4)
            return json_file
    
    def save_data(self, data:dict) -> str:
        if 'key' in data:
            key_value = data['key']
            json_file = self.index_meta[key_value]['json_file']
            if os.path.isfile(json_file):
                with open(json_file, 'w') as f:
                    json.dump(data, f, indent=4)
                return json_file
        # add data
        return self.add_data(data)
                
    def get_data(self, key_value:str) -> dict:
        if key_value in self.index_meta:
            json_file = self.index_meta[key_value]['json_file']
            if os.path.isfile(json_file):
                with open(json_file, 'r') as f:
                    data = json.load(f)
                    return data
        return {}