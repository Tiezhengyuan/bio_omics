'''
'''
import os
import json
import math
from typing import Iterable
from pathlib import Path

class IntegrateData:
    def __init__(self, meta):
        '''
        args: meta: either the entity directory (str or Path) or a dict
              with 'entity_path' and optional 'index_file' keys
        NB: the old version only accepted a dict; every caller that passed
        a plain path (genBank, refseq, integrate_ncbi_protein,
        iedb_antigen...) crashed with TypeError
        '''
        if isinstance(meta, dict):
            self.entity_path = str(meta['entity_path'])
            self.index_meta_file = meta.get(
                'index_file', os.path.join(self.entity_path, 'index.json'))
        else:
            self.entity_path = str(meta)
            self.index_meta_file = os.path.join(self.entity_path, 'index.json')
        self.index_meta = self.get_index_meta()

    def get_index_meta(self) -> dict:
        '''
        self.index_meta
        key: a certain accession
        value: dictionary
        '''
        if os.path.isfile(self.index_meta_file):
            with open(self.index_meta_file, 'r') as f:
                index_meta = json.load(f)
                return index_meta
        return {}

    def get_meta(self, input:dict=None) -> dict:
        '''
        merge input into self.index_meta and return it
        '''
        if isinstance(input, dict):
            self.index_meta.update(input)
        return self.index_meta

    def save_index_meta(self, input:dict=None) -> bool:
        if isinstance(input, dict):
            self.index_meta.update(input)
        if self.index_meta:
            with open(self.index_meta_file, 'w') as f:
                json.dump(self.index_meta, f, indent=4)
            return True
        return False
        
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
        # NB: non-numeric keys (e.g. accessions) must not break max()
        ids = [int(i['key']) for i in self.index_meta.values()
               if str(i.get('key', '')).isdigit()]
        if ids:
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
        Path(path).mkdir(parents=True, exist_ok=True)
        json_file = os.path.join(path, f'{new_id}.json')
        return json_file

    def key_json_path(self, key_value:str):
        '''
        slice key_value per 3-char as sub-dir
        '''
        id_prefix = str(key_value)[:-2]
        sub_dirs = [id_prefix[i:i+3] for i in range(0, len(id_prefix), 3)]
        path = os.path.join(self.entity_path, *sub_dirs)
        Path(path).mkdir(parents=True, exist_ok=True)
        json_file = os.path.join(path, f'{key_value}.json')
        return json_file
    
    def add_data(self, data:dict, key_value:str=None, source:str=None):
        '''
        'key' is added into new data
        key is unique id for identification of data
        key could be new_id or accession
        '''
        try:
            json_file = None
            if key_value is None:
                key_value = self.next_id()
                json_file = self.new_json_path(key_value)
            else:
                json_file = self.key_json_path(key_value)
                        
            new_data = {'key': key_value}
            new_data.update(data)
            # update index_meta
            self.index_meta[key_value] = {
                'key': key_value,
                'json_file': json_file,
                'source': [source if source else "UNKNOWN",],
            }
            with open(json_file, 'w') as f:
                json.dump(new_data, f, indent=4)
                return json_file
        except Exception as e:
            print(f"Error: can't save json file. error={e}")

    
    def save_data(self, data:dict, source:str=None) -> str:
        '''
        new data is added or data is updated
        '''
        # NB: .get so that the fallback never hits an unbound key_value
        key_value = data.get('key')
        try:
            if source and source not in self.index_meta[key_value]['source']:
                self.index_meta[key_value]['source'].append(source)
            
            json_file = self.index_meta[key_value]['json_file']
            if os.path.isfile(json_file):
                with open(json_file, 'w') as f:
                    json.dump(data, f, indent=4)
                return json_file
        except Exception as e:
            print(f"Error: can't update json file. error={e}")
        # fallback: key known but json file missing - re-add under the same key
        return self.add_data(data, key_value=key_value, source=source)
                
    def get_data(self, key_value:str) -> dict:
        if key_value in self.index_meta:
            json_file = self.index_meta[key_value]['json_file']
            if os.path.isfile(json_file):
                with open(json_file, 'r') as f:
                    data = json.load(f)
                    return data
        return {}