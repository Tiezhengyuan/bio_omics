'''
'''
from biosequtils import Dir
import os
import json
import math
from typing import Iterable
from datetime import datetime

class IntegrateData:
    def __init__(self, entity_path:str):
        '''
        args: entity_path: store integrated data.
        args: key_index: use key as index or build new index if
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
        print(f"Try to get path of json files from {self.index_meta_file}")
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
        file_name = f"{updated_meta.get('source', '')}_meta.json"
        meta_file = os.path.join(self.entity_path, file_name)
        if os.path.isfile(meta_file):
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
        meta['start_time'] = datetime.now()
        return meta

    def save_meta(self, meta:dict):
        end_time = datetime.now()
        delta = end_time - meta['start_time']
        meta['duration'] = delta.seconds
        meta['start_time'] = meta['start_time'].strftime("%m/%d/Y, %H:%M:%S")
        meta['end_time'] = end_time.strftime("%m/%d/Y, %H:%M:%S")
        # save
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
            ids = [int(i['key']) for i in self.index_meta.values()]
            return str(max(ids) + 1)
        return '1'
    
    def new_json_path(self, new_id:str) -> str:
        '''
        id = '1234'
        path: ./12/34/1234.json
        '''
        id_prefix = str(math.floor(int(new_id)/1000))
        sub_dirs = [id_prefix[i:i+2] for i in range(0, len(id_prefix), 2)]
        path = os.path.join(self.entity_path, 'data', *sub_dirs)
        Dir(path).init_dir()
        json_file = os.path.join(path, f'{new_id}.json')
        return json_file

    def key_json_path(self, key_value:str):
        id_prefix = str(key_value)[:-2]
        sub_dirs = [id_prefix[i:i+3] for i in range(0, len(id_prefix), 3)]
        path = os.path.join(self.entity_path, 'data',  *sub_dirs)
        Dir(path).init_dir()
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
                json_file = self.new_json_path(self.next_id())
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
        try:
            key_value = data['key']
            if source and source not in self.index_meta[key_value]['source']:
                    self.index_meta[key_value]['source'].append(source)
            json_file = self.index_meta[key_value]['json_file']
            if os.path.isfile(json_file):
                with open(json_file, 'w') as f:
                    json.dump(data, f, indent=4)
                return json_file
        except exception as e:
            print(f"Error: can't update json file. error={e}")
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