from datetime import datetime
import json
import os
from biosequtils import Dir

class ProteinMeta:
    
    def __init__(self, result_dir:str, entity:str, source:str):
        self.result_dir = result_dir
        self.entity = entity
        self.source = source

    def entity_path(self):
        entity_path = os.path.join(self.result_dir, self.entity, self.source)
        Dir(entity_path).init_dir()
        return entity_path

    def meta_path(self):
        meta_path = os.path.join(self.result_dir, self.entity, 'meta')
        Dir(meta_path).init_dir()
        return meta_path

    def meta_file(self):
        meta_path = self.meta_path()
        name = f'{self.source}_meta.json'
        return os.path.join(meta_path, name)

    def index_file(self):
        meta_path = self.meta_path()
        name = f'{self.source}_index.json'
        return os.path.join(meta_path, name)

    def get_meta(self, default_meta:dict=None):
        '''
        return meta mapping
        '''
        # override default meta if IEDB_meta.json exists
        meta_file = self.meta_file()
        if os.path.isfile(meta_file):
            with open(meta_file, 'r') as f:
                meta = json.load(f)
        else:
            # create new meta
            meta = {} if default_meta is None else default_meta
            meta.update({
                'entity': self.entity,
                'entity_path': self.entity_path(),
                # meta path store IEDB_meta.json and IEBD_index.json
                'meta_path': self.meta_path(),
                'meta_file': meta_file,
                'index_file': self.index_file(),
            })
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

    def get_index_meta(self, index_meta_file:str=None):
        index_file = self.index_file()
        if os.path.isfile(index_file):
            with open(index_file, 'r') as f:
                return json.load(f)
        return {}