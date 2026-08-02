import os
import json
from copy import deepcopy

class Utils:

    @staticmethod
    def scan_json(dir_path:str):
        for root, _dir, files in os.walk(dir_path):
            for file_name in files:
                path = os.path.join(root, file_name)
                if file_name.endswith('json'):
                    with open(path, 'r') as f:
                        data = json.load(f)
                        yield path, data
    
    @staticmethod
    def update_json(json_path:str, new_data:dict, keys:list=None):
        data = {}
        if os.path.isfile(json_path):
            with open(json_path, 'r') as f:
                data = json.load(f)
        
        if keys:
            _d = data
            for k in  keys[:-1]:
                if k not in _d:
                    _d[k] = {}
                _d = _d[k]
            _d[keys[-1]] = new_data
        else:
            data.update(new_data)
        print(f'update {json_path}.')
        with open(json_path, 'w') as f:
            json.dump(data, f, indent=4)
            return True

