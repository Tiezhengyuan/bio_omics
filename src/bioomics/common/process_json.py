'''
'''
import os
import json
from pathlib import Path
class ProcessJson:
            
    @staticmethod
    def load_json(infile, default=None):
        if default is None:
            default = {}
        try:
            with open(infile, 'r') as f:
                return json.load(f)
        except Exception as e:
            print(f"Failure: can't read {infile}, error={e}")
        return default

    @staticmethod
    def save_json(data, outfile) -> bool:
        try:
            with open(outfile, 'w') as f:
                json.dump(data, f, indent=4, sort_keys=True, default=str)
            return True
        except Exception as e:
            print(f"Failure: can't write {outfile}, error={e}")
        return False

    @staticmethod
    def iter_to_json(_iter, outdir:Path, chunk_size:int=10_000):
        pool, file_list = [], []
        n, ix = 0, 0
        for record in _iter:
            pool.append(record)
            n += 1
            if len(pool) >= chunk_size:
                outfile = outdir / f"{ix}.json"
                with open(outfile, 'w') as f:
                    json.dump(pool, f, indent=4, default=str)
                pool = []
                ix += 1
                file_list.append(outfile)
        else:
            if pool:
                outfile = outdir / f"{ix}.json"
                with open(outfile, 'w') as f:
                    json.dump(pool, f, indent=4, default=str)
                file_list.append(outfile)
        print('total', n)
        return file_list
    
    @staticmethod
    def get_value(infile:Path, keys:list):
        if infile.isfile() and keys:
            data = ProcessJson.load_json(infile)
            for key in keys[:-1]:
                if key in data:
                    data = data[key]
                else:
                    return None
            return data.get(keys[-1])
        return None