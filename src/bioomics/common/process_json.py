'''
'''
import os
import json

class ProcessJson:
            
    @staticmethod
    def load_json(infile):
        try:
            with open(infile, 'r') as f:
                return json.load(f)
        except Exception as e:
            pass
        return {}

    @staticmethod
    def save_json(data, outfile):
        try:
            with open(outfile, 'w') as f:
                json.dump(data, f, indent=4, sort_keys=True)
            return True
        except Exception as e:
            pass
        return False

    @staticmethod
    def iter_to_json(_iter, outdir:str, chunk_size:int=10_000):
        pool, file_list = [], []
        n, ix = 0, 0
        for record in _iter:
            pool.append(record)
            n += 1
            if len(pool) >= chunk_size:
                outfile = os.path.join(outdir, f"{ix}.json")
                with open(outfile, 'w') as f:
                    json.dump(pool, f, indent=4)
                pool = []
                ix += 1
                file_list.append(outfile)
        else:
            if pool:
                outfile = os.path.join(outdir, f"{ix}.json")
                with open(outfile, 'w') as f:
                    json.dump(pool, f, indent=4)
                file_list.append(outfile)
        print('total', n)
        return file_list
    
    @staticmethod
    def get_value(infile:str, keys:list):
        if os.path.isfile(infile) and keys:
            data = ProcessJson.load_json(infile)
            for key in keys[:-1]:
                if key in data:
                    data = data[key]
                else:
                    return None
            return data.get(keys[-1])
        return None