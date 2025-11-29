'''
'''
import os
from .dir import Dir


class ProcessFile:

    @staticmethod
    def get_meta(outputs_dir):
        ready = {}
        file_iter = Dir(outputs_dir).recrusive_files()
        for path in file_iter:
            if path.endswith('meta.json'):
                pdb_id = os.path.basename(os.path.dirname(path))
                ready[pdb_id] = path
        return ready

    @staticmethod
    def iter_files(outputs_dir, file_name):
        file_iter = Dir(outputs_dir).recrusive_files()
        for path in file_iter:
            if path.endswith(file_name):
                pdb_id = os.path.basename(os.path.dirname(path))
                yield pdb_id, path

    @staticmethod
    def to_label_csv(df, outfile):
        df = df.sample(frac=1, replace=False, random_state=1)
        df.to_csv(outfile, index=False, header=False)
        print(df.shape, outfile)

    @staticmethod
    def iter_to_txt(_iter, outprefix:str, chunk_size:int=10_000):
        pool, file_list = [], []
        n, ix = 0, 0
        for pdb_id, path in _iter:
            pool.append(path)
            n += 1
            if len(pool) >= chunk_size:
                outfile = f"{outprefix}_{ix}.txt"
                with open(outfile, 'w') as f:
                    f.write('\n'.join(pool))
                pool = []
                ix += 1
                file_list.append(outfile)
        else:
            if pool:
                outfile = f"{outprefix}_{ix}.txt"
                with open(outfile, 'w') as f:
                    f.write('\n'.join(pool))
                file_list.append(outfile)
        print('total', n)
        return file_list
