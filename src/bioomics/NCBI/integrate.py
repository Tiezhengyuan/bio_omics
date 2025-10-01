

class Integrate:
    source = 'NCBI'
    meta_file_name = 'ncbi_meta.json'
    default_count = {
        'proteins': 0,
        'updated_proteins': 0,
        'skipped_proteins': 0,
        'new_proteins': 0,
    }

    def __init__(self, local_path:str, overwrite:bool=None, nworks:int=None):
        self.local_path = local_path
        self.overwrite = True if overwrite else False
        self.nworks = nworks

    def _join_meta(self, results, count_key:str):
        if count_key:
            if count_key not in self.meta:
                self.meta[count_key] = dict(self.default_count)
            for count, index_meta in results:
                for k,v in count.items():
                    self.meta[count_key][k] += v
                self.integrate.index_meta.update(index_meta)
        else:
            for count, index_meta in results:
                self.meta.update(count)
                self.integrate.index_meta.update(index_meta)