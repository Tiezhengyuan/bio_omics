
'''
retrieve table and compare with local files 
return filtered data
'''
import os
from typing import Iterable

from .query_complex import QueryComplex
from ..common.dir import Dir

class FilterDb:
    def __init__(self, data_dir:str):
        self.data_dir = data_dir

    def uniprot_align(self) -> Iterable:
        # collect aln files of which data not in DB
        pdb_ids = QueryComplex().unique_values('uniprot_align', 'pdb_id')
        file_iter = Dir(self.data_dir).recursive_files()
        for path in file_iter:
            pdb_id = os.path.basename(os.path.dirname(path))
            if pdb_id not in pdb_ids:
                file_name = os.path.basename(path)
                # outfmt=5
                if file_name == 'chain_sequences.uniprot.blastp5':
                    yield path