
'''
retrieve table and compare with local files 
return filtered data
'''
import os
from typing import Iterable

from .query_complex import QueryComplex
from ..common.dir import Dir

class FilterDb:
    def __init__(self, data_dir:str, verbose:bool=False):
        self.data_dir = data_dir
        self.verbose = verbose

    def uniprot_align_iter(self) -> Iterable:
        file_iter = Dir(self.data_dir).recursive_files()
        for path in file_iter:
            # outfmt=5
            if path.endswith('chain_sequences.uniprot.blastp5'):
                yield path

    def uniprot_align_check_iter(self, table_name:str) -> Iterable:
        # collect aln files of which data not in DB
        pdb_ids = QueryComplex().get_unique_values(table_name, 'pdb_id')
        align_iter = self.uniprot_align_iter()
        for path in align_iter:
            pdb_id = os.path.basename(os.path.dirname(path))
            if pdb_id not in pdb_ids:
                if self.verbose:
                    print(f"Try to parse {pdb_id}:{path} for table {table_name}...")
                yield path

    def igblastp_align_iter(self) -> Iterable:
        file_iter = Dir(self.data_dir).recursive_files()
        for path in file_iter:
            if path.endswith('chain_sequences.vregion.igblastp3'):
                _dir = os.path.dirname(path)
                fa_file = os.path.join(_dir, 'chain_sequences.faa')
                if os.path.isfile(fa_file):
                    yield path, fa_file

    def igblastp_align_check_iter(self, table_name:str) -> Iterable:
        pdb_ids = QueryComplex().get_unique_values(table_name, 'pdb_id')
        align_iter = self.igblastp_align_iter()
        for aln_file, fa_file in align_iter:
            pdb_id = os.path.basename(os.path.dirname(aln_file))
            if pdb_id not in pdb_ids:
                if self.verbose:
                    print(f"Try to parse {pdb_id}:{aln_file} and {fa_file} for table {table_name}...")
                yield aln_file, fa_file
      
