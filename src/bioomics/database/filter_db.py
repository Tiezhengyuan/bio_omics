
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
        file_iter = Dir(self.data_dir).recursive_files()
        for path in file_iter:
            pdb_id = os.path.basename(os.path.dirname(path))
            if pdb_id not in pdb_ids:
                file_name = os.path.basename(path)
                # outfmt=5
                if file_name == 'chain_sequences.uniprot.blastp5':
                    if self.verbose:
                        print(f"Try to parse {pdb_id}:{path} for the table {table_name}...")
                    yield path

    def igblastp_align_iter(self) -> Iterable:
        dir_iter = Dir(self.data_dir).recursive_dirs()
        for dir_path in dir_iter:
            aln_file, fa_file = None, None
            for file_name in os.listdir(dir_path):
                if file_name == 'chain_sequences.vregion.igblastp3':
                    aln_file = os.path.join(dir_path, file_name)
                elif file_name == 'chain_sequences.faa':
                    fa_file = os.path.join(dir_path, file_name)
            if aln_file and fa_file:
                yield aln_file, fa_file

    def igblastp_align_check_iter(self, table_name:str) -> Iterable:
        pdb_ids = QueryComplex().get_unique_values(table_name, 'pdb_id')
        dir_iter = Dir(self.data_dir).recursive_dirs()
        for dir_path in dir_iter:
            pdb_id = os.path.basename(dir_path)
            if pdb_id not in pdb_ids:
                aln_file, fa_file = None, None
                for file_name in os.listdir(dir_path):
                    if file_name == 'chain_sequences.vregion.igblastp3':
                        aln_file = os.path.join(dir_path, file_name)
                    elif file_name == 'chain_sequences.faa':
                        fa_file = os.path.join(dir_path, file_name)
                if aln_file and fa_file:
                    yield aln_file, fa_file