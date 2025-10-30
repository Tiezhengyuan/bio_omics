
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
      
    def chain_iter(self, table_name:str=None, pdb_ids:list=None):
        qc = QueryComplex()
        ready_chain_ids = qc.unique_values(table_name, 'chain_id')
        chain_ids_pool = qc.get_chain_ids(pdb_ids)
        chain_ids = set(chain_ids_pool).difference(set(ready_chain_ids))
        if self.verbose:
            print(f"{len(chain_ids)} chains are not in {table_name}.")
        for chain_id in list(chain_ids):
            pdb_id, model_no, chain_no = chain_id.split('_')
            chains_dir = os.path.join(self.data_dir, pdb_id[:2], pdb_id, 'chains')
            faa_file = os.path.join(chains_dir, f'{model_no}_{chain_no}.faa')
            pdb_file = os.path.join(chains_dir, f'{pdb_id}_{chain_no}.pdb')
            if os.path.isfile(faa_file) and os.path.isfile(pdb_file):
                yield pdb_id, faa_file, pdb_file