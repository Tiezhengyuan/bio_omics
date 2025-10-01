import pandas as pd
from typing import Iterable

from .ncbi import NCBI

class ParseID:
    def __init__(self, nrows:int=None, nworks:int=None, chunk_size:int=None):
        self.nrows = nrows
        self.chunk_size = chunk_size if isinstance(chunk_size, int) \
            and chunk_size > 100 else 100000
        self.nworks = nworks if isinstance(nworks, int) \
            and nworks > 1 else 8

    def gene_map0(self, gz_file:str, index_key:str) -> Iterable:
        '''
        no chunk
        '''
        print(f"Try to read {gz_file}...")
        df = pd.read_csv(gz_file, compression='gzip', sep='\t', \
            header=0, nrows=self.nrows)
        df = df.set_index(df[index_key]).drop(index_key, axis=1)

        print('Aggregate and slice data into chunk parts...')
        data_pool, chunk_data, n = [], [], 0
        groups = df.groupby(df.index, observed=True)
        return groups

    def gene_map(self, gz_file:str, index_key:str) -> Iterable:
        '''
        example:
            file name: gene_refseq_uniprotkb_collab.gz
            index by refseq or unirpotkb accession
        '''
        print(f"Try to read {gz_file}...")
        df = pd.read_csv(gz_file, compression='gzip', sep='\t', \
            header=0, nrows=self.nrows)
        df = df.set_index(df[index_key]).drop(index_key, axis=1)

        print('Aggregate and slice data into chunk parts...')
        data_pool, chunk_data, n = [], [], 0
        groups = df.groupby(df.index, observed=True)
        for acc, group in groups:
            n += 1
            # print(acc, end='\t')
            chunk_data.append((acc, group))
            if len(chunk_data) >= self.chunk_size:
                print(f"\tPush chunk data {n}...")
                data_pool.append(chunk_data)
                chunk_data = []
            if len(data_pool) >= self.nworks:
                yield data_pool
                data_pool = []
        else:
            if chunk_data:
                data_pool.append(chunk_data)
        yield data_pool

    def gene_map2(self, gz_file:str, index_key:str) -> Iterable:
        '''
        used for ThreadPoolExecutor()
        example:
            file name: gene_refseq_uniprotkb_collab.gz
            index by refseq or unirpotkb accession
        '''
        print(f"Try to read {gz_file}...")
        df = pd.read_csv(gz_file, compression='gzip', sep='\t', \
            header=0, nrows=self.nrows)
        df = df.set_index(df[index_key]).drop(index_key, axis=1)

        print('Aggregate and slice data into chunk parts...')
        chunk_data, n = [], 0
        groups = df.groupby(df.index, observed=True)
        for acc, group in groups:
            n += 1
            # print(acc, end='\t')
            chunk_data.append((acc, group))
            if len(chunk_data) >= self.chunk_size:
                print(f"\tPush chunk data {n}...")
                yield chunk_data
                chunk_data = []
        else:
            if chunk_data:
                print(f"\tPush chunk data {n}...")
                yield chunk_data
