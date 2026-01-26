'''
example: 
    - python src/pipelines/ab_abybank.py
functions:
    - download data from abYbank
    - put data to table ab_abybank
'''

from helper import *
from bioomics import Abybank


if __name__ == "__main__":
    meta = defaultdict(int)
    args = {
        'verbose': True,
        'chunk_size': 20,
        'table_name': 'ab_abybank'
    }
    # pull records from abYbank
    abYbank_pdb = Abybank().pull_pdb_records()
    meta['rows'] = len(abYbank_pdb)

    # insertion
    bc = BuildComplex(abYbank_pdb, args['verbose'], args['chunk_size'])
    bc.empty_table(args['table_name'])
    bc.insert_batch_records(args['table_name'])
    
    print(meta)

