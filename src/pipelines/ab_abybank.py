'''
example: 
    - python src/pipelines/ab_abybank.py
functions:
    - download data from abYbank
    - put data to table ab_abybank
'''
from collections import defaultdict
import os
import sys
src_dir = os.path.dirname(os.path.dirname(__file__))
if src_dir not in sys.path:
    sys.path.append(src_dir)

from bioomics import BuildComplex, Abybank


if __name__ == "__main__":
    meta = defaultdict(int)
    args = {
        'verbose': True,
        'chunk_size': 20,
    }
    # pull records from abYbank
    abYbank_pdb = Abybank().pull_pdb_records()
    meta['rows'] = len(abYbank_pdb)

    # insertion
    table_name = 'ab_abybank'
    bc = BuildComplex(abYbank_pdb, args['verbose'], args['chunk_size'])
    bc.empty_table(table_name)
    bc.insert_ab_abybank()
    
    print(meta)

