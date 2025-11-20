'''
example: 
    - python src/pipelines/imgt_inn_chain.py
functions:
    - download data from IMGT
    - put data to table imgt_inn_chain
'''
from collections import defaultdict
from dotenv import load_dotenv
import os
import sys
src_dir = os.path.dirname(os.path.dirname(__file__))
if src_dir not in sys.path:
    sys.path.append(src_dir)

from bioomics import BuildComplex, Imgt, Dir, ParseImgtAnnot

def pull_data(args, meta):
    indir = os.path.join(args['imgt_dir'], '3Dstructure-DB', 'IMGT3DFlatFiles')
    for gz_file in Dir(indir).recursive_files():
        if gz_file.endswith('.inn.gz'):
            try:
                inn_data = ParseImgtAnnot(gz_file)()
                meta['inn_data'] += 1
                yield inn_data
            except Exception as e:
                print(f"{gz_file}, error={e}")
                meta['invalid_inn'] += 1

if __name__ == "__main__":
    meta = defaultdict(int)
    load_dotenv()
    args = {
        'verbose': True,
        'chunk_size': 50,
        'imgt_dir': os.getenv('imgt_dir'),
    }
    # pull records from IMGT
    inn = Imgt(args['imgt_dir'], args['verbose'])
    data_iter = pull_data(args, meta)

    # insertion
    table_name = 'imgt_inn_chain'
    bc = BuildComplex(data_iter, args['verbose'], args['chunk_size'])
    bc.empty_table(table_name)
    bc.insert_imgt_inn_chain()
    
    print(meta)

