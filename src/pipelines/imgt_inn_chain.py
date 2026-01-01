'''
example: 
    - python src/pipelines/imgt_inn_chain.py
functions:
    - download data from IMGT
    - put data to table imgt_inn_chain
'''

from helper import *
from bioomics import Imgt, ParseImgtAnnot

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
        'table_name': 'imgt_inn_chain',
    }
    # pull records from IMGT
    inn = Imgt(args['imgt_dir'], args['verbose'])
    data_iter = pull_data(args, meta)

    # insertion
    bc = BuildComplex(data_iter, args['verbose'], args['chunk_size'])
    bc.empty_table(ars['table_name'])
    bc.insert_imgt_inn_chain()
    
    print(meta)

