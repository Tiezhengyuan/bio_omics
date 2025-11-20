'''
example: 
    - parallel -j4 python src/pipelines/absd_seq.py
    - python src/pipelines/absd_seq.py
functions:
    - download data from ABSD
    - put data to table absd_seq
'''
from dotenv import load_dotenv
from collections import defaultdict
import os
import sys
src_dir = os.path.dirname(os.path.dirname(__file__))
if src_dir not in sys.path:
    sys.path.append(src_dir)

from bioomics import BuildComplex, Absd, QueryComplex

def update_pool(seq_pool:dict, seq:str):
    prefix = seq[:4]
    if prefix in seq_pool:
        if seq not in seq_pool[prefix]:
            seq_pool[prefix].append(seq)
            return True
    else:
        seq_pool[prefix] = [seq,]
        return True

def get_data(seq_pool:dict, args):
    ver = args['absd_version']
    data_dir = os.path.join(args['absd_dir'], ver, 'data')
    a = Absd(data_dir, args['verbose'])
    for specie, path in a.scan_files('fasta'):
        for record in a.scan_data(path):
            seq = str(record.seq)
            is_new = update_pool(seq_pool, seq)
            if is_new:
                absd_id = record.id.split("|||")[0]
                yield absd_id, specie, ver, len(seq), seq, record.id


if __name__ == "__main__":
    meta = defaultdict(int)
    load_dotenv()
    args = {
        'verbose': True,
        'chunk_size': 100,
        'absd_dir': os.getenv('absd_dir'),
        'absd_version': os.getenv('absd_version'),
    }
    table_name = 'absd_source'
    seq_pool = {}
    seqs = QueryComplex(args['verbose']).unique_values(table_name, 'pro_seq')
    for seq in seqs:
        update_pool(seq_pool, seq)

    # filter
    record_iter = get_data(seq_pool, args)

    # insertion
    bc = BuildComplex(record_iter, args['verbose'], args['chunk_size'])
    bc.insert_absd_seq()
    
    print(meta)

