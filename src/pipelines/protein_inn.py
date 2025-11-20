#!/usr/bin/python
'''
requirements:
    - view_pdb_chain
example:
    - python src/pipelines/protein_inn.py
function:
    - put INN sequences into table pro_seq and pro_source
'''
from collections import defaultdict
import os
import sys
src_dir = os.path.dirname(os.path.dirname(__file__))
if src_dir not in sys.path:
    sys.path.append(src_dir)

from bioomics import BuildComplex, QueryComplex

def pull_seq(args, meta):
    table_name = args['table_name']
    qc = QueryComplex(args['verbose'])
    query = f'''
        SELECT chain_seq AS seq
        FROM {table_name}
        EXCEPT
        SELECT pro_seq AS seq
        FROM protein_seq
    ;
    '''
    df = qc.list_data(query, True)
    meta['num_seq'] = len(df)
    if len(df) > 0:
        df['prefix'] = df['seq'].map(lambda x: x[:3])
        next_id = qc.next_id('protein_seq', 'pro_id')
        df['id'] = range(next_id, next_id+len(df))
        for i, row in df.iterrows():
            meta['records'] += 1
            pro_len = len(row['seq'])
            rec = (row['id'], row['prefix'], pro_len, row['seq'],)
            yield rec

def pull_id(args, meta):
    table_name = args['table_name']
    query = f'''
        WITH tmp_chain AS (
            SELECT chain_id, chain_seq
            FROM {table_name}
            WHERE chain_id NOT IN (
                SELECT source_id
                FROM protein_source
                WHERE source_name = '{table_name}'
            )
        )
        SELECT C.chain_id, S.pro_id
        FROM tmp_chain C
        LEFT JOIN protein_seq S
        ON C.chain_seq = S.pro_seq;
    '''
    rows = QueryComplex(args['verbose']).list_data(query)
    meta['num_records'] = len(rows)
    for row in rows:
        rec = (row['pro_id'], table_name, row['chain_id'])
        yield rec
        
if __name__ == "__main__":
    meta = defaultdict(int)
    args = {
        'verbose': True,
        'chunk_size': 100,
        'table_name': 'imgt_inn_chain',
    }

    # pull chain_seq and update pro_seq
    record_iter = pull_seq(args, meta)
    bc = BuildComplex(record_iter, args['verbose'], args['chunk_size'])
    bc.insert_protein_seq()

    # pull chain_id and update pro_source
    record_iter = pull_id(args, meta)
    bc = BuildComplex(record_iter, args['verbose'], args['chunk_size'])
    bc.insert_protein_source()

    print(meta)