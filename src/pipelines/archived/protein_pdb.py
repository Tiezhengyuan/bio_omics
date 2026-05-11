#!/usr/bin/python
'''
requirements:
    - view_pdb_chain
example:
    - python app.py protein_pdb
function:
    - put pdb sequences into table pro_seq and pro_source
'''
from src.helper import *

def pull_seq(min_len, max_len, params, meta):
    table_name = params['table_name']
    qc = QueryComplex(params['verbose'])
    query = f'''
        SELECT chain_seq AS seq
        FROM {table_name}
        WHERE pro_len >= {min_len}
        AND pro_len < {max_len}
        EXCEPT
        SELECT pro_seq AS seq
        FROM protein_seq
        WHERE pro_len >= {min_len}
            AND pro_len < {max_len}
    ;
    '''
    df = qc.list_data(query, True)
    meta['num_seq'] = len(df)
    if len(df) > 0:
        df['prefix'] = df['seq'].map(lambda x: x[:3])
        next_id = qc.next_id('protein_seq', 'pro_id')
        df['id'] = range(next_id, next_id+len(df))
        for i, row in df.iterrows():
            meta['num_seq'] += 1
            pro_len = len(row['seq'])
            rec = (row['id'], row['prefix'], pro_len, row['seq'],)
            yield rec

def pull_id(min_len, max_len, params, meta):
    table_name = params['table_name']
    query = f'''
        WITH tmp_chain AS (
            SELECT chain_id, chain_seq
            FROM {table_name}
            WHERE pro_len >= {min_len}
                AND pro_len < {max_len}
                AND chain_id NOT IN (
                SELECT source_id
                FROM protein_source
                WHERE source_name = '{table_name}'
                    AND pro_len >= {min_len}
                    AND pro_len < {max_len}
            )
        )
        SELECT C.chain_id, S.pro_id
        FROM tmp_chain C
        LEFT JOIN protein_seq S
        ON C.chain_seq = S.pro_seq;
    '''
    rows = QueryComplex(params['verbose']).list_data(query)
    meta['num_id'] = len(rows)
    for row in rows:
        rec = (row['pro_id'], table_name, row['chain_id'])
        yield rec
        
if __name__ == "__main__":
    params.update({
        'chunk_size': 500,
        'table_name': 'pdb_chain_seq',
    })
    for min_len in range(8, 5000, 50):
        max_len = min_len + 50
        meta['min_len'] = min_len
        meta['max_len'] = max_len

        # pull chain_seq and update pro_seq
        seq_iter = pull_seq(min_len, max_len, params, meta)
        bc = BuildComplex(seq_iter, params['verbose'], params['chunk_size'])
        bc.insert_protein_seq()

        # pull chain_id and update pro_source
        id_iter = pull_id(min_len, max_len, params, meta)
        bc = BuildComplex(id_iter, params['verbose'], params['chunk_size'])
        bc.insert_protein_source()
    
        footer(meta)