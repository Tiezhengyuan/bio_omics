#!/usr/bin/python
'''
requirements:
    - view_pdb_chain
example:
    - python src/pipelines/protein_absd.py
function:
    - put ABSD sequences into table pro_seq and pro_source
'''

from helper import *

def pull_seq(min_len, max_len, args, meta):
    table_name = args['table_name']
    qc = QueryComplex(args['verbose'])
    query = f'''
        SELECT pro_seq AS seq
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

def pull_id(min_len, max_len, args, meta):
    table_name = args['table_name']
    query = f'''
        WITH tmp AS (
            SELECT seq_id, pro_seq
            FROM {table_name}
            WHERE pro_len >= {min_len}
                AND pro_len < {max_len}
                AND seq_id NOT IN (
                SELECT source_id
                FROM protein_source
                WHERE source_name = '{table_name}'
                    AND pro_len >= {min_len}
                    AND pro_len < {max_len}
            )
        )
        SELECT C.seq_id, S.pro_id
        FROM tmp C
        LEFT JOIN protein_seq S
        ON C.pro_seq = S.pro_seq
        WHERE S.pro_id IS NOT NULL;
    '''
    rows = QueryComplex(args['verbose']).list_data(query)
    meta['num_id'] = len(rows)
    for row in rows:
        rec = (row['pro_id'], table_name, row['seq_id'])
        yield rec
        
if __name__ == "__main__":
    args = {
        'verbose': True,
        'chunk_size': 500,
        'table_name': 'absd_seq',
    }
    step = 10
    for min_len in range(8, 200, step):
        max_len = min_len + step
        meta['min_len'] = min_len
        meta['max_len'] = max_len

        # # pull chain_seq and update pro_seq
        # seq_iter = pull_seq(min_len, max_len, args, meta)
        # bc = BuildComplex(seq_iter, args['verbose'], args['chunk_size'])
        # bc.insert_protein_seq()

        # pull chain_id and update pro_source
        id_iter = pull_id(min_len, max_len, args, meta)
        bc = BuildComplex(id_iter, args['verbose'], args['chunk_size'])
        bc.insert_protein_source()
        
        print(meta)
