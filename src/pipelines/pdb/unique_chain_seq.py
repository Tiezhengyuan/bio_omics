'''
requirements:
    - The table pdb_chain_seq is ready.
example:
    - python app.py unique_chain_seq
function:
    - assign unique id to chain protein sequences
'''

from src.helper import *
from bioomics import UpdateComplex

def unique_seqs(params, meta):
    table_name = params['table_name']
    query = f"""
        SELECT first_chain_id, chain_seq
        FROM {table_name}
        WHERE chain_seq IS NOT NULL
            AND first_chain_id IS NOT NULL
    ;"""
    rows = QueryComplex(params['verbose']).list_data(query)
    meta['unique_chains'] = len(rows)
    res = defaultdict()
    for row in rows:
        seq, chain_id = row['chain_seq'], row['first_chain_id']
        prefix = seq[:3]
        if prefix not in res:
            res[prefix] = defaultdict()
        res[prefix][seq] = chain_id
    return res

def parse_seq(seq_dict, params, meta):
    table_name = params['table_name']
    query = f"""
        SELECT pdb_id, chain_id, chain_seq
        FROM {table_name}
        WHERE first_chain_id IS NULL
            AND chain_seq IS NOT NULL
        ORDER BY pdb_id, chain_id
    ;"""
    rows = QueryComplex(params['verbose']).list_data(query)
    meta['rows'] = len(rows)
    for row in rows:
        seq = row['chain_seq']
        prefix = seq[:3]
        rec = None
        if prefix in seq_dict and seq in seq_dict[prefix]:
            meta['dup_seq'] += 1
            rec = (seq_dict[prefix][seq], row['pdb_id'], row['chain_id'])
        else:
            meta['new_seq'] += 1
            if prefix not in seq_dict:
                seq_dict[prefix] = defaultdict()
            seq_dict[prefix][seq] = row['chain_id']
            rec = row['chain_id'], row['pdb_id'], row['chain_id']
        if params['verbose']:
            print(rec)
        yield rec

##################################
if __name__ == "__main__":
    params.update({
        'table_name': 'pdb_chain_seq',
    })

    # get unique sequences with first_chain_id
    seq_dict = unique_seqs(params, meta)

    # detect unique sequences and parse sequences with first_chain_id
    data = parse_seq(seq_dict, params, meta)

    print(f"Try to update table {params['table_name']}...")
    c = UpdateComplex(params['verbose'])
    c.update_chain_seq(params['table_name'], data)

    footer(meta)