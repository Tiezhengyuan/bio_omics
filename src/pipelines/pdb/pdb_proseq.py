'''
example: 
    - bomics pdb_proseq
functions:
    - add new PDB sequences into proseq_* tables
    - parse pdb_chain_seq.chain_id with view_proseq.seq_id
'''
from src.helper import *
from bioomics import TABLE_COLUMNS

def retrieve(params, seq_len):
    table_name = params['table_name']
    # column name must be 'seq'
    query = f"""
        SELECT pdb_id, chain_id, chain_seq AS seq
        FROM pdb_chain_seq
        WHERE pro_len = {seq_len}
            AND chain_seq IS NOT NULL
            AND chain_id NOT IN (
                SELECT chain_id FROM {table_name}
                WHERE chain_id IS NOT NULL
            )
    """
    qc = QueryComplex(params['verbose'])
    rows = qc.list_data(query)
    return qc.pack_proseq(rows)

def build(data_iter, params, meta):
    cols = TABLE_COLUMNS[params['table_name']]
    qc = QueryComplex(params['verbose'])
    bc = BuildComplex(params['verbose'])
    
    for (table_name, seq_len), rows in data_iter:
        # detect known seqs
        parsed, unparsed = qc.parse_seqid(table_name, seq_len, rows)
        
        # insert new seqs
        num_insertion, fail = bc.insert_proseq(params['verbose'], table_name, unparsed)
        meta['new_insertion'] += num_insertion
        meta['fail_insert'] += len(fail)
        # parse new seqs
        new_parsed, fail = qc.parse_seqid(table_name, seq_len, unparsed)
        meta['fail_insert'] += len(fail)

        #build record
        for row in parsed + new_parsed:
            yield tuple([row[i] for i in cols])
            meta['records'] += 1    

if __name__ == "__main__":
    params.update({
        'source': 'PDB',
        'table_name': 'pdb_proseq',
    })

    pool = QueryComplex(params['verbose']).unique_values('pdb_chain_seq', 'pro_len')
    for seq_len in pool:
        meta['seq_len'] = seq_len
        data_iter = retrieve(params, seq_len)
        record_iter = build(data_iter, params, meta)
        bc = BuildComplex(params['verbose'])
        bc.insert_batch_records(record_iter, params['table_name'])
        print(meta)
