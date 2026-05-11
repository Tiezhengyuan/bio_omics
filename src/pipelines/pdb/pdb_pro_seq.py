'''
example: 
    - python app.py pdb_pro_seq
functions:
    - add new PDB sequences into proseq_* tables
    - parse pdb_chain_seq.chain_id with view_proseq.seq_id
'''
from src.helper import *
from bioomics import TABLE_COLUMNS

def retrieve(params, pdb_group):
    table_name = params['table_name']
    # column name must be 'seq'
    query = f"""
        SELECT s.pdb_id, s.chain_id, s.chain_seq AS seq
        FROM pdb_chain_seq s
        LEFT JOIN meta_pdb m ON s.pdb_id = m.pdb_id
        WHERE m.pdb_group = {pdb_group}
            AND s.chain_seq IS NOT NULL
            AND s.chain_id NOT IN (
                SELECT chain_id FROM {table_name}
                WHERE chain_id IS NOT NULL
            )
    """
    qc = QueryComplex(params['verbose'])
    return qc.pack_proseq(query)

def build(data_iter, params, meta):
    cols = TABLE_COLUMNS[params['table_name']]
    qc = QueryComplex(params['verbose'])
    bc = BuildComplex(params['verbose'])
    
    for conditions, rows in data_iter:
        # detect known seqs
        parsed, unparsed = qc.parse_seqid(conditions, rows)
        # insert new seqs
        num_insertion, fail = bc.insert_proseq(unparsed, params['source'])
        meta['new_insertion'] += num_insertion
        meta['fail_insert'] += len(fail)
        # parse new seqs
        new_parsed, fail = qc.parse_seqid(conditions, unparsed)
        meta['fail_insert'] += len(fail)

        #build record
        for row in parsed + new_parsed:
            yield tuple([row[i] for i in cols])
            meta['records'] += 1    

if __name__ == "__main__":
    params.update({
        'source': 'PDB',
        'table_name': 'pdb_pro_seq',
    })
    
    qc = QueryComplex(params['verbose'])
    max_group = qc.max_value('meta_pdb', 'pdb_group')
    for pdb_group in range(int(max_group)+1):
        data_iter = retrieve(params, pdb_group)
        record_iter = build(data_iter, params, meta)
        bc = BuildComplex(params['verbose'])
        bc.insert_batch_records(record_iter, params['table_name'])
        print(meta)
