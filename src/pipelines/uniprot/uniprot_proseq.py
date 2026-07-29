'''
example: 
    - bomics uniprot_proseq
functions:
    - input pro_seq into table uniprot_protseq
    - add new UniProt sequences into proseq_* tables
    - parse uniprot.uniprot_acc with view_proseq.seq_id
'''
from src.helper import *
from bioomics import TABLE_COLUMNS

def retrieve(params):
    table_name = params['table_name']
    chunk_size = params['chunk_size']
    # column name must be 'seq'
    query = f"""
        SELECT uniprot_acc, pro_seq AS seq
        FROM uniprot
        WHERE pro_seq IS NOT NULL
            AND uniprot_acc NOT IN (
                SELECT uniprot_acc FROM {table_name}
            )
        LIMIT {chunk_size}
    """
    qc = QueryComplex(params['verbose'])
    return qc.pack_proseq(query)

def build(data_iter, params, meta):
    uniprot_cols = TABLE_COLUMNS[params['table_name']]
    qc = QueryComplex(params['verbose'])
    bc = BuildComplex(params['verbose'])
    
    tag = 0
    for (table_name, seq_len), rows in data_iter:
        tag = 1
        # detect known seqs
        parsed, unparsed = qc.parse_seqid(table_name, seq_len, rows)
        
        # insert new seqs into proseq_* tables
        proseq_cols = ['source', 'seq']
        new_records = [(params['source'], i['seq']) for i in unparsed]
        num_insertion, fail = bc.insert_batch_records(new_records, table_name, proseq_cols)
        meta['new_insertion'] += num_insertion
        meta['fail_insert'] += len(fail)
        # parse new seqs with seqid
        new_parsed, fail = qc.parse_seqid(table_name, seq_len, unparsed)
        meta['fail_insert'] += len(fail)

        #build record for inertion into uniprot_proseq
        for row in parsed + new_parsed:
            yield tuple([row[i] for i in uniprot_cols])
            meta['records'] += 1
    if tag == 0:
        params['retrieve'] = False

if __name__ == "__main__":
    params.update({
        'chunk_size': 100_000,
        'source': 'UniProt',
        'table_name': 'uniprot_proseq',
    })

    params['retrieve'] = True
    while params['retrieve']:
        data_iter = retrieve(params)
        record_iter = build(data_iter, params, meta)
        bc = BuildComplex(params['verbose'])
        bc.insert_batch_records(record_iter, params['table_name'])
        print(meta)
