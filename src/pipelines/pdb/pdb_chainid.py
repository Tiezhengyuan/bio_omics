'''
note:
    - PDB files should be synchronized to the local.
    - The table pdb is ready.
example:
    - parallel -j8 python app.py pdb_chainid ::: {0..240}
    - python app.py pdb_chainid 0
function:
    - parse chain by pdb-model-chain and assign with chain_id.
'''
from src.helper import *
from bioomics import Pdb

@header
def retrieve_pdb(params, meta):
    qc = QueryComplex(params['verbose'])
    rows = qc.meta_pdb(params['pdb_group'], params['table_name'])
    print(f"New PDB: {len(rows)}")
    meta['rows'] = len(rows)
    for row in rows:
        row['pdb_file'] = os.path.join(params['rawdata_pdb_dir'], row['relative_pdb'])
        yield row

def build(row_iter, params, meta):
    for row in row_iter:
        '''
        chain_id is unique and could be null
        '''
        p = Pdb(row['pdb_file'], None, params['verbose'])
        p.load_structure()

        tag = 0
        for item in p.parse_chains():
            rec = [item[i] for i in params['table_cols']]
            yield tuple(rec)
            tag = 1
            meta['records'] += 1
        if tag == 0:
            if params['verbose']:
                print(f"Warning: no chains detected in {row['pdb_file']}.")
            yield (row['pdb_id'],) + (None,) * 4,
            meta['null'] += 1
    
##################################
if __name__ == "__main__":
    params.update({
        'chunk_size': 50,
        'table_name': 'pdb_chainid',
        'table_cols': ['chain_id', 'pdb_id', 'model_no', 'chain_no', 'idx'],
    })
    row_iter = retrieve_pdb(params, meta)
    record_iter = build(row_iter, params, meta)

    print(f"Try to insert chain id data into table {params['table_name']}...")
    bc = BuildComplex(params['verbose'], params['chunk_size'])
    bc.insert_batch_records(record_iter, params['table_name'], params['table_cols'])

    footer(meta)
