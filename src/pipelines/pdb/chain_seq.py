'''
requirements:
    - The table pdb_chainid is ready.
example:
    - parallel -j4 python app.py chain_seq ::: {0..240}
    - python app.py chain_seq 0
function:
    - assign unique id to chain protein sequences
'''

from src.helper import *
from bioomics import Pdb

@header
def retrieve(args, meta):
    pdb_group = args['pdb_group']
    table_name = args['table_name']
    query = f"""
        SELECT DISTINCT C.pdb_id, P.relative_path
        FROM pdb_chainid   C
        LEFT JOIN meta_pdb P ON C.pdb_id = P.pdb_id
        WHERE P.pdb_group = {pdb_group}
            AND C.chain_id NOT IN (
                SELECT chain_id FROM {table_name}
            );
    """
    rows = QueryComplex(args['verbose']).list_data(query)
    meta['rows'] = len(rows)
    for row in rows:
        row['pdb_file'] = os.path.join(args['rawdata_pdb_dir'], row['relative_path'])
        yield row

def build(row_iter, params, meta):
    for row in row_iter:
        pdb_id = row['pdb_id']
        pdb_file = row['pdb_file']
        p = Pdb(pdb_file, None, params['verbose'])
        p.load_structure()

        tag = 0 
        for item in p.parse_chains():
            # the column pro_seq is not null
            chain_seq = p._chain_seq(item['chain'])
            rec = (pdb_id, item['chain_id'], chain_seq)
            yield rec
            meta['records'] += 1
            tag = 1
        if tag == 0:
            if params['verbose']:
                print(f"Warning: no chains {pdb_id} are detectin in {pdb_file}")
            yield (pdb_id, None, None,)
            meta['null'] += 1

##################################
if __name__ == "__main__":
    params.update({
        'table_name': 'pdb_chain_seq',
    })
    row_iter = retrieve(params, meta)
    record_iter = build(row_iter, params, meta)

    print(f"Try to insert PDB data into table {params['table_name']}...")
    bc = BuildComplex(params['verbose'])
    bc.insert_batch_records(record_iter, params['table_name'])
   
    footer(meta)