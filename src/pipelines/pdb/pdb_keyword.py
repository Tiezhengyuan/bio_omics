'''
note:
    - PDB files should be synchronized to the local.
    - The table pdb is ready.
example:
    - parallel -j4 python app.py pdb_keyword ::: {0..240}
    - python app.py pdb_keyword 0
function:
    synchronize PDB/UniProt/alignments with data in database COMPLEX.
'''

from src.helper import *
from bioomics import Pdb

@header
def retrieve(params, meta):
    pdb_group = params['pdb_group']
    table_name = params['table_name']
    query = f"""
        SELECT pdb_id, relative_path FROM meta_pdb
        WHERE pdb_group = {pdb_group}
            AND pdb_id NOT IN (
                SELECT DISTINCT pdb_id FROM {table_name}
            );
    """
    rows = QueryComplex(params['verbose']).list_data(query)
    meta['rows'] = len(rows)
    for row in rows:
        row['pdb_file'] = os.path.join(params['rawdata_pdb_dir'], row['relative_path'])
        yield row

def build(row_iter, params, meta):
    for row in row_iter:
        p = Pdb(row['pdb_file'], None, params['verbose'])
        p.load_structure()
        for keyword in p.keywords:
            meta['records'] += 1
            yield (row['pdb_id'], keyword)

##################################
if __name__ == "__main__":
    params.update({
        'chunk_size': 500,
        'table_name': 'pdb_keyword',
    })
    
    row_iter = retrieve(params, meta)
    record_iter = build(row_iter, params, meta)

    print(f"Try to insert PDB data into table {params['table_name']}...")
    bc = BuildComplex(params['verbose'], params['chunk_size'])
    bc.insert_batch_records(record_iter, params['table_name'])

    footer(meta)