'''
requirements:
    - PDB files should be synchronized to the local.
    - The table pdb is ready.
example:
    - parallel -j8 python app.py pdb_compound ::: {0..237}
    - python app.py pdb_compound 0
function:
    - retrieve compound in head of PDB
'''

from src.helper import *

def retrieve_pdb(params, meta):
    pdb_group = params['pdb_group']
    table_name = params['table_name']
    query = f"""
        SELECT pdb_id, relative_path
        FROM meta_pdb
        WHERE pdb_group = {pdb_group}
            AND pdb_id NOT IN (
                SELECT DISTINCT pdb_id FROM {table_name}
            );
    """
    rows = QueryComplex(params['verbose']).list_data(query)
    meta['pdb_group'] = int(pdb_group)
    meta['rows'] += len(rows)
    for row in rows:
        pdb_file = os.path.join(params['pdb_dir'], row['relative_path'])
        yield {
            'pdb_id': row['pdb_id'],
            'pdb_file': pdb_file,
        }

##################################
if __name__ == "__main__":
    params.update({
        'chunk_size': 5,
        'table_name': 'pdb_compound',
    })

    print(f"Try to insert PDB data into table {params['table_name']}...")
    record_iter = retrieve_pdb(params, meta)
    bc = BuildComplex(record_iter, params['verbose'], params['chunk_size'])
    bc.insert_pdb_compound()
    print(f"Synchronization of PDB with {params['table_name']} is good.")

    footer(meta)

