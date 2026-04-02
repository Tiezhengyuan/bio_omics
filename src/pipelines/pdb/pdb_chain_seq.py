'''
requirements:
    - The table pdb_chainid is ready.
example:
    - parallel -j4 python app.py chain_seq ::: {0..237}
    - python app.py chain_seq 0
function:
    - assign unique id to chain protein sequences
'''

from src.helper import *

def retrieve_data(args, meta):
    pdb_group = args['pdb_group']
    table_name = args['table_name']
    query = f"""
        SELECT DISTINCT C.pdb_id, P.relative_path
        FROM pdb_chainid C
        LEFT JOIN meta_pdb P
            ON C.pdb_id = P.pdb_id
        WHERE P.pdb_group = {pdb_group}
            AND C.chain_id NOT IN (
                SELECT chain_id FROM {table_name}
            );
    """
    rows = QueryComplex(args['verbose']).list_data(query)
    meta['pdb_group'] = int(pdb_group)
    meta['rows'] = len(rows)
    for row in rows:
        pdb_file = os.path.join(args['pdb_dir'], row['relative_path'])
        yield {
            'pdb_id': row['pdb_id'],
            'pdb_file': pdb_file,
        }

##################################
if __name__ == "__main__":
    args = {
        'pdb_group': sys.argv[1],
        'verbose': True,
        'chunk_size': 100,
        'pdb_dir': os.getenv("pdb_dir"),
        'table_name': 'pdb_chain_seq',
    }
    rows = retrieve_data(args, meta)

    print(f"Try to insert PDB data into table {args['table_name']}...")
    bc = BuildComplex(rows, args['verbose'], args['chunk_size'])
    bc.insert_pdb_chain_seq()
    print(f"Synchronization of PDB with {args['table_name']} is good.")
    
    footer(meta)