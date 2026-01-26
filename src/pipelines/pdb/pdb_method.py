'''
note:
    - PDB files should be synchronized to the local.
    - The table pdb is ready.
example:
    - parallel -j4 python src/pipelines/pdb_method.py ::: {0..117}
    - python src/pipelines/pdb_method.py 0
function:
    synchronize PDB/UniProt/alignments with data in database COMPLEX.
'''

from helper import *

def retrieve_pdb(pdb_group, args, meta):
    table_name = args['table_name']
    query = f"""
        SELECT pdb_id, relative_path
        FROM meta_pdb
        WHERE pdb_group = {pdb_group}
            AND pdb_id NOT IN (
                SELECT DISTINCT pdb_id FROM {table_name}
            );
    """
    rows = QueryComplex(args['verbose']).list_data(query)
    meta['pdb_group'] = pdb_group
    meta['rows'] += len(rows)
    for row in rows:
        pdb_file = os.path.join(args['pdb_dir'], row['relative_path'])
        yield {
            'pdb_id': row['pdb_id'],
            'pdb_file': pdb_file,
        }

##################################
if __name__ == "__main__":
    args = {
        'verbose': True,
        'pdb_dir': os.getenv('pdb_dir'),
        'chunk_size': 100,
        'table_name': 'pdb_method',
    }
    pdb_group = sys.argv[1]
    args_iter = retrieve_pdb(pdb_group, args, meta)

    table_name = args['table_name']
    print(f"Try to insert PDB data into table {table_name}...")
    bc = BuildComplex(args_iter, args['verbose'], args['chunk_size'])
    bc.insert_pdb_method()
    print(f'Synchronization of PDB with {table_name} is good.')

    print(meta)
