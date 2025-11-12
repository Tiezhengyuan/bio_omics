'''
note:
    - PDB files should be synchronized to the local.
    - The table pdb is ready.
example:
    - python src/pipelines/put_keyword.py 46
function:
    synchronize PDB/UniProt/alignments with data in database COMPLEX.
'''
from dotenv import load_dotenv
from collections import defaultdict
import os
import sys
src_dir = os.path.dirname(os.path.dirname(__file__))
if src_dir not in sys.path:
    sys.path.append(src_dir)

from bioomics import QueryComplex,  BuildComplex

def retrieve_pdb(pdb_group, args, meta):
    query = f"""
        SELECT pdb_id, pdb_file
        FROM meta_pdb
        WHERE pdb_group = {pdb_group}
            AND pdb_id NOT IN (
                SELECT DISTINCT pdb_id
                FROM pdb_keyword
            );
    """
    rows = QueryComplex(args['verbose']).list_data(query)
    meta['pdb_group'] = int(pdb_group)
    meta['rows'] += len(rows)
    return rows

##################################
if __name__ == "__main__":
    meta =defaultdict(int)
    load_dotenv()
    args = {
        'verbose': True,
        'pdb_dir': os.getenv("pdb_dir"),
        'chunk_size': 100,
    }
    pdb_group = sys.argv[1]
    args_iter = retrieve_pdb(pdb_group, args, meta)

    table_name = 'pdb_keyword'
    print(f"Try to insert PDB data into table {table_name}...")
    bc = BuildComplex(args_iter, args['verbose'], args['chunk_size'])
    bc.insert_pdb_keyword()
    print(f'Synchronization of PDB with {table_name} is good.')

