'''
requirements:
    - The table pdb_chain is ready.
example:
    - parallel -j4 python src/pipelines/put_chain_seq.py ::: {0..46}
    - python src/pipelines/put_chain_seq.py 46
function:
    - assign unique id to chain protein sequences
'''
from collections import defaultdict
import os
import sys
src_dir = os.path.dirname(os.path.dirname(__file__))
if src_dir not in sys.path:
    sys.path.append(src_dir)

from bioomics import QueryComplex,  BuildComplex


def retrieve_data(pdb_group, args, meta):
    query = f"""
        SELECT DISTINCT C.pdb_id, P.pdb_file
        FROM pdb_chainid C
        LEFT JOIN meta_pdb P
            ON C.pdb_id = P.pdb_id
        WHERE P.pdb_group = {pdb_group}
            AND C.chain_id NOT IN (
                SELECT chain_id
                FROM pdb_chain_seq
            );
    """
    rows = QueryComplex(args['verbose']).list_data(query)
    meta['pdb_group'] = int(pdb_group)
    meta['rows'] = len(rows)
    return rows
##################################
if __name__ == "__main__":
    pdb_group = sys.argv[1]
    meta = defaultdict(int)
    args = {
        'verbose': True,
        'chunk_size': 100,
    }
    rows = retrieve_data(pdb_group, args, meta)

    table_name = 'pdb_proseq'
    print(f"{pdb_group}. Try to insert PDB data into table {table_name}...")
    bc = BuildComplex(rows, args['verbose'], args['chunk_size'])
    bc.insert_pdb_chain_seq()
    print(meta)
    print(f'Synchronization of PDB with {table_name} is good.')
