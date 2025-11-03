'''
note:
    - PDB files should be synchronized to the local.
    - The table pdb is ready.
example:
    - parallel -j4 python src/pipelines/put_chain.py ::: {0..46}
    - python src/pipelines/put_chain.py 46
function:
    synchronize PDB/UniProt/alignments with data in database COMPLEX.
'''
from collections import defaultdict
import os
import sys
src_dir = os.path.dirname(os.path.dirname(__file__))
if src_dir not in sys.path:
    sys.path.append(src_dir)

from bioomics import QueryComplex,  BuildComplex


def retrieve_pdb(pdb_group, args):
    query = f"""
        SELECT pdb_id, pdb_file
        FROM meta_pdb
        WHERE pdb_group = {pdb_group}
            AND pdb_id NOT IN (
                SELECT DISTINCT pdb_id
                FROM pdb_chain
            );
    """
    rows = QueryComplex(args['verbose']).list_data(query)
    rows = [i['pdb_file'] for i in rows]
    return rows
##################################
if __name__ == "__main__":
    pdb_group = sys.argv[1]
    args = {
        'verbose': True,
        'pdb_dir': '/home/yuan/data/pdb/rsync',
        'chunk_size': 100,
    }
    meta = defaultdict(int)
    rows = retrieve_pdb(pdb_group, args)
    meta['rows'] = len(rows)

    if len(rows) > 0:
        print(f"Try to insert {len(rows)} PDB data into table pdb_chain...")
        bc = BuildComplex(rows, args['verbose'], args['chunk_size'])
        bc.insert_pdb_chain()
    else:
        print('Synchronization of PDB with pdb_chain is good.')
    print(pdb_group, meta)
