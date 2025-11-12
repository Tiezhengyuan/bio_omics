'''
requirements:
    - bfactor.pkl
    - table meta_bfactor
example: 
    - parallel -j4 python src/pipelines/put_bfactor.py ::: {0..46}
    - python src/pipelines/put_bfactor.py 46
functions:
    - put average bfactor to table pdb_bfactor
'''
from collections import defaultdict
import os
import sys
src_dir = os.path.dirname(os.path.dirname(__file__))
if src_dir not in sys.path:
    sys.path.append(src_dir)

from bioomics import QueryComplex, BuildComplex

def retrieve_pdb(pdb_group, args, meta):
    query = f"""
        SELECT B.pdb_id, B.pickle_file
        FROM meta_bfactor B, meta_pdb P
        WHERE B.pdb_id=P.pdb_id
            AND P.pdb_group = {pdb_group}
            AND B.pickle_file IS NOT NULL
            AND B.pdb_id NOT IN (
                SELECT DISTINCT pdb_id
                FROM pdb_bfactor
            );
    """
    # each record must have pdb_id and pdb_file
    rows = QueryComplex(args['verbose']).list_data(query)
    print(f"{len(rows)} PDB files ready...")
    meta['pdb_group'] = int(pdb_group)
    meta['rows'] += len(rows)
    return rows


if __name__ == "__main__":
    meta = defaultdict(int)
    args = {
        'verbose': True,
        'chunk_size': 100,
    }

    # retrieve pdb_ids
    pdb_group = sys.argv[1]
    args_iter = retrieve_pdb(pdb_group, args, meta)

    # put bfactor into table meta_bfactor
    print(f"Try put meta data to meta_bfactor...")
    bc = BuildComplex(args_iter, args['verbose'], args['chunk_size'])
    bc.insert_pdb_bfactor()

    print(meta)

