'''
example: 
    - parallel -j4 python src/pipelines/export_chain_pdb.py ::: {0..46}
    - python src/pipelines/export_chain_pdb.py 46
functions:
    - retrieve pdb to pdb by chain
    - data are required by B-factor, DSSP analysis
'''
from dotenv import load_dotenv
from itertools import islice
from collections import defaultdict
import os
import sys
src_dir = os.path.dirname(os.path.dirname(__file__))
if src_dir not in sys.path:
    sys.path.append(src_dir)

from bioomics import Pdb, QueryComplex, BuildComplex

def chunk_iterable(iterable, chunk_size):
    it = iter(iterable)
    while True:
        chunk = list(islice(it, chunk_size))
        if not chunk:
            return None
        yield chunk

def retrieve_pdb(pdb_group, args, meta):
    query = f"""
        SELECT DISTINCT P.pdb_id, M.pdb_file
        FROM pdb_chain P
        LEFT JOIN meta_pdb M
        ON P.pdb_id = M.pdb_id
        WHERE M.pdb_group = {pdb_group}
            AND P.chain_seq IS NOT NULL
            AND P.model_no <= 10
            AND P.chain_id NOT IN (
                SELECT chain_id
                FROM meta_chain_pdb
            );
    """
    rows = QueryComplex(args['verbose']).list_data(query)
    meta['pdb_group'] = int(pdb_group)
    meta['rows'] += len(rows)
    print(f"{len(rows)} PDB files ready for processing...")
    for chunk_rows in chunk_iterable(rows, 20):
        yield chunk_rows

def pull_chain_pdb(rows, args, meta):
    records = []
    for row in rows:
        print(row)
        p = Pdb(row['pdb_file'], args['data_dir'], args['verbose'])
        p.load_structure()
        # export chain sequences in fasta to data/ in default
        _records = p.export_chain_pdb(meta, args['overwrite'])
        records += _records
    meta['db_records'] += len(records)
    return records

if __name__ == "__main__":
    meta = defaultdict(int)
    load_dotenv()
    args = {
        'verbose': True,
        'overwrite': True,
        'data_dir': os.getenv("data_dir"),
    }
    pdb_group = sys.argv[1]
    rows_iter = retrieve_pdb(pdb_group, args, meta)

    for rows in rows_iter:
        records = pull_chain_pdb(rows, args, meta)
        # insert
        bc = BuildComplex(records, args['verbose'])
        bc.insert_meta_chain_pdb()

    print(meta)
