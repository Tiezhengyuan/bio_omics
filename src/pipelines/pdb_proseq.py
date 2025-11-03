'''
requirements:
    - The table pdb_chain is ready.
example:
    - python src/pipelines/pdb_proseq.py
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


def retrieve_proseq(args):
    qc = QueryComplex(args['verbose'])

    # set chunk_no
    chunk_no, curr_num = 0, 0
    query = """
        SELECT chunk_no, COUNT(*) AS count
        FROM pdb_proseq
        GROUP BY chunk_no ORDER BY chunk_no DESC;
    """
    rows = qc.list_data(query)
    if rows:
        chunk_no = int(rows[0]['chunk_no'])
        curr_num = int(rows[0]['count'])
    print(chunk_no, curr_num)

    query = f"""
        SELECT DISTINCT chain_seq
        FROM pdb_chain
        WHERE chain_seq NOT IN (
                SELECT pro_seq
                FROM pdb_proseq
            );
    """
    rows = qc.list_data(query)
    meta['rows'] = len(rows)

    # build records
    records = []
    for row in rows:
        if curr_num < args['chunk_size']:
            curr_num += 1
        else:
            chunk_no += 1
            curr_num = 1
        rec = (chunk_no, row['chain_seq'])
        records.append(rec)
    return records
##################################
if __name__ == "__main__":
    args = {
        'verbose': True,
        'chunk_size': 100,
    }
    meta = defaultdict(int)
    rows = retrieve_proseq(args)

    table_name = 'pdb_proseq'
    if len(rows) > 0:
        print(f"Try to insert {len(rows)} PDB data into table {table_name}...")
        bc = BuildComplex(rows, args['verbose'])
        bc.insert_pdb_proseq()
    else:
        print(f'Synchronization of PDB with {table_name} is good.')
    print(meta)
