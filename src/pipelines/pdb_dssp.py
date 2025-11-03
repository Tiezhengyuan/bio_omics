'''
requirements:
    - PDB raw data
    - table pdb, meta_chain_pdb
example: 
    - parallel -j4 python /home/yuan/bio/bio_omics/src/pipelines/pdb_dssp.py ::: {0..46}
    - python /home/yuan/bio/bio_omics/src/pipelines/pdb_dssp.py 46
functions:
    - Analyze secondary structure by chains
    - put dssp sequence to table pdb_dssp
'''
from collections import defaultdict
import os
import sys
src_dir = os.path.dirname(os.path.dirname(__file__))
if src_dir not in sys.path:
    sys.path.append(src_dir)

from bioomics import QueryComplex, BuildComplex


def search_pdb(pdb_group, meta):
    # 1. ==pdb_group, 2. pdb_id not in pdb_secondary 3. chain is protein
    query = f"""
        SELECT M.pdb_id, M.model_no, M.chain_combo as chain_no, M.pdb_file
        FROM meta_pdb P, meta_chain_pdb M
        WHERE P.pdb_group = {pdb_group}
            AND P.pdb_id = M.pdb_id
            AND LENGTH(M.chain_combo) = 1
            AND P.pdb_id NOT IN (
                SELECT DISTINCT pdb_id
                FROM pdb_secondary)
            AND P.pdb_id IN (
                SELECT DISTINCT B.pdb_id
                FROM meta_chain_faa A, pdb_chain B
                WHERE A.chain_id=B.chain_id
            );
    """
    rows = QueryComplex().list_data(query)
    meta['rows'] = len(rows)
    return rows


if __name__ == "__main__":
    args = {
        'verbose': True,
        'data_dir': '/home/yuan/data/pdb/data',
        'chunk_size': 20,
    }
    meta = defaultdict(int)
    pdb_group = sys.argv[1]

    print(f"Group {pdb_group}: Collect *.pdb of chains ...")
    rows = search_pdb(pdb_group, meta)
    # put bfactor into table pdb_bfactor
    if rows:
        print(f"Try to insert {len(rows)} data into pdb_secondary...")
        bc = BuildComplex(rows, args['verbose'], args['chunk_size'])
        bc.insert_pdb_secondary()
