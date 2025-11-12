'''
example: 
    - parallel -j4 python src/pipelines/export_pair_pdb.py ::: {0..46}
    - python src/pipelines/export_pair_pdb.py 46
functions:
    - retrieve pdb to pdb by chains
    - for B-factor, DSSP analysis
'''
from dotenv import load_dotenv
from collections import defaultdict
import os
import sys
src_dir = os.path.dirname(os.path.dirname(__file__))
if src_dir not in sys.path:
    sys.path.append(src_dir)

from bioomics import Pdb, QueryComplex, BuildComplex

def retrieve_pdb(pdb_group, args, meta):
    query = f"""
        SELECT DISTINCT P.pdb_id, M.pdb_file
        FROM pdb_chain P 
        LEFT JOIN meta_pdb M
        ON P.pdb_id = M.pdb_id
        WHERE M.pdb_group = {pdb_group}
            AND P.chain_seq IS NOT NULL
            AND P.pdb_id NOT IN (
                SELECT pdb_id
                FROM meta_combo_pdb
            );
    """
    rows = QueryComplex(args['verbose']).list_data(query)
    print(f"{len(rows)} PDB files ready for processing...")
    meta['pdb_group'] = int(pdb_group)
    meta['rows'] += len(rows)
    return rows

if __name__ == "__main__":
    meta = defaultdict(int)
    load_dotenv()
    args = {
        'verbose': True,
        'overwrite': False,
        'outdir': os.path.join(os.getenv("data_dir"), 'data'),
        'chunk_size': 10,
    }

    pdb_group = sys.argv[1]
    rows = retrieve_pdb(pdb_group, args, meta)

    for row in rows:
        p = Pdb(row['pdb_file'], args['outdir'], args['verbose'])
        p.load_structure()
        # export chain sequences in fasta to data/ in default
        db_records = p.export_combo_pdb(combo=2, meta=meta, overwrite=args['overwrite'])
        if db_records:
            # insert
            meta['db_records'] += len(db_records)
            bc = BuildComplex(db_records, args['verbose'], args['chunk_size'])
            bc.insert_meta_combo_pdb()

    print(meta)
