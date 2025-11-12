'''
example: 
    - parallel -j4 python src/pipelines/export_chain_faa.py ::: {0..46}
    - python src/pipelines/export_chain_faa.py 46
functions:
    - retrieve pdb to fasta by chains
    - blastp against uniprot_sprot
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
        SELECT pdb_id, pdb_file
        FROM meta_pdb
        WHERE pdb_group = {pdb_group}
            AND pdb_id NOT IN (
                SELECT DISTINCT P.pdb_id
                FROM pdb_chain P
                LEFT JOIN meta_chain_faa M
                ON M.chain_id = P.chain_id
                WHERE M.faa_file IS NOT NULL
                    OR P.chain_seq is NULL
            );
    """
    rows = QueryComplex(args['verbose']).list_data(query)
    print(f"{len(rows)} PDB files ready...")
    meta['pdb_group'] = int(pdb_group)
    meta['rows'] += len(rows)
    return rows

if __name__ == "__main__":
    meta = defaultdict(int)
    load_dotenv()
    args = {
        'verbose': True,
        'overwrite': False,
        'outdir': os.path.join(os.getenv('data_dir'), 'data'),
    }
    pdb_group = sys.argv[1]
    rows = retrieve_pdb(pdb_group, args, meta)

    for row in rows:
        p = Pdb(row['pdb_file'], args['outdir'], args['verbose'])
        p.load_structure()
        # export chain sequences in fasta to data/ in default
        db_records = p.export_chain_faa(meta, args['overwrite'])
        # insert
        meta['db_records'] += len(db_records)
        bc = BuildComplex(db_records, args['verbose'])
        bc.insert_meta_chain_faa()

    print(meta)
