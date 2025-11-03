'''
example: 
    - parallel -j4 python /home/yuan/bio/bio_omics/src/pipelines/export_chain_faa.py ::: {0..46}
    - python /home/yuan/bio/bio_omics/src/pipelines/export_chain_faa.py 46
functions:
    - retrieve pdb to fasta by chains
    - blastp against uniprot_sprot
'''
from collections import defaultdict
import os
import sys
src_dir = os.path.dirname(os.path.dirname(__file__))
if src_dir not in sys.path:
    sys.path.append(src_dir)

from bioomics import Pdb, QueryComplex, BuildComplex

if __name__ == "__main__":
    meta = defaultdict(int)
    args = {
        'verbose': True,
        'overwrite': False,
        'outdir': '/home/yuan/data/pdb/data',
    }

    pdb_group = sys.argv[1]
    query = f"""
        SELECT pdb_id, pdb_file
        FROM meta_pdb
        WHERE pdb_group = {pdb_group}
            AND pdb_id NOT IN (
                SELECT DISTINCT P.pdb_id
                FROM meta_chain_faa M, pdb_chain P
                WHERE M.chain_id = P.chain_id
            );
    """
    rows = QueryComplex().list_data(query)

    print(f"{len(rows)} PDB files ready...")
    meta['rows'] += len(rows)
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
