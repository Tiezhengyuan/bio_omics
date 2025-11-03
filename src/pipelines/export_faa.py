'''
requirements:
    - PDB data
    - table meta_pdb
example: 
    - parallel -j4 python /home/yuan/bio/bio_omics/src/pipelines/export_faa.py ::: {0..46}
    - python /home/yuan/bio/bio_omics/src/pipelines/export_faa.py 46
functions:
    - retrieve pdb to fasta
    - blastp against uniprot_sprot
'''
from collections import defaultdict
import os
import sys
src_dir = os.path.dirname(os.path.dirname(__file__))
if src_dir not in sys.path:
    sys.path.append(src_dir)

from bioomics import Pdb, QueryComplex, BuildComplex

def build_faa(args, rows, meta):
    for row in rows:
        p = Pdb(row['pdb_file'], args['outdir'], args['verbose'])
        p.load_structure()
        # export chain sequences in fasta to data/ in default
        db_record = p.export_fasta(meta, args['overwrite'])
        if db_record:
            meta['db_records'] += 1
            yield db_record

if __name__ == "__main__":
    meta = defaultdict(int)
    args = {
        'verbose': True,
        'overwrite': True,
        'outdir': '/home/yuan/data/pdb/data',
    }

    # retrieve pdb_ids
    pdb_group = sys.argv[1]
    query = f"""
        SELECT pdb_id, pdb_file
        FROM meta_pdb
        WHERE pdb_group = {pdb_group}
            AND pdb_id NOT IN (
                SELECT DISTINCT pdb_id
                FROM meta_faa
            );
    """
    rows = QueryComplex().list_data(query)
    meta['rows'] = len(rows)

    print(f"{len(rows)} PDB files ready for processing...")
    record_iter = build_faa(args, rows, meta)

    # insert
    bc = BuildComplex(record_iter, args['verbose'])
    bc.insert_meta_faa()

    print(meta)

