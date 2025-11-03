#!/usr/bin/python
'''
requirements:
    - table meta_faa
    - igblastp ready
example:
    - parallel -j4 python src/pipelines/put_align_imgt.py ::: {0..46}
    - python src/pipelines/put_align_imgt.py 46
function:
    - build table pdb and insert data into table align_imgt
'''
from collections import defaultdict
import os
import sys
src_dir = os.path.dirname(os.path.dirname(__file__))
if src_dir not in sys.path:
    sys.path.append(src_dir)

from bioomics import BuildComplex, QueryComplex

def retrieve_data(pdb_group, args, meta):
    query = f"""
        SELECT F.pdb_id, F.faa_file, I.aln_file 
        FROM meta_pdb P, meta_faa F, meta_pdb_igblastp I
        WHERE P.pdb_group = {pdb_group}
            AND P.pdb_id = F.pdb_id
            AND P.pdb_id = I.pdb_id
            AND F.pdb_id NOT IN(
                SELECT DISTINCT pdb_id
                FROM align_imgt
            );
    """
    rows = QueryComplex(args['verbose']).list_data(query)
    meta['rows'] = len(rows)
    return rows

if __name__ == "__main__":
    pdb_group = sys.argv[1]
    meta = defaultdict(int)
    args = {
        'verbose': True,
    }
    rows = retrieve_data(pdb_group, args, meta)

    table_name = 'align_imgt'
    print(f"Try to insert data into {table_name}")
    bc = BuildComplex(rows, args['verbose'])
    bc.insert_align_imgt()

    print(f"PDB group: {pdb_group}")
    print(meta)

