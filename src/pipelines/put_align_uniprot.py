'''
requirements:
    - blastp alignment *.blastp5 files
    - table meta_blastp5
example: 
    - parallel -j4 python src/pipelines/put_align_uniprot.py ::: {0..46}
    - python src/pipelines/put_align_uniprot.py 46
functions:
    - pust alignment results of blastp into align_uniprot
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
        SELECT B.pdb_id, B.aln_file
        FROM meta_blastp5 B, meta_pdb P
        WHERE P.pdb_group = {pdb_group}
            AND B.pdb_id = P.pdb_id
            AND B.pdb_id NOT IN (
            SELECT DISTINCT pdb_id
            FROM align_uniprot
        );
    """
    rows = QueryComplex(args['verbose']).list_data(query)
    meta['pdb_group'] = int(pdb_group)
    meta['align_uniprot'] = len(rows)
    return rows

if __name__ == "__main__":
    args = {
        'verbose': True,
    }
    meta = defaultdict(int)

    pdb_group = sys.argv[1]
    rows = retrieve_pdb(pdb_group, args, meta)

    # insert alignment into table align_uniprot
    BuildComplex(rows, args['verbose']).insert_align_uniprot()

    print(meta)
