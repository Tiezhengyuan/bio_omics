'''
example:
    - python src/pipelines/ab_blastp_db.py
functions:
    - collect antibody chain sequences
    - build igblastp database as reference
'''
from collections import defaultdict
from dotenv import load_dotenv
import os
import sys
src_dir = os.path.dirname(os.path.dirname(__file__))
if src_dir not in sys.path:
    sys.path.append(src_dir)

from bioomics import Dir, QueryComplex, ProcessFasta, ProcessBlast

def pull_data(args, meta):
    query = '''
        SELECT DISTINCT pro_id AS id, chain_seq AS seq
        FROM view_antibody
        where pro_id IS NOT NULL
            AND chain_seq IS NOT NULL;
    '''
    rows = QueryComplex(args['verbose']).list_data(query)
    meta['num_rows'] = len(rows)
    return rows

if __name__ == "__main__":
    meta = defaultdict(int)
    load_dotenv()
    args = {
        'verbose': True,
        'overwrite': True,
        'data_dir': os.getenv('data_dir'),
        'blast_bin_dir': os.getenv('blast_bin_dir'),
    }
    args['outdir'] = os.path.join(args['data_dir'], 'blastp_db')
    Dir(args['outdir']).init_dir()

    # prepare fasta
    name = 'pdb_chain_antibody'
    fa_file = os.path.join(args['outdir'], f'{name}.faa')
    if args['overwrite'] or not os.path.isfile(fa_file):
        rows = pull_data(args, meta)
        ProcessFasta(fa_file, args['verbose']).row_seq(rows)

    # build igblastp database
    pb = ProcessBlast(args['outdir'], args['verbose'])
    pb.build_blastp_db(fa_file)

    print(meta)

