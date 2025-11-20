'''
example:
    - python src/pipelines/imgt_inn_blastp.py
functions:
    - collect antibody chain sequences of INN
    - blastp against pdb antibody
    - put alignment to table
'''
from collections import defaultdict
from dotenv import load_dotenv
import os
import sys
src_dir = os.path.dirname(os.path.dirname(__file__))
if src_dir not in sys.path:
    sys.path.append(src_dir)

from bioomics import Dir, QueryComplex, ProcessFasta, ProcessBlast, BuildComplex

def pull_data(args, meta):
    query = '''
        SELECT DISTINCT pro_id AS id, chain_seq AS seq
        FROM imgt_inn_chain C
        LEFT JOIN protein_source S
        ON C.chain_id = S.source_id
        WHERE S.source_name = 'imgt_inn_chain';
    '''
    rows = QueryComplex(args['verbose']).list_data(query)
    meta['num_records'] = len(rows)
    return rows

if __name__ == "__main__":
    meta = defaultdict(int)
    load_dotenv()
    args = {
        'verbose': True,
        'overwrite': True,
        'data_dir': os.getenv('data_dir'),
        'db_pdb_antibody': os.getenv("db_pdb_antibody"),
    }
    args['outdir'] = os.path.join(args['data_dir'], 'collect')
    Dir(args['outdir']).init_dir()

    # prepare query fasta
    name = 'inn_mabs'
    fa_file = os.path.join(args['outdir'], f'{name}.faa')
    if args['overwrite'] or not os.path.isfile(fa_file):
        rows = pull_data(args, meta)
        ProcessFasta(fa_file).row_seq(rows)

    # blastp, reference is pdb antibody
    outfmt = 5
    aln_file = os.path.join(args['outdir'], f'{name}_pdb_antibody.blastp{outfmt}')
    if args['overwrite'] or not os.path.isfile(aln_file):
        pb = ProcessBlast(args['outdir'], args['verbose'])
        pb.blastp(fa_file, args['db_pdb_antibody'], aln_file, outfmt)
    

    # insert aignment
    table_name = 'imgt_inn_blastp'
    bc = BuildComplex([aln_file,], args['verbose'])
    bc.empty_table(table_name)
    bc.insert_imgt_inn_blastp()

    print(meta)

