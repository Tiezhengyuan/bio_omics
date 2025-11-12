'''
requirements:
    - local */chain_sequences.faa are ready
    - table meta_faa, meta_pdb
    - blastp database
example: 
    - parallel -j4 python src/pipelines/pdb_blastp_uniprot.py ::: {0..46}
    - python src/pipelines/pdb_blastp_uniprot.py 46
functions:
    - blastp against uniprot_sprot
    - put aln_file into meta_blastp5
'''
from collections import defaultdict
from dotenv import load_dotenv
import os
import sys
src_dir = os.path.dirname(os.path.dirname(__file__))
if src_dir not in sys.path:
    sys.path.append(src_dir)

from bioomics import ProcessBlast, QueryComplex, BuildComplex


def query_faa_iter(pdb_group, args, meta):
    query = f"""
        SELECT F.pdb_id, F.faa_file
        FROM meta_faa F, meta_pdb P
        WHERE P.pdb_group = {pdb_group}
            AND F.pdb_id = P.pdb_id
            AND F.faa_file IS NOT NULL
            AND F.pdb_id not in (
                SELECT distinct pdb_id
                FROM meta_blastp5
            );
    """
    rows = QueryComplex(args['verbose']).list_data(query)
    meta['pdb_group'] = int(pdb_group)
    meta['rows'] = len(rows)
    for row in rows:
        row.update(args)
        yield row

# blastp of chain sequences against uniprot_sprot
def uniprot_blastp(args_pool, meta, outfmt=5):
    for args in args_pool:
        query_fa = args['faa_file']
        cascade_dir = os.path.dirname(query_fa)
        aln_file = os.path.join(cascade_dir, f'chain_sequences.uniprot.blastp{outfmt}')
        if os.path.isfile(query_fa):
            if args['overwrite'] or not os.path.isfile(aln_file):
                pb = ProcessBlast(verbose=args['verbose'])
                status = pb.blastp(query_fa, args['db_uniprot'], aln_file, outfmt)
                meta[status] += 1
            else:
                meta['skip_alignment'] += 1
        else:
            meta['no_query_faa'] += 1
        if os.path.isfile(aln_file):
            meta['aln_file'] += 1
            yield args['pdb_id'], aln_file

if __name__ == "__main__":
    meta = defaultdict(int)
    load_dotenv()
    args = {
        'verbose': True,
        'overwrite': False,
        'db_uniprot': os.getenv('db_uniprot'),
        'chunk_size': 10,
    }
    # collect query
    pdb_group = sys.argv[1]
    args_pool = query_faa_iter(pdb_group, args, meta)

    # blastp alignment
    aln5_pool = uniprot_blastp(args_pool, meta, outfmt=5)

    #insert aln_file to table meta_blastp5
    bc = BuildComplex(aln5_pool, args['verbose'], args['chunk_size'])
    bc.insert_meta_blastp5()

    print(meta)
