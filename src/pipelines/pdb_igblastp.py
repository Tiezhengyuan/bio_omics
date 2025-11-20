'''
requirments:
    - igblast is ready. run script in that directory
        cd /home/yuan/data/ncbi-igblast-1.22.0
    - igblastp database is ready
example: 
    - parallel -j4 python /home/yuan/bio/bio_omics/src/pipelines/pdb_igblastp.py ::: {0..46}
    - python /home/yuan/bio/bio_omics/src/pipelines/pdb_igblastp.py 46
functions:
    - igblastp against IMGT v-region
    - put meta data to table meta_pdb_igblastp
'''
from collections import defaultdict
from dotenv import load_dotenv
import os
import sys
src_dir = os.path.dirname(os.path.dirname(__file__))
if src_dir not in sys.path:
    sys.path.append(src_dir)

from bioomics import QueryComplex, ProcessIgblast, BuildComplex


def retrieve_query(pdb_group, args, meta):
    query = f"""
        SELECT F.pdb_id, F.faa_file AS query_fa
        FROM meta_faa F, meta_pdb P
        WHERE P.pdb_group = {pdb_group}
            AND P.pdb_id = F.pdb_id
            AND F.faa_file IS NOT NULL
            AND F.pdb_id NOT IN(
                SELECT pdb_id
                FROM meta_pdb_igblastp
            );
    """
    rows = QueryComplex(args['verbose']).list_data(query)
    meta['rows'] = len(rows)
    meta['pdb_group'] = int(pdb_group)
    return rows

def vregion_igblastp(rows, args, meta, outfmt:int=3):
    for row in rows:
        query_fa = row['query_fa']
        if os.path.isfile(query_fa):
            outdir = os.path.dirname(query_fa)
            aln_file = os.path.join(outdir, f'chain_sequences.vregion.igblastp{outfmt}')
            if args['overwrite'] or not os.path.isfile(aln_file):
                pb = ProcessIgblast(args['igblast_bin'], None, args['verbose'])
                status = pb.igblastp(query_fa, args['db_vregion'], aln_file, outfmt)
                meta[status] += 1
            else:
                meta['skip_igblastp'] += 1
        else:
            meta['no_query_files'] += 1
        if os.path.isfile(aln_file):
            yield row['pdb_id'], aln_file

if __name__ == "__main__":
    meta = defaultdict(int)
    load_dotenv()
    args = {
        'verbose': True,
        'overwrite': True,
        'chunk_size': 50,
        'db_vregion': os.getenv("db_vregion"),
        'igblast_bin': os.getenv("igblast_bin"),
    }
    pdb_group = sys.argv[1]
    rows = retrieve_query(pdb_group, args, meta)
    # igblastp agaist IMGT V-region
    aln_iter = vregion_igblastp(rows, args, meta)

    # put meta data to table meta_pdb_igblastp
    bc = BuildComplex(aln_iter, args['verbose'], args['chunk_size'])
    bc.insert_meta_pdb_igblastp()

    print(meta)
