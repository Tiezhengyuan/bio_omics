'''
example: 
    - parallel -j8 python /home/yuan/bio/bio_omics/src/pipelines/absd_igblastp.py ::: {0..400}
    - python /home/yuan/bio/bio_omics/src/pipelines/absd_igblastp.py 0
functions:
    - igblastp against v-region
'''
from collections import defaultdict
from dotenv import load_dotenv
import os
import sys
src_dir = os.path.dirname(os.path.dirname(__file__))
if src_dir not in sys.path:
    sys.path.append(src_dir)

from bioomics import BuildComplex, ProcessIgblast, QueryComplex

def get_data(idx, args, meta):
    end = (idx + 1 ) * 10_000
    start = end - 10_000
    query = f"""
        SELECT seq_id , faa_file
        FROM absd_faa
        WHERE seq_id >= {start}
            AND seq_id < {end}
            AND seq_id NOT IN (
                SELECT seq_id
                FROM absd_igblastp
            );
    """
    rows= QueryComplex(args['verbose']).list_data(query)
    meta['idx'] = idx
    meta['rows'] = len(rows)
    return rows 

def align(rows, args, meta):
    for row in rows:
        query_fa = row['faa_file']
        # aln file
        file_name = os.path.splitext(os.path.basename(query_fa))[0]
        outdir = os.path.dirname(query_fa)
        outfmt = 3
        aln_file = os.path.join(outdir, f"{file_name}.igblastp{outfmt}")

        # alignment
        if args['overwrite'] or not os.path.isfile(aln_file):
            pb = ProcessIgblast(outdir, args['verbose'])
            status = pb.igblastp(query_fa, args['db_vregion'], aln_file, outfmt=outfmt)
            meta[status] += 1
        else:
            meta['skip'] += 1
        if os.path.isfile(aln_file):
            meta['record'] += 1
            yield row['seq_id'], aln_file

if __name__ == "__main__":
    meta= defaultdict(int)
    load_dotenv()
    args = {
        'verbose': True,
        'overwrite': False,
        'absd_dir': os.getenv('absd_dir'),
        'db_vregion': os.getenv("db_vregion"),
        'chunk_size': 50,
    }
    idx = int(sys.argv[1])
    rows = get_data(idx, args, meta)

    if len(rows) > 0:
        record_iter = align(rows, args, meta)
        # insertion
        table_name = 'absd_igblastp'
        bc = BuildComplex(record_iter, args['verbose'], args['chunk_size'])
        bc.insert_absd_igblastp()
    
    print(meta)