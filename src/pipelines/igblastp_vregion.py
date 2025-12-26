'''
requirments:
    - igblast is ready. run script in that directory
        cd /home/yuan/data/ncbi-igblast-1.22.0
    - igblastp database is ready
example: 
    - parallel -j4 python /home/yuan/bio/bio_omics/src/pipelines/igblastp_vregion.py ::: {0..116}
    - python /home/yuan/bio/bio_omics/src/pipelines/igblastp_vregion.py 0
functions:
    - igblastp against IMGT v-region
    - put meta data to table meta_igblastp
note:
    - Execution: Step into igblastp bin, and use absolute path
    - only consider model.id==0
'''

from helper import *
from bioomics import ProcessIgblast

def retrieve(pdb_group, args, meta):
    table_name = args['table_name']
    query = f"""
        SELECT F.pdb_id, F.model_no, F.relative_faa
        FROM pdb_faa F
        LEFT JOIN meta_pdb P ON F.pdb_id = P.pdb_id
        WHERE P.pdb_group = {pdb_group}
            AND F.relative_faa IS NOT NULL
            AND F.model_no = 0
            AND F.pdb_id NOT IN(
                SELECT pdb_id FROM {table_name}
            )
    ;"""
    rows = QueryComplex(args['verbose']).list_data(query)
    meta['rows'] = len(rows)
    meta['pdb_group'] = int(pdb_group)
    return rows

def run(rows, args, meta):
    outfmt = args['outfmt']
    for row in rows:
        pdb_id = row['pdb_id']
        model_no = row['model_no']
        endpoint = os.path.join(args['table_name'], pdb_id[:2], pdb_id)
        outfile_name = f"{pdb_id}_{model_no}.vregion.igblastp{outfmt}"
        
        query_fa = os.path.join(args['data_dir'], row['relative_faa'])
        if os.path.isfile(query_fa):
            outdir = os.path.join(args['data_dir'], endpoint)
            Dir(outdir).init_dir()

            aln_file = os.path.join(outdir, outfile_name)
            if args['overwrite'] or not os.path.isfile(aln_file):
                pb = ProcessIgblast(args['igblast_bin'], None, args['verbose'])
                status = pb.igblastp(query_fa, args['db_vregion'], aln_file, outfmt)
                meta[status] += 1
            else:
                meta['skip_igblastp'] += 1
        else:
            meta['no_query_files'] += 1

        # build record
        if os.path.isfile(aln_file):
            relative_aln = os.path.join(endpoint, outfile_name)
            yield pdb_id, model_no, relative_aln
        else:
            yield pdb_id, model_no, None

if __name__ == "__main__":
    args = {
        'verbose': True,
        'overwrite': True,
        'chunk_size': 50,
        'data_dir': os.getenv('data_dir'),
        'db_vregion': os.getenv("db_vregion"),
        'outfmt': 3,
        'igblast_bin': os.getenv("igblast_bin"),
        'table_name': 'meta_igblastp',
    }
    pdb_group = sys.argv[1]
    rows = retrieve(pdb_group, args, meta)
    # igblastp agaist IMGT V-region
    record_iter = run(rows, args, meta)

    # put meta data to table meta_pdb_igblastp
    bc = BuildComplex(record_iter, args['verbose'], args['chunk_size'])
    bc.insert_batch_records(args['table_name'])

    print(meta)
