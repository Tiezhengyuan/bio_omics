'''
requirments:
    - aligner igblast is ready. run script in that directory
        /home/yuan/bio/ncbi-igblast-1.22.0
    - aligner blastp is ready.
    - database for igblastp and blastp alignment is ready
example: 
    - parallel -j4 python app.py meta_vdj ::: {0..117}
    - python app.py meta_vdj 0
functions:
    - igblastp against IMGT v-region
    - blastp against IMGT d/j-region
    - put meta data to table meta_vdj
note:
    - Execution: Step into igblastp bin, and use absolute path
    - only consider model.id==0
'''

from src.helper import *
from bioomics import ProcessIgblast, ProcessBlast

@header
def retrieve(args, meta):
    pdb_group = args['pdb_group']
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
    return rows

def run(rows, args, meta):
    for row in rows:
        pdb_id = row['pdb_id']
        query_fa = os.path.join(args['data_dir'], row['relative_faa'])
        if not os.path.isfile(query_fa):
            meta['no_query_files'] += 1
            continue

        model_no = row['model_no']
        rec = [pdb_id, model_no, None, None, None,]
        endpoint = os.path.join(args['table_name'], pdb_id[:2], pdb_id)
        outdir = os.path.join(args['data_dir'], endpoint)

        # v-region using igblastp
        outfmt = 3
        outfile_name = f"{pdb_id}_{model_no}.vregion.igblastp{outfmt}"
        aln_file = os.path.join(outdir, outfile_name)
        if args['overwrite'] or not os.path.isfile(aln_file):
            c = ProcessIgblast(args['igblast_bin'], None, args['verbose'])
            status = c.igblastp(query_fa, args['db_vregion'], aln_file, outfmt)
            meta[f'vregion_{status}'] += 1
        else:
            meta['vregion_skip'] += 1
        if os.path.isfile(aln_file):
            rec[2] = os.path.join(endpoint, outfile_name)

        # d-region using blastp
        outfmt = 5
        outfile_name = f"{pdb_id}_{model_no}.dregion.blastp{outfmt}"
        aln_file = os.path.join(outdir, outfile_name)
        if args['overwrite'] or not os.path.isfile(aln_file):
            c = ProcessBlast(None, args['verbose'])
            status = c.blastp(query_fa, args['db_dregion'], aln_file, outfmt)
            meta[f'dregion_{status}'] += 1
        else:
            meta['dregion_skip'] += 1
        if os.path.isfile(aln_file):
            rec[3] = os.path.join(endpoint, outfile_name)

        # j-region using blastp
        outfmt = 5
        outfile_name = f"{pdb_id}_{model_no}.jregion.blastp{outfmt}"
        aln_file = os.path.join(outdir, outfile_name)
        if args['overwrite'] or not os.path.isfile(aln_file):
            c = ProcessBlast(None, args['verbose'])
            status = c.blastp(query_fa, args['db_jregion'], aln_file, outfmt)
            meta[f'jregion_{status}'] += 1
        else:
            meta['jregion_skip'] += 1
        if os.path.isfile(aln_file):
            rec[4] = os.path.join(endpoint, outfile_name)

        yield tuple(rec)

if __name__ == "__main__":
    args = {
        'pdb_group': sys.argv[1],
        'verbose': True,
        'overwrite': False,
        'chunk_size': 50,
        'data_dir': os.getenv('data_dir'),
        'igblast_bin': os.getenv("igblast_bin"),
        'db_vregion': os.getenv("db_vregion"),
        'db_dregion': os.getenv("db_dregion"),
        'db_jregion': os.getenv("db_jregion"),
        'table_name': 'meta_vdj',
    }
    rows = retrieve(args, meta)
    # igblastp agaist IMGT V/D/J-region
    record_iter = run(rows, args, meta)

    # put meta data to table meta_pdb_igblastp
    bc = BuildComplex(record_iter, args['verbose'], args['chunk_size'])
    bc.insert_batch_records(args['table_name'])

    print(meta)
