'''
requirements:
    - local */chain_sequences.faa are ready
    - table meta_faa, meta_pdb
    - blastp database
example: 
    - parallel -j4 python app.py blastp_uniprot ::: {0..237}
    - python app.py blastp_uniprot 0
functions:
    - blastp against uniprot_sprot
    - put aln_file into meta_blastp5
note:
    - only consider model.id=0
'''

from src.helper import *
from bioomics import ProcessBlast

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
            AND F.pdb_id not in (
                SELECT distinct pdb_id FROM {table_name}
            );
    """
    rows = QueryComplex(args['verbose']).list_data(query)
    meta['rows'] = len(rows)
    return rows

# blastp of chain sequences against uniprot_sprot
def run(rows, args, meta):
    outfmt = args['outfmt']
    for row in rows:
        pdb_id = row['pdb_id']
        model_no = row['model_no']
        endpoint = os.path.join(args['table_name'], pdb_id[:2], pdb_id)
        outfile_name = f"{pdb_id}_{model_no}.uniprot.blastp{outfmt}"
        
        query_fa = os.path.join(args['data_dir'], row['relative_faa'])
        if os.path.isfile(query_fa):
            outdir = os.path.join(args['data_dir'], endpoint)
            aln_file = os.path.join(outdir, outfile_name)
            if args['overwrite'] or not os.path.isfile(aln_file):
                Dir(outdir).init_dir()
                pb = ProcessBlast(verbose=args['verbose'])
                status = pb.blastp(query_fa, args['db_uniprot'], aln_file, outfmt)
                meta[status] += 1
            else:
                meta['skip_alignment'] += 1
        else:
            meta['no_query_faa'] += 1
        
        # build record
        if os.path.isfile(aln_file):
            meta['aln_file'] += 1
            relative_aln = os.path.join(endpoint, outfile_name)
            yield pdb_id, model_no, relative_aln
        else:
            yield pdb_id, model_no, None

if __name__ == "__main__":
    args = {
        'pdb_group': sys.argv[1],
        'verbose': True,
        'overwrite': True,
        'data_dir': os.getenv('data_dir'),
        'db_uniprot': os.getenv('db_uniprot'),
        'chunk_size': 10,
        'table_name': 'meta_blastp5',
        'outfmt': 5,
    }
    # collect query
    rows = retrieve(args, meta)

    # blastp alignment
    record_iter = run(rows, args, meta)

    #insert aln_file to table meta_blastp5
    bc = BuildComplex(record_iter, args['verbose'], args['chunk_size'])
    bc.insert_batch_records(args['table_name'])

    footer(meta)
