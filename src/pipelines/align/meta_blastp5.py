'''
requirements:
    - local */chain_sequences.faa are ready
    - table meta_faa, meta_pdb
    - blastp database
example: 
    - parallel -j8 bomics meta_blastp5 ::: {0..241}
    - bomics meta_blastp5 0
functions:
    - blastp against uniprot_sprot
    - put aln_file into meta_blastp5
note:
    - only consider model.id=0
'''

from src.helper import *
from bioomics import ProcessBlast

@header
def retrieve(params, meta):
    pdb_group = params['pdb_group']
    table_name = params['table_name']
    query = f"""
        SELECT F.pdb_id, F.model_no, F.relative_faa
        FROM pdb_faa F
        LEFT JOIN meta_pdb P ON F.pdb_id = P.pdb_id
        WHERE P.pdb_group = {pdb_group}
            AND F.relative_faa IS NOT NULL
            AND F.model_no = 0
            AND F.pdb_id not in (
                SELECT pdb_id FROM {table_name}
            );
    """
    rows = QueryComplex(params['verbose']).list_data(query)
    meta['rows'] = len(rows)
    return rows

# blastp of chain sequences against uniprot_sprot
def run(rows, params, meta):
    outfmt = params['outfmt']
    for row in rows:
        pdb_id = row['pdb_id']
        model_no = row['model_no']
        endpoint = Path(params['table_name']) / pdb_id[:2] / pdb_id
        outfile_name = f"{pdb_id}_{model_no}.uniprot.blastp{outfmt}"
        
        query_fa = params['output_pdb_dir'] / row['relative_faa']
        if query_fa.is_file():
            outdir = params['output_pdb_dir'] / endpoint
            aln_file = outdir / outfile_name
            if params['overwrite'] or not aln_file.is_file():
                outdir.mkdir(parents=True, exist_ok=True)
                pb = ProcessBlast(verbose=params['verbose'])
                status = pb.blastp(query_fa, params['db_uniprot'], aln_file, outfmt)
                meta[status] += 1
            else:
                meta['skip_alignment'] += 1
        else:
            meta['no_query_faa'] += 1
        
        # build record
        relative_aln = None
        if aln_file.is_file():
            meta['aln_file'] += 1
            relative_aln = str(endpoint / outfile_name)
        yield pdb_id, model_no, relative_aln

if __name__ == "__main__":
    params.update({
        'overwrite': False,
        'db_uniprot': os.getenv('db_uniprot'),
        'chunk_size': 10,
        'table_name': 'meta_blastp5',
        'outfmt': 5,
    })

    # collect query
    rows = retrieve(params, meta)
    # blastp alignment
    record_iter = run(rows, params, meta)

    #insert aln_file to table meta_blastp5
    bc = BuildComplex(params['verbose'], params['chunk_size'])
    bc.insert_batch_records(record_iter, params['table_name'])

    footer(meta)
