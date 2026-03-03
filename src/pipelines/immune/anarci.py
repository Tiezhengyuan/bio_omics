'''
requirments:
    - chain_id
    - ANARCI is installed
example: 
    - parallel -j8 python app.py anarci ::: {0..117}
    - python app.py anarci 0
functions:
    - antibody numbering
    - put numbering output to ab_anarci
'''

from src.helper import *
from bioomics import Anarci

@header
def retrieve(args, meta):
    pdb_group = args['pdb_group']
    table_name = args['table_name']
    query = f"""
        SELECT F.pdb_id, F.chain_id, F.relative_path
        FROM meta_chain_faa F
        LEFT JOIN meta_pdb P ON F.pdb_id = P.pdb_id
        WHERE P.pdb_group = {pdb_group}
            AND F.relative_path IS NOT NULL
            AND F.chain_id NOT IN(
                SELECT chain_id FROM {table_name}
            )
    ;"""
    df = QueryComplex(args['verbose']).list_data(query, True)
    meta['rows'] = len(df)
    if len(df) > 0:
        g = df.groupby('pdb_id')
        return g

def run(g, args, meta):
    for pdb_id, sub in g:
        endpoint = os.path.join(args['table_name'], pdb_id[:2], pdb_id)
        outdir = os.path.join(args['data_dir'], endpoint)
        Dir(outdir).init_dir()
        a = Anarci(outdir, args['verbose'])

        # antibody numbering
        for _, row in sub.iterrows():        
            query_fa = os.path.join(args['data_dir'], row['relative_path'])
            if os.path.isfile(query_fa):
                status = a.run_anarci(query_fa, row['chain_id'])
                meta[status] += 1
            else:
                meta['no_query_files'] += 1

        # retrieve numbering results
        csv_files = [os.path.join(outdir, i) for i in os.listdir(outdir)]
        results = a.parse_numbering(csv_files)
        for rec in results:
            meta['record'] += 1
            yield (rec['Id'], rec['domain_no'], rec['hmm_species'], rec['chain_type'],
                rec['e-value'], rec['score'], rec['seqstart_index'], rec['seqend_index'],
                rec['identity_species'], rec['v_gene'], rec['v_identity'],
                rec['j_gene'], rec['j_identity'],)
        meta['pdb'] += 1

if __name__ == "__main__":
    args = {
        'pdb_group': sys.argv[1],
        'verbose': True,
        'overwrite': True,
        'chunk_size': 20,
        'data_dir': os.getenv('data_dir'),
        'table_name': 'ab_anarci',
    }
    g = retrieve(args, meta)
    record_iter = run(g, args, meta)
    bc = BuildComplex(record_iter, args['verbose'], args['chunk_size'])
    bc.insert_batch_records(args['table_name'])

    footer(meta)
