#!/usr/bin/python
'''
requirements:
    - pdb fixed: table fixed_chain
example:
    - parallel -j8 python app.py chain_af ::: {0..237}
    - python app.py chain_af 0
function:
    - prefix monomer
    - insert data into table chain_alphafold
'''

from src.helper import *
from bioomics import Predictor

@header
def retrieve_pdb(params, meta):
    pdb_group = params['pdb_group']
    table_name = params['table_name']
    query = f"""
        WITH tmp AS(
            SELECT DISTINCT I.pdb_id, S.first_chain_id, F.relative_path
                FROM pdb_chain_seq       S
                LEFT JOIN pdb_chainid    I ON S.first_chain_id = I.chain_id
                LEFT JOIN meta_chain_faa F ON S.first_chain_id = F.chain_id
                WHERE F.relative_path IS NOT NULL
                    AND S.first_chain_id IS NOT NULL
                    AND S.first_chain_id IN (
                        SELECT chain_id from view_antibody
                    )
                    AND S.first_chain_id NOT IN (
                        SELECT chain_id FROM {table_name}
                    )  
        )
        SELECT T.pdb_id, T.first_chain_id AS chain_id, T.relative_path
        FROM tmp T
        LEFT JOIN meta_pdb P ON T.pdb_id=P.pdb_id
        WHERE P.pdb_group = {pdb_group}
            AND P.relative_path IS NOT NULL
    ;"""
    rows = QueryComplex(params['verbose']).list_data(query)
    print(f"{len(rows)} PDB files ready...")
    meta['rows'] += len(rows)
    return rows

if __name__ == "__main__":
    params.update({
        'predict_dir': os.getenv('predict_dir'),
        'table_name': 'chain_alphafold',
        'num_recycle': '3',
        'relax': True,
        'max_iterations': '2000',
    })
    add_params(params)

    rows = retrieve_pdb(params, meta)
    host = socket.gethostname()
    for row in rows:
        t0 = datetime.now()
        pdb_id, chain_id = row['pdb_id'], row['chain_id']
        endpoint = os.path.join(params['table_name'], pdb_id[:2], pdb_id, chain_id)
        outdir = os.path.join(params['predict_dir'], endpoint)
        Dir(outdir).init_dir()

        # predict
        _params = {
            'faa_file': os.path.join(params['data_dir'], row['relative_path']),
            'output_dir': outdir,
        }
        _params.update(params)
        if params['verbose']:
            print('Arguments for prediction: ', _params)
        Predictor(_params).colabfold()

        # build record
        done_file = os.path.join(outdir, f"{chain_id}.done.txt")
        if os.path.isfile(done_file):
            meta['succeed'] += 1
        else:
            meta['fail'] += 1
            endpoint = None
        rec = [pdb_id, chain_id, endpoint, host, t0]

        # insertion
        if params['verbose']:
            print('### Try to insert ', rec)
        bc = BuildComplex([tuple(rec),], params['verbose'])
        bc.insert_batch_records(params['table_name'])

        # delete *_env for saving space
        if endpoint:
            _env_dir = os.path.join(outdir, f"{chain_id}_env")
            if os.path.isdir(_env_dir):
                shutil.rmtree(_env_dir)
            #copy results to local server if that is remote
            copy_data(outdir)

    footer(meta, params)